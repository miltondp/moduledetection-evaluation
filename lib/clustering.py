from collections import defaultdict

import numpy as np
import pandas as pd
import sklearn.cluster

import os
import subprocess as sp

from modulecontainers import Module, Modules
from simdist import simdist

# rpy2 imports for R-based methods (CLAMP)
try:
    import rpy2.robjects as ro
    import rpy2.robjects.pandas2ri
    import rpy2.robjects.numpy2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.conversion import localconverter

    # Try to activate converters (method may vary by rpy2 version)
    try:
        pandas2ri.activate()
    except Exception:
        pass  # May be deprecated, will use localconverter instead
    try:
        numpy2ri.activate()
    except Exception:
        pass
except ImportError:
    pass  # rpy2 not available

def standardize(X):
    return (X - X.mean())/(X.std())

def dummy(E, n=10, **kwargs):
    labels = np.random.randint(0, n, len(E.columns))
    modules = convert_labels2modules(labels, E.columns)
    return modules

def agglom(E, k=100, linkage="complete", simdist_function="pearson_correlation", **kwargs):
    distances = simdist(E, simdist_function, similarity = False)
    agglom = sklearn.cluster.AgglomerativeClustering(n_clusters=int(k), affinity = "precomputed", linkage = linkage)
    agglom.fit(distances)
    modules = convert_labels2modules(agglom.labels_, E.columns)
    return modules

def agglom_pearson_abs(E, k=100, linkage="complete", simdist_function="pearson_correlation_absolute", **kwargs):
    distances = simdist(E, simdist_function, similarity = False)
    agglom = sklearn.cluster.AgglomerativeClustering(n_clusters=int(k), affinity = "precomputed", linkage = linkage)
    agglom.fit(distances)
    modules = convert_labels2modules(agglom.labels_, E.columns)
    return modules

# Commented out clustermatch methods - not needed for CLAMP testing
# def agglom_clustermatch(E, k=100, linkage="complete", simdist_function="clustermatch", **kwargs):
#     distances = simdist(E, simdist_function, similarity = False)
#     agglom = sklearn.cluster.AgglomerativeClustering(n_clusters=int(k), affinity = "precomputed", linkage = linkage)
#     agglom.fit(distances)
#     modules = convert_labels2modules(agglom.labels_, E.columns)
#     return modules
#
# def agglom_clustermatch_linear(E, k=100, linkage="complete", simdist_function="clustermatch_linear", **kwargs):
#     distances = simdist(E, simdist_function, similarity = False)
#     agglom = sklearn.cluster.AgglomerativeClustering(n_clusters=int(k), affinity = "precomputed", linkage = linkage)
#     agglom.fit(distances)
#     modules = convert_labels2modules(agglom.labels_, E.columns)
#     return modules

def ica_zscore(E, k=200, stdcutoff=1e-3, seed=None, **kwargs):
    source = _ica_fastica(E, k, seed)
    modules = _ica_zscore(E, source, stdcutoff)

    return modules

def _ica_fastica(E, k, seed=None):
    ica = sklearn.decomposition.FastICA(n_components=int(k), random_state=seed)
    source = ica.fit_transform(standardize(E).T)

    return source

def _ica_zscore(E, source, stdcutoff):
    modules = []
    for source_row in source.T:
        genes = E.columns[source_row < -source_row.std() * stdcutoff].tolist() + E.columns[source_row > +source_row.std() * stdcutoff].tolist()

        modules.append(Module(genes))
    return modules

def meanshift(E, bandwidth=None, cluster_all=True, **kwargs):
    if bandwidth is None or bandwidth == "auto":
        meanshift = sklearn.cluster.MeanShift(cluster_all=cluster_all)
    else:
        meanshift = sklearn.cluster.MeanShift(bandwidth=bandwidth, cluster_all=cluster_all)

    meanshift.fit(standardize(E).T)
    meanshift.labels_

    modules = convert_labels2modules(meanshift.labels_, E.columns)

    return modules

def baseline_permuted(modules, **kwargs):
    modules = Modules(modules)
    modules = modules.shuffle()
    return modules

## CLAMP-based methods
def clamp_base(E, k=100, adaptive_p=0.05, qvalcutoff=1e-3, **kwargs):
    """
    CLAMP base method: SVD-based matrix factorization without pathway priors.

    Parameters:
    -----------
    E : pandas.DataFrame
        Expression matrix (samples × genes)
    k : int
        Number of latent variables to extract
    adaptive_p : float
        Percentile for adaptive sparsity threshold (default: 0.05)
    qvalcutoff : float
        FDR q-value cutoff for gene selection (default: 1e-3)

    Returns:
    --------
    list of Module objects
    """
    importr("CLAMP")

    # Standardize expression data (samples × genes)
    E_std = standardize(E)

    # Run CLAMPbase (expects genes as rows, samples as columns - transpose)
    # Use localconverter to handle pandas DataFrame conversion
    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv["E_t"] = E_std.T

    # Ensure it's a numeric matrix in R and remove any NA values
    ro.r("E_t = as.matrix(E_t)")
    ro.r("E_t[is.na(E_t)] = 0")  # Replace NA with 0

    ro.r(f"set.seed(1)")
    ro.r(f"clamp_result = CLAMPbase(E_t, k={int(k)}, adaptive.p={adaptive_p})")

    # Extract Z matrix (genes × latent variables)
    Z_matrix = np.array(ro.r["clamp_result"].rx2("Z"))

    # Convert to modules using FDR thresholding
    modules = _ica_fdrtool(E, Z_matrix, qvalcutoff)

    return modules


def clamp_full(E, k=100, adaptive_p=0.05, qvalcutoff=1e-3, pathway_source="CellMarker_2024", **kwargs):
    """
    CLAMP full method: Matrix factorization with pathway-guided refinement.

    Parameters:
    -----------
    E : pandas.DataFrame
        Expression matrix (samples × genes)
    k : int
        Number of latent variables to extract
    adaptive_p : float
        Percentile for adaptive sparsity threshold (default: 0.05)
    qvalcutoff : float
        FDR q-value cutoff for gene selection (default: 1e-3)
    pathway_source : str
        Pathway library to use: "CellMarker_2024", "KEGG_2021_Human", or "GO_BP_2025"

    Returns:
    --------
    list of Module objects
    """
    importr("CLAMP")

    # Standardize expression data
    E_std = standardize(E)

    # Run CLAMPbase first
    # Use localconverter to handle pandas DataFrame conversion
    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv["E_t"] = E_std.T

    # Ensure it's a numeric matrix in R and remove any NA values
    ro.r("E_t = as.matrix(E_t)")
    ro.r("E_t[is.na(E_t)] = 0")  # Replace NA with 0

    ro.r(f"set.seed(1)")
    ro.r(f"base_result = CLAMPbase(E_t, k={int(k)}, adaptive.p={adaptive_p})")

    # Prepare pathway priors - download from Enrichr
    pathway_urls = {
        "CellMarker_2024": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=CellMarker_2024",
        "KEGG_2021_Human": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human",
        "GO_BP_2025": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2025"
    }

    url = pathway_urls.get(pathway_source, pathway_urls["CellMarker_2024"])

    ro.r(f"""
    # Get pathway annotations
    gmt <- getGMT('{url}', '{pathway_source}')
    pathMat <- gmtListToSparseMat(list(pathways = gmt))
    matchedPaths <- getMatchedPathwayMat(pathMat, rownames(E_t), min.genes=2)
    """)

    # Run CLAMPfull with pathway priors
    ro.r(f"""
    full_result = CLAMPfull(
        E_t,
        priorMat = matchedPaths,
        clamp.base.result = base_result,
        use_cpp = TRUE
    )
    """)

    # Extract Z matrix (genes × latent variables)
    Z_matrix = np.array(ro.r["full_result"].rx2("Z"))

    # Convert to modules using FDR thresholding
    modules = _ica_fdrtool(E, Z_matrix, qvalcutoff)

    return modules


def _ica_fdrtool(E, source, qvalcutoff):
    """
    Convert continuous gene loadings to discrete modules using FDR thresholding.

    Parameters:
    -----------
    E : pandas.DataFrame
        Expression matrix (samples × genes) - used to get gene names
    source : numpy.ndarray
        Gene loadings matrix (genes × components)
    qvalcutoff : float
        FDR q-value threshold for significance

    Returns:
    --------
    list of Module objects
    """
    importr("fdrtool")
    rfdrtool = ro.r["fdrtool"]

    modules = []

    # Iterate through each component (column of source matrix)
    for source_row in source.T:
        # Compute q-values for gene loadings using fdrtool
        rresults = rfdrtool(ro.FloatVector(source_row), plot=False, cutoff_method="fndr", verbose=False)
        qvals = np.array(rresults.rx2("qval"))

        # Select significant genes
        genes = E.columns[qvals < qvalcutoff]

        modules.append(Module(genes))

    return modules


## utility functions
def convert_labels2modules(labels, G, ignore_label=None):
    modules = defaultdict(Module)
    for label, gene in zip(labels, G):
        if label != ignore_label:
            modules[label].add(gene)
    return list(modules.values())

def convert_modules2labels(modules, G):
    labels = {}
    for i, module in enumerate(modules):
        for g in module:
            labels[g] = i
    return labels
