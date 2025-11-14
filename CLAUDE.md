# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a module detection evaluation framework forked from https://github.com/saeyslab/moduledetection-evaluation. It benchmarks various gene module detection methods on multiple gene expression datasets.

## Environment Setup

**IMPORTANT:** The conda environment for this project is named `clamp-module-eval` and conda is installed at `~/software/miniforge3/`.

Create and activate the conda environment:
```bash
conda env create -n clamp-module-eval -f environment.yml
conda activate clamp-module-eval
```

Or activate the existing environment:
```bash
source ~/software/miniforge3/etc/profile.d/conda.sh
conda activate clamp-module-eval
```

Compile Cython evaluation metrics:
```bash
cd lib
python setup.py build_ext --inplace
cd ..
```

Set up clustermatch integration (required for clustermatch-based methods):
```bash
# Clone https://github.com/greenelab/clustermatch-gene-expr first
export PYTHONPATH=[CLUSTERMATCH_REPO_DIR]/libs:$(pwd)/lib:$PYTHONPATH
export NUMBA_NUM_THREADS=3  # cores for clustermatch
```

Verify installation:
```bash
python -c "from clustermatch.coef import cm"
```

## Running Module Detection Methods

### Workflow Overview

1. **Generate jobs** - Creates parameter grid for a method across datasets
2. **Run jobs in parallel** - Executes all parameter combinations
3. **Evaluate** - Computes performance scores
4. **Generate plots** - Visualizes results

### Detailed Steps

For a specific method (e.g., `agglom_pearson_abs`):

```bash
export METHOD=agglom_pearson_abs
cd notebooks/

# Step 1: Generate job commands
papermill --log-output generate_jobs.ipynb ${METHOD}-generate_jobs.ipynb -p method_name ${METHOD}

# Step 2: Run jobs in parallel
cd ..
# For non-clustermatch methods (use 3 cores):
parallel -j 3 -a tmp/paramexplo/${METHOD}.txt

# For clustermatch* methods (already parallelized, use 1 core):
parallel -j 1 -a tmp/paramexplo/${METHOD}.txt

# Step 3: Evaluate performance
cd notebooks/
papermill --log-output evaluate.ipynb ${METHOD}-evaluate.ipynb -p method_name ${METHOD}
# Optional: specify cores with -p n_jobs 1

# Step 4: Generate plots
# Open performance_plots.ipynb in browser, set METHOD variable at top, and run
```

## Architecture

### Core Components

**lib/** - Core library code
- `modulecontainers.py` - Module/Bicluster data structures (`Module`, `Modules`, `Bicluster` classes)
- `modulescomparison.py` - Module comparison and scoring logic using Cython-accelerated metrics
- `methods/clustering.py` - Clustering method implementations (60+ methods including kmeans, agglomerative, spectral, WGCNA, PCA, ICA, CLAMP, etc.)
- `methods/biclustering.py` - Biclustering method implementations
- `methods/moduleni.py` - Network inference-based module detection
- `methods/directni.py` - Direct network inference methods (GENIE3, TIGRESS, CLR)
- `simdist.py` - Similarity/distance functions (Pearson, clustermatch, etc.)
- `clustervalidityindices.py` - Cluster validation metrics
- `*.pyx` - Cython implementations for performance-critical metrics (ebcubed, jaccard, cfisher)

**scripts/** - Entry point scripts
- `moduledetection.py` - Main script that runs a method+dataset combination. Takes 3 args: method JSON, dataset JSON, output folder
- `moduledetection_baseline.py` - Baseline method evaluation

**conf/** - Configuration files
- `paramexplo_blueprints.py` - Parameter grids for all methods (blueprints dict defines staticparams and dynparams for each method)
- `datasets/*.json` - Dataset configurations
- `paramexplo/{method}/*.json` - Generated parameter combinations

**notebooks/** - Jupyter notebooks for orchestration
- `generate_jobs.ipynb` - Generates parameter sweep commands for a method
- `evaluate.ipynb` - Computes performance scores for completed runs
- `performance_plots.ipynb` - Creates final visualizations and analysis

**data/** and **results/** - Downloaded separately from Zenodo (https://zenodo.org/record/5532578)

### Method Implementation Pattern

Methods are Python functions in `lib/methods/*.py` that:
1. Accept an expression matrix `E` (pandas DataFrame with genes as columns)
2. Accept method-specific parameters
3. Return a list of `Module` objects (or `Bicluster` objects)

Example signature:
```python
def agglom_pearson_abs(E, k=100, linkage="complete", simdist_function="pearson_correlation_absolute", **kwargs):
    # Implementation
    return modules  # List of Module objects
```

### Parameter Exploration

The `conf/paramexplo_blueprints.py` blueprints dict defines parameter spaces:
- `staticparams` - Fixed parameters for the method
- `dynparams` - Parameters to explore via grid search
- `type` - Method type ("moduledetection", "moduleni", "directni")

Methods with `_auto` suffix use cluster validity indices (CVI) to auto-select module numbers.

### Evaluation Metrics

Two evaluation approaches:
1. **Known modules** (non-human datasets) - Compare against ground truth modules using precision/recall
2. **Regulator coverage** (human datasets) - Evaluate biological meaningfulness via TF enrichment

Cython modules accelerate metric computation (must be compiled before use).

## Adding New Methods

1. Implement method function in `lib/methods/clustering.py` (or appropriate file)
2. Add import in `scripts/moduledetection.py`
3. Define parameter blueprint in `conf/paramexplo_blueprints.py`
4. Run the standard workflow with your method name

## Key Dependencies

- Python 3.9
- Scientific: numpy, pandas, scipy, scikit-learn, statsmodels
- R integration: rpy2 (many methods use R packages)
- Clustering: sklearn, rpy2 interfaces to R packages (WGCNA, cluster, etc.)
- Notebooks: jupyter, jupytext, papermill
- Performance: Cython, numba
- Custom: clustermatch (external repo), cython-munkres-wrapper

## Important Environment Variables

- `PERSOFTWARELOCATION` - Location of external software (FLAME, Click, TransClust, etc.)
- `PYTHONPATH` - Must include clustermatch libs and local lib/ directory
- `NUMBA_NUM_THREADS` - Controls clustermatch parallelization

## Cython Compilation

The Cython evaluation metrics MUST be compiled before running evaluations:
```bash
cd lib
python setup.py build_ext --inplace
```

This creates .so files for ebcubed.pyx, jaccard.pyx, and cfisher.pyx.

## Dataset Structure

Datasets are JSON configs pointing to expression matrices and ground truth modules. Expression files can be:
- Tab-separated text (.tsv/.txt)
- HDF5 (.hdf)
- Pickle (.pkl)

Format: Genes in columns, samples in rows.
