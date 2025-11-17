# CLAMP Integration into Module Detection Evaluation Framework

## Summary

Successfully integrated the CLAMP R package into the gene module detection evaluation framework. Two CLAMP variants have been implemented:
- **clamp_base**: SVD-based matrix factorization without pathway priors
- **clamp_full**: Matrix factorization with pathway-guided refinement

## Quick Start: Complete Workflow

```bash
# 1. Activate conda environment and setup PYTHONPATH (from root directory)
source ~/software/miniforge3/etc/profile.d/conda.sh
conda activate clamp-module-eval
export PYTHONPATH=`realpath lib`

# 2. Generate jobs for clamp_base
cd notebooks/
papermill --log-output generate_jobs.ipynb clamp_base-generate_jobs.ipynb -p method_name clamp_base

# 3. Run jobs in parallel
cd ..
parallel -j 3 -a tmp/paramexplo/clamp_base.txt

# 4. Evaluate performance
cd notebooks/
# Re-export if in new shell
# export PYTHONPATH=`realpath ../lib`
papermill --log-output evaluate.ipynb clamp_base-evaluate.ipynb -p method_name clamp_base -p n_jobs 3

# 5. (Optional) Repeat for clamp_full
cd ..
#export PYTHONPATH=`realpath lib`
cd notebooks/
papermill --log-output generate_jobs.ipynb clamp_full-generate_jobs.ipynb -p method_name clamp_full
cd ..
parallel -j 1 -a tmp/paramexplo/clamp_full.txt
cd notebooks/
#export PYTHONPATH=`realpath ../lib`
papermill --log-output evaluate.ipynb clamp_full-evaluate.ipynb -p method_name clamp_full -p n_jobs 3

# 6. Generate plots (open in Jupyter/browser)
# Open notebooks/performance_plots.ipynb and add "clamp_base" and "clamp_full" to methods list
```

## Files Modified

### 1. `lib/methods/clustering.py`
**Added:**
- `clamp_base(E, k, adaptive_p, qvalcutoff, **kwargs)` - CLAMP base method (line ~398)
- `clamp_full(E, k, adaptive_p, qvalcutoff, pathway_source, **kwargs)` - CLAMP full method (line ~443)
- Both methods placed in the decomposition/factorization section alongside `pca`, `ica_*`, and `nmf_*` methods

**Modified:**
- `_ica_fdrtool(E, source, qvalcutoff)` - Updated with comprehensive docstring (line ~566)
- rpy2 imports updated with proper pandas2ri/numpy2ri conversion and deprecation handling

**Note:** The file `lib/clustering.py` previously created during development has been deprecated. All CLAMP methods are now properly integrated into `lib/methods/clustering.py` following the project architecture.

### 2. `lib/simdist.py`
**Modified:**
- Commented out clustermatch imports and functions (not needed for CLAMP)

### 3. `conf/paramexplo_blueprints.py`
**Added:**
- `clamp_base` blueprint with parameter grid:
  - `k`: 25-300 in steps of 25 (12 values)
  - `qvalcutoff`: 10^(-1 to -10) (10 values)
  - `adaptive_p`: 0.05 (fixed, not explored)
  - Total: 120 parameter combinations per dataset

- `clamp_full` blueprint with same parameter grid plus:
  - `pathway_source`: "CellMarker_2024" (static)
  - Total: 120 parameter combinations per dataset

**Note:** `adaptive_p` was removed from parameter exploration (see git commit 7edaec6) to match PCA's parameter grid. The fixed value of 0.05 is used for all runs.

- Added entries to `methodparamsoi` and `methodparams_modulenumber` dicts

## Implementation Details

### CLAMP Base Workflow
1. Standardize expression data using framework's `standardize()` function
2. Transpose to genes×samples format (CLAMP expects genes as rows)
3. Convert to R matrix and handle NA values
4. Run `CLAMPbase()` to extract latent variables
5. Apply FDR thresholding using `fdrtool` to convert continuous loadings to discrete modules

### CLAMP Full Workflow
1. Same preprocessing as clamp_base
2. Run `CLAMPbase()` for initialization
3. Download pathway annotations from Enrichr (CellMarker_2024, KEGG, or GO_BP)
4. Match pathways to gene space
5. Run `CLAMPfull()` with pathway priors
6. Apply FDR thresholding to extract modules

### Key Technical Decisions

1. **Module Extraction**: Using FDR-based thresholding (`_ica_fdrtool`) similar to PCA implementation, rather than top-k or z-score approaches

2. **Preprocessing**: Using framework's `standardize()` function for consistency with other methods, rather than CLAMP's built-in preprocessing

3. **Pandas/R Conversion**: Using `localconverter` with `pandas2ri.converter` for proper DataFrame→R matrix conversion (newer rpy2 pattern)

4. **NA Handling**: Replacing NA values with 0 after standardization to prevent CLAMP errors

5. **Parameter Ranges**:
   - `k` (25-300): Standard range used by PCA, ICA, and other factorization methods
   - `adaptive_p`: Fixed at 0.05 (not explored in parameter sweep)
   - `qvalcutoff` (10^-10 to 10^-1): FDR threshold for significance

## Dependencies Installed

- R package: **fdrtool** (for FDR-based thresholding)
- Python package: **rpy2** (already in environment)

## Next Steps for Testing

### 0. Setup Environment (Required)
From the root directory of the project:
```bash
source ~/software/miniforge3/etc/profile.d/conda.sh
conda activate clamp-module-eval
export PYTHONPATH=`realpath lib`
```

This must be run before any papermill commands to ensure the conda environment is active and custom modules in `lib/` are importable.

### 1. Generate Parameter Sweep Jobs
```bash
cd notebooks/
papermill --log-output generate_jobs.ipynb clamp_base-generate_jobs.ipynb -p method_name clamp_base
```

This will create ~1,200 jobs (120 param combinations × 10 datasets)

### 2. Run Jobs in Parallel
```bash
cd ..
parallel -j 3 -a tmp/paramexplo/clamp_base.txt
```

Note: Use `-j 3` or lower since CLAMP is R-based and may use multiple cores

### 3. Evaluate Performance
```bash
cd notebooks/
# Make sure PYTHONPATH is still set
export PYTHONPATH=`realpath ../lib`
papermill --log-output evaluate.ipynb clamp_base-evaluate.ipynb -p method_name clamp_base -p n_jobs 3
```

### 4. Repeat for clamp_full
```bash
# From root directory, ensure environment is active and PYTHONPATH is set
source ~/software/miniforge3/etc/profile.d/conda.sh
conda activate clamp-module-eval
export PYTHONPATH=`realpath lib`

# Generate jobs
cd notebooks/
papermill --log-output generate_jobs.ipynb clamp_full-generate_jobs.ipynb -p method_name clamp_full

# Run jobs (may be slower due to pathway downloads)
cd ..
parallel -j 1 -a tmp/paramexplo/clamp_full.txt  # Use -j 1 since pathway downloads can be network-intensive

# Evaluate
cd notebooks/
export PYTHONPATH=`realpath ../lib`  # Ensure PYTHONPATH is set from notebooks directory
papermill --log-output evaluate.ipynb clamp_full-evaluate.ipynb -p method_name clamp_full -p n_jobs 3
```

### 5. Generate Performance Plots
Open `notebooks/performance_plots.ipynb` and add "clamp_base" and "clamp_full" to the methods list.

## Expected Outputs

### Scores Metrics
- **Known modules comparison** (non-human datasets):
  - Recovery, Relevance, F1rr
  - Recall, Precision, F1rp
  - F1rprr (combined metric)
  - Comparison vs. permuted baseline

- **Regulator coverage** (human datasets):
  - AUCODDS (area under curve of odds)
  - Comparison vs. permuted baseline

### Performance Comparisons
- CLAMP vs. PCA (similar matrix factorization)
- CLAMP vs. ICA (independent component analysis)
- CLAMP vs. clustering methods (agglomerative, k-means, etc.)
- clamp_base vs. clamp_full (with/without pathway priors)

## Known Limitations

1. **Clustermatch methods disabled**: The `agglom_clustermatch` methods are commented out since clustermatch is not currently installed

2. **Synthetic data testing**: The test script with random synthetic data doesn't work well because:
   - Gene names don't match pathway databases (for clamp_full)
   - fdrtool requires realistic data distributions
   - Real datasets from the benchmark should be used for testing

3. **Pathway availability**: clamp_full downloads pathways from Enrichr on each run. For large-scale evaluations, consider caching or pre-downloading.

4. **Gene name matching**: clamp_full requires gene symbols that match the pathway databases. Ensembl IDs or other identifiers may need mapping.

## Comparison to PCA

CLAMP is most similar to PCA in this benchmark:
- Both use matrix factorization
- Both extract latent variables/components
- Both use FDR thresholding to convert loadings to modules

**Key differences:**
- CLAMP uses adaptive sparsity (sets weak loadings to zero)
- CLAMP can incorporate pathway priors (clamp_full variant)
- CLAMP alternates between updating gene loadings and sample scores

Expected performance: CLAMP should perform comparably or better than PCA, especially clamp_full on datasets where pathway priors are informative.

## Files Created/Modified

- `lib/methods/clustering.py` - CLAMP methods integrated into main clustering module
- `test_clamp.py` - Simple test script (works with real data only)
- `CLAMP_INTEGRATION.md` - This documentation
- `lib/clustering.py` - **DEPRECATED** - Previously used during development, now superseded by proper integration

## Important: fdrtool Error Fix (Nov 2024)

**If you see errors like:**
```
Error in optimize(nlogL, lower = lo, upper = up) : 'xmin' not less than 'xmax'
```

This has been **FIXED** as of November 2024. The fix handles CLAMP's sparse gene loadings robustly by:
1. Only applying fdrtool to non-zero gene loadings
2. Using percentile-based fallbacks when fdrtool fails
3. Adapting to different levels of sparsity

**See `CLAMP_FDRTOOL_FIX.md` for detailed explanation.**

The error messages may still appear in R console output but are harmless - they're caught and handled gracefully with fallback thresholds. Jobs will complete successfully.

## Troubleshooting

### Issue: "ModuleNotFoundError: No module named 'clustering'" or similar import errors
**Solution:** Set PYTHONPATH before running papermill:
```bash
export PYTHONPATH=`realpath lib`  # From root directory
# OR
export PYTHONPATH=`realpath ../lib`  # From notebooks/ directory
```

### Issue: "fdrtool not installed"
**Solution:**
```bash
R -e "install.packages('fdrtool', repos='https://cloud.r-project.org')"
```

### Issue: "Conversion 'py2rpy' not defined for DataFrame"
**Solution:** Already handled using `localconverter` with `pandas2ri.converter`

### Issue: "requires numeric/complex matrix/vector arguments"
**Solution:** Already handled by converting to R matrix and replacing NA values

### Issue: clamp_full finds "0 genes in intersection"
**Solution:** This happens when gene names don't match pathway databases. Use real gene symbols (not synthetic IDs).
