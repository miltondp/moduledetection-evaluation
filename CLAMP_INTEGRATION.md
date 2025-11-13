# CLAMP Integration into Module Detection Evaluation Framework

## Summary

Successfully integrated the CLAMP R package into the gene module detection evaluation framework. Two CLAMP variants have been implemented:
- **clamp_base**: SVD-based matrix factorization without pathway priors
- **clamp_full**: Matrix factorization with pathway-guided refinement

## Files Modified

### 1. `lib/clustering.py`
**Added:**
- `clamp_base(E, k, adaptive_p, qvalcutoff, **kwargs)` - CLAMP base method
- `clamp_full(E, k, adaptive_p, qvalcutoff, pathway_source, **kwargs)` - CLAMP full method
- `_ica_fdrtool(E, source, qvalcutoff)` - Helper function for FDR-based module extraction
- rpy2 imports with proper pandas2ri/numpy2ri conversion setup

**Modified:**
- Commented out `agglom_clustermatch` and `agglom_clustermatch_linear` (clustermatch dependency not available)

### 2. `lib/simdist.py`
**Modified:**
- Commented out clustermatch imports and functions (not needed for CLAMP)

### 3. `conf/paramexplo_blueprints.py`
**Added:**
- `clamp_base` blueprint with parameter grid:
  - `k`: 25-300 in steps of 25 (12 values)
  - `adaptive_p`: [0.01, 0.05, 0.1] (3 values)
  - `qvalcutoff`: 10^(-1 to -10) (10 values)
  - Total: 360 parameter combinations per dataset

- `clamp_full` blueprint with same parameter grid plus:
  - `pathway_source`: "CellMarker_2024" (static)
  - Total: 360 parameter combinations per dataset

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
   - `adaptive_p` (0.01-0.1): Controls sparsity in gene loadings
   - `qvalcutoff` (10^-10 to 10^-1): FDR threshold for significance

## Dependencies Installed

- R package: **fdrtool** (for FDR-based thresholding)
- Python package: **rpy2** (already in environment)

## Next Steps for Testing

### 1. Generate Parameter Sweep Jobs
```bash
cd notebooks/
papermill --log-output generate_jobs.ipynb clamp_base-generate_jobs.ipynb -p method_name clamp_base
```

This will create ~3,600 jobs (360 param combinations × 10 datasets)

### 2. Run Jobs in Parallel
```bash
cd ..
parallel -j 3 -a tmp/paramexplo/clamp_base.txt
```

Note: Use `-j 3` or lower since CLAMP is R-based and may use multiple cores

### 3. Evaluate Performance
```bash
cd notebooks/
papermill --log-output evaluate.ipynb clamp_base-evaluate.ipynb -p method_name clamp_base -p n_jobs 3
```

### 4. Repeat for clamp_full
```bash
# Generate jobs
papermill --log-output generate_jobs.ipynb clamp_full-generate_jobs.ipynb -p method_name clamp_full

# Run jobs (may be slower due to pathway downloads)
cd ..
parallel -j 1 -a tmp/paramexplo/clamp_full.txt  # Use -j 1 since pathway downloads can be network-intensive

# Evaluate
cd notebooks/
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

## Files Created

- `lib/clustering.py` - CLAMP method implementations
- `test_clamp.py` - Simple test script (works with real data only)
- `CLAMP_INTEGRATION.md` - This documentation

## Troubleshooting

### Issue: "fdrtool not installed"
**Solution:** `R -e "install.packages('fdrtool', repos='https://cloud.r-project.org')"`

### Issue: "Conversion 'py2rpy' not defined for DataFrame"
**Solution:** Already handled using `localconverter` with `pandas2ri.converter`

### Issue: "requires numeric/complex matrix/vector arguments"
**Solution:** Already handled by converting to R matrix and replacing NA values

### Issue: clamp_full finds "0 genes in intersection"
**Solution:** This happens when gene names don't match pathway databases. Use real gene symbols (not synthetic IDs).
