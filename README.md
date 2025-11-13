# Module detection evaluation framework

This repository is a fork of https://github.com/saeyslab/moduledetection-evaluation

# Dependencies
This section was modified from the [original repository](https://github.com/saeyslab/moduledetection-evaluation).
Please, check the README.md file there in case you need more instructions, but prioritize the ones here.

Data can be downloaded from zenodo: https://zenodo.org/record/5532578 (new version since 2021/10/04).
Unzip the files and you should get a data and results folder in the root directory of the repo.

After installing [Anaconda or Miniconda](https://www.continuum.io/downloads), please, run the command below instead to create the Conda environment (since other dependencies were added):
```bash
conda env create -n moduledetection_evaluation -f environment.yml
conda activate moduledetection_evaluation
```

The Cython code for the evaluation metrics should also be compiled. For this, change directory to /lib and run
```bash
cd lib
python setup.py build_ext --inplace
```

## Install CLAMP

```bash
conda install -c conda-forge -c bioconda r-base=4.4 r-essentials r-bigstatsr r-data.table r-dplyr r-rsvd r-glmnet r-matrix r-knitr r-here r-tibble r-tidyr r-devtools r-irkernel  r-knitr r-pkgdown hdf5 r-hdf5r bioconductor-rhdf5 bioconductor-qvalue bioconductor-biomart bioconductor-recount3 bioconductor-preprocesscore bioconductor-fgsea bioconductor-org.Hs.eg.db bioconductor-ensdb.hsapiens.v86 bioconductor-experimenthub bioconductor-genomeinfodbdata bioconductor-bioccheck bioconductor-complexheatmap bioconductor-biocstyle bioconductor-experimenthubdata r-ggrepel r-irlba rpy2 r-fdrtool

#REPO_PATH=~/path/to/CLAMP
REPO_PATH=/home/miltondp/projects/clamp/CLAMP

Rscript -e "devtools::install_local('$REPO_PATH', force=TRUE, dependencies=FALSE)"

Rscript -e "library(CLAMP); cat('CLAMP version:', packageVersion('CLAMP'), '\n')"
```

# Running

For a specific clustering method `${METHOD}`:

```bash
# select one method
export METHOD=agglom_pearson_abs

# generate the jobs
cd notebooks/
papermill --log-output generate_jobs.ipynb ${METHOD}-generate_jobs.ipynb -p method_name ${METHOD}

# next, choose the following two commands depending on whether the method is clustermatch* or not
# not clustermatch*:
cd ..
parallel -j 3 -a tmp/paramexplo/${METHOD}.txt

# clustermatch* (use only one core in parallel since the method already parallelizes
cd ..
parallel -j 1 -a tmp/paramexplo/${METHOD}.txt
```

Then compute the scores for each method (you can specify the number of cores here with `n_jobs` as shown below):
```bash
cd notebooks/
papermill --log-output evaluate.ipynb ${METHOD}-evaluate.ipynb -p method_name ${METHOD} # -p n_jobs 1
```

Finally, open `performance_plots.ipynb` from the browser, specify your method name (the one you specified with `${METHOD}`) at the
top of the notebook, and run the notebook to get the final scores and plots.
