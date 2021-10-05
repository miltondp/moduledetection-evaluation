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

You also need to clone the [clustermatch](https://github.com/greenelab/clustermatch-gene-expr) repository and export this variable:

```bash
export PYTHONPATH=[CLUSTERMATCH_REPO_DIR]/libs:`pwd`/lib:$PYTHONPATH
export NUMBA_NUM_THREADS=3  # number of cores that will be used by clustermatch
```

If everything was installed properly, this command should not raise an error:
```bash
python -c "from clustermatch.coef import cm"
```

# Running

For a specific clustering method `${METHOD}`:

```bash
# select one method
export METHOD=agglom_pearson_abs

# generate the jobs
cd notebooks/
papermill --log-output ${METHOD}-generate_jobs.ipynb ${METHOD}-generate_jobs.ipynb

# from the previous run you will see a message indicating which command you have to run next
# for example:
cd ..
parallel -j 3 -a tmp/paramexplo/{METHOD}.txt
```
