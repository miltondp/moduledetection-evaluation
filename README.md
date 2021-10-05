# Module detection evaluation framework

This repository is a fork of https://github.com/saeyslab/moduledetection-evaluation

# Dependencies
This section was modified from the [original repository](https://github.com/saeyslab/moduledetection-evaluation).
Please, check the README.md file there in case you need more instructions, but prioritize the ones here.

Data can be downloaded from zenodo: https://zenodo.org/record/5532578 (Note: we added a new version since 2021/10/04).
Unzip the files and you should get a data and results folder in the root directory of the repo.

After installing [Anaconda or Miniconda](https://www.continuum.io/downloads), please, run the command below instead to create the Conda environment (since other dependencies were added):
```bash
conda env create -n moduledetection_evaluation -f environment.yml
```

The Cython code for the evaluation metrics should also be compiled. For this, change directory to /lib and run `python setup.py build_ext --inplace`

You also need to clone the [clustermatch](https://github.com/greenelab/clustermatch-gene-expr) repository and export this variable:

```bash
export PYTHONPATH=[CLUSTERMATCH_REPO_DIR]/libs:$PYTHONPATH
export NUMBA_NUM_THREADS=3  # number of cores that will be used
```

If everything was installed properly, this command should not raise an error:
```bash
python -c "from clustermatch.coef import cm"
```

Code was tested in Python 3.9.
