# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill,-trusted
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Overall performance of module detection methods

# %% scrolled=true
import sys
import os
# sys.path.insert(0,os.path.abspath("../lib/"))

import json

from util import JSONExtendedEncoder

# %load_ext autoreload
# %autoreload 2

# %matplotlib inline
from matplotlib.pyplot import *

import pandas as pd
import numpy as np

import multiprocessing as mp

from itertools import product

import itertools
import shutil

import os

conf_folder = "conf/"

# %% [markdown] tags=[]
# # Settings

# %%
N_JOBS = 1
# N_JOBS = mp.cpu_count()-1

# %% tags=["parameters"]
method_name = None
assert method_name is not None, "You have to specify a method_name"

# %% [markdown] tags=[]
# # Running a method on different parameter settings and datasets

# %% [markdown]
# Note: If you downloaded the results from zenodo, you don't need to rerun this for "dummy", "agglom", "ica_zscore", "spectral_biclust" and "meanshift"

# %% [markdown]
# The following code will explore the parameters of a module detection method on every dataset using a grid-search approach.

# %% [markdown]
# If you want to run your own method, you should wrap it into a python function and add its parameters to `conf/paramexplo_blueprints.py`. We will show the whole workflow here for a "dummy"  (but fast) clustering method, which will simply group genes randomly.

# %% [markdown]
# Every module detection method is wrapped in a python function (see `scripts/moduledetection.py`)
#
# Because module detection methods usually take a while to run, we generate the files necessary to run a method on the several parameter settings and datasets here. These can then be easily called from the commandline, for example on a computer cluster or locally using GNU `parallel`.
#
# This function will be called by scripts/moduledetection.py , which will save the modules in the correct format along with additional run information (such as running times).

# %%
# datasets to run
datasetnames = [
    "ecoli_colombos",
    "ecoli_dream5",
    "yeast_gpl2529",
    "yeast_dream5",
    "synth_ecoli_regulondb",
    "synth_yeast_macisaac",
    "human_tcga",
    "human_gtex",
    "human_seek_gpl5175",
    "ecoli_precise2"
]

# choose the method to evaluate
method_name = "agglom_pearson_abs" # use the dummy method to check if everything works correctly
# method_name = "agglom" # this method runs very fast, and has the best performance among clustering methods
# method_name = "ica_zscore" # this method runs very slow, but has approx. the highest performance in the benchmark
# method_name = "spectral_biclust" # top biclustering method
# method_name = "meanshift"

# %% [markdown]
# To add your own method, create a function with "your_method_name" in the `lib/clustering.py` file (or any other file as long as it's imported in `scripts/moduledetection.py`.
# This function should accept an `E` object (which is a dataframe with genes in columns) and any additional parameters
# Then add reasonable parameter setting of your method to `conf/paramexplo_blueprints.py`.

# %% [markdown]
# method_name = "your_method_name"

# %%
# paramexplo_blueprints.py stores for every method the parameters which will be varied using a grid-search approach.
# %run ../conf/paramexplo_blueprints.py
methodblueprint = blueprints[method_name]

# %%
methodblueprint

# %% [markdown]
# Generate different parameter settings using a grid-search.

# %%
params_folder = "conf/paramexplo/" + method_name + "/"
if os.path.exists("../" + params_folder):
    shutil.rmtree("../" + params_folder)
os.makedirs("../" + params_folder)

methodsettings = []
method_locations = []
i = 0
for dynparam_combination in list(itertools.product(*[methodblueprint["dynparams"][param] for param in sorted(methodblueprint["dynparams"].keys())])):
    method = {"params":{}}
    method["params"] = methodblueprint["staticparams"].copy()
    method["params"].update(dict(zip(sorted(methodblueprint["dynparams"].keys()), dynparam_combination)))
    method["location"] = params_folder + str(i) + ".json"
    method["seed"] = 0

    methodsettings.append(method)

    json.dump(method, open("../" + method["location"], "w"), cls=JSONExtendedEncoder)

    method_locations.append(method["location"])

    i+=1

# %% [markdown]
# Now combine the different parameter settings and datasets. Then generate the different python commands to run every parameter setting and dataset in parallel.

# %%
settings_name = "paramexplo/{method_name}".format(method_name = method_name)
settings = []
for datasetname in datasetnames:
    for setting_ix, methodsetting in enumerate(methodsettings):
        settingid = datasetname + "_" + str(setting_ix)
        settings.append({
            "dataset_location":"conf/datasets/" + datasetname + ".json",
            "dataset_name":datasetname,
            "method_location":methodsetting["location"],
            "output_folder":"results/" + methodblueprint["type"] + "/{settings_name}/{settingid}/".format(settings_name=settings_name, settingid=settingid),
            "settingid":settingid
        })
json.dump(settings, open("../conf/settings/{settings_name}.json".format(settings_name=settings_name), "w"))

# %%
settings_dataset = pd.DataFrame([dict(settingid=setting["settingid"], **json.load(open("../" + setting["dataset_location"]))["params"]) for setting in settings])
settings_method = pd.DataFrame([dict(settingid=setting["settingid"], **json.load(open("../" + setting["method_location"]))["params"]) for setting in settings])

# %%
commands = ""
for i, setting in enumerate(settings):
    #commands += "python scripts/moduledetection.py {method_location} {dataset_location} {output_folder} 0 test\n".format(**setting)
    commands += "python3 scripts/" + methodblueprint["type"] + ".py {method_location} {dataset_location} {output_folder}\n".format(**setting)

commands_location = "tmp/{settings_name}.txt".format(**locals())
os.makedirs("../" + os.path.dirname(commands_location), exist_ok=True)
with open("../" + commands_location, "w") as outfile:
    outfile.write(commands)
commands_location = "tmp/{settings_name}.txt".format(**locals())
os.makedirs(os.path.dirname("../tmp/" + commands_location), exist_ok=True)
with open("../tmp/" + commands_location, "w") as outfile:
    outfile.write(commands)
    
#script_location = generate_batchcode(commands_location, settings_name, len(settings), {"memory":"10G", "numcores":1}, "biclust_comp2")

# this command can be used on most linux computers to run the different parameter settings in parallel
print("Run the following command:\nparallel -j 4 -a " + commands_location)

# %%
