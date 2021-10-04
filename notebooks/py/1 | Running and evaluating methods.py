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
sys.path.insert(0,os.path.abspath("../lib/"))

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
method_name = "dummy" # use the dummy method to check if everything works correctly
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
print("parallel -j 4 -a " + commands_location)

# %% [markdown]
# # Evaluating the method

# %%
from modulescomparison import ModevalKnownmodules, ModevalCoverage

# %% [markdown]
# Note: If you downloaded the results from zenodo, you don't need to rerun this for "dummy", "agglom", "ica_zscore", "spectral_biclust" and "meanshift"

# %% [markdown]
# ## By comparing with known modules

# %% [markdown]
# Evaluate by comparing with known modules

# %%
# create pool of processors
if "pool" in locals().keys():
    pool.close()
pool = mp.Pool(N_JOBS)

# %%
settings_filtered = [setting for setting in settings if not setting["dataset_name"].startswith("human")] # only evaluate non-human datasets
modeval = ModevalKnownmodules(settings_filtered, baseline = True)

# %%
modeval.run(pool)
modeval.save(settings_name)

# %%
modeval.load(settings_name)

# %%
modeval.scores

# %% [markdown]
# ## Using the coverage of regulators

# %%
# create pool of processors
if "pool" in locals().keys():
    pool.close()
pool = mp.Pool(N_JOBS)

# %%
settings_filtered = [setting for setting in settings if setting["dataset_name"].startswith("human")] # only evaluate human datasets
modeval = ModevalCoverage(settings_filtered, baseline = True)

# %%
modeval.run(pool)
modeval.save(settings_name)

# %%
modeval.load(settings_name)

# %%
modeval.scores


# %% [markdown]
# ## Comparing with other methods

# %% [markdown]
# This compares all methods as was done in the paper. Essentially, we will calculate test scores by choosing optimal parameters from one dataset and check how they performed on another dataset. We only compare between [ecoli, yeast], [synthetic] and [human] datasets.

# %%
def score_method(scores):
    methodscores = []
    for ((datasetoi, goldstandardoi), scoresoi), ((datasetor, goldstandardor), scoresor) in product(scores.groupby(["datasetname", "goldstandard"]), scores.groupby(["datasetname", "goldstandard"])):
        if (datasetor.split("_")[0]=="synth" and datasetoi.split("_")[0]!="synth") or (datasetor.split("_")[0]!="synth" and datasetoi.split("_")[0]=="synth"):
            continue
                
        if (goldstandardoi.split("#")[-1] != goldstandardor.split("#")[-1]):
            if (datasetoi.startswith("human") != datasetor.startswith("human")):
                ""
            else:
                continue

        # find the most optimal method parameters in the reference dataset (test dataset)
        bestparams = scoresor[paramsoi].loc[scoresor["score"].idxmax()]
        
        try:
            rowids = scoresoi.index[np.where(np.all([scoresoi[param] == paramvalue for param, paramvalue in bestparams.items()], 0))[0]]
        except:
            print(scoresoi)

        # now find these parameters in the dataset of interest (training dataset)
        rowids = scoresoi.index[np.where(np.all([scoresoi[param] == paramvalue for param, paramvalue in bestparams.items()], 0))[0]]
            
        if len(rowids) == 0:
            print("parameters could not be matched!!", datasetoi, datasetor)
            print(bestparams)
            print([scoresoi[param] == paramvalue for param, paramvalue in bestparams.items()])
        if len(rowids) > 1:
            print(datasetoi)
            print("multiple matched parameters")
            print(scoresoi.loc[rowids][paramsoi])

        methodscores.append({
            "datasetoi":datasetoi,
            "datasetor":datasetor,
            "score":scoresoi.loc[rowids,"score"].max(),
            "method":methodname,
            "goldstandardoi":goldstandardoi,
            "goldstandardor":goldstandardor,
            "ofinterest":datasetoi + "#" + goldstandardoi,
            "ofreference":datasetor + "#" + goldstandardor,
            "runningtime":scoresoi.loc[rowids, "runningtime"].mean() if "runningtime" in scoresoi.columns else 0,
            "moduledef":scoresoi.loc[rowids, "moduledef"].tolist()[0],
            "organismoi":scoresoi.loc[rowids, "organism"].tolist()[0],  
        })
    
    return pd.DataFrame(methodscores)


# %%
methodnames = ["dummy", "agglom", "ica_zscore", "spectral_biclust", "meanshift"]

# %%
finalscores = []
for methodname in methodnames:
    settings_name = "paramexplo/" + methodname
    settings = json.load(open("../conf/settings/{}.json".format(settings_name)))
    settings_dataset = pd.DataFrame([dict(settingid=setting["settingid"], **json.load(open("../" + setting["dataset_location"]))["params"]) for setting in settings])
    settings_method = pd.DataFrame([dict(settingid=setting["settingid"], **json.load(open("../" + setting["method_location"]))["params"]) for setting in settings])
    
    print(methodname)
    paramsoi = methodparamsoi[methodname]

    scores = pd.DataFrame()
    
    modeval = ModevalKnownmodules(settings_name)
    modeval.load(settings_name)
    modeval.scores["score"] = modeval.scores["F1rprr_permuted"]
    modeval.scores["moduledef"] = [modulesname if modulesname in ["minimal", "strict"] else "interconnected" for modulesname in modeval.scores["knownmodules_name"]]
    modeval.scores = modeval.scores.merge(settings_dataset, on="settingid").merge(settings_method, on="settingid")
    scores = scores.append(modeval.scores, ignore_index=True)
    
    modeval = ModevalCoverage(settings_name)
    modeval.load(settings_name)
    modeval.scores["score"] = modeval.scores["aucodds_permuted"]
    modeval.scores = modeval.scores.merge(settings_dataset, on="settingid").merge(settings_method, on="settingid")
    scores = scores.append(modeval.scores, ignore_index=True)
    
    methodscores = score_method(scores)
    
    methodscores["organismnetoi"] = [dataset.split("_")[0] for dataset in methodscores["goldstandardoi"]]
    methodscores["organismnetor"] = [dataset.split("_")[0] for dataset in methodscores["goldstandardor"]]

    finalscores.append(methodscores)
finalscores = pd.concat(finalscores, ignore_index=True)

# %% [markdown]
# The final scores contains all the comparisons we made, together with a final score in the score column:

# %%
finalscores

# %%
finalscores.query("method == 'ica_zscore'")


# %% [markdown]
# We add weights to the test scores, because e.g. E. coli datasets will have many more test scores as there are more "reference" datasets available.

# %%
def add_weights(scores):
    weights = []
    scores["moduledef"] = scores["moduledef"].fillna("")
    for organismoi, subscores in scores.groupby("organismoi"):
        moduledef_weights = 1/subscores.groupby("moduledef")["score"].count()
        for moduledef, weight in moduledef_weights.items():
            weights.append({
                    "organism":organismoi,
                    "moduledef":moduledef,
                    "weight":weight / len(moduledef_weights)
                })
    weights = pd.DataFrame(weights).set_index(["organism", "moduledef"])["weight"]
    
    scores["weight"] = weights.loc[pd.Index(scores[["organismoi", "moduledef"]])].tolist()
    
    return scores


# %%
trainingscores_ = add_weights(finalscores.loc[(finalscores["ofinterest"] == finalscores["ofreference"])])
testscores_ = add_weights(finalscores.loc[(finalscores["ofinterest"] != finalscores["ofreference"]) & (finalscores["organismnetoi"] != finalscores["organismnetor"])])

# %% [markdown]
# Do a weighted mean:

# %%
trainingscores = trainingscores_.groupby("method").apply(lambda x: np.average(x.score, weights=x.weight))
testscores = testscores_.groupby("method").apply(lambda x: np.average(x.score, weights=x.weight))

# %%
testscores_.to_csv("../results/testscores_.tsv", sep="\t")
trainingscores_.to_csv("../results/trainingscores_.tsv", sep="\t")

# %%
trainingscores

# %%
testscores

# %% [markdown]
# Visualization of overall training and test scores:

# %%
import matplotlib as mpl

# %%
# A bar chart is actually not the ideal representation here, given that we're working with ratios. 
# This way of plotting is kept here because it most closely resembles that of the paper.

fig, ax = subplots(figsize=(len(trainingscores)/2, 5))

methodorder = testscores.sort_values(ascending=False).index

ax.axhline(1, color = "black")
ax.bar(range(len(methodorder)), trainingscores[methodorder], color="grey")
ax.bar(range(len(methodorder)), testscores[methodorder], color="#333333")
ax.set_xticks(np.arange(len(methodorder)))
ax.set_xticklabels(methodorder, rotation=45, ha="right", va="top")
ax.tick_params(labelsize=14)
""

# %%
# A better way to visualize the data would be dotplot

fig, ax = subplots(figsize=(5, len(trainingscores)/2))

methodorder = testscores.sort_values(ascending=True).index

ax.axvline(1, color = "black")
for y, method in enumerate(methodorder):
    ax.plot([trainingscores[method], testscores[method]], [y, y], zorder = 0, color = "grey")
ax.scatter(trainingscores[methodorder], range(len(methodorder)), color="grey", s = 20)
ax.scatter(testscores[methodorder], range(len(methodorder)), color="#333333", s = 100)
ax.set_yticks(np.arange(len(methodorder)))
ax.set_yticklabels(methodorder)
ax.tick_params(labelsize=14)
ax.set_xlim([0, ax.get_xlim()[1]])
""

# %% [markdown]
# You can also calculate scores for a particular organism, ...:

# %%
trainingscores = trainingscores_.query("organismoi == 'ecoli'").groupby("method").apply(lambda x: np.average(x.score, weights=x.weight))
testscores = testscores_.query("organismoi == 'ecoli'").groupby("method").apply(lambda x: np.average(x.score, weights=x.weight))

# %%
# A better way to visualize the data would be dotplot

fig, ax = subplots(figsize=(5, len(trainingscores)/2))

methodorder = testscores.sort_values(ascending=True).index

ax.axvline(1, color = "black")
for y, method in enumerate(methodorder):
    ax.plot([trainingscores[method], testscores[method]], [y, y], zorder = 0, color = "grey")
ax.scatter(trainingscores[methodorder], range(len(methodorder)), color="grey", s = 20)
ax.scatter(testscores[methodorder], range(len(methodorder)), color="#333333", s = 100)
ax.set_yticks(np.arange(len(methodorder)))
ax.set_yticklabels(methodorder)
ax.tick_params(labelsize=14)
ax.set_xlim([0, ax.get_xlim()[1]])
""

# %%
