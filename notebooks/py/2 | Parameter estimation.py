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

# %% jupyter={"outputs_hidden": false}
import sys
import os
sys.path.insert(0,os.path.abspath("../lib/"))

import json

from util import JSONExtendedEncoder

from modulecontainers import Modules

# %load_ext autoreload
# %autoreload 2

# %matplotlib inline
from matplotlib.pyplot import *

import pandas as pd
import numpy as np

import multiprocessing as mp

from collections import OrderedDict

from itertools import product
from collections import defaultdict

import itertools
import shutil

conf_folder = "conf/"

import os

# %% [markdown]
# This runs parameter estimation using either functional enrichment as in Figure 3 of our paper.

# %% [markdown]
# ## Run a method

# %% [markdown]
# To run individual methods:

# %% tags=[]
from modulescomparison import ModevalFunctional, ModevalKnownmodules, ModevalCoverage

# %% tags=[]
methodname = "agglom"
settings_name = "paramexplo/" + methodname
settings = json.load(open("../conf/settings/{}.json".format(settings_name)))

# %% jupyter={"outputs_hidden": false}
# create pool of processors
if "pool" in locals().keys():
    pool.close()
pool = mp.Pool(1)

# %% jupyter={"outputs_hidden": false}
modeval_functional = ModevalFunctional([setting for setting in settings if not setting["dataset_name"].startswith("synth")])

# %% jupyter={"outputs_hidden": false}
modeval_functional.run(pool)

# %% jupyter={"outputs_hidden": false}
modeval_functional.save(settings_name)

# %% [markdown]
# ## Combine output from different methods

# %% [markdown]
# Combine different methods to compare:

# %%
method_names = ["dummy", "agglom"]

# %% jupyter={"outputs_hidden": false}
scores = {}
for methodname in method_names:
    settings_name = "paramexplo/" + methodname
    settings = json.load(open("../conf/settings/{}.json".format(settings_name)))
    settings_dataset = pd.DataFrame([dict(settingid=setting["settingid"], **json.load(open("../" + setting["dataset_location"]))["params"]) for setting in settings])
    settings_method = pd.DataFrame([dict(settingid=setting["settingid"], **json.load(open("../" + setting["method_location"]))["params"]) for setting in settings])
    
    modeval_functional.load(settings_name)

    modeval_knownmodules = ModevalKnownmodules(settings)
    modeval_knownmodules.load(settings_name)
    modeval_coverage = ModevalCoverage(settings)
    modeval_coverage.load(settings_name)

    modeval_coverage.scores["refscore"] = modeval_coverage.scores["aucodds_permuted"]
    modeval_knownmodules.scores["refscore"] = modeval_knownmodules.scores["F1rprr_permuted"]
    modeval_knownmodules.scores["moduledef"] = [modulesname if modulesname in ["minimal", "strict"] else "interconnected" for modulesname in modeval_knownmodules.scores["knownmodules_name"]]

    refscores = pd.concat([modeval_knownmodules.scores, modeval_coverage.scores], ignore_index=True)[["settingid", "refscore", "goldstandard", "moduledef"]]
    refscores = refscores.merge(settings_dataset, on="settingid")
    
    proxyscores = modeval_functional.scores.rename(columns={"score":"proxyscore", "scorename":"proxyname"})
    
    combined = refscores.merge(proxyscores, on="settingid")
    
    scores[methodname] = combined

# %% tags=[]
optimalproxy_funcs = {"asw":np.nanargmax, "ch":np.nanargmax, "dbindex":np.nanargmin, "dbstarindex":np.nanargmin, "bhi#go":np.nanargmax, "bhi#kegg":np.nanargmax, "bhi#reactome":np.nanargmax, 
                      "faucodds#go":np.nanargmax, "faucodds#kegg":np.nanargmax, "faucodds#reactome":np.nanargmax}

# %%
methodscores

# %% jupyter={"outputs_hidden": false}
finalscores = []

for method, methodscores in scores.items():
    methodscores["proxyscore"][np.isinf(methodscores["proxyscore"])] = np.nan
    for (dataset, goldstandard, proxyname), subscores in methodscores.groupby(["datasetname", "goldstandard", "proxyname"]):
        if proxyname not in optimalproxy_funcs.keys():
            continue
        if proxyname.endswith("reactome"):
            continue

        optimalproxy_func = optimalproxy_funcs[proxyname]
        if np.all(np.isnan(subscores["proxyscore"])):
            print("random?", method)
            print(dataset, refgrouper,proxygrouper)
            optimalproxy = np.random.choice(subscores.index)
        else:
            optimalproxy = subscores.index[optimalproxy_func(subscores["proxyscore"])]
            
        row = subscores.loc[optimalproxy].copy()
        row["optimalrefscore"] = subscores["refscore"].max()
        row["method"] = method
        finalscores.append(row)
finalscores = pd.DataFrame(finalscores).rename(columns={"organism":"organismoi", "goldstandard":"goldstandardoi", "refscore":"score"})

# %% tags=[]
from collections import Counter


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


# %% jupyter={"outputs_hidden": false}
proxyscores_ = add_weights(finalscores)
proxyscores = proxyscores_.groupby(["proxyname", "method"]).apply(lambda x: np.average(x.score, weights=x.weight)).unstack()

proxyscores_.groupby("method").apply(lambda x: np.average(x.optimalrefscore, weights=x.weight))

# %% tags=[]
testscores_original = pd.read_table("../results/testscores_.tsv", index_col=0)
testscores_original = testscores_original.query("organismoi != 'synth'").groupby("method").apply(lambda x: np.average(x.score, weights=x.weight))

trainingscores_original = pd.read_table("../results/trainingscores_.tsv", index_col=0)
trainingscores_original = trainingscores_original.query("organismoi != 'synth'").groupby("method").apply(lambda x: np.average(x.score, weights=x.weight))

proxyscores.loc["training"] = trainingscores_original.loc[proxyscores.columns]
proxyscores.loc["test"] = testscores_original.loc[proxyscores.columns]

# %% jupyter={"outputs_hidden": false} scrolled=true
xs = pd.Series(np.arange(0, len(proxyscores.columns)), index=proxyscores.columns)
fig = figure(figsize=(len(xs)*0.5,4))
ax = fig.add_subplot(111)

markersettings = {
    "training":["#999999", "_", 0],
    "test":["#333333", "_", 0],
    "bhi#go":["#00A6BF", "o", 50],
    "faucodds#go":["#0067bd", "o", 50],
    "bhi#kegg":["#00A6BF", "x", 50],
    "faucodds#kegg":["#0067bd", "x", 50]
}

proxyorder = ["test", "bhi#go", "faucodds#go", "bhi#kegg", "faucodds#kegg", "training"]
proxyscores = proxyscores.loc[proxyorder]
w = 0.8
for method, methodscores in proxyscores.iteritems():
    x = xs[method]
    ax.add_patch(
        Rectangle(
            (x-w/2, methodscores["test"]),
            w,
            methodscores["training"]-methodscores["test"],
            fc="#AAAAAA", zorder=0, lw=0
        )
    )
    for i, (proxygroup, score) in enumerate(methodscores.items()):
        if proxygroup in ["training"]:#, "test"]:
            ""
        else:
            markersetting = markersettings[proxygroup]
            scatter(x, score, color=markersetting[0], marker=markersetting[1], s=markersetting[2], zorder=i+5)
            #scatter(x, score, color=markersetting[0], marker="_", s=250)
            line =  Line2D((x-w/2,x+w/2),(score, score),color=markersetting[0], zorder=i+5, solid_capstyle="butt", lw=1)
            ax.add_artist(line)
        
ax.set_ylim(0, 6.5)
ax.set_xlim(-.6, xs.max()+0.6)
#axis("off")
ax.set_xticks(range(len(proxyscores.columns)))
ax.set_xticklabels(proxyscores.columns, rotation=45, ha="right")
ax.set_yticks(np.arange(0, 6.01, 2))

axhline(1, linestyle="--", color="#555555", lw=1)

ax.tick_params(labelsize=14)
