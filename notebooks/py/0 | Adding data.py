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
# Note: if you have downloaded the results and data folder from zenodo, you do not have to run this notebook.

# %% [markdown]
# # Adding a new dataset

# %%
import pandas as pd

# %% [markdown]
# ## Gene expression data

# %%
# ! mkdir -p ../data/ecoli_precise2

# %%
# ! wget https://raw.githack.com/SBRG/precise2/master/data/precise2/log_tpm_norm_qc.csv -O ../data/ecoli_precise2/E_raw.csv

# %%
# ! wget https://raw.githack.com/SBRG/precise2/master/data/precise2/gene_info.csv  -O ../data/ecoli_precise2/gene_info.csv

# %%
# ! head ../data/ecoli_precise2/gene_info.csv

# %%
E = pd.read_csv("../data/ecoli_precise2/E_raw.csv", index_col = 0)

# %%
gene_info = pd.read_csv("../data/ecoli_precise2/gene_info.csv", index_col = 0)

# %%
gene_info.index

# %%
E.index = gene_info["gene_name"][E.index]

# %%
E.T.to_csv("../data/ecoli_precise2/E.tsv", sep = "\t")

# %% [markdown]
# ## Download and process most recent regulondb network

# %%
# ! mkdir -p ../data/ecoli_precise2/knownmodules

# %%
# ! wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/network_tf_gene.txt -O ../data/ecoli_precise2/knownmodules/regulondb_raw.txt

# %%
# ! head -n 40 ../data/ecoli_precise2/knownmodules/regulondb_raw.txt

# %%
regnet = pd.read_table("../data/ecoli_precise2/knownmodules/regulondb_raw.txt", comment = "#", names = ["tf_id", "tf_name", "target_id", "target_name", "effect", "evidence", "evidence_type", "-"])
regnet = regnet.query("target_name in @E.index")

# %%
modules = []
for tf_name, tf_regnet in regnet.groupby("tf_name"):
    if tf_regnet.shape[0] > 5:
        modules.append(tf_regnet["target_name"].tolist())

# %%
len(modules)

# %%
import json

# %%
# !mkdir -p ../data/ecoli_precise2/knownmodules/ecoli_regulondb

# %%
json.dump(modules, open("../data/ecoli_precise2/knownmodules/ecoli_regulondb/minimal.json", "w"))

# %%
tfs = list(set([tf[0].lower() + tf[1:] for tfs in regnet["tf_name"].unique() for tf in tfs.split("-")]) & set(E.index))

json.dump(tfs, open("../data/ecoli_precise2/regulators.json", "w"))

# %% [markdown]
# ## Write configuration file

# %%
config = {
    "location": "conf/datasets/ecoli_precise2.json",
    "knownmodules": {
        "ecoli_regulondb": {
            "minimal": "data/ecoli_precise2/knownmodules/ecoli_regulondb/minimal.json"
        }
    },
    "gsets": {
    },
    "regulators": "data/ecoli_colombos/regulators.json",
    "binding": {
    },
    "datasetname": "ecoli_precise2",
    "regnets": [
        "ecoli_regulondb"
    ],
    "params": {
        "datasetname": "ecoli_precise2",
        "organism": "ecoli"
    },
    "baselinename": "ecoli_precise2",
    "expression": "data/ecoli_precise2/E.tsv"
}
json.dump(config, open("../conf/datasets/ecoli_precise2.json", "w"))
