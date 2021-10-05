import sys
from pathlib import Path
from datetime import datetime
import hashlib

import numpy as np
import pandas as pd

from clustermatch.corr import clustermatch


DATA_HASH_DIR = Path(__file__).parent / "data_hash"
DATA_HASH_FILENAME_TEMPLATE = "{hash}.pkl"


def standardize(X):
    return (X - X.mean())/(X.std())

###### Distance/similarity functions #######

def cal_triangular(E, func):
    """
    Calculates a given function for every gene pair in E
    """
    Emat = E.T.values
    correlations = np.identity(Emat.shape[0])
    todo = Emat.shape[0] * Emat.shape[0] / 2 - Emat.shape[0]
    for i, x in enumerate(Emat):
        row = []
        start = datetime.now()
        sys.stdout.write(str(i) + "/"+str(Emat.shape[0]))
        for y in Emat[i+1:]:
            sys.stdout.write(">")
            row.append(func(x, y))

        correlations[i, 1+i:] = row
        correlations[1+i:, i] = row

        todo -= len(row)
        if len(row) > 20:
            sys.stdout.write(
                "\nETA: " + str(todo * (datetime.now() - start)/len(row))+
                "\nPer: " + str((datetime.now() - start)/len(row)))

        sys.stdout.write("\n")
    return correlations


def get_func_output(E, func):
    """
    TODO
    """
    DATA_HASH_DIR.mkdir(exist_ok=True)

    data_hash = hashlib.sha256(pd.util.hash_pandas_object(E, index=True).to_numpy()).hexdigest()
    filepath = DATA_HASH_DIR / DATA_HASH_FILENAME_TEMPLATE.format(hash=hash)

    if filepath.exists():
        return pd.read_pickle(filepath)

    # compute the similarity matrix
    simdist_matrix = func(E)
    simdist_matrix.to_pickle(filepath)

    return simdist_matrix



def simdist(E, simdist_function, similarity=True, **kwargs):
    choices = {
        "pearson_correlation":[True, lambda E: np.corrcoef(E.T)],
        "pearson_correlation_absolute":[True, lambda E: np.abs(np.corrcoef(E.T))],
        "clustermatch": [True, lambda E: clustermatch(E.T)],
        "clustermatch_linear": [True, lambda E: clustermatch(E.T, internal_n_clusters=list(range(2, 2 + 1)))],
    }

    measure_similarity, func = choices[simdist_function]

    if simdist_function in ("clustermatch", "clustermatch_linear"):
        simdist_matrix = get_func_output(E, func)
    else:
        simdist_matrix = func(E)
        simdist_matrix = pd.DataFrame(simdist_matrix, columns=E.columns, index=E.columns)

    if (measure_similarity and similarity) or (not measure_similarity and not similarity):
        ""
    else:
        simdist_matrix =  (-simdist_matrix) + simdist_matrix.max().max()
        print(simdist_matrix)
    return simdist_matrix

def _cal_dcor(x, y):
    A = _cal_A(x)
    Avar = _cal_dvar(A)
    B = _cal_A(y)
    Bvar = _cal_dvar(B)

    return _cal_dcov(A, B)/(np.sqrt(Avar * Bvar))

def _cal_dcor_row_star(args):
    return _cal_dcor_row(*args)
def _cal_dcor_row(x, Y, i):
    A = _cal_A(x)
    Avar = _cal_dvar(A)

    start = datetime.now()
    dcors = []
    for y in Y[i:]:
        B = _cal_A(y)
        Bvar = _cal_dvar(B)

        dcors.append(_cal_dcov(A, B)/(np.sqrt(Avar * Bvar)))

    if len(Y[i:]) > 0:
        print(len(Y[i:]))
        print((datetime.now() - start).total_seconds() / len(Y[i:]))

    return dcors

def _cal_A(x):
    d = np.abs(x[:, None] - x)
    return d - d.mean(0) - d.mean(1)[:,None] + d.mean()
def _cal_dvar(A):
    return np.sqrt((A**2).sum() / len(A)**2)
def _cal_dcov(A,B):
    return np.sqrt((A * B).sum() / len(A)**2)
