"""
Script for running the autoML pipeline TPOT, given a labeled feature matrix.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pickle
import numpy as np
from datetime import datetime as dt
import time
from tpot import TPOTClassifier
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import PrecisionRecallDisplay
import matplotlib.pyplot as plt

def fmt_data(fmat_file):
    
    # read in data
    with open(fmat_file, 'rb') as handle:
        fmat = pickle.load(handle)
    
    # get train/test data
    fmat = fmat[fmat['label'].isin([-1,1])]
    fmat.sort_values('label', inplace=True, ascending=False)
    fmat.reset_index(inplace=True, drop=True)
    
    # define cols
    label_cols = ['ID', 'label', 'super_group']
    data_cols = [c for c in fmat.columns.values.tolist() if c not in label_cols]

    # make data, target, and group arrays
    X = fmat[data_cols].to_numpy()
    y = fmat[label_cols[1]].to_numpy()
    groups = fmat[label_cols[2]].to_numpy()
    
    return(X, y, groups)
    

def def_grp_split(method='GroupShuffleSplit', num_splits=5, train_size=0.7, seed=None):
    
    if method == 'GroupShuffleSplit':
        gs = GroupShuffleSplit(n_splits = num_splits, train_size=train_size, random_state=seed)
    elif method == 'GroupKFold':
        gs = GroupKFold(n_splits = num_splits)
    elif method == 'StratifiedGroupKFold':
        gs = StratifiedGroupKFold(n_splits = num_splits)
    else:
        print('Invalid group split strategy specified; please choose one of: "GroupShuffleSplit", "GroupKFold", or "StratifiedGroupKFold"')
        print('Review documentation for more information on each method: https://scikit-learn.org/stable/modules/classes.html#module-sklearn.model_selection')
    
    return(gs)

""" Main """
def main():
    
    t0 = time.time()
    
    print(f'[{dt.now()}] Extracting train/test data from {args.featmat}...')
    X, y, groups = fmt_data(args.featmat)
    
    print(f'[{dt.now()}] Loading complex group split method parameters ...')
    gs = def_grp_split(args.group_split_method, args.num_splits, args.train_size, args.seed)
    
    print(f'[{dt.now()}] Loading TPOTClassifier() pipeline ...')
    pipeline_opt = TPOTClassifier()
    pipeline_opt = TPOTClassifier(generations=args.generations, population_size=args.pop_size, random_state=args.seed, cv=5, verbosity=2)
    
    outfile = args.outdir+'tpot_pipeline'
    model_out = args.outdir+'tpot_model'

    print(f'[{dt.now()}] Running TPOT for {args.num_splits} total {args.group_split_method} splits ...')
    for i, (test_idx, train_idx) in enumerate(gs.split(X, y, groups)):

        X_train = X[train_idx]
        y_train = y[train_idx]
        X_test = X[test_idx]
        y_test = y[test_idx]

        print(f'[{dt.now()}] Running TPOT for split #{i+1} ...')
        print(f"--> # train PPIs = {len(X[train_idx])}")
        print(f"--> # test PPIs = {len(X[test_idx])}")

        pipeline_opt.fit(X_train, y_train) 
        print(f'[{dt.now()}] Test set score for split #{i+1}: {pipeline_opt.score(X_test, y_test)}')
            
        print(f"[{dt.now()}] Writing TPOT results to {outfile+'_'+str(i+1)} ...")
        pipeline_opt.export(outfile+'_'+str(i+1))
        
        print(f"[{dt.now()}] Writing optimized TPOT model object to {model_out+'_'+str(i+1)}.pkl ...")
        model = pipeline_opt.fitted_pipeline_
        with open(f'{model_out}_{i+1}.pkl', 'wb') as f:
            pickle.dump(model,f)
    
    print(f"[{dt.now()}] ---------------------------------------------------------")
    print(f"[{dt.now()}] Total run time: {round((time.time()-t0)/60, 2)} minutes.")
    print(f"[{dt.now()}] ---------------------------------------------------------")

""" When executed from the command line: """
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify feature matrix
    parser.add_argument("-f", "--featmat", help="(Required) Path to labeled, grouped, and pickled feature matrix. PPI ID column is a column named 'ID' and protein names are separated by a space, positive/negative labels are given as 1/-1 in a column named 'label', groups are given in a column named 'super_group'. ")
    
    # specify output directory
    parser.add_argument("-o", "--outdir", action="store", help="(Required) Path to directory to write results.")
    
    # specify group split strategy
    parser.add_argument("--group_split_method", action="store", default="GroupShuffleSplit", help="(Optional) Specify method for splitting protein complex groups; one of: 'GroupShuffleSplit', 'GroupKFold', or 'StratifiedGroupKFold' (default=GroupShuffleSplit)")
    
    # specify number of group splits to test for cross-validation
    parser.add_argument("--num_splits", action="store", default=5, type=int, help="(Optional) Specify number of grouped cross-validation splits to run recursive feature elimination on (default=5). Set to 1 if you only want to run TPOT one time on a single split of the data.")
    
    # specify proportion of data to train on for every group split
    parser.add_argument("--train_size", action="store", type=float, default=0.7, help="(Optional) Specify fraction of data to train on for every group split iteration (default=0.7).")
    
    # specify # of TPOT generations
    parser.add_argument("--generations", action="store", type=int, default=10, help="(Optional) Specify number of TPOT generations, or the number of iterations to run the pipeline optimization process (default=10).")
    
    # specify TPOT population size
    parser.add_argument("--pop_size", action="store", type=int, default=20, help="(Optional) Specify TPOT population size, or the number of individuals to retain in the GP population every generation. Generally, TPOT will work better when you give it more individuals (and therefore time) to optimize the pipeline. (default=20).")
    
    # specify seed to make results consistent
    parser.add_argument("--seed", action="store", type=int, default=None, help="(Optional) Specify seed to make group and RFECV train/test splits reproducible.")

    args = parser.parse_args()
    main()