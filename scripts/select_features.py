"""
Script for performing recursive feature elimination for a pickled scikit-learn model, given a labeled feature matrix.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import pandas  as pd
import numpy as np
import pickle
import argparse
import random
import time
from datetime import datetime as dt
from sklearn.ensemble import *
from sklearn.model_selection import *
from sklearn.metrics import *
from sklearn.preprocessing import *
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import RFECV
from functools import reduce

""" Functions """
def load_data(fmat_file):
    # load data
    with open(fmat_file, 'rb') as handle:
        fmat = pickle.load(handle) 
    # print # +/- labels
    print(f" ► Total # positive labels: {len(fmat[fmat['label']==1])}")
    print(f" ► Total # negative labels: {len(fmat[fmat['label']==-1])}")
    # define cols
    label_cols = ['ID', 'label', 'super_group']
    data_cols = [c for c in fmat.columns.values.tolist() if c not in label_cols]
    return(fmat, label_cols, data_cols)

def load_model(model_file, seed=None):
    with open(model_file, 'rb') as f:
        model = pickle.load(f)
    try:
        n_steps = len(model.named_steps)
        idx = n_steps-1
        model_name = type(model[idx]).__name__
        if seed:
            setattr(model[idx], 'random_state', seed)
    except:
        model_name = type(model).__name__
        if seed:
            try:
                setattr(model, 'random_state', seed)
            except:
                print(f'WARNING: failed to set random state for {model} to {seed}!')
    return(model, model_name)

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

def fmt_outdir(outdir):
    if outdir.endswith('/'):
        return(outdir)
    else:
        return(outdir+'/')

def fmt_arrays(fmat, label_cols, data_cols):
    # convert df cols to arrays
    X = fmat[data_cols].to_numpy()
    y = fmat[label_cols[1]].to_numpy()
    groups = fmat[label_cols[2]].to_numpy()
    return(X, y, groups)

def split_data(X, y, train_idx, test_idx):
    # get splits
    X_train = X[train_idx]
    y_train = y[train_idx]
    X_test = X[test_idx]
    y_test = y[test_idx]
    # check split & label balance
    label, counts = np.unique(y_train, return_counts=True)
    label_counts_train = dict(zip(label, counts))
    label, counts = np.unique(y_test, return_counts=True)
    label_counts_test = dict(zip(label, counts))
    print(f' ► # train PPIs = {len(X_train)}')
    print(f' ► train +/- label counts: {label_counts_train}')
    print(f' ► # test PPIs = {len(X_test)}')
    print(f' ► test +/- label counts: {label_counts_test}')
    return(X_train, y_train, X_test, y_test)

def set_rfe_params(model, step=7, min_feats=1, num_threads=4, seed=None):
    # extract estimator if TPOT pipeline passed to script
    try:
        n_steps = len(model.named_steps)
        idx = n_steps-1
        est = model[idx]
    except AttributeError:
        est = model
    print(f'Selected method: {est}')
    # feature selector (rfe w/ cross-validation)
    cv = StratifiedKFold(5)
    rfecv = RFECV(
        estimator=est,
        step=step,
        cv=cv,
        scoring="accuracy",
        min_features_to_select=min_feats,
        n_jobs=num_threads)
    
    print(f'Selected RFECV parameters:')
    print(rfecv)
    return(rfecv)

def fit_rfe(model, X_train, y_train):
    # run recursive feature elimination
    model.fit(X_train, y_train)
    print(f"Optimal number of features: {model.n_features_}")
    return(model)
        
def plot_results(rfecv_fit, fold, X_test, y_test, outname, outdir):
    # plot results
    print(f' ► Generating RFECV evaluation plots ...')
    min_features_to_select=1
    # mean test accuracy vs number of feats selected
    import matplotlib.pyplot as plt
    n_scores = len(rfecv_fit.cv_results_["mean_test_score"])
    plt.figure()
    plt.xlabel(f"Number of features selected\n(optimal={rfecv_fit.n_features_} features)")
    plt.ylabel("Mean CV test accuracy")
    plt.errorbar(
        range(min_features_to_select, n_scores + min_features_to_select),
        rfecv_fit.cv_results_["mean_test_score"],
        yerr=rfecv_fit.cv_results_["std_test_score"],
    )
    plt.title(f"Recursive feature elimination\nwith correlated features (fold #{fold+1})")
    print(f' ► Saving mean accuracy plot to {outdir+outname}_nfeats-vs-acc_internal_cvtest.png ...')
    plt.savefig(f'{outdir+outname}_nfeats-vs-acc_internal_cvtest.png', dpi=300, transparent=True)
    # PR curve
    from sklearn.metrics import PrecisionRecallDisplay
    PrecisionRecallDisplay.from_estimator(rfecv_fit, X_test, y_test)
    print(f' ► Saving PR plot to {outdir+outname}_prcurve_holdout_test.png ...')
    plt.savefig(f'{outdir+outname}_prcurve_holdout_test.png', dpi=300, transparent=True)
    
def get_importances(rfecv_fit, data_cols, fold):
    # get result table
    results = pd.DataFrame({'feature':data_cols, 'rank':rfecv_fit.ranking_, 'support':rfecv_fit.support_})
    sel_feats = results[results['support'] == True]
    sel_feats_scored = sel_feats.head(rfecv_fit.n_features_)
    sel_feats_scored = sel_feats_scored.drop(['support'], axis=1)
    sel_feats_scored = sel_feats_scored.drop(['rank'], axis=1)
    try:
        sel_feats_scored['mdi'] = rfecv_fit.estimator_.feature_importances_
        sel_feats_scored = sel_feats_scored.sort_values('mdi')
    except AttributeError:
        coef_array = rfecv_fit.estimator_.coef_
        sel_feats_scored['coef'] = coef_array.flatten()
        sel_feats_scored = sel_feats_scored.sort_values(by='coef', key=abs, ascending=False)
    sel_feats_scored['fold'] = int(fold+1)
    print(sel_feats_scored)
    return(sel_feats_scored)

""" Main """
def main():
    
    t0 = time.time()
    # load & format data
    print(f'[{dt.now()}] Loading feature matrix ...')
    fmat, label_cols, data_cols = load_data(args.featmat)
    X, y, groups = fmt_arrays(fmat, label_cols, data_cols)
    
    print(f'[{dt.now()}] Loading complex group split method parameters ...')
    gs = def_grp_split(args.group_split_method, args.num_splits, args.train_size, args.seed)
    
    print(f'[{dt.now()}] Loading model parameters ...')
    model, model_name = load_model(args.model, args.seed)

    # define feature selector
    print(f'[{dt.now()}] Loading recursive feature selection parameters ...')
    rfecv_params = set_rfe_params(model, step=args.remove_per_step,
                                  min_feats=1, num_threads=args.threads)
    
    # get gss splits for each iteration
    print(f"[{dt.now()}] ---------------------------------------------------------")
    print(f'[{dt.now()}] Running recursive feature elimination for ')
    print(f'[{dt.now()}] {args.num_splits} {args.group_split_method} train/test sets.')
    print(f"[{dt.now()}] ---------------------------------------------------------")
    
    fold_df_lst = []
    optimal_n_dict = dict()
    
    # TODO: organize blocks into functions to make main() more readable
    for i, (train_idx, test_idx) in enumerate(gs.split(X, y, groups)):
        
        print(f'[{dt.now()}] Executing RFE for train/test split #{i+1} ...')
        
        # get train/test splits
        print(f"Getting test/train splits ({i+1}/{args.num_splits}) ...")
        X_train, y_train, X_test, y_test = split_data(X, y, train_idx, test_idx)
          
        # run rfe
        try:
            rfecv_fit = fit_rfe(rfecv_params, X_train, y_train)
        except:
            print(f'ERROR: {model_name} does not provide logic for feature selection.')
            print('Please provide a model with a feature importance attribute (examples: ExtraTreesClassifier, LinearSVC, SGDClassifier).')
            return

        # define output vars
        suffix = 'fold'+str(i+1)
        outname = f'top_feats_{model_name}_RFECV_{suffix}'
        outdir = fmt_outdir(args.outdir)
        
        # extract & write results
        plot_results(rfecv_fit, i, X_test, y_test, outname=outname, outdir=outdir)
        
        print(f'[{dt.now()}] Writing feature importances for train/test split #{i+1} to {outdir+outname}.csv ...')
        featsel_df = get_importances(rfecv_fit, data_cols, i)
        featsel_df.to_csv(f'{outdir+outname}.csv', index=False)
        
        fold_n = len(featsel_df)
        optimal_n_dict.update({int(i+1): int(fold_n)})
        fold_df_lst.append(featsel_df)
        print()

    print(f'[{dt.now()}] Getting selected feature importances & summary stats across all folds ...')
    all_res = pd.concat(fold_df_lst)
    gb = all_res.groupby(['feature'])
    counts = gb.size().to_frame(name='counts')
    
    if 'mdi' in all_res.columns:
        score = 'mdi'
    else:
        score = 'coef'
        
    agg_res = (counts
               .join(gb.agg({'fold': lambda x: ', '.join(set(x.astype(str).dropna()))}))
               .join(gb.agg({f'{score}': 'mean'}).rename(columns={f'{score}': f'mean_{score}'}))
               .join(gb.agg({f'{score}': 'min'}).rename(columns={f'{score}': f'min_{score}'}))
               .join(gb.agg({f'{score}': 'max'}).rename(columns={f'{score}': f'max_{score}'}))
               .join(gb.agg({f'{score}': 'std'}).rename(columns={f'{score}': f'stdev'}))
               #.sort_values(['counts', 'stdev', f'mean_{score}'], key=abs, ascending=[False, True, False])
               #.reset_index()
              )
    
    agg_res['rsd'] = agg_res.stdev/agg_res[f'mean_{score}']
    agg_res = (agg_res
               .assign(rsd_sort=agg_res['rsd'].abs(), score_sort=agg_res[f'mean_{score}'].abs())
               .sort_values(['counts', 'rsd_sort', 'score_sort'],ascending=[False, True, False])
               .drop(['rsd_sort', 'score_sort'], 1)
               .reset_index()
              )
    
    feat_intxn = agg_res[agg_res['counts'] == args.num_splits]
    
    # get & format results for optimal # of features
    all_n_feats = list(optimal_n_dict.values())
    avg_n_feats = sum(all_n_feats)/len(all_n_feats)
    optimal_n_dict.update({'avg': avg_n_feats})
    opt_n_df = pd.DataFrame(list(optimal_n_dict.items()), columns=['fold', 'num_optimal_features'])
    
    print(f'[{dt.now()}] # common features: {len(feat_intxn)}')
    print(f'[{dt.now()}] # optimal features averaged across all folds: {avg_n_feats}')
    print(f'[{dt.now()}] Writing summary results to {outdir} ...')
    agg_res.to_csv(outdir+f'top_feats_{model_name}_RFECV_all.csv', index=False)
    feat_intxn.to_csv(outdir+f'top_feats_{model_name}_RFECV_intersection.csv', index=False)
    opt_n_df.to_csv(outdir+f'top_feats_{model_name}_RFECV_optimalnum.csv', index=False)
    
    print(f"[{dt.now()}] ---------------------------------------------------------")
    print(f"[{dt.now()}] Total run time: {round((time.time()-t0)/60, 2)} minutes.")
    print(f"[{dt.now()}] ---------------------------------------------------------")

""" When executed from the command line: """
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify feature matrix
    parser.add_argument("-f", "--featmat", help="(Required) Path to labeled, grouped, and pickled feature matrix. PPI ID column is a column named 'ID' and protein names are separated by a space, positive/negative labels are given as 1/-1 in a column named 'label', groups are given in a column named 'super_group'. ")
    
    # specify model
    parser.add_argument("-m", "--model", action="store", help="(Required) Path to (pickle) file containing a scikit-learn model object with pre-optimized parameters (e.g., as output by run_tpot.py).")
    
    # specify output directory
    parser.add_argument("-o", "--outdir", action="store", help="(Required) Path to directory to write results.")
    
    # specify group split strategy
    parser.add_argument("--group_split_method", action="store", default="GroupShuffleSplit", help="(Optional) Specify method for splitting protein complex groups; one of: 'GroupShuffleSplit', 'GroupKFold', or 'StratifiedGroupKFold' (default=GroupShuffleSplit)")
    
    # specify number splits to test
    parser.add_argument("--num_splits", action="store", default=5, type=int, help="(Optional) Specify number of train/test splits to run recursive feature elimination on (default=5).")
    
    # specify proportion of data to train on for every GSS
    parser.add_argument("--train_size", action="store", type=float, default=0.7, help="(Optional) Specify fraction of data to train on for every group-shuffle-split iteration (default=0.7).")
    
    # specify number of features to eliminate per round of RFE
    parser.add_argument("--remove_per_step", action="store", type=int, default=1, help="(Optional) Specify number of features to remove at every step of recursive feature elimination (default=1).")
    
    # specify number threads to give RFE
    parser.add_argument("--threads", action="store", type=int, default=4, help="(Optional) Specify number of threads to use (default=4).")
    
    # specify seed to make results consistent
    parser.add_argument("--seed", action="store", type=int, default=None, help="(Optional) Specify seed to make GSS and RFECV train/test splits reproducible.")

    args = parser.parse_args()
    main()