"""
Script for performing recursive feature elimination for an Extra Trees Classifier, given a labeled feature matrix.
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
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import RFECV
from functools import reduce

def load_data(fmat_file):
    
    # load data
    print(f'Reading in {fmat_file} ...')
    with open(fmat_file, 'rb') as handle:
        fmat = pickle.load(handle)
        
    # print # +/- labels
    print(f"Total # positive labels: {len(fmat[fmat['label']==1])}")
    print(f"Total # negative labels: {len(fmat[fmat['label']==-1])}")
    
    # define cols
    label_cols = ['ID', 'label', 'super_group']
    data_cols = [c for c in fmat.columns.values.tolist() if c not in label_cols]
    
    return(fmat, label_cols, data_cols)

def fmt_arrays(fmat, label_cols, data_cols):

    # convert df cols to arrays
    X = fmat[data_cols].to_numpy()
    y = fmat[label_cols[1]].to_numpy()
    groups = fmat[label_cols[2]].to_numpy()
    
    return(X, y, groups)

def set_rfe_params(step=7, min_feats=1, num_threads=4, seed=None):
    
    # define model to train on
    clf = ExtraTreesClassifier(n_estimators=100, random_state=seed)
    print(f'Selected model: {clf}')
    
    # feature selector (rfe w/ cross-validation)
    cv = StratifiedKFold(5)
    rfecv = RFECV(
        estimator=clf,
        step=step,
        cv=cv,
        scoring="accuracy",
        min_features_to_select=min_feats,
        n_jobs=num_threads)
    
    return(rfecv)


def fit_rfe(model, fold, X_train, y_train):
    
    # run recursive feature elimination
    model.fit(X_train, y_train)
    print(f"Optimal number of features: {model.n_features_}")
    
    return(model)

def plot_results(rfecv_fit, fold, X_test, y_test, outname, outdir):
    
    # plot results
    print(f'Generating RFECV evaluation plots ...')
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
    print(f'Saving mean accuracy plot to {outdir+outname}_nfeats-vs-acc_internal_test.png ...')
    plt.savefig(f'{outdir+outname}_cv-test_nfeats-vs-acc.png', dpi=300, transparent=True)
    
    # PR curve
    from sklearn.metrics import PrecisionRecallDisplay
    PrecisionRecallDisplay.from_estimator(rfecv_fit, X_test, y_test)
    plt.savefig(f'{outdir+outname}_holdout-test_prcurve.png', dpi=300, transparent=True)
    print(f'Saving PR plot to {outdir+outname}_prcurve_holdout_test.png ...')
    
def get_importances(rfecv_fit, data_cols, fold):
    
    # get result table
    print(f'Extracting RFECV selected features & scores ...')
    results = pd.DataFrame({'feature':data_cols, 'rank':rfecv_fit.ranking_, 'support':rfecv_fit.support_})
    sel_feats = results[results['support'] == True]
    sel_feats_scored = sel_feats.head(rfecv_fit.n_features_)
    sel_feats_scored.drop(['support'], axis=1, inplace=True)
    sel_feats_scored.drop(['rank'], axis=1, inplace=True)
    sel_feats_scored['mdi'] = rfecv_fit.estimator_.feature_importances_
    sel_feats_scored = sel_feats_scored.sort_values('mdi')
    sel_feats_scored['fold'] = int(fold+1)

    return(sel_feats_scored)

def main():
    
    # define seed, if any
    if args.seed:
        seed = args.seed
    
    # define split strategy
<<<<<<< HEAD
    gss = GroupShuffleSplit(n_splits=args.num_splits, train_size=float(args.train_size), random_state=args.seed)
=======
    gss = GroupShuffleSplit(n_splits = args.num_splits, train_size=args.train_size, random_state=args.seed)
>>>>>>> 34a20383574397ee9f1d795b4fe573e4242463a1
    # define feature selector
    rfecv_params = set_rfe_params(step=args.remove_per_step, min_feats=1, num_threads=args.threads)
    print(rfecv_params)
    
    # load & format data
    fmat, label_cols, data_cols = load_data(args.featmat)
    X, y, groups = fmt_arrays(fmat, label_cols, data_cols)
    
    # get gss splits for each iteration
    print(f'----- Running recursive feature elimination for {args.num_splits} GSS splits -----')
    fold_df_lst = []
    for i, (train_idx, test_idx) in enumerate(gss.split(X, y, groups)):

        # define test/train splits
        X_train = X[train_idx]
        y_train = y[train_idx]
        X_test = X[test_idx]
        y_test = y[test_idx]
        
        # check split & label balance
        label, counts = np.unique(y_train, return_counts=True)
        label_counts_train = dict(zip(label, counts))
        label, counts = np.unique(y_test, return_counts=True)
        label_counts_test = dict(zip(label, counts))
        print(f'# train PPIs = {len(X_train)}')
        print(f' --> +/- label balance: {label_counts_train}')
        print(f'# test PPIs = {len(X_test)}')
        print(f' --> +/- label balance: {label_counts_test}')
        
        print(f'Executing RFE for GSS split #{i+1} ...')
        # run rfe
        rfecv_fit = fit_rfe(rfecv_params, i, X_train, y_train)
        
        # define output vars
        suffix = 'gssfold'+str(i+1)
        outname = f'featsel_xtrees_{suffix}'
        outdir = args.outdir
        
        # extract & write results
<<<<<<< HEAD
        plot_results(rfecv_fit, i, X_test, y_test, outname, outdir)
=======
        
        plot_metrics(rfecv_fit, i, X_test, y_test, outname=outname, outdir=outdir)
>>>>>>> 34a20383574397ee9f1d795b4fe573e4242463a1
        
        print(f'Writing feature importances for GSS split #{i+1} to {outdir+outname}.csv ...')
        featsel_df = get_importances(rfecv_fit, data_cols, i)
        featsel_df.to_csv(f'{outdir+outname}.csv', index=False)
        fold_df_lst.append(featsel_df)

    print(f'Getting selected feature importances & summary stats across all folds ...')
    all_res = pd.concat(fold_df_lst)
    gb = all_res.groupby(['feature'])
    counts = gb.size().to_frame(name='counts')
    agg_res = (counts
               .join(gb.agg({'fold': lambda x: ', '.join(set(x.astype(str).dropna()))}).rename(columns={'fold': 'gss_folds'}))
               .join(gb.agg({'mdi': 'mean'}).rename(columns={'mdi': 'mean_mdi'}))
               .join(gb.agg({'mdi': 'min'}).rename(columns={'mdi': 'min_mdi'}))
               .join(gb.agg({'mdi': 'max'}).rename(columns={'mdi': 'max_mdi'}))
               .sort_values(['counts', 'mean_mdi'], ascending=[False, False])
               .reset_index()
              )
    feat_intxn = agg_res[agg_res['counts'] == num_splits]
    print(f'# common features: {len(feat_intxn)}')
    print(f'Writing results to {outdir} ...')
    agg_res.to_csv(f'{out_dir}featsel_xtrees_allres.csv', index=False)
    feat_intxn.to_csv(f'{out_dir}featsel_xtrees_intxn.csv', index=False)
    
    

""" When executed from the command line: """
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify feature matrix
    parser.add_argument("-f", "--featmat", help="(Required) Path to labeled, grouped, and pickled feature matrix. PPI ID column is a column named 'ID' and protein names are separated by a space, positive/negative labels are given as 1/-1 in a column named 'label', groups are given in a column named 'super_group'. ")
    
    # specify output directory
    parser.add_argument("-o", "--outdir", action="store", help="(Required) Path to directory to write results.")
    
    # specify number of group-shuffle-splits to test
    parser.add_argument("--num_splits", action="store", default=5, type=int, help="(Optional) Specify number of group-shuffled train/test splits to run recursive feature elimination on (default=5).")
    
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