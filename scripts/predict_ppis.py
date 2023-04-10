"""
Script for predicting protein interactions, given a labeled feature matrix and a pre-optimized model object.
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
import pandas as pd
from functools import reduce
from beautifultable import BeautifulTable
from sklearn.ensemble import *
from sklearn.model_selection import *
from sklearn.metrics import *
from sklearn.preprocessing import *
from sklearn.pipeline import make_pipeline

""" Functions for model set up """

def get_features(fmat, label_cols, feature_file=None, n_feats2sel=None):
    all_cols = fmat.columns.values.tolist()
    if not feature_file:
        data_cols = [c for c in all_cols if c not in label_cols]
    else:
        fsel = pd.read_csv(feature_file)
        select_feats = fsel['feature'].head(int(n_feats2sel)).tolist()
        data_cols = [c for c in all_cols if c not in label_cols and c in select_feats]
    return(data_cols)
    
def fmt_data(fmat, subset, label_cols, data_cols, keep_groups=True):
    # get desired subset
    fmat_fmt = fmat[fmat['label'].isin(subset)]
    fmat_fmt.reset_index(inplace=True, drop=True)
    # format data for sklearn input
    X = fmat_fmt[data_cols].to_numpy()
    y = fmat_fmt[label_cols[1]].to_numpy()
    ids = fmat_fmt[label_cols[0]]
    # keep/drop group col option
    if keep_groups:
        groups = fmat_fmt[label_cols[2]].to_numpy()
        return(X, y, ids, groups)
    else:
        return(X, y, ids)
    
def fmt_outdir(outdir):
    if outdir.endswith('/'):
        return(outdir)
    else:
        return(outdir+'/')

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

def write_results(all_res, test_scores, model_name, outdir, fdr_cutoff):
    # format path, file names
    outdir_fmt = fmt_outdir(outdir)
    fdr_fmt = int(fdr_cutoff*100)
    bases = ['scored_interactions_all', f'scored_interactions_fdr{fdr_fmt}', 'precision_recall' ]
    all_out = outdir_fmt+f'{bases[0]}_{model_name}.csv'
    high_conf_out = outdir_fmt+f'{bases[1]}_{model_name}.csv'
    pr_out = outdir_fmt+f'{bases[2]}_{model_name}.csv'
    # format all scored PPIs
    all_res.sort_values('ppi_score', inplace=True, ascending=False)
    all_res.reset_index(inplace=True, drop=True)
    print(f' ► Writing all scored PPIs to {all_out} ...')
    all_res.to_csv(all_out, index=False)
    # get PR curve & threshold PPIs at given FDR
    pr_curve, high_conf_ppis = threshold_ppis(test_scores, all_res, fdr_cutoff)
    print(f' ► Writing precision-recall results to {pr_out} ...')
    pr_curve.to_csv(pr_out, index=False
                   )
    print(f' ► Writing results to {high_conf_out} ...')
    high_conf_ppis.to_csv(high_conf_out, index=False)

""" Functions for model fitting & evaluation """

def calc_pr(df):
    # compute precision/recall
    print(f" ► Computing precision/recall ...")
    tp_count = 0
    fp_count = 0
    p_list = []
    r_list = []
    f_list = []
    all_pos = len(df[df['label'] == 1])
    for i in range(len(df)):
        if df['label'][i] == 1:
            tp_count += 1
        else:
            fp_count += 1
        tps = tp_count
        fps = fp_count
        fns = all_pos - tps
        precision = tps/(tps+fps)
        recall = tps/(tps+fns)
        fdr = 1 - precision
        p_list.append(float(precision))
        r_list.append(float(recall))
        f_list.append(float(fdr))
    df['precision'] = p_list
    df['recall'] = r_list
    df['fdr'] = f_list
    return(df)

def get_scores(model, array, labels, ids):
    # get scores
    try:
        scores = model.predict_proba(array)
    except:
        scores = model.predict_proba_lr(array)
    probabilities = np.split(scores, 2, axis=1)
    neg_prob = probabilities[0]
    pos_prob = probabilities[1]
    # format into df
    df = pd.DataFrame()
    df['ID'] = ids
    df['label'] = labels
    df['ppi_score'] = pos_prob
    df.sort_values('ppi_score', inplace=True, ascending=False)
    df.reset_index(inplace=True, drop=True)
    return(df)

def threshold_ppis(test_scores, all_res, fdr_cutoff):
    # compute precision/recall
    test_scores_pr = calc_pr(test_scores)
    thres_df = test_scores_pr[test_scores_pr['fdr'] <= fdr_cutoff]
    prob_cutoff = min(thres_df['ppi_score'])
    print(f' ► PPI score cutoff for {int(fdr_cutoff*100)}% FDR: {prob_cutoff}')
    # theshold results
    thres_df = all_res[all_res['ppi_score'] >= prob_cutoff]
    ids = thres_df['ID'].str.split(' ', expand=True)
    print(f' ► # total PPIs evaluated: {len(all_res)}')
    print(f' ► # PPIs above threshold: {len(thres_df)}')
    try: # only works if there are no self-self PPIs
        uniq_prots = np.unique(ids[[0, 1]].values)
        print(f' ► # unique proteins above threshold: {len(uniq_prots)}')
    except:
        print('WARNING: Problem with IDs detected.')
        print(ids[0:51])
    df_out = thres_df[['ID','ppi_score','set']]
    return(test_scores_pr, df_out)
    

""" Main """
def main():
    
    t0 = time.time()
    ## read in data & define model params
    print(f'[{dt.now()}] Loading feature matrix ...')
    with open(args.featmat, 'rb') as handle:
        fmat = pickle.load(handle)
    
    label_cols = ['ID', 'label', 'super_group']
    data_cols = get_features(fmat, label_cols, args.feature_selection, args.num_feats)
    X, y, ids, groups = fmt_data(fmat, [1,-1], label_cols, data_cols, keep_groups=True)
    X_pred, y_pred, ids_pred = fmt_data(fmat, [0], label_cols, data_cols, keep_groups=False)
    
    print(f'[{dt.now()}] Loading complex group split method parameters ...')
    gs = def_grp_split(args.group_split_method, args.num_splits, args.train_size, args.seed)
    
    print(f'[{dt.now()}] Loading model parameters ...')
    model, model_name = load_model(args.model, args.seed)
    
    table = BeautifulTable()
    table.columns.header = ["PARAMETER", "SETTING"]
    table.rows.append(["Protein complex split method", args.group_split_method])
    table.rows.append(["# train/test (CV) splits", args.num_splits])
    table.rows.append(["Machine learning method", model_name])
    table.rows.append(["% data for training (approx*):", f'{int(args.train_size*100)}%'])
    table.rows.append(["% data for testing (approx*):", f'{int((1-args.train_size)*100)}%'])
    table.rows.append(["% FDR threshold", f'{int(args.fdr_cutoff*100)}%'])
    table.rows.append(["# features to be used", len(data_cols)])
    
    print()
    print(f'[{dt.now()}] Pipeline settings detected:')
    print(table)
    print('*Note: Group split method will attempt to get as close to these settings as possible, but results may vary depending on the seed, group sizes, and number of cross-validation splits specified.')
    
    ## fit model, compute precision/recall, and output results
    for i, (train_idx, test_idx) in enumerate(gs.split(X, y, groups)):
        
        # get train/test splits
        print()
        print(f"[{dt.now()}] Getting test/train splits ({i+1}/{args.num_splits}) ...")
        X_train, y_train, X_test, y_test = split_data(X, y, train_idx, test_idx)
          
        # fit model
        print(f"[{dt.now()}] Fitting {model_name} ...")
        
        # check if model has predict probabilities functionality
        if 'predict_proba' in dir(model[-1]):
            model.fit(X_train, y_train)
        # if not, add it and refit (required for SVMs)
        else:
            from sklearn.calibration import CalibratedClassifierCV
            model = CalibratedClassifierCV(model) 
            model.fit(X_train, y_train)
          
        # extract results
        test_ids = ids[test_idx]
        train_ids = ids[train_idx]
        pred_ids = ids_pred
        
        if len(test_ids) > len(train_ids) and args.num_splits == 1:
            print("ERROR: Imbalanced train/test split. Try a different seed using the --seed argument.")
            return
          
        print(f"[{dt.now()}] Getting probability scores ...")
        test_scores = get_scores(model, X_test, y_test, test_ids)
        train_scores = get_scores(model, X_train, y_train, train_ids)
        pred_scores = get_scores(model, X_pred, y_pred, pred_ids)
        
        # format & write results
        test_scores['set'] = 'test'
        train_scores['set'] = 'train'
        pred_scores['set'] = 'predict'
        
        if args.num_splits > 1:
            out_name = f'{model_name}{i+1}'
        else:
            out_name = model_name
          
        all_res = pd.concat([test_scores, train_scores, pred_scores])
        write_results(all_res, test_scores, out_name, args.outdir, args.fdr_cutoff)
        
    print(f"[{dt.now()}] ---------------------------------------------------------")
    print(f"[{dt.now()}] Total run time: {round((time.time()-t0)/60, 2)} minutes.")
    print(f"[{dt.now()}] ---------------------------------------------------------")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify feature matrix
    parser.add_argument("-f", "--featmat", help="(Required) Path to labeled, grouped, and pickled feature matrix. PPI ID column is a column named 'ID' and protein names are separated by a space, positive/negative labels are given as 1/-1 in a column named 'label', groups are given in a column named 'super_group' (as output by label_featmat.py).")
    
    # specify model
    parser.add_argument("-m", "--model", action="store", help="(Required) Path to (pickle) file containing a scikit-learn model object with pre-optimized parameters (e.g., as output by run_tpot.py).")
    
    # specify output directory
    parser.add_argument("-o", "--outdir", action="store", help="(Required) Path to directory to write results.")
    
    # specify feature selection file
    parser.add_argument("--feature_selection", action="store", default=None, help="(Optional) Path to CSV file output by select_features.py where the features have been arranged from most important to least important.")
    
    # specify number features to select
    parser.add_argument("--num_feats", action="store", default=None, help="(Optional**) Required if feature_selection file specified. Number of features to select from the feature selection file.")
    
    # specify group split strategy
    parser.add_argument("--group_split_method", action="store", default="GroupShuffleSplit", help="(Optional) Specify method for splitting protein complex groups; one of: 'GroupShuffleSplit', 'GroupKFold', or 'StratifiedGroupKFold' (default=GroupShuffleSplit)")
    
    # specify number of group splits to test for cross-validation
    parser.add_argument("--num_splits", action="store", default=1, type=int, help="(Optional) Specify number of group-based splits and models to generate (default=1). If group_split_method is set to GroupKFold or StratifiedGroupKFold, can be used for Kfold cross-validation where K=num_splits.")
    
    # specify proportion of data to train on for every group split
    parser.add_argument("--train_size", action="store", type=float, default=0.7, help="(Optional) Specify fraction of data to train on for every group split iteration (default=0.7).")
    
    # specify false discovery rate threshold
    parser.add_argument("--fdr_cutoff", action="store", type=float, default=0.1, help="(Optional) Specify FDR threshold for PPI predictions. Accepts a number between 0 and 1, e.g., 0.001, 0.01, and 0.1 translate to a 0.1%, 1% and 10% FDR, respectively (default=0.1).")
    
    # specify seed to make results consistent
    parser.add_argument("--seed", action="store", type=int, default=None, help="(Optional) Specify seed to make train/test splits reproducible.")

    args = parser.parse_args()
    main()