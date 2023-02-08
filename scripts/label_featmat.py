"""
Script for labeling a PPI feature matrix given a gold standard PPI file (one protein complex per line). Sorts complexes into non-redundant super groups for cross-validation.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pandas as pd
import pickle
import numpy as np
import re
import random
import time
import datetime as dt
from itertools import combinations
from itertools import accumulate
from collections import defaultdict

"""
Functions for iterative row-wise labeling:
"""

# label positive & negative PPIs
def match_label(pair, pos_pairs, neg_ppi_dict):
    if pair in pos_pairs:
        return(1)
    elif pair in neg_ppi_dict.keys():
        return(-1)
    else:
        return(0)

# label protein complex groups
def match_group(pair, pos_ppi_dict, neg_ppi_dict):
    group_list = []
    for k, v in pos_ppi_dict.items():
        if pair in v:
            group_list.append(k)
    if len(group_list) > 1:
        return(group_list)
    elif len(group_list) == 1:
        return(group_list[0])
    elif pair in neg_ppi_dict.keys():
        return(neg_ppi_dict.get(pair))
    else:
        return(0)

# label non-redunant protein complex groups
def match_spr_grp(group, super_grp_dict):
    if type(group) == list:
        return(super_grp_dict.get(group[0]))
    else:
        return(super_grp_dict.get(group, group))
    

"""
Functions for making positive, negative, & group label dictionaries:
"""

# get positive labels
def make_pos_dict(gs_file):
    print(f'[{dt.datetime.now()}] Generating grouped positive PPI labels from gold standard complexes ...')
    pos_ppi_dict = dict()
    group_no = 1
    dupes = []
    with open(gs_file, 'r') as f:
        ppis = f.read().splitlines() 
        for p in ppis:
            ogs = p.split(' ')
            fsets = [frozenset({i, j}) for i,j in list(combinations(ogs, 2))]
            pos_ppi_dict.update({int(group_no): fsets})
            group_no += 1
    print(f'[{dt.datetime.now()}] Finished generating positive PPI labels!')
    return(pos_ppi_dict)

# make negative labels
# --> TODO: rework to get negative PPIs from observed PPIs in input featmat
#def make_neg_dict(pos_dict, fmat_file, target=20000):
def make_neg_dict(pos_dict):
    # make negative PPIs from positive PPIs
    print(f'[{dt.datetime.now()}] Getting random proteins from positive PPIs to generate negative PPIs ...')
    random_prots = set()
    for group_no, fsets in pos_dict.items():
        prot_set = set()
        if len(fsets) > 1:
            for pair in fsets:
                prot_set.add(list(pair)[0])
                prot_set.add(list(pair)[1])
            neg_prots = random.sample(list(prot_set), 3)
            for p in neg_prots:
                random_prots.add(p)
    print(f'[{dt.datetime.now()}] Generating negative PPIs ...')
    neg_ppis = [frozenset({i, j}) for i,j in list(combinations(random_prots, 2))]

    # remove any overlap between neg & pos PPIs
    print(f'[{dt.datetime.now()}] Locating overlap between negative PPIs & positive PPIs ...')
    t0 = time.time()
    all_pos_ppis = list(pos_dict.values())
    flat_pos_ppis = [pair for pair_list in all_pos_ppis for pair in pair_list]
    overlap = set(neg_ppis).intersection(set(flat_pos_ppis))
    overlap_count = len(overlap)
    print(f'[{dt.datetime.now()}] Removing overlap between negative PPIs & positive PPIs ...')
    for i in overlap:
        neg_ppis.remove(i)
    print(f'[{dt.datetime.now()}] # overlapping negative PPIs found & removed = {overlap_count}; total time: {time.time() - t0} seconds)')
    
    
    
    # get negative PPI splits
    num_groups = len(pos_dict)
    print(f'[{dt.datetime.now()}] Randomly splitting negative PPIs into {num_groups} groups ...')
    neg_cmplx_sizes = [random.randint(2, 30) for x in range(num_groups)]
    neg_ppi_grouped = [neg_ppis[x - y: x] for x, y in zip(
            accumulate(neg_cmplx_sizes), neg_cmplx_sizes)]
    print(f'[{dt.datetime.now()}] # of negative PPI groups =', len(neg_ppi_grouped))
    
    # get negative PPI groups
    print(f'[{dt.datetime.now()}] Generating grouped negative PPI labels ...')
    neg_ppi_dict = dict()
    group_sizes = []
    group_no = 1
    for group in neg_ppi_grouped:
        group_sizes.append(len(group))
        for pair in group:
            neg_ppi_dict.update({pair: int(group_no)})
        group_no += 1
    
    print(f'[{dt.datetime.now()}] Finished generating negative PPI labels!')
    return(neg_ppi_dict)

# merges overlapping PPIs in different complexes into unique super groups
def merge_groups(lsts):
    sets = [set(lst) for lst in lsts if lst]
    merged = True
    while merged:
        merged = False
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = True
                    common |= x
            results.append(common)
        sets = results
    return(sets)

# get non-redundant protein complex super groups
def make_sprgrp_dict(labeled_fmat):
    mgroup_list = []
    for i in range(len(labeled_fmat)):
        group = labeled_fmat['group'][i]
        if type(group) == list:
            mgroup_list.append(group)
    merged = merge_groups(mgroup_list)
    super_grp_num = 1
    super_grp_dict = dict()
    for group in merged:
        for old_gnum in group:
            super_grp_dict[int(old_gnum)] = int(super_grp_num)
        super_grp_num += 1
    return(super_grp_dict)


"""
Functions for reading in, labeling, & writing out PPI feature matrices for input into ML pipeline:
"""

def label_fmat(fmat_file, pos_dict, neg_dict):
    print(f'[{dt.datetime.now()}] Loading features from {fmat_file}...')
    # TODO: add pickle or CSV option
    with open(fmat_file, 'rb') as handle:
        fmat = pickle.load(handle)

    print(f'[{dt.datetime.now()}] Formatting feature matrix ID columns & rows ...')
    fmat[['ID1','ID2']] = fmat['ID'].str.split(' ',expand=True)
    fmat = fmat[fmat['ID2'].notna()]
    
    t0 = time.time()
    # maybe TODO: long step; potentially optimize ..?
    print(f'[{dt.datetime.now()}] Labeling feature matrix (takes awhile) ...')
    print(fmat.dtypes)
    # get all positive pairs
    pos_pairs = [item for sublist in list(pos_dict.values()) for item in sublist]
    # label pairs
    fmat['label'] = [match_label(frozenset({i, j}), pos_pairs, neg_dict) for i, j in zip(fmat['ID1'], fmat['ID2'])]
    fmat['group'] = [match_group(frozenset({i, j}), pos_dict, neg_dict) for i, j in zip(fmat['ID1'], fmat['ID2'])]
    
    # debug group matching
    print('before:')
    print(fmat.dtypes['group'])
    fmat['group'] = fmat['group'].astype(int)
    print()
    print('after:')
    print(fmat.dtypes['group'])
    
    print(f'[{dt.datetime.now()}] Total time to label {len(fmat)} rows: {time.time() - t0} seconds')
    
    num_pos = len(fmat[(fmat['label'] == 1)])
    num_neg = len(fmat[(fmat['label'] == -1)])
    print(f'[{dt.datetime.now()}] Total # positive PPIs = {num_pos}')
    print(f'[{dt.datetime.now()}] Total # negative PPIs = {num_neg}')
    return(fmat)

def label_fmat_supergrps(labeled_fmat):
    print(f'[{dt.datetime.now()}] Generating merged complex groups ...')
    sdict = make_sprgrp_dict(labeled_fmat)
    print(f'[{dt.datetime.now()}] Labeling non-redundant complex groups ...')
    labeled_fmat['super_group'] = [match_spr_grp(i, sdict) for i in labeled_fmat['group']]
    return(labeled_fmat)

def format_fmat(labeled_fmat, keep_overlap_groups=False, shuffle_feats=False, shuffle_rows=False):
    # drop ID split cols
    labeled_fmat.drop(['ID1', 'ID2'], axis=1, inplace=True)
    # get col names for labels, features
    print(f'[{dt.datetime.now()}] Reformatting columns ...')
    label_cols = ['ID', 'group', 'super_group', 'label']
    feature_cols = [c for c in labeled_fmat.columns.values.tolist() if c not in label_cols]
    # optionally shuffle feature order
    if shuffle_feats:
        print(f'[{dt.datetime.now()}] Shuffling feature columns ...')
        random.shuffle(feature_cols)
    # reorder columns
    fmat_fmt = labeled_fmat[label_cols + feature_cols]
    # optionally drop group_col with redundant PPIs
    # probably always want to drop it tho tbh
    if not keep_overlap_groups:
        print(f'[{dt.datetime.now()}] Dropping redundant complex groups ...')
        fmat_fmt = fmat_fmt.drop(['group'], axis=1)
    # optionally shuffle row order;
    # --> technically will be done later w/ sklearn.model_selection.GroupShuffleSplit
    # --> but it's here if you want to shuffle at this step for some reason
    if shuffle_rows:
        print(f'[{dt.datetime.now()}] Shuffling non-redundant protein complex super groups ...')
        grps = fmat_fmt['super_group'].unique()
        random.shuffle(grps)
        fmat_fmt = fmat_fmt.set_index('super_group').loc[grps].reset_index()
    print('Final feature matrix:')
    print(fmat_fmt.head())
    return(fmat_fmt)

def write_fmat_files(labeled_fmat, fmat_file, outfile=None):
    # format outfile path/name if none specified
    if not outfile:
        plist = fmat_file.split('/', 1)
        outpath = '/'.join(plist[:-1])+'/'
        outfile = outpath+'featmat_labeled'
    
    # write out matrices
    print(f'[{dt.datetime.now()}] Writing out matrices:')
    print(f"[{dt.datetime.now()}] \t► Full matrix (labeled + unlabeled) --> {outfile}")
    print(f"[{dt.datetime.now()}] \t► Positive & negative PPIs --> {outfile+'_traintest'}")
    print(f"[{dt.datetime.now()}] \t► Gold standard (positive) PPIs only --> {outfile+'_goldstd'}")
    
    # gold standard (positive/known) PPIs only
    t1 = time.time()
    print(f'[{dt.datetime.now()}] Extracting gold standard PPIs ...')
    goldstd = labeled_fmat[(labeled_fmat['label'] == 1)]
    goldstd.reset_index(drop=True, inplace=True)
    print(f"[{dt.datetime.now()}] Writing serialized gold standard matrix to {outfile+'_goldstd.pkl'} ... ")
    goldstd.to_pickle(outfile+'_goldstd.pkl')
    print(f"[{dt.datetime.now()}] Writing comma-separated gold standard matrix to {outfile+'_goldstd'} ... ")
    goldstd.to_csv(outfile+'_goldstd', index=False)
    print(f"[{dt.datetime.now()}] Total time to write gold standard feature matrix of shape {goldstd.shape}: {time.time() - t1} seconds")
    
    # positive + negative PPIs only
    t2 = time.time()
    print(f'[{dt.datetime.now()}] Extracting train/test rows ...')
    traintest = labeled_fmat[(labeled_fmat['label'] == 1) | (labeled_fmat['label'] == -1)]
    traintest.reset_index(drop=True, inplace=True)
    print(f"[{dt.datetime.now()}] Writing serialized train/test matrix to {outfile+'_traintest.pkl'} ... ")
    traintest.to_pickle(outfile+'_traintest.pkl')
    print(f"[{dt.datetime.now()}] Writing comma-separated train/test matrix to {outfile+'_traintest'} ... ")
    traintest.to_csv(outfile+'_traintest', index=False)
    print(f"[{dt.datetime.now()}] Total time to write train/test feature matrix of shape {traintest.shape}: {time.time() - t2} seconds")
    
    # all data, labeled & unlabeled
    t3 = time.time()
    print(f"[{dt.datetime.now()}] Writing full serialized matrix to {outfile+'.pkl'} ... ")
    labeled_fmat.to_pickle(outfile+'.pkl')
    print(f"[{dt.datetime.now()}] Writing full comma-separated matrix to {outfile} ... ")
    labeled_fmat.to_csv(outfile, index=False)
    print(f"[{dt.datetime.now()}] Total time to write full feature matrix of shape {labeled_fmat.shape}: {time.time() - t3} seconds")

""" All wrapped together: """
def main():
    # if specified, set seed
    t0 = time.time()
    if args.seed:
        random.seed(args.seed)
    # make dicts for +/- labels
    pos_dict = make_pos_dict(args.gold_std)
    neg_dict = make_neg_dict(pos_dict)
    # label feature matrix
    labeled_fmat = label_fmat(args.featmat, pos_dict, neg_dict)
    labeled_fmat_final = label_fmat_supergrps(labeled_fmat)
    # format feature matrix
    fmat_out = format_fmat(labeled_fmat_final, args.keep_cmplx_overlap, args.shuffle_feats, args.shuffle_rows)
    # write results
    write_fmat_files(fmat_out, args.featmat, args.outfile_name)
    print(f"[{dt.datetime.now()}] ---------------------------------------------------------")
    print(f"[{dt.datetime.now()}] Total run time: {time.time() - t0} seconds")
    print(f"[{dt.datetime.now()}] ---------------------------------------------------------")
    

""" When executed from the command line: """
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify feature matrix
    parser.add_argument("-f", "--featmat", help="(Required) Path to feature matrix you want to label. PPI ID column is a column named 'ID' and protein names are separated by a space.")

    # specify positive labels
    parser.add_argument("-g", "--gold_std", help="(Required) Path to file containing gold standard PPIs; 1 complex per line, subunits are space separated.")
    
    # specify outfile name (default='featmat', written to the given data directory)
    parser.add_argument("-o", "--outfile_name", action="store", default=None, help="(Optional) Specify the outfile path/name. Defaults are 'featmat_labeled' and 'featmat_labeled.pkl', written to the same directory containing the input.")
    
    # specify outfile name (default='featmat', written to the given data directory)
    parser.add_argument("-s", "--seed", action="store", default=None, help="(Optional) Specify seed to make randomized negative PPIs reproducible.")
    
    # specify if you want to keep the redundant group col
    parser.add_argument("--keep_cmplx_overlap", action="store_true", default=False, help="(Optional) Specify if you want to keep protein complex group numbers with redundant PPIs. Default=False (recommended).")
    
    # specify if you want to shuffle feature columns
    parser.add_argument("--shuffle_feats", action="store_true", default=False, help="(Optional) Specify if you want to shuffle feature columns. Default=False.")
    
    # specify if you want to shuffle rows
    parser.add_argument("--shuffle_rows", action="store_true", default=False, help="(Optional) Specify if you want to shuffle rows. Default=False. Note: Can be performed later with sklearn.")

    # specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main()