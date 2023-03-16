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
import os
from datetime import datetime as dt
from tqdm import tqdm
from itertools import combinations
from itertools import accumulate
from collections import defaultdict

"""
Functions for iterative row-wise labeling:
"""

# convert space-separated IDs to iterable sets
def make_fset(x, drop=True):
    if len(set(x.split(' '))) < 2:
        if drop == False:
            print(f"WARNING: Features for '{x}' (self-self PPI) detected ...")
            x1 = x.split(' ')[0]
            fset = frozenset({x1})
            return(fset)
        else:
            return()   
    else:
        x1 = x.split(' ')[0]
        x2 = x.split(' ')[1]
        fset = frozenset({x1,x2})
        return(fset)

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

# get gold standard ppis
def make_gs_dict(gs_file):
    gs_dict = dict()
    group_no = 1
    dupes = []
    with open(gs_file, 'r') as f:
        ppis = f.read().splitlines() 
        for p in ppis:
            ogs = p.split(' ')
            fsets = [frozenset({i, j}) for i,j in list(combinations(ogs, 2))]
            
            # make sure PPI combinations are behaving as expected
            num_prots = len(ogs)
            expected_ppi_number = (num_prots*(num_prots-1))/2
            assert expected_ppi_number == len(fsets), "ERROR: Problem with gold standard complexes; make sure each complex contains unique protein IDs (no repeated subunits)."
            
            gs_dict.update({int(group_no): fsets})
            group_no += 1
    return(gs_dict)

# get all possible negative ppis
def get_neg_ppis(gs_dict):

    # get proteins from gold standard PPIs to generate negative PPIs
    random_prots = set()
    flat_gs_ppis = [pair for pair_list in list(gs_dict.values()) for pair in pair_list]
    print(f' ► # total possible gold standard PPIs = {len(set(flat_gs_ppis))}')
    all_gs_prots = [p for pair in flat_gs_ppis for p in list(pair)]
    uniq_gs_prots = set(all_gs_prots)
    print(f' ► # unique gold standard prots = {len(set(uniq_gs_prots))}')
    
    print(f' ► Getting all protein combinations ...')
    fsets = [frozenset({i, j}) for i,j in list(combinations(list(uniq_gs_prots), 2))]
    print(f' ► # total pairwise PPIs = {len(fsets)}')
    
    print(f' ► Removing known gold standard PPIs from all possible PPI combinations ...')
    neg_ppis = set(fsets).difference(set(flat_gs_ppis))
    
    # make sure there is no overlap between positive & negative PPI pairs
    assert len(set(neg_ppis) & set(flat_gs_ppis)) == 0, "Overlap between positive and negative PPIs detected."
    
    print(f' ► # total possible negative PPIs = {len(neg_ppis)}')
    return(neg_ppis)

# get ppis actually osberved in data
def find_obs_labels(fmat_file, all_neg_ppis, gs_dict):
    
    # TODO: add pickle or CSV option
    with open(fmat_file, 'rb') as handle:
        fmat = pickle.load(handle)
    fmat_ppis = [make_fset(i) for i in fmat['ID']]
    print(f' ► # total PPIs observed in data = {len(fmat_ppis)}')
    
    print(f' ► Finding overlap between observed PPIs and gold standard PPIs ...')
    flat_gs_ppis = [pair for pair_list in list(gs_dict.values()) for pair in pair_list]
    neg_overlap = set(fmat_ppis).intersection(set(all_neg_ppis))
    pos_overlap = set(fmat_ppis).intersection(set(flat_gs_ppis))
    print(f' ► # total negative PPIs observed in data = {len(neg_overlap)}')
    print(f' ► # total positive PPIs observed in data = {len(pos_overlap)}')
    return(neg_overlap, pos_overlap)

# make label dictionaries
def make_label_dicts(obs_neg, obs_pos, gs_dict, num_neg_labels=None):
    
    # specify number of negative labels to include
    if not num_neg_labels:
        num_neg_labels = 3*len(obs_pos)
    
    final_grp_nums = []
    pos_ppi_dict = dict()
    for grp, cmplx in gs_dict.items():
        if any(ppi in obs_pos for ppi in cmplx):
            pos_ppi_dict.update({grp: cmplx})
            for ppi in cmplx:
                final_grp_nums.append(grp)
    
    print(f' ► (# total GS groups in data)/(# total possible GS groups) = {len(pos_ppi_dict)}/{len(gs_dict)}')
    
    print(f' ► Randomly sampling {num_neg_labels} negative PPIs from {len(obs_neg)} total observed negative PPIs ...')
    random.shuffle(list(obs_neg))
    neg_ppis = random.sample(list(obs_neg), num_neg_labels)
    neg_grp_nums = random.choices(final_grp_nums, k=num_neg_labels)
    
    print(f' ► Assigning group numbers ...')
    neg_ppi_dict = dict()
    for i in range(len(neg_grp_nums)):
        neg_ppi_dict.update({neg_ppis[i]: int(neg_grp_nums[i])})
    
    print(f' ► Finished generating positive and negative PPI labels!')
    return(neg_ppi_dict, pos_ppi_dict)

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

def format_outdir(featmat, outfile_name=None):
    # format outfile paths/names
    if outfile_name:
        path, filename = os.path.split(os.path.realpath(outfile_name))
        outpath = path+'/'
        fmat_out = outpath+'/'+outfile_name
    else:
        path, filename = os.path.split(os.path.realpath(featmat))
        outpath = path+'/'
        fmat_out = outpath+'/'+'featmat_labeled'
    return(outpath, fmat_out)

def label_fmat(fmat_file, pos_dict, neg_dict):
    # TODO: add pickle or CSV option, eventually
    with open(fmat_file, 'rb') as handle:
        fmat = pickle.load(handle)

    fmat['frozen_pair'] = [make_fset(i, drop=True) for i in fmat['ID']]
    
    t0 = time.time()
    
    # get all positive pairs
    pos_pairs = [pair for cmplx in list(pos_dict.values()) for pair in cmplx]
    # label pairs
    fmat['label'] = [match_label(i, pos_pairs, neg_dict) for i in tqdm(fmat['frozen_pair'])]
    print(f' ► Labeling complex groups ...')
    fmat['group'] = [match_group(i, pos_dict, neg_dict) for i in tqdm(fmat['frozen_pair'])]
    print(f' ► Total time to label {len(fmat)} rows: {time.time() - t0} seconds')
    
    num_pos = len(fmat[(fmat['label'] == 1)])
    num_neg = len(fmat[(fmat['label'] == -1)])
    print(f' ► Total # positive PPIs = {num_pos}')
    print(f' ► Total # negative PPIs = {num_neg}')
    
    print(f' ► Generating merged complex groups ...')
    sdict = make_sprgrp_dict(fmat)
    print(f' ► Labeling non-redundant complex groups ...')
    fmat['super_group'] = [match_spr_grp(i, sdict) for i in tqdm(fmat['group'])]
    return(fmat)

def format_fmat(labeled_fmat, keep_overlap_groups=False, shuffle_feats=False, shuffle_rows=False):
    # drop ID split cols
    #labeled_fmat.drop(['ID1', 'ID2'], axis=1, inplace=True)
    # get col names for labels, features
    labeled_fmat.drop(['frozen_pair'], axis=1, inplace=True)
    label_cols = ['ID', 'group', 'super_group', 'label']
    feature_cols = [c for c in labeled_fmat.columns.values.tolist() if c not in label_cols]
    
    # optionally shuffle feature order
    if shuffle_feats:
        print(f' ► Shuffling feature columns ...')
        random.shuffle(feature_cols)
    
    # reorder columns
    fmat_fmt = labeled_fmat[label_cols + feature_cols]
    
    # optionally drop group_col with redundant PPIs
    # probably always want to drop it tho tbh
    if not keep_overlap_groups:
        print(f' ► Dropping redundant complex groups ...')
        fmat_fmt = fmat_fmt.drop(['group'], axis=1)
    
    # optionally shuffle row order;
    # --> technically can be done later w/ sklearn.model_selection.GroupShuffleSplit
    # --> but it's here if you want to shuffle at this step for some reason
    if shuffle_rows:
        print(f' ► Shuffling non-redundant protein complex super groups ...')
        grps = fmat_fmt['super_group'].unique()
        random.shuffle(grps)
        fmat_fmt = fmat_fmt.set_index('super_group').loc[grps].reset_index()
    
    print('Final feature matrix:')
    print(fmat_fmt)
    return(fmat_fmt)

def write_fmat_files(labeled_fmat, fmat_file, outfile=None):
    
    # gold standard (positive/known) PPIs only
    t1 = time.time()
    print(f' ► Writing {outfile}_goldstd ...')
    goldstd = labeled_fmat[(labeled_fmat['label'] == 1)]
    goldstd.reset_index(drop=True, inplace=True)
    goldstd.to_pickle(outfile+'_goldstd.pkl')
    goldstd.to_csv(outfile+'_goldstd', index=False)
    print(f" ► Total time to write gold standard feature matrix of shape {goldstd.shape}: {time.time() - t1} seconds")
    
    # positive + negative PPIs only
    t2 = time.time()
    print(f' ► Writing {outfile}_traintest ...')
    traintest = labeled_fmat[(labeled_fmat['label'] == 1) | (labeled_fmat['label'] == -1)]
    traintest.reset_index(drop=True, inplace=True)
    traintest.to_pickle(outfile+'_traintest.pkl')
    traintest.to_csv(outfile+'_traintest', index=False)
    print(f" ► Total time to write train/test feature matrix of shape {traintest.shape}: {time.time() - t2} seconds")
    
    # all data, labeled & unlabeled
    t3 = time.time()
    print(f' ► Writing {outfile} ...')
    labeled_fmat.to_pickle(outfile+'.pkl')
    labeled_fmat.to_csv(outfile, index=False)
    print(f" ► Total time to write full feature matrix of shape {labeled_fmat.shape}: {time.time() - t3} seconds")

""" All wrapped together: """
def main():
    
    t0 = time.time()
    
    # if specified, set seed
    if args.seed:
        random.seed(args.seed)
        
    # format outfile paths/names
    outpath, fmat_outfile = format_outdir(args.featmat, args.outfile_name)
    print(f'[{dt.now()}] Output directory: {outpath}')
    pos_ppi_outfile = outpath+'positive_ppi_dict.pkl'
    neg_ppi_outfile = outpath+'negative_ppi_dict.pkl'
        
    # make dicts for +/- labels
    print(f'[{dt.now()}] Generating grouped positive PPI labels from gold standard complexes ...')
    gs_dict = make_gs_dict(args.gold_std_file)
    print(f'[{dt.now()}] Getting proteins from gold standard PPIs to generate negative PPIs ...')
    all_neg_ppis = get_neg_ppis(gs_dict)
    print(f'[{dt.now()}] Loading PPIs from {args.featmat}...')
    obs_neg, obs_pos = find_obs_labels(args.featmat, all_neg_ppis, gs_dict)
    print(f'[{dt.now()}] Getting complex group labels ...')
    neg_dict, pos_dict = make_label_dicts(obs_neg, obs_pos, gs_dict, num_neg_labels=args.num_negatives)
    
    # write out positive/negative ppis
    with open(pos_ppi_outfile, 'wb') as handle:
        pickle.dump(pos_dict, handle)
    with open(neg_ppi_outfile, 'wb') as handle:
        pickle.dump(neg_dict, handle)
        
    # label feature matrix
    print(f'[{dt.now()}] Labeling feature matrix (takes awhile) ...')
    labeled_fmat = label_fmat(args.featmat, pos_dict, neg_dict)
    
    # format feature matrix
    print(f'[{dt.now()}] Reformatting feature matrix columns ...')
    fmat_out = format_fmat(labeled_fmat, args.keep_cmplx_overlap, args.shuffle_feats, args.shuffle_rows)
        
    # write final feature matrix results
    # write out matrices
    print(f'[{dt.now()}] Writing out matrices:')
    print(f"[{dt.now()}] Full matrix (labeled + unlabeled) --> {outfile}")
    print(f"[{dt.now()}] Positive & negative PPIs --> {outfile+'_traintest'}")
    print(f"[{dt.now()}] Gold standard (positive) PPIs only --> {outfile+'_goldstd'}")
    write_fmat_files(fmat_out, args.featmat, fmat_outfile)
        
    print(f"[{dt.now()}] ---------------------------------------------------------")
    print(f"[{dt.now()}] Total run time: {round((time.time()-t0)/60, 2)} minutes.")
    print(f"[{dt.now()}] ---------------------------------------------------------")
    

""" When executed from the command line: """
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify feature matrix
    parser.add_argument("-f", "--featmat", help="(Required) Path to feature matrix you want to label. PPI ID column is a column named 'ID' and protein names are separated by a space.")

    # specify positive labels
    parser.add_argument("-g", "--gold_std_file", action="store", help="(Required) Path to file containing gold standard PPIs; 1 complex per line, subunits are space separated.")
    
    # specify outfile name (default='featmat', written to the given data directory)
    parser.add_argument("-o", "--outfile_name", action="store", default=None, help="(Optional) Specify the outfile path/name. Defaults are 'featmat_labeled' and 'featmat_labeled.pkl', written to the same directory containing the input.")
    
    # specify outfile name (default='featmat', written to the given data directory)
    parser.add_argument("-n", "--num_negatives", action="store", default=None, help="(Optional) Specify max number of negative labels (default = 3X the number of positive labels).")
    
    # specify outfile name (default='featmat', written to the given data directory)
    parser.add_argument("-s", "--seed", action="store", default=None, help="(Optional) Specify seed to make randomized negative PPIs reproducible.")
    
    # specify if you want to keep the redundant group col
    parser.add_argument("--keep_cmplx_overlap", action="store_true", default=False, help="(Optional) Specify if you want to keep protein complex group numbers with redundant PPIs. Default=False (only recommended to keep for running through sanity_checks.py, which will evaluate merge behavior).")
    
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