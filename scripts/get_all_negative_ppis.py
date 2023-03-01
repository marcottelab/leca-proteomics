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
ct=dt.datetime.now()

def make_fset(x):
    if len(set(x.split(' '))) < 2:
        print(f"WARNING: Correlation metrics for '{x}' (self-self PPI) detected; make sure you mean for this to be included ...")
        x1 = x.split(' ')[0]
        fset = frozenset({x1})
        return(fset)
    else:
        x1 = x.split(' ')[0]
        x2 = x.split(' ')[1]
        fset = frozenset({x1,x2})
        return(fset)

def make_gs_dict(gs_file):
    print(f'[{ct}] Generating grouped positive PPI labels from gold standard complexes ...')
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
    return(pos_ppi_dict)

# make negative labels
# reworking to get negative PPIs from observed PPIs in input featmat
def get_neg_ppis(gs_dict):

    print(f'[{ct}] Getting proteins from gold standard PPIs to generate negative PPIs ...')
    random_prots = set()
    flat_gs_ppis = [pair for pair_list in list(gs_dict.values()) for pair in pair_list]
    print(f'--> # gold standard pairwise PPIs = {len(set(flat_gs_ppis))}')
    all_gs_prots = [p for pair in flat_gs_ppis for p in list(pair)]
    uniq_gs_prots = set(all_gs_prots)
    print(f'--> # unique gold standard prots = {len(set(uniq_gs_prots))}')
    
    print(f'[{ct}] Getting all protein combinations ...')
    fsets = [frozenset({i, j}) for i,j in list(combinations(list(uniq_gs_prots), 2))]
    print(f'--> # total pairwise PPIs = {len(fsets)}')
    
    print(f'[{ct}] Removing known gold standard PPIs from all possible PPI combinations ...')
    neg_ppis = set(fsets).difference(set(flat_gs_ppis))
    
    print(f'--> # total possible negative PPIs = {len(neg_ppis)}')
    return(neg_ppis)

def find_obs_labels(fmat_file, all_neg_ppis, gs_dict):
    
    print(f'[{ct}] Loading features from {fmat_file}...')
    # TODO: add pickle or CSV option
    with open(fmat_file, 'rb') as handle:
        fmat = pickle.load(handle)
    fmat_ppis = [make_fset(i) for i in fmat['ID']]
    print(f'--> # total PPIs observed in data = {len(fmat_ppis)}')
    
    print(f'[{ct}] Finding overlap between observed PPIs and gold standard PPIs ...')
    flat_gs_ppis = [pair for pair_list in list(gs_dict.values()) for pair in pair_list]
    neg_overlap = set(fmat_ppis).intersection(set(all_neg_ppis))
    pos_overlap = set(fmat_ppis).intersection(set(flat_gs_ppis))
    print(f'--> # total negative PPIs observed in data = {len(neg_overlap)}')
    print(f'--> # total positive PPIs observed in data = {len(pos_overlap)}')
    return(neg_overlap, pos_overlap)

def make_label_dicts(obs_neg, obs_pos, gs_dict, num_neg_labels=None):
    
    # specify number of negative labels to include
    if not num_neg_labels:
        num_neg_labels = 3*len(obs_pos)
    
    print(f'[{ct}] Getting complex group labels ...')
    final_grp_nums = []
    pos_ppi_dict = dict()
    for grp, cmplx in gs_dict.items():
        if any(ppi in obs_pos for ppi in cmplx):
            pos_ppi_dict.update({grp: cmplx})
            for ppi in cmplx:
                final_grp_nums.append(grp)
    print(f'--> (# total GS groups in data)/(# total possible GS groups) = {len(pos_ppi_dict)}/{len(gs_dict)}')
    
    print(f'[{ct}] Randomly sampling {num_neg_labels} negative PPIs from {len(obs_neg)} total observed negative PPIs ...')
    random.shuffle(list(obs_neg))
    neg_ppis = random.sample(list(obs_neg), num_neg_labels)
    
    print(f'[{ct}] Assigning group numbers ...')
    neg_ppi_dict = dict()
    for i in range(len(final_grp_nums)):
        neg_ppi_dict.update({neg_ppis[i]: int(final_grp_nums[i])})
    
    print(f'[{ct}] Finished generating positive and negative PPI labels!')
    return(neg_ppi_dict, pos_ppi_dict)
    
gs_file = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/gold_stds/all.gold.cmplx.noRibos.merged.txt'
fmat_file = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/featmats/featmat.pkl'
outfile = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/gold_stds/all.gold.cmplx.noRibos.merged_allnegatives.txt'

t0 = time.time()
gs_dict = make_gs_dict(gs_file)
all_neg_ppis = get_neg_ppis(gs_dict)
obs_neg, obs_pos = find_obs_labels(fmat_file, all_neg_ppis, gs_dict)
neg_dict, pos_dict = make_label_dicts(obs_neg, obs_pos, gs_dict)
print(f'[{ct}] Total run time: {time.time() - t0} seconds')