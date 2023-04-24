"""
Script for calculating precision-recall for each cut of a dendogram representing clustered protein-protein interactions.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pickle
import pandas as pd
from itertools import combinations
import time
from datetime import datetime as dt

""" Functions """

def get_ppis(df, col, all_pos, all_neg):
    df = df[['ID', col]]
    grouped = df.groupby(col)
    pos_set = set()
    neg_set = set()
    for cluster, data in grouped:
        prot_lst = [p for p in data['ID']]
        ppis = [frozenset({i, j}) for i,j in list(combinations(prot_lst, 2))]
        pos = set(ppis).intersection(set(all_pos))
        if len(pos) > 0:
            pos_set.update(set(pos))
        neg = set(ppis).intersection(set(all_neg))
        if len(neg) > 0:
            neg_set.update(set(neg))
    return(pos_set, neg_set)

def calc_cut_pr(pos, neg, total_pos):
    tp = len(pos)
    fp = len(neg)
    fn = total_pos - tp
    precision = tp/(tp+fp)
    recall = tp/(tp+fn)
    return(precision, recall)


""" Main """
def main():
    
    t0 = time.time()
    
    print(f'[{dt.now()}] Reading in data ...')
    with open(args.positive_ppi_dict, 'rb') as f:
        pos_dict = pickle.load(f)
    with open(args.negative_ppi_dict, 'rb') as f:
        neg_dict = pickle.load(f)    
    neg_ppis = list(neg_dict.keys())
    pos_ppis = [item for sublist in list(pos_dict.values()) for item in sublist]

    df = pd.read_csv(args.cluster_file)
    ogs = df['ID'].tolist()
    cuts = [c for c in df.columns.values.tolist() if 'cut' in c]
    all_pairs = [frozenset({i, j}) for i,j in list(combinations(ogs, 2))]
    print(f'[{dt.now()}] Total number of pairs in results: {len(all_pairs)}')

    print(f'[{dt.now()}] Getting all observed positive PPIs in results ...')
    all_obs_pos = set(all_pairs).intersection(set(pos_ppis))
    total_obs_pos = len(all_obs_pos)
    print(f'[{dt.now()}] Total possible positive PPIs:', total_obs_pos)

    print(f'[{dt.now()}] Getting all observed negative PPIs in results ...')
    all_obs_neg = set(all_pairs).intersection(set(neg_ppis))
    total_obs_neg = len(all_obs_neg)
    print(f'[{dt.now()}] Total possible negative PPIs:', total_obs_neg)
    
    print(f'[{dt.now()}] Calculating precision, recall for {len(cuts)} cuts ...')
    res_dict = dict()
    for c in cuts:
        n_clst = c.replace('cut_','') 
        pos, neg = get_ppis(df, c, all_obs_pos, all_obs_neg)
        precision, recall = calc_cut_pr(pos, neg, total_obs_pos)
        res_dict[n_clst] = [float(precision), float(recall)]
        print(f'[{dt.now()}] {n_clst} clusters: P={precision}, R={recall}')
    
    res_df = pd.DataFrame.from_dict(res_dict, orient='index', columns=['precision','recall'])
    res_df = res_df.reset_index().rename(columns={'index': f'n_clusters(n_prots={len(ogs)})'})
    print(f'[{dt.now()}] Writing results to {args.outfile} ...')
    res_df.to_csv(args.outfile, index=False)


    print(f"[{dt.now()}] ---------------------------------------------------------")
    print(f"[{dt.now()}] Total run time: {round((time.time()-t0)/60, 2)} minutes.")
    print(f"[{dt.now()}] ---------------------------------------------------------")
    
if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # specify infile
    parser.add_argument("-c", "--cluster_file", help="(Required) Path to file containing protein clusters output from a community detection algorithm (e.g., walktrap). Must contain one row per protein in a column named 'ID' and cluster assignments must be given in columns that contain the word 'cut', e.g., a column named 'cut_222' might correspond to a dendogram cut into 222 clusters and each row contains a cluster assignment for each protein.")
    
    # specify positive PPI dict file
    parser.add_argument("-p", "--positive_ppi_dict", action="store", help="(Required) Specify path to pickled dictionary containing all possible positive PPI pairs.")
    
    # specify negative PPI dict file
    parser.add_argument("-n", "--negative_ppi_dict", action="store", help="(Required) Specify path to pickled dictionary containing all possible negative PPI pairs.")
    
    # specify outfile
    parser.add_argument("-o", "--outfile", action="store", help="(Required) Specify the outfile path/name.")

    # specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main()