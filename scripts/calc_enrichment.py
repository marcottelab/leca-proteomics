"""
Script for calculating external PPI enrichment relative to CFMS PPIs binned by ML score.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import time
import os
from datetime import datetime as dt
from tqdm import tqdm

import pandas as pd
import numpy as np
import math
from itertools import combinations

""" Functions """

def fmt_external(ext_file, int_filter_file=None):
    
    ext_df = pd.read_csv(ext_file)
    if int_filter_file:
        int_filter = set([line.strip() for line in open(int_filter_file, 'r')])
        ext_fs = [frozenset({i, j}) for i, j in zip(ext_df['ID1'], ext_df['ID2']) if i in int_filter and j in int_filter]
    else:
        ext_fs = [frozenset({i, j}) for i, j in zip(ext_df['ID1'], ext_df['ID2'])]
    
    return(ext_fs)

def map_internal(int_file, ext_fs):
    
    int_df = pd.read_csv(int_file)
    int_df[['ID1', 'ID2']] = int_df['ID'].str.split(' ', n=1, expand=True)
    int_df['fs'] = [frozenset({i, j}) for i, j in zip(int_df['ID1'], int_df['ID2'])]
    
    int_fs = [frozenset({i, j}) for i, j in zip(int_df['ID1'], int_df['ID2'])]
    intersection = set(int_fs) & set(ext_fs)
    int_df['ext'] = [1 if i in intersection else 0 for i in int_df['fs']]
    
    return(int_df)

def fmt_outdir(int_file, outfile_name=None):
    
    if outfile_name:
        path, filename = os.path.split(os.path.realpath(outfile_name))
        outname = path+'/'+filename
    else:
        path, filename = os.path.split(os.path.realpath(int_file))
        outname = path+'/'+'ppi_enrichment.csv'
        
    return(outname)

def calc_binned_enrichment(df, ext_set):
    
    all_leca_ppis = set([frozenset({i, j}) for i,j in zip(df['ID1'], df['ID2'])])
    all_prots = set([p for pair in all_leca_ppis for p in list(pair)])
    all_ppis = set([frozenset({i, j}) for i,j in list(combinations(list(all_prots), 2))])
    all_assayable_ppis = all_ppis.difference(all_leca_ppis)
    
    leca_only = len(all_leca_ppis.difference(ext_set))
    ext_only = len(all_assayable_ppis.intersection(ext_set))
    neither = len(all_assayable_ppis.difference(ext_set))
    both = len(df[df['ext'] == 1])
    
    odds_ratio = (both/ext_only)/(leca_only/neither)
    prob = (both/(both+ext_only))/(leca_only/(leca_only+neither))
    
    try:
        
        log_odds = math.log10(odds_ratio)
        log_prob = math.log10(prob)
        
    except ValueError: 
        
        # 0.5 is added to avoid indeterminate ratios (Haldane-Anscombe correction)
        leca_only = leca_only + 0.5
        ext_only = ext_only + 0.5
        neither = neither + 0.5
        both = both + 0.5
        
        odds_ratio = (both/ext_only)/(leca_only/neither)
        prob = (both/(both+ext_only))/(leca_only/(leca_only+neither))
        log_odds = math.log10(odds_ratio)
        log_prob = math.log10(prob)
                                   
    return(odds_ratio, log_odds, prob, log_prob)

""" Main """

def main():
    
    t0 = time.time()
    
    print(f'[{dt.now()}] Getting external PPIs ...')
    ext_fs = fmt_external(args.ext_file, args.int_filter_file)
    
    print(f'[{dt.now()}] Mapping external PPIs onto internal PPIs ...')
    int_df = map_internal(args.int_file, ext_fs)
    
    print(f'[{dt.now()}] Binning internal PPIs ...')
    int_df['bin'] = np.divmod(np.arange(len(int_df)), args.bin_width)[0]+1
    int_binned = int_df.groupby('bin')
    avg = int_binned.agg({'ppi_score':'mean'}).rename(columns={'ppi_score':'avg_ppi_score'})
    
    print(f'[{dt.now()}] Calculating external PPI enrichment per bin ...')
    res_dict = dict()
    cols = ['odds_ratio', 'log_odds', 'prob', 'log_prob']
    for bin_num, group in tqdm(int_binned):
        res_dict[bin_num] = list(calc_binned_enrichment(group, ext_fs))
    res_df = pd.DataFrame.from_dict(res_dict, orient='index', columns=cols)
    res_df.index.name = 'bin'
    out_df = avg.join(res_df)
    print(out_df)
    
    outname = fmt_outdir(args.int_file, args.outfile_name)
    print(f'[{dt.now()}] Writing results to {outname} ...')
    out_df.to_csv(outname)
    
    print(f"[{dt.now()}] ---------------------------------------------------------")
    print(f"[{dt.now()}] Total run time: {round((time.time()-t0)/60, 2)} minutes.")
    print(f"[{dt.now()}] ---------------------------------------------------------")

""" When executed from the command line: """
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify internal PPI file
    parser.add_argument("-i", "--int_file", help="(Required) Path to scored interactions file with a column called 'ID' containing space-separated PPIs and a column called 'ppi_score' containing PPI scores in descending order (i.e., as output by ppi_predict.py).")

    # specify external PPI file
    parser.add_argument("-e", "--ext_file", action="store", help="(Required) Path to external PPI file with columns called 'ID1' and 'ID2' for each external positive PPI.")
    
    # specify ID filter file (optional)
    parser.add_argument("-f", "--int_filter_file", action="store", default=None, help="(Optional) Specify a file containing a list of IDs to filter the external PPIs on, with one ID per line.")
    
    # specify outfile name (optional)
    parser.add_argument("-o", "--outfile_name", action="store", default=None, help="(Optional) Specify the outfile path/name. Default is 'ppi_enrichment.csv' in the same directory as the input internal PPI file.")
    
    # specify bin width (optional)
    parser.add_argument("-bw", "--bin_width", action="store", default=1000, type=int, help="(Optional) Specify number the number of PPIs per bin (bin width) to group internal PPIs by; default = 1000.")
    
    args = parser.parse_args()
    main()