"""
Script for extracting unique peptide assignments from CFMS data.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pandas as pd
import re
import os
from functools import reduce
import numpy as np

def fmt_data(file):
    
    # read in .prot_count file
    print(f'Processing {file} ...')
    df = pd.read_csv(file, sep='\t', header=None)

    # extract & format peptides 
    pep_list = []
    for i in df.iloc[:, 0]:
        pep = str.split(i, '.')[-1]
        pep_list.append(pep)

    # extract & format protein IDs
    id_list = []
    for i in df.iloc[:, 1]:
        ids = str.split(i, '(')[0]
        id_list.append(ids)
        
    # make new df
    df_fmt = pd.DataFrame()
    df_fmt['peptide'] = pep_list
    df_fmt['protein'] = id_list
    
    return(df_fmt)

def unique_peps(df):
    
    # get protein IDs with unique peptides
    df['num_matches'] = df['protein'].apply(lambda x: len(str.split(x, ',')))
    df = df[df['num_matches'] == 1]
    df_uniq = df.drop(['num_matches'], axis=1)
    
    return(df_uniq)

def count_peps(df, frac):
    
    # get peptide counts
    counts = df.groupby(['peptide']).size().sort_values(ascending=False)
    count_dict = dict()
    for i in counts.items():
        pep = i[0]
        count = i[1]
        count_dict[pep] = count
        
    # join count info back onto df
    count_col = 'frac_count'+str(frac)
    df[count_col] = [count_dict[i] for i in df['peptide']]
    df_counts = df.drop_duplicates()
    
    return(df_counts)

def process_fracs(data_dir, file_list):
    
    df_list = []
    frac_count = 0
    for f in file_list:
        frac_count += 1
        df_fmt = fmt_data(data_dir+f)
        df_uniq = unique_peps(df_fmt)
        df_counts = count_peps(df_uniq, frac_count)
        df_list.append(df_counts)

    if len(df_list) > 1:
        df_joined = reduce(lambda x, y: pd.merge(x, y, on=['peptide', 'protein'], how='outer'), df_list)
        df_joined.fillna(0, inplace=True)
    else:
        df_joined = df_list[0]

    count_cols = []
    for c in df_joined.columns:
        if df_joined[c].dtype == float:
            df_joined[c] = df_joined[c].astype(int)
            count_cols.append(c)

    tcol = 'total_counts'
    df_joined[tcol] = df_joined[count_cols].sum(axis=1)
    
    return(df_joined)

def clean_peps(df):

    # get experiment-wide unique peptide assignments
    print('Getting unique peptide assignments ...')
    counts = df.groupby(['peptide']).size().sort_values(ascending=False)
    bad_peps = []
    for p in counts.items():
        if p[1] > 1:
            bad_peps.append(p[0])
    df_out = df[df['peptide'].apply(lambda x: x not in bad_peps)]
    
    # check for errors
    bad_pep_sum = 0
    for c in counts.items():
        if c[1] > 1:
            bad_pep_sum += c[1]

    actual_sum = len(df) - len(df_out)
    
    # check for error
    #assert(actual_sum == bad_pep_sum, 'Something went wrong ...')
    if actual_sum != bad_pep_sum:
        print("Something went wrong ...")
    
    return(df_out)

def main(args):
    
    data_dir = args.data_dir
    pep_files = [f for f in os.listdir(data_dir) if re.match('.*prot_list$', f)]
    df_joined = process_fracs(data_dir, pep_files)
    df_clean = clean_peps(df_joined)
    df_write = df_clean[['peptide', 'protein', 'total_counts']]
    print(f'Writing peptide assignments & total counts to {args.outfile} ...')
    df_write.to_csv(args.outfile, index=False)
    print('Done!')

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # path to directory containing data
    parser.add_argument("-d", "--data_dir", action="store", dest="data_dir", help="Path to directorying containing CFMS files ending in '.prot_list'")

    # outfile name
    parser.add_argument("-o", "--outfile", action="store", dest="outfile", help="Name for outfile")

    # optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbosity (-v, -vv, etc)")

    # specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)