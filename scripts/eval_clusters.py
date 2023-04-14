"""
Script for evaluating PPI clustering results based on a gold standard protein complex file.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
from datetime import datetime as dt
import pandas as pd
from functools import reduce

""" Functions """

def get_obs_gold(gold_std_file, results_df):
    with open(gold_std_file, 'r') as f:
        ppis = f.read().splitlines()
        prot_list = []
        cmplx_list = []
        for i, p in enumerate(ppis):
            ogs = p.split(' ')
            for prot in ogs:
                prot_list.append(prot)
                cmplx_list.append(i+1)

    df_gold = pd.DataFrame()
    df_gold['ID'] = prot_list
    df_gold['cmplx'] = cmplx_list
    df_gold = df_gold[df_gold['ID'].isin(results_df['ID'])]
    return(df_gold)

def calc_max_proportion(df, level, expected_gold_sizes):
    cmplx_counts = df.groupby([level, 'cmplx']).size().reset_index(name='counts')
    prop_list = []
    cmplx_list = []
    for i in range(len(cmplx_counts)):
        cmplx_num = cmplx_counts['cmplx'][i]
        obs = cmplx_counts['counts'][i]
        actual = expected_gold_sizes[cmplx_num]
        # calculate proportion of obs gold std complex in same cluster
        prop = round(obs/actual, 2)
        # record results
        cmplx_list.append(cmplx_num)
        prop_list.append(prop)
    cut_res = pd.DataFrame()
    cut_res['cmplx'] = cmplx_list
    cut_res[level] = prop_list
    max_obs = cut_res.loc[cut_res.groupby(['cmplx'])[level].idxmax()]
    return(max_obs)

""" Main """
def main():
    
    print(f'[{dt.now()}] Reading in cluster file ...')
    res = pd.read_csv(args.cluster_file)
    cuts = [c for c in res.columns.values.tolist() if 'cut' in c]
    df_clst = res[['ID'] + cuts]
    
    print(f'[{dt.now()}] Reading in gold standard file ...')
    df_gold = get_obs_gold(args.gold_std_file, res)
    df_clst_gold = df_clst.merge(df_gold, how='inner', left_on=['ID'], right_on=['ID']).sort_values('cmplx')
    
    print(f'[{dt.now()}] Calculating expected gold standard complex sizes ...')
    expected_cmplx_sizes = df_gold.groupby('cmplx').size()
    
    print(f'[{dt.now()}] Calculating max % of each gold standard complex clustered together for each cut ...')
    df_list = []
    for c in cuts:
        print(f'[{dt.now()}] Evaluating {c} ...')
        cut_df = df_clst_gold[['ID', c, 'cmplx']]
        df = calc_max_proportion(cut_df, c, expected_cmplx_sizes)
        df_list.append(df)
        
    for df in df_list:
        df.set_index(['cmplx'], inplace=True)
    final_df = reduce(lambda x, y: x.join(y, how='outer'), df_list)
    final_df = final_df.fillna(0)
    
    print(f'[{dt.now()}] Reordering rows from strongest to weakest complex ...')
    final_df['sum'] = final_df[cuts].sum(axis=1)
    final_df = final_df.sort_values('sum', ascending=False)
    final_df.drop(labels='sum', inplace=True, axis=1)
    #x = final_df.to_numpy()
    #clst = linkage(x, metric = "braycurtis", method = "average")
    #clst_order = list(leaves_list(clst))
    #final_df.reset_index(inplace=True)
    #ord_df = final_df.reindex(clst_order, copy=False)
    print(f'[{dt.now()}] Writing results to {args.outfile} ...')
    print(final_df)
    final_df.to_csv(args.outfile, index=True)
    print(f'[{dt.now()}] Done!')

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # specify infile
    parser.add_argument("-c", "--cluster_file", help="(Required) Path to file containing protein clusters output from a community detection algorithm (e.g., walktrap). Must contain one row per protein in a column named 'ID' and cluster assignments must be given in columns that contain the word 'cut', e.g., a column named 'cut_222' might correspond to a dendogram cut into 222 clusters and each row contains a cluster assignment for each protein.")

    # specify gold standard file
    parser.add_argument("-g", "--gold_std_file", action="store", help="(Required) Path to file containing gold standard PPIs; 1 complex per line, subunits are space separated.")

    # specify outfile name
    parser.add_argument("-o", "--outfile", action="store", help="(Required) Specify the outfile path/name.")

    # specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main()