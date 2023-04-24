"""
Script for detecting protein complex communities given a graph of protein-protein edges weighed by their interaction score.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pandas as pd
import igraph as ig
import random
import numpy as np
import time
from functools import reduce
from datetime import datetime as dt

""" Functions """
def make_graph(score_file):
    # read in data
    scores = pd.read_csv(score_file)
    # format graph data frame
    gdf = pd.DataFrame()
    gdf[['ID1','ID2']] = scores['ID'].str.split(' ', expand=True)
    gdf['weight'] = scores['ppi_score']
    graph = ig.Graph.TupleList(gdf.itertuples(index=False), directed=False, weights=True)
    return(graph)

def walktrap(graph, n_steps=4, n_clusters=None):
    # run walktrap & get clusters
    clusters = graph.community_walktrap(weights='weight', steps=n_steps).as_clustering(n_clusters)
    # write clusters & IDs to dict
    clst_dict = dict()
    for cluster, id_list in enumerate(clusters):
        clst_dict.update({cluster: id_list})
    # get node ids & number of complexes
    nodes = graph.get_vertex_dataframe()
    n_cmplx = len(clusters)
    # format & return data frame
    clst_df = (pd.DataFrame.from_dict(clst_dict, orient='index').T
           .melt(var_name='id', value_name='value')
           .dropna(subset=['value'])
           .astype(int)
           .rename(columns={'id':f'cut_{n_cmplx}', 'value':'ID'}))
    clst_df['ID'].replace(nodes['name'], inplace=True)
    out_df = clst_df[['ID', f'cut_{n_cmplx}']].reset_index(drop=True)
    return(out_df)

""" Main """
def main():

    t0 = time.time()
    # set seed if specified
    if args.seed:
        random.seed(args.seed)

    # format data into igraph object
    print(f'[{dt.now()}] Loading PPI scores into graph ...')
    ppi_graph = make_graph(args.scores)
    total_prots = ppi_graph.vcount()

    # get dendrogram w/ optimal number of clusters
    # optimal number maximizes modularity of the network
    print(f'[{dt.now()}] Initializing walktrap with default parameters to find optimal number of clusters ...')
    df_opt = walktrap(ppi_graph)
    n_opt = len(df_opt.iloc[:,1].drop_duplicates())
    print(f'[{dt.now()}] Optimal # of clusters that maximize network modularity:', n_opt)
    # get range of cuts 
    print(f'[{dt.now()}] Computing dendogram cuts for more inclusive clusters ...')
    front_cuts = np.linspace(n_opt/2, n_opt*2, 6, endpoint=True)
    front_cuts = np.insert(front_cuts, 1, n_opt)
    print(f'[{dt.now()}] Computing dendogram cuts for more exclusive clusters ...')
    back_cuts = np.linspace(n_opt*2, 0.75*total_prots, 4, endpoint=False)
    back_cuts = np.delete(back_cuts, 0)
    
    # format cuts
    cuts = np.concatenate((front_cuts, back_cuts), axis=None)
    cuts = np.floor(cuts)
    cuts = np.unique(cuts)
    print("Final # of clusters per cut:", cuts)
    
    # cut walktrap dendrogram for each cut
    df_list = []
    print(f'[{dt.now()}] Executing random walks with {args.steps} steps per vertice ... ')
    for i in cuts:
        print(f'[{dt.now()}] Cutting dendrogram into {int(i)} clusters ...')
        clst = walktrap(ppi_graph, n_steps=args.steps, n_clusters=int(i))
        clst.set_index(['ID'], inplace=True)
        df_list.append(clst)

    # merge all cuts
    print(f'[{dt.now()}] Merging all walktrap cuts ...')
    df = reduce(lambda x, y: x.join(y, how='outer'), df_list)

    # sort clusters
    sort_cols = df.columns.values[0:].tolist()
    print("Sort columns:", sort_cols)
    df_out = df.sort_values(sort_cols)

    # join annotations if specified
    if args.annotations:
        print(f'[{dt.now()}] Joining annotations ...')
        annot_df = pd.read_csv(args.annotations)
        df_out = df_out.merge(annot_df, how='left', on=['ID'])

    # write results
    outname = args.outfile.split('.csv', 1)[0]
    print(f'[{dt.now()}] Writing results to {outname}...')
    df_out.to_csv(outname+'.csv', index=False)
    df_out.to_excel(outname+'.xlsx', index=False)
    
    print(f"[{dt.now()}] ---------------------------------------------------------")
    print(f"[{dt.now()}] Total run time: {round((time.time()-t0)/60, 2)} minutes.")
    print(f"[{dt.now()}] ---------------------------------------------------------")
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify PPI score file
    parser.add_argument("-i", "--scores", help="(Required) Path to scored comma-separated PPI file with two columns called 'ID' and 'ppi_score'. Column 'ID' contains space separated protein IDs.")
    
    # specify outfile name
    parser.add_argument("-o", "--outfile", action="store", help="(Required) Path/name for outfile.")
    
    # specify number of walktrap steps
    parser.add_argument("--steps", action="store", default=4, type=int, help="(Optional) Specify the number of random walks the algorithm should take. Recommended values are between 3 and 6 (default=4).")
    
    # specify seed to make results consistent
    parser.add_argument("--seed", action="store", type=int, default=None, help="(Optional) Specify seed to make community clusters reproducible.")
    
    # specify annotation file
    parser.add_argument("--annotations", action="store", help="(Optional) Specify the number of random walks the walktrap clustering algorithm should take. Recommended values are between 3 and 6 (default=4).")
    
    
    args = parser.parse_args()
    main()