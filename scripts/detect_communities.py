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
import datetime as dt

""" Functions """
def make_graph(score_file):
    # read in data
    scores = pd.read_csv(score_file)
    # format graph dataframe
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
    
    ct = dt.datetime.now()
    t0 = time.time()
    
    # set seed if specified
    if args.seed:
        random.seed(args.seed)

    # format data into igraph object
    print(f'[{ct}] Loading scores into graph ...')
    ppi_graph = make_graph(args.scores)
    total_prots = ppi_graph.vcount()

    # get dendrogram w/ optimal number of clusters
    # optimal number maximizes modularity of the network
    print(f'[{ct}] Initializing walktrap with default parameters to find optimal number of clusters ...')
    df_opt = walktrap(ppi_graph)
    n_opt = len(df_opt.iloc[:,1].drop_duplicates())
    print(f'[{ct}] Optimal # of clusters that maximize network modularity: ', n_opt)

    # get range of cuts 
    print(f'[{ct}] Computing dendogram cuts for more exclusive clusters ...')
    cuts = np.linspace(n_opt, total_prots, 8, endpoint=False)
    cuts = np.floor(cuts)
    cuts = np.delete(cuts, 0)
    
    # cut walktrap dendrogram for each cut
    df_list = []
    print(f'[{ct}] Executing random walks with {args.steps} steps per vertice ... ')
    for i in cuts:
        print(f'[{ct}] Cutting dendrogram into {int(i)} clusters ...')
        clst = walktrap(ppi_graph, n_steps=args.steps, n_clusters=int(i))
        df_list.append(clst)

    # merge all cuts
    print(f'[{ct}] Merging all walktrap cuts ...')
    for df in df_list:
        df_opt = df_opt.merge(df, how='left', on='ID')

    # sort clusters
    sort_cols = df_opt.columns.values[1:].tolist()
    df_out = df_opt.sort_values(sort_cols)

    # join annotations if specified
    if args.annotations:
        print(f'[{ct}] Joining annotations ...')
        annot_df = pd.read_csv(args.annotations)
        df_out = df_out.merge(annot_df, how='left', on=['ID'])

    # write results
    outname = args.outfile.split('.csv', 1)[0]
    print(f'[{ct}] Writing results to {outname}...')
    df_out.to_csv(outname+'.csv', index=False)
    df_out.to_excel(outname+'.xlsx', index=False)
    rt = round((time.time()-t0)/60, 2)
    print(f'[{ct}] Done! Total run time = {rt} minutes.')
    
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