"""
Script for hierarchically clustering elution data given various methods and metrics.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pickle
import time
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

def read_input(infile):
    
    if infile.endswith('pkl'):
        print(f"Reading in pickled file {infile} ...")
        df = pd.read_pickle(infile)
        df.set_index(df.iloc[:, 0].name, inplace=True)
        df.head()
    else:
        print(f"Reading in CSV file {infile} ...")
        df = pd.read_csv(infile)
        df.set_index(df.iloc[:, 0].name, inplace=True)
        df.head()
    
    return df

def fmt_outfile(infile, outdir, metric, method):
    
    if '/' in infile:
        path = infile.rsplit('/', 1)[0]
        filename = infile.rsplit('/', 1)[1]
    else:
        filename = infile
    
    if filename.endswith('pkl'):
        base = filename.rsplit('.', 2)[0]
        ext = filename.rsplit('.', 2)[1]
    else:
        base = filename.rsplit('.', 1)[0]
        ext = filename.rsplit('.', 1)[1]
    
    if not outdir.endswith('/'):
        outfile = outdir+'/'+base+'.clst.'+metric+'.'+method+'.elut'
    else:
        outfile = outdir+base+'.clst.'+metric+'.'+method+'.elut'
    
    return outfile

def cluster_data(df, metric, method):
    
    x = df.to_numpy()
    np.unique(x)
    t0 = time.time()
    
    print(f"Clustering matrix of shape {x.shape} using a '{metric}' distance metric and '{method}' clustering algorithm ...")
    clst = linkage(x, metric = metric, method = method)
    print(f'Time to cluster {x.shape} matrix = {time.time() - t0} seconds')
    
    return clst
                   
def reorder_input(df, clst, outfile, pickle=False):
    
    print("Re-ordering input elution matrix ...")
    clst_order = list(leaves_list(clst))
    df.reset_index(inplace=True)
    clst_df = df.reindex(clst_order, copy=False)
    clst_df.head(10)
    
    t0 = time.time()
    print(f"Writing clustered elution matrix to {outfile} ...")
    clst_df.to_csv(outfile, index=False)
    
    if pickle:
        clst_df.to_pickle(outfile+'.pkl')
    print(f"Time to write file(s) = {time.time() - t0}")
    
    return clst_df
                   
def plot_heatmap(clst_df, p_outfile):
    
    import seaborn as sns
    import matplotlib.pyplot as plt
    cscale = sns.cm.rocket_r
    
    t0 = time.time()
    clst_df.set_index(clst_df.iloc[:, 0].name, inplace=True)
    x = clst_df.to_numpy()
    print(f"Generating heat map for matrix of shape {x.shape}...")
    p = sns.heatmap(x, xticklabels=False, yticklabels=False, cmap = cscale)
    p.tick_params(left=False, bottom=False)
    print(f"Saving heat map to {p_outfile} ...")
    plt.savefig(p_outfile, dpi=300)
    print(f"Total time to make & save heat map plot = {time.time() - t0}")

def main():
    
    """ Main entry point """
    
    t0_wp = time.time()
    
    metric = args.distance_metric
    method = args.cluster_method
    df = read_input(args.infile)
    clst = cluster_data(df, metric, method)
    outfile = fmt_outfile(args.infile, args.output_dir, metric, method)
    clstered_df = reorder_input(df, clst, outfile, args.pickle_output)
    
    if args.plot:
        p_outfile = outfile+'.png'
        plot_heatmap(clstered_df, p_outfile)
        
    print(f"Total time to completion: {time.time() - t0_wp} seconds")

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # path to input file
    parser.add_argument("-i", "--infile", action="store", help="Elution matrix (can be given as a CSV or pickle file; if pickled, make sure extension for the input file is '.pkl').")

    # distance metric
    parser.add_argument("-d", "--distance_metric", action="store", default='euclidean', help = "Default = euclidean. Other distance metrics you can specify: braycurtis, canberra, chebyshev, cityblock, correlation, cosine, dice, euclidean, hamming, jaccard, jensenshannon, kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto, russellrao, seuclidean, sokalmichener, sokalsneath, sqeuclidean, yule. Descriptions: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html)")

    # clustering algorithm
    parser.add_argument("-m", "--cluster_method", action="store", default='single', help="Default = single. Other clustering algorithms you can specify: single, complete, average, weighted, centroid, median, ward. Descriptions: https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.hierarchy.linkage.html")
                   
    # output directory
    parser.add_argument("-o", "--output_dir", action="store", default='single', help="Path to directory where output should be written; outfile name will be the same as the input file name with '.clst' + metric + method appended. It behaves like this because I'm lazy & it fits my use case.")
    
    # generate clustered output as .pkl
    parser.add_argument("-c", "--pickle_output", action="store_true", default=False, help="Optionally specify this argument to generate clustered matrix as compressed .pkl file (in addition to a CSV file).")               
                   
    # plot heat map (optional)
    parser.add_argument("-p", "--plot", action="store_true", default=False, help="Optionally specify this argument to output a clustered heat map as a .png. File name will will be the same as the input file name with '.clst' + metric + method + '.png' appended. It behaves like this because I'm lazy & it fits my use case.")

    # optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbosity (-v, -vv, etc)")

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main()