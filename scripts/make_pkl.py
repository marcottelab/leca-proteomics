"""
Script for pickling CSV or TSV files.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pickle
import pandas as pd
import time

def main():
    t0 = time.time()
    print(f'Reading in {args.infile} ...')
    df = pd.read_csv(args.infile)
    if not args.outfile_name:
        outfile_name = args.infile+'.pkl'
    print(f'Writing serialized output to {outfile_name} ...')
    df.to_pickle(outfile_name)
    print(f'Done!')
    print(f'Total time to pickle file = {time.time() - t0} seconds.')

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # specify infile (required)
    parser.add_argument("-i", "--infile", help="Required: path to file you want to pickle.")
    
    # specify outfile name (default is infile name + 'pkl' suffix)
    parser.add_argument("-o", "--outfile_name", action="store", default=None, help="Optional; specify the outfile path/name (defaults to the infile name/directory with a '.pkl' suffix).")
    
    # specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))
    
    args = parser.parse_args()
    main()