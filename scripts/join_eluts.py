"""
Sort alphabetically & concatenate all .elut files in a given directory with a given extension.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

import argparse
import os
import pandas as pd
import time
from tqdm import tqdm
from functools import reduce

def get_files(input_dir, ext):
    
    print("Finding files ...")
    files = [fn for fn in os.listdir(input_dir) if fn.endswith(ext)]
    files.sort()
    
    return files

def read_csv_pgbar(csv_path, chunksize, dtype=object):
    
    rows = sum(1 for _ in open(csv_path, 'r')) - 1 # minus the header
    chunk_list = []
    with tqdm(total=rows, desc='Rows read: ') as bar:
        for chunk in pd.read_csv(csv_path, chunksize=chunksize,
                                 dtype=dtype, index_col=0):
            chunk_list.append(chunk)
            bar.update(len(chunk))
    df = pd.concat((f for f in chunk_list), axis=0)
    print('Done!')
    
    return df

def parse_files(input_dir, files):
    
    df_list = []
    for i in range(len(files)):
        f = files[i]
        file_path = input_dir+f
        print(f'Parsing {f} (exp #{i+1}) ...')
        df = read_csv_pgbar(file_path, 10**4)
        
        # // TODO: make better file type check //
        cols = [c for c in df]
        if len(cols) < 2:
            try:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
            except:
                print(f'{f} in an invalid file type.')
        
        df.index.names = ['orthogroup']
        df.reset_index(inplace=True)
        df_list.append(df)
    
    return df_list

def main():
    
    t0 = time.time()
    
    file_list = get_files(args.input_dir, args.ext)
    df_list = parse_files(args.input_dir, file_list)
    
    print('Merging all experiments ...')
    # // TODO: progress indicator //
    df_merged = reduce(lambda x, y: pd.merge(x, y, on='orthogroup', how='outer'), df_list)
    df_merged.fillna(0, inplace=True)
    print(df_merged.head())
    print(f'Writing output to {args.outfile} ...')
    df_merged.to_csv(args.outfile, index=False)
    
    print(f"Total time to completion: {time.time() - t0} seconds")
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Concatenate all .elut files in a given directory with a given extension.")
    parser.add_argument("--input_dir", action="store", required=True, \
                        help="Directory containing elution profiles")
    parser.add_argument("--ext", action="store", required=True, \
                        help="Extension for files being targeted; e.g., 'mFDRpsm001.unique.norm.elut'")
    parser.add_argument("--outfile", action="store", required=True, \
                        help="Path/filename for output file")
    
    # verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbosity (-v, -vv, etc)")
    
    args = parser.parse_args()
    main()