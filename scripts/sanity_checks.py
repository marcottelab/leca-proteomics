"""
Script for sanity checking the output of various stages of the PPI machine learning pipeline.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pandas as pd
import pickle

"""
Sanity check functions for label_featmay.py
"""

def make_fset(x, drop=True):
    if len(set(x.split(' '))) < 2:
        if drop == False:
            print(f"WARNING: Features for '{x}' (self-self PPI) detected ...")
            x1 = x.split(' ')[0]
            fset = frozenset({x1})
            return(fset)
        else:
            return()   
    else:
        x1 = x.split(' ')[0]
        x2 = x.split(' ')[1]
        fset = frozenset({x1,x2})
        return(fset)

def check_uniq_ppis(fmat):
    ppi_counts = fmat.groupby(['ID']).size().sort_values(ascending=False)
    assert any(x > 1 for x in ppi_counts) == False, "Non-unique PPI labels detected."
    
    fmat['fsets'] = [make_fset(i, drop=True) for i in fmat['ID']]
    fset_counts = fmat.groupby(['fsets']).size().sort_values(ascending=False)
    assert any(x > 1 for x in fset_counts) == False, "Non-unique PPI labels detected."
    
def check_group_merge(fmat):

"""
Sanity check functions for predict_ppis.py
"""

def check_traintest_overlap(res):
    
    all_ppi_counts = res.groupby(['ID']).size().sort_values(ascending=False)
    assert any(x > 1 for x in all_ppi_counts) == False, "Non-unique PPI labels detected."
    
    tt_df = res[(res.set == 'test') | (res.set == 'train')]
    tt_ppi_counts = tt_df.groupby(['ID']).size().sort_values(ascending=False)
    assert any(x > 1 for x in tt_ppi_counts) == False, "Non-unique PPI labels detected."
    
    train_df = tt_df[(tt_df.set == 'train')]
    train_ppis = [make_fset(i, drop=True) for i in train_df['ID']]
    test_df = tt_df[(tt_df.set == 'test')]
    test_ppis = [make_fset(i, drop=True) for i in test_df['ID']]
    assert len(set(test_ppis) & set(train_ppis)) == 0, "Overlap between train and test PPIs detected."

""" Main """

def main():
    
    if args.featmat_file:
        with open(fmat_file, 'rb') as handle:
            fmat = pickle.load(handle)
        
        # assert no ppis are repeated in labeled feature matrix
        check_uniq_ppis(fmat)
        
    if args.results_file:
        res = pd.read_csv(results_file)
        
        # assert that there is no overlap between train and test ppis
        check_traintest_overlap(res)
        
""" When executed from the command line: """
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify feature matrix
    parser.add_argument("-f", "--featmat_file", help="(Optional) Path to pickled feature matrix that was output by `label_featmat.py`. PPI ID column is a column named 'ID' with space separated protein ID pairs. If both 'group' and 'super_group' columns are present, removal of complex overlap will be evaluated.")

    # specify positive labels
    parser.add_argument("-r", "--results_file", action="store", help="(Optional) Path to results file output by `predict_ppis.py`. Column names are 'ID', 'label', 'ppi_score' and 'set'. The program will assert that there is only 1 row per unique PPI pair, and there is no overlap in the train/test PPI pairs (e.g., no circularity due to PPIs belonging to multiple complexes).")
    