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
from beautifultable import BeautifulTable

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

def count_grps(x):
    if type(x) == list:
        return(len(x))
    else:
        return(1)

def count_multi(df, col_name):
    n = 0
    for i in range(len(df)):
        if type(df[col_name][i]) == list:
            n += 1
    return(n)

def check_uniq_ppis(fmat):
    table = BeautifulTable()
    table.columns.header = ["TEST DESCRIPTION", "STATUS"]
    ppi_counts = fmat.groupby(['ID']).size().sort_values(ascending=False)
    assert any(x != 1 for x in ppi_counts) == False, "Non-unique PPI labels detected."
    #print("Assert all labeled rows are unique (group by & count IDs, all n==1)\tPASS")
    table.rows.append(["Assert exactly one label and complex group for each PPI pair\n(group by IDs & count rows, pass if all n==1)", "PASS"])
    
    fmat['fsets'] = [make_fset(i, drop=True) for i in fmat['ID']]
    fset_counts = fmat.groupby(['fsets']).size().sort_values(ascending=False)
    assert any(x != 1 for x in fset_counts) == False, "Non-unique PPI labels detected."
    #print("Assert all labeled PPI pairs are unique (group by & count frozen sets, all n==1)\tPASS")
    table.rows.append(["Assert each PPI pair is unique & appears exactly once\n(group by frozen set & count rows, pass if all n==1)", "PASS"])
    print(table,'\n')

def check_group_merge(fmat):
    table = BeautifulTable()
    table.columns.header = ["", "BEFORE MERGE", "AFTER MERGE"]
    fmat_tt = fmat[(fmat.label==1)|(fmat.label==-1)]
    fmat_tt = fmat_tt.reset_index(drop=True)
    counts = (fmat_tt.explode('group')).nunique()
    n_uniq_pairs, n_pre, n_post = counts[0:3]
    n_ppi_pre = count_multi(fmat_tt, 'group')
    n_ppi_post = count_multi(fmat_tt, 'super_group')
    table.rows.append(["# unique PPI pairs labeled", n_uniq_pairs, n_uniq_pairs])
    table.rows.append(["# protein complex groups", n_pre, n_post])
    table.rows.append(["# PPIs assigned to >1 complex group", n_ppi_pre, n_ppi_post])
    print(table,'\n')

"""
Sanity check functions for predict_ppis.py
"""

def check_traintest_overlap(res):
    table = BeautifulTable()
    table.columns.header = ["TEST DESCRIPTION", "STATUS"]
    all_ppi_counts = res.groupby(['ID']).size().sort_values(ascending=False)
    assert any(x != 1 for x in all_ppi_counts) == False, "Non-unique PPI labels detected."
    table.rows.append(["Assert all unlabeled PPI pairs are unique\n(group by IDs & count rows, pass if all n==1)", "PASS"])
    
    tt_df = res[(res.set == 'test') | (res.set == 'train')]
    tt_ppi_counts = tt_df.groupby(['ID']).size().sort_values(ascending=False)
    assert any(x != 1 for x in tt_ppi_counts) == False, "Non-unique PPI labels detected."
    table.rows.append(["Assert all train/test PPI pairs are unique\n(group by IDs & count rows, pass if all n==1)", "PASS"])
    
    train_df = tt_df[(tt_df.set == 'train')]
    train_ppis = [make_fset(i, drop=True) for i in train_df['ID']]
    test_df = tt_df[(tt_df.set == 'test')]
    test_ppis = [make_fset(i, drop=True) for i in test_df['ID']]
    assert len(set(test_ppis) & set(train_ppis)) == 0, "Overlap between train and test PPIs detected."
    table.rows.append(["Assert no overlap between train/test PPI frozen sets\n(pass if len(set(test_ppis) & set(train_ppis)) == 0)", "PASS"])
    print(table,'\n')

""" Main """

def main():
    
    if args.featmat_file:
        print(f'Reading in {args.featmat_file} ...')
        with open(args.featmat_file, 'rb') as handle:
            fmat = pickle.load(handle)
        print('Evaluating labeled feature matrix ...')
        # assert no ppis are repeated in labeled feature matrix
        check_uniq_ppis(fmat)
        
        # evaluate complex merge behavior
        if 'group' in fmat.columns and 'super_group' in fmat.columns:
            print('Evaluating protein complex group merge ...')
            check_group_merge(fmat)
        
    if args.results_file:
        print(f'Reading in {args.results_file} ...')
        res = pd.read_csv(args.results_file)
        print('Evaluating machine learning data/results ...')
        # assert that there is no overlap between train and test ppis
        check_traintest_overlap(res)
        
""" When executed from the command line: """
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify feature matrix
    parser.add_argument("-f", "--featmat_file", help="(Optional) Path to pickled feature matrix that was output by `label_featmat.py`. PPI ID column is a column named 'ID' with space separated protein ID pairs. If both 'group' and 'super_group' columns are present, removal of complex overlap will be evaluated.")

    # specify positive labels
    parser.add_argument("-r", "--results_file", action="store", help="(Optional) Path to results file output by `predict_ppis.py`. Column names are 'ID', 'label', 'ppi_score' and 'set'. The program will assert that there is only 1 row per unique PPI pair, and there is no overlap in the train/test PPI pairs (e.g., no circularity due to PPIs belonging to multiple complexes).")
    
    args = parser.parse_args()
    main()
    