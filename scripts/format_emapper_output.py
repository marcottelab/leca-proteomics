"""
Script for formatting orthology mapping results as output by emapper 2.0.5
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import re
import pandas as pd

""" Functions """
def format_results(ogs):
    # format eggnog results
    res_fmt = []
    results = ogs.split(',')
    for i in results:
        og = re.search(".*(?=@)", i)[0]
        if og and og[0].isdigit():
            og = 'ENOG50'+og
        lvl = re.search("(?<=@).*(?=\|)", i)[0]
        taxa = re.search("(?<=\|).*", i)[0]
        res_fmt.append([og, lvl, taxa])
    return(res_fmt)


""" Main """
def main():
    
    enog = pd.read_csv(args.infile, sep='\t')
    enog.columns.values[0] = 'ProteinID'
    enog.columns.values[4] = 'ID'
    enog = enog[['ProteinID','ID']]

    results_dict = dict()
    for i in range(len(enog)):
        prot_id = enog['ProteinID'][i]
        if args.format_uniprot_ids:
            prot_id = re.search("(?<=\|).*(?=\|)", prot_id)[0]
        ogs = enog['ID'][i]
        res_fmt = format_results(ogs)
        results_dict[prot_id] = res_fmt

    df = pd.DataFrame.from_dict(results_dict, orient='index')
    df = (df
          .reset_index()
          .melt(id_vars='index')
          .drop(labels='variable',axis=1)
          .dropna()
          .rename(columns={'index':'ProteinID'})
         )

    df[['ID','level','taxonomy']] = pd.DataFrame(df.value.tolist(), index=df.index)
    df['level'] = df['level'].astype(int)
    df.drop(labels='value', axis=1, inplace=True)

    level_df = df[df.level==args.level]
    out_df = level_df[['ProteinID', 'ID']]
        
if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # specify infile
    parser.add_argument("-i", "--infile", help="(Required) path to '.annotations' file output by emapper (v2.0.5)")

    # specify desired taxonomic level
    parser.add_argument("-l", "--level", action="store", type=int, help="(Required) specify desired taxonomic mapping to extract; full list of options: http://eggnog5.embl.de/#/app/downloads")

    # specify outfile name
    parser.add_argument("-o", "--outfile", action="store", help="(Required) specify the outfile path/name")
    
    # specify if you want to format your uniprot IDs (e.g., if using for annotation reasons)
    parser.add_argument("--format_uniprot_ids", action="store_true", default=False, help="(Optional) specify if you want to format your uniprot IDs, typically for downstream annotation reasons; e.g., 'sp|P17442|PHO81_YEAST' becomes 'P17442'")

    # specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main()