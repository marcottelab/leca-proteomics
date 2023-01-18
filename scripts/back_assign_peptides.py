import csv
import argparse
import pandas as pd
from Bio import SeqIO
from time import sleep
from tqdm import tqdm

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Take a comma-separated .csv file containing peptide assigments to \
    	eggNOG-collapsed protein groups and relate those peptides back to their original, uncollapsed FASTA entries. \
    	Requires unique peptide matches. Outputs a CSV file mapping the peptides assigned to each orthogroup \
    	back to protein FASTA entries assigned to that orthogroup.")
    parser.add_argument("-p", "--peptides", action="store", required=True,
                                        help="Comma-separated file containing eggNOG groups and peptide matches")
    parser.add_argument("-f", "--fasta", action="store", required=True,
                                        help="Original FASTA file (that was mapped to eggNOG groups) containing \
                                        uncollapsed entries")
    parser.add_argument("-m", "--mapping", action="store", required=True,
                                        help="Tab-delimited file mapping eggNOG groups to FASTA entries")
    parser.add_argument("-o", "--outfile", action="store", required=False,
                                        help="Name for outfile (optional)")
    args = parser.parse_args()


fasta = args.fasta

if args.outfile == None:
	writefile = args.peptides+'.back_assignments'
else:
	writefile = args.outfile

print('Reading in data ...')
print(f'\tCFMS peptide file:\n\t{args.peptides}')
print(f'\teggNOG mapping file:\n\t{args.mapping}')
print(f'\tFASTA file:\n\t{args.fasta}')
# read in peptide assignments
pep_df = pd.read_csv(args.peptides)
count_col = pep_df.columns[-1]
pep_df.head()

# read in eggnog mapping file
map_df = pd.read_csv(args.mapping, sep='\t')
map_df.head()

# define orthogroup lookup
og_lookup = [o for o in pep_df['protein'].tolist() if o.startswith(('KOG', 'ENOG'))]

print(f'Performing look up for {len(set(og_lookup))} orthogroups ...')
# perform lookup for each orthogroup
df_list = []
fasta = args.fasta
for og in tqdm(set(og_lookup)):
    family = map_df[map_df['ID'] == og]['ProteinID'].tolist()
    pep_lookup = pep_df[pep_df['protein'] == og]['peptide'].tolist()
    lookup_res = []
    for record in SeqIO.parse(open(fasta,"r"), "fasta"):
        prot_id = record.id
        seq = str(record.seq.upper())
        if prot_id in family:
            for pep in pep_lookup:
                if pep in seq:
                    count = pep_df.loc[pep_df['peptide'] == pep][count_col].tolist()
                    lookup_res.append([og, pep, prot_id, count[0]])
    df = pd.DataFrame(lookup_res, columns=['orthogroup', 'peptide', 'protein_match', 'count'])
    df_list.append(df)

# merge lookups for each orthogroup
print('Merging look up results for each orthogroup ...')
final_df = pd.concat(df_list, ignore_index=True, sort=False)

# write result
print(f'Writing result to {writefile} ...')
final_df.to_csv(writefile, index=False)
print(final_df)
print('Done!')