import csv
import argparse
import pandas as pd
from Bio import SeqIO

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Take a comma-separated .csv file containing peptide assigments to \
    	eggNOG-collapsed protein groups and relate those peptides back to their original, uncollapsed FASTA entries. \
    	Requires unique peptide matches. Outputs a tab-delimited file mapping the peptides assigned to each orthogroup \
    	back to protein FASTA entries assigned to that orthogroup.")
    parser.add_argument("-p", "--peptides", action="store", required=True,
                                        help="Comma-separated file containing eggNOG groups and peptide matches")
    parser.add_argument("-f", "--fasta", action="store", required=True,
                                        help="Original FASTA file (that was mapped to eggNOG groups) containing \
                                        uncollapsed entries")
    parser.add_argument("-o", "--out", action="store", required=False,
                                        help="Name for outfile (optional)")
    parser.add_argument("-c", "--collapse", action="store_true", required=False,
                                        help="Specify this argument to group proteins as a comma-separated list if more \
                                        than one is assigned to a peptide (optional)")
    args = parser.parse_args()


fasta = args.fasta

if args.out == None:
	writefile = args.peptides+'.back_assignments'
else:
	writefile = args.out+'.back_assignments'
	
pep_df = pd.read_csv(args.peptides, header=None)
pep_df.columns = ['eggNOG_ID','peptide_match']

pep_list = []
for pep in pep_df.iloc[:,1]:
	pep_fmt = pep.replace(".","")
	pep_list.append(pep_fmt)

pep_df['peptide_fmt'] = pep_list
print(pep_df.head())

pep_lookup = set(pep_list)

lookup_df = pd.DataFrame(columns=['peptide_fmt', 'protein_match'])

for record in SeqIO.parse(open(fasta,"r"), "fasta"):
	sequence = str(record.seq.upper())
	for pep in pep_lookup:
		num_occurrences = sequence.count(pep)
		if num_occurrences > 0:
			lookup_df = lookup_df.append({'peptide_fmt':pep, 'protein_match':record.id}, 
				ignore_index=True)

final_df = pd.merge(pep_df, lookup_df, how='left', on=['peptide_fmt'])

if args.collapse:
	collapsed_df = final_df.groupby(['eggNOG_ID','peptide_match','peptide_fmt'],sort=False).agg(lambda x: set(x)).reset_index()
	collapsed_df['protein_match'] = collapsed_df['protein_match'].apply(list)
	collapsed_df['protein_match'] = [','.join(map(str,x)) for x in collapsed_df['protein_match']]
	collapsed_df.to_csv(writefile+'.collapsed', index=False, sep='\t')
	print(collapsed_df)
else:
	final_df.to_csv(writefile, index=False, sep='\t')
	print(final_df)




