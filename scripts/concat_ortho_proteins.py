import csv
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Collapses proteins in a FASTA file based on shared orthology, where the protein sequences are concatenated with a triple lysine separator.")
    parser.add_argument("-m", "--orthomap", action="store", required=True,
                                        help="Tab-separated file mapping protein groups to eggNOG groups where \
                                        the columns are labeled 'ID' for the eggNOG groups and 'ProteinID' for the protein. \
                                        This is the output from the format_emapper script.")
    parser.add_argument("-f", "--fasta", action="store", required=True,
                                        help="FASTA file to be collapsed")
    parser.add_argument("-o", "--out", action="store", required=False,
                                        help="Name for outfile (optional)")
    parser.add_argument("-l", "--limit_length", action="store_true", required=False,
                                        help="Limit sequence size of collapsed orthogroups; groups that grow too large \
                                        (>=100,000aa) are broken back up into the constituent proteins.")

    args = parser.parse_args()


if args.out == None:
    new_ext = '.collapsed.fasta'
    writefile = new_ext.join(args.fasta.rsplit('.fasta', 1))
    print("collapsed fasta outfile name:", writefile)
else:
    writefile = args.out
    print("collapsed fasta outfile name:", writefile)

def write_fasta_entry(ID, seq, file):
    file.write('>'+ID+'\n')
    file.write(str(seq)+'\n')
    print('>'+ID)
    print(seq)

# write to-be-collapsed FASTA to dicitonary
fasta_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))

# write orthology mapping to dictionary
ortho_dict = {}
with open(args.orthomap, 'r') as f:
    for line in f:
        proteins = line.strip().split('\t')[0]
        ortho_dict[proteins] = line.strip().split('\t')[1]

# write dictionary that concatenates FASTA entries based on shared orthology
concat_dict = {}
for protein in fasta_dict:
    if protein not in ortho_dict:
        concat_dict[protein] = [str(fasta_dict[protein].seq)]
    elif ortho_dict[protein] not in concat_dict:
        concat_dict[ortho_dict[protein]] = [str(fasta_dict[protein].seq)]
    else:
        concat_dict[ortho_dict[protein]] += ['KKK'+str(fasta_dict[protein].seq)]

for ids in concat_dict:
    concat_dict[ids] = ''.join(map(str,concat_dict[ids]))

# if max entry length is NOT enabled:
# write concatenated dictonary to file
if args.limit_length == False:
    with open(writefile, 'w') as f:
        for i in concat_dict:
            write_fasta_entry(i, concat_dict[i], f)

# if max entry length is enabled: 
# find the ortho IDs with the biggest entries and break them into their component protein sequences
else:
    big_seq_OG_ids = []
    for ids in concat_dict:
        seq_len = len(concat_dict[ids])
        if seq_len >= 100000:
            big_seq_OG_ids.append(ids)

    big_seq_prot_ids = []
    for prots, ogs in ortho_dict.items():
        if ogs in big_seq_OG_ids:
            if prots in fasta_dict:
                big_seq_prot_ids.append(prots)

    limited_dict = {}
    for ids in concat_dict:
        if ids not in big_seq_OG_ids:
            limited_dict[ids] = concat_dict[ids]
    
    for ids in big_seq_prot_ids:
        limited_dict[ids] = str(fasta_dict[ids].seq)

    with open(writefile, 'w') as f:
        for i in limited_dict:
            write_fasta_entry(i, limited_dict[i], f)

################################
# same thing but with pandas
################################

#ortho_df = pd.read_csv(args.orthomap, sep='\t')

#print(ortho_df.head())

#identifiers = []
#seqs = []
#for record in SeqIO.parse(open(fasta, 'r'), 'fasta'):
#    identifiers.append(record.id)
#    seqs.append(''.join(record.seq))

#fasta_df = pd.DataFrame()
#fasta_df['ProteinID'] = identifiers
#fasta_df['Sequences'] = seqs

#print(fasta_df.head())

#all_df = pd.merge(fasta_df, ortho_df, how='left', on=['ProteinID'])
#all_df['ProteinID'] = all_df['ProteinID'].str.split('|').str[1]
#all_df['ID'] = all_df['ID'].fillna(all_df['ProteinID'])

#print(all_df.head())

#grouped_df = all_df[['ID','Sequences']]
#grouped_df = grouped_df.groupby(['ID'])
#concatseq_df = grouped_df['Sequences'].agg(lambda x: 'KKK'.join(x)).reset_index()
#print(concatseq_df)

#concatseq_df.set_index(['ID'])
#concatseqs_dict = dict(zip(concatseq_df.ID, concatseq_df.Sequences))

#new_entries = []
#with open(writefile+'.test', 'w') as f:
#    for ids in concatseqs_dict:
#        record = SeqRecord(Seq(concatseqs_dict[ids]), ids, '', '')
#        new_entries.append(record)
#        f.write('>{}\n{}\n'.format(ids, concatseqs_dict[ids]))

#SeqIO.write(new_entries, writefile, "fasta")

    # a different way to write SeqRecord objects to file:
    #new_entries = []
    #for ids in pd_dict:
    #    record = SeqRecord(Seq(pd_dict[ids]), ids, '', '')
    #    new_entries.append(record)

    #with open(writefile+'.pd', 'w') as f:
    #    SeqIO.write(new_entries, f, "fasta")







