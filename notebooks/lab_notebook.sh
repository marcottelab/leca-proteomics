#######################################################################################
#######################################################################################
# PROJECT DIRECTORY: /project/rmcox/LECA/
# CODE FOR DATA ANALYSIS AND FIGURES: /project/rmcox/LECA/LECA.Rmd
#######################################################################################
#######################################################################################

#######################################################################################
# Mapping positive controls (qualifying exam, February 2020)
#######################################################################################

## eggnog map positive controls derived from Koumandou et al. 2013
git clone https://github.com/jhcepas/eggnog-mapper.git

# newest version (v2), diamond
nice -5 python /project/rmcox/programs/eggnog-mapper/emapper.py -i ${exp}.fasta --output ${exp}.mapping --output_dir /project/rmcox/LECA/positive_controls/eggnog_mapping --scratch_dir /project/rmcox/eggnog_mapping_scratch -m diamond --predict_ortho --keep_mapping_files --cpu 18 --override

while read exp; do echo "nice -5 python /project/rmcox/programs/eggnog-mapper/emapper.py -i ${exp}.fasta --output ${exp}.diamond --output_dir /project/rmcox/LECA/positive_controls/eggnog_mapping --scratch_dir /project/rmcox/eggnog_mapping_scratch -m diamond --predict_ortho --target_taxa Eukaryotes --keep_mapping_files --cpu 32 --override"; done < explist.txt > map_LECA_controls_diamond.sh

cat map_LECA_controls_diamond.sh | parallel -j32

# older version (v1), hmmer
#################### couldn't get this to work
nice -5 python /project/rmcox/programs/eggnog-mapper/emapper.py -i ${exp}.fasta --output ${exp}.mapping --output_dir /project/rmcox/LECA/positive_controls/eggnog_mapping --scratch_dir /project/rmcox/eggnog_mapping_scratch -m hmmer -d euNOG --predict_ortho --keep_mapping_files --cpu 18 --override

while read exp; do echo "nice -5 python /project/rmcox/programs/eggnog-mapper/emapper.py -i ${exp}.fasta --output ${exp}.hmmer --output_dir /project/rmcox/LECA/positive_controls/eggnog_mapping --scratch_dir /project/rmcox/eggnog_mapping_scratch -m hmmer -d euNOG --predict_ortho --keep_mapping_files --cpu 18 --override"; done < explist.txt > map_LECA_controls_hmmer.sh

cat map_LECA_controls_hmmer.sh | parallel -j18
#################### couldn't get this to work

python /project/rmcox/scripts/format_emapper_output.py -f clathrin_complex_GO71439.diamond.emapper.annotations -o clathrin_complex_GO71439.diamond.mapping --diamond

#######################################################################################
# Building gene trees (qualifying exam, February 2020)
#######################################################################################

## ETE toolkit

#################### couldn't get this to work
# use miniconda to install
wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/anaconda_ete/
export PATH=~/anaconda_ete/bin:$PATH;

conda install -c etetoolkit ete3 ete_toolchain # install ETE

ete3 build check # check all external apps are available

# run this line before using ETE:
export PATH=~/anaconda_ete/bin:$PATH; # add to .bashrc?


## getting genome set
while read path; do echo "rsync --copy-links --recursive --times --verbose ${path} /project/rmcox/LECA/notung/genomes/"; done < genome_ftp_paths.txt > get_genomes.sh
#################### couldn't get this to work

## FastTree (this looks promising)

## general info:
# takes in MSA file to make a tree; uses JTT, WAG or LG models of amino acid evolution
# to account for varying rates of evolution across sites, uses a single rate for each site ("CAT" approximation); can override this with Gamma option (-gamma). let's do this since it allows statistical comparisons of likelihoods, and it's only 5% slower

# installation:
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm

# in ~/bin:
ln -s /stor/work/Marcotte/project/rmcox/programs/FastTree

# general use:
FastTree alignment.file > tree_file

# test on AP1B1 family (KOG1061)
FastTree -gamma KOG1061_alignment_trimmed.fasta > KOG1061_tree.nwk

# consider using IQ-tree instead

# going to try using the UniProt reference proteomes to build trees instead
# Claire already mapped the eukaryotic ref proteomes. Now I'm going to do the same with prokaryotes

while read proteome; do echo "wget -nc ftp://ftp.uniprot.org:21/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/${proteome}"; done < uniprot_proteome_list_bact.txt > record_ftp_commands_reference_proteomes_bact.sh

# for some reason this isn't working
while read proteome; do echo "wget -nc ftp://ftp.uniprot.org:21/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Archaea/${proteome}"; done < uniprot_proteome_list_arch.txt > record_ftp_commands_reference_proteomes_arch.sh

# going to try this
wget -nc ftp://ftp.uniprot.org:21/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Archaea/*fasta.gz

# weirdly I only got 295 archaea proteomes... there isn't that many on the ftp site
# uniprot database says ~2100 proteomes, ftp does not agree... perhaps only includes reference-status proteomes?
# in any case, these are the commands for eggnog mapping

# archaea:
python /project/rmcox/programs/eggnog-mapper/emapper.py -i /project/rmcox/LECA/notung/proteomes_arch/${proteome} --output /project/rmcox/LECA/notung/eggnog_mapping/${proteome}.euNOG.diamond -d euNOG --override --resume --scratch_dir /project/rmcox/LECA/notung/eggnog_mapping/scratch --no_file_comments --keep_mapping_files -m diamond

while read proteome; do echo "python /project/rmcox/programs/eggnog-mapper/emapper.py -i /project/rmcox/LECA/notung/proteomes_arch/${proteome} --output /project/rmcox/LECA/notung/eggnog_mapping/${proteome}.euNOG.diamond -d euNOG --override --resume --scratch_dir /project/rmcox/LECA/notung/eggnog_mapping/scratch --no_file_comments --keep_mapping_files -m diamond"; done < uniprot_proteome_list_arch.txt > record_eggnog_mapping_arch.sh

cat record_eggnog_mapping_arch.sh | parallel -j25

# bacteria:
ls | sort -R | tail -100 > uniprot_proteome_list_bact.txt

while read proteome; do echo "nice -5 python /project/rmcox/programs/eggnog-mapper/emapper.py -i /project/rmcox/LECA/notung/proteomes_bact/${proteome} --output /project/rmcox/LECA/notung/eggnog_mapping/${proteome}.euNOG.diamond -d euNOG --override --resume --scratch_dir /project/rmcox/LECA/notung/eggnog_mapping/scratch --no_file_comments --keep_mapping_files -m diamond"; done < uniprot_proteome_list_bact.txt > record_eggnog_mapping_bact.sh

cat record_eggnog_mapping_bact.sh | parallel -j25

#######################################################################################
# RefSeq genome stats (qualifying exam, February 2020)
#######################################################################################

# format stat files
for f in *.txt; do sed -i 's/,//' $f; done
for f in *.txt; do sed -i '/^$/d' $f; done

for f in *.txt; do awk '{print $4, $5}' $f > ${f}_fmt.txt; done
for f in *_fmt.txt; do sed -i '/^$/d' $f; done

#######################################################################################
# Processing Halobacterium
#######################################################################################

# data: /MS/processed/Fusion_data/Ophelia/Halobacterium_032020/
cp -r /MS/processed/Fusion_data/Ophelia/Halobacterium_032020/* halobacterium/ 

# run MSblender on hopper
vim exp_list.txt  # names of folders containing DB, mzXML, working and output directories

cat uniprot_HALSA_proteome_UP000000554.fasta contam.fasta > uniprot_HALSA_proteome_UP000000554_contam.combined.fasta

/stor/work/Marcotte/project/rmcox/programs/run_msblender$ bash scripts/create_commands.sh exp_list.txt ../MSblender/runMS2.sh ../SearchGUI-3.3.5/

/stor/work/Marcotte/project/rmcox/programs/run_msblender$ for f in halobacterium_COMMANDS.txt; do cat $f | parallel -j20; done

python /stor/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/msblender2elution.py --prot_count_files /project/rmcox/programs/run_msblender/halobacterium/output/*.group --output_filename halobacterium.prot_count_mFDRpsm001.unique.elut --fraction_name_from_filename --parse_uniprot_id --remove_zero_unique

python /project/cmcwhite/github/protein_complex_maps/protein_complex_maps/preprocessing_util/elutionTidy.py --input_elution /stor/work/Marcotte/project/rmcox/programs/run_msblender/halobacterium/halobacterium.prot_count_mFDRpsm001.unique.elut --outfile halobacterium.unique.protcount.tidy --firstcol ProteinID --valuecol SpectralCounts

awk -F',' '{print $2}' halobacterium.unique.protcount.tidy | sort -u | wc -l  # 2337 unique protein IDs

grep -c '^>' uniprot_HALSA_proteome_UP000000554_contam.combined.fasta  # 2674 possible protein IDs

## edit for archaea
python /stor/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/elutionCluster.py --input_elutions /project/rmcox/programs/run_msblender/halobacterium/halobacterium.prot_count_mFDRpsm001.unique.elut --outfile halobacterium.unique.protcount.mFDRpsm001.xlsx --hclust_metric pearson --hclust_method average --path_to_annotations /project/rmcox/proteomes/halsa/uniprot-proteome_UP000000554.tab

#######################################################################################
# LECA Genome Derivation
# Dollo Parsimony of QoF/SwissProt Reference Species Orthogroups (May 2020)
#######################################################################################

#### download data:

# 78 reference organisms (RELEASE = 22-Apr-2020; RETRIEVED 4-May-2020)
# taxonomic division folders (Archaea, Bacteria and Eukaryota) containing
#  *.gene2acc, *.fasta (canonical sequences), *_DNA.fasta, *.xml and *.idmapping
#  files for 78 species, plus 42 *_additional.fasta (isoform sequences) and 50
#  *_DNA.miss files
wget ftp://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/QfO_release_2020_02.tar.gz
tar -xvzf QfO_release_2020_02.tar.gz

# extracted the table statistics from the README & cleaned up in R to make this file: 
cat /project/rmcox/LECA/qfo_tree/speciesinfo.txt
# columns (8): proteome_id, tax_id, species_code, tax_group, num_fasta_entries, num_isoforms, num_gene2acc, species_name

ls qfo_tree/Archaea/*[0-9].fasta | wc -l   # 7 archaea
realpath Archaea/*[0-9].fasta > archaea_proteomes.txt

ls qfo_tree/Bacteria/*[0-9].fasta | wc -l   # 23 bacteria
realpath Bacteria/*[0-9].fasta > bacteria_proteomes.txt

ls qfo_tree/Eukaryota/*[0-9].fasta | wc -l   # 48 eukarya
realpath Eukaryota/*[0-9].fasta > eukarya_proteomes.txt

cat *proteomes.txt > all_proteome_paths.txt


#### map eggnog groups:

# ***NOTE:*** this will exclude isoforms

# create commands
cat all_proteome_paths.txt | { while read path; do proteome=`basename $path .fasta`; echo "nice -5 python /project/rmcox/programs/eggnog-mapper/emapper.py -i ${path} --output /project/rmcox/LECA/qfo_tree/eggnog_mapping/${proteome}.euNOG.diamond -d euNOG --override --resume --scratch_dir /project/rmcox/LECA/qfo_tree/eggnog_mapping/scratch --no_file_comments --keep_mapping_files -m diamond"; done } > map_qfo_organisms.sh

# execute eggnog mapper
cat map_qfo_organisms.sh | parallel -j32


# format eggnog output for eukarya-level taxonomy (2759)
cat all_proteome_paths.txt | { while read path; do proteome=`basename $path .fasta`; echo "Rscript ../scripts/format_emapper_output.R -f eggnog_mapping/${proteome}.euNOG.diamond.emapper.annotations -o eggnog_mapping/${proteome}.euNOG.diamond.mapping.2759 -s diamond -l 2759 -v 2 -p ENOG50"; done } > format_emapper_output_2759.sh

# format eggnog output for highest level taxonomy (1)
cat all_proteome_paths.txt | { while read path; do proteome=`basename $path .fasta`; echo "Rscript ../scripts/format_emapper_output.R -f eggnog_mapping/${proteome}.euNOG.diamond.emapper.annotations -o eggnog_mapping/${proteome}.euNOG.diamond.mapping.1 -s diamond -l 1 -v 2 -p ENOG50"; done } > format_emapper_output_1.sh

cat records/format_emapper_output_* | parallel -j32

cat eggnog_mapping/*mapping.1 > allspecies_1.mapping
wc -l allspecies_1.mapping # 842063 genes mapped

cat eggnog_mapping/*mapping.2759 > allspecies_2759.mapping
wc -l allspecies_2759.mapping # 745713 genes mapped

# format uniprot IDs for future purposes:
for x in *mapping.2759; do Rscript /project/rmcox/LECA/scripts/format_uniprot_ids.R -f $x; done

# create matrix of counts for each species and ENOG
# actually I'm going to do this in R

# so the species tree contains 80 species that are not included in what I mapped to OGs
# need to retrieve those proteomes from UniProt by their taxonomy IDs...

# perl access to the uniprot API is written into a script

# to get full proteomes:
my $query = "https://www.uniprot.org/uniprot/?query=organism:$taxon&format=fasta&include=no";

# to get reference proteomes:
my $query = "https://www.uniprot.org/uniprot/?query=organism:$taxon&format=fasta&fil=reference&include=no"

# to get annotations for those proteoms
my $query = "https://www.uniprot.org/uniprot/?query=organism:$taxon&columns=id,entry name,reviewed,protein names,genes,organism,length&format=tab";

# execute for all tax ids:
while read id; do echo "perl ../../scripts/retrieve_proteomes.pl ${id}"; done < taxid_list.txt > retrieve_proteomes.sh

readlink -f *fasta > proteome_paths.txt

cat proteome_paths.txt | { while read path; do proteome=`basename $path .fasta`; echo "nice -5 python /project/rmcox/programs/eggnog-mapper/emapper.py -i ${path} --output /project/rmcox/LECA/qfo_tree/eggnog_mapping/${proteome}.euNOG.diamond -d euNOG --override --resume --scratch_dir /project/rmcox/LECA/qfo_tree/eggnog_mapping/scratch --no_file_comments --keep_mapping_files -m diamond"; done } > map_additional_organisms.sh

cat map_additional_organisms.sh | parallel -j40

# formatting additional outputs
cat all_proteome_paths.txt | { while read path; do proteome=`basename $path .fasta`; echo "Rscript ../scripts/format_emapper_output.R -f eggnog_mapping/${proteome}.euNOG.diamond.emapper.annotations -o eggnog_mapping/${proteome}.euNOG.diamond.mapping.1 -s diamond -l 1 -v 2 -p ENOG50"; done } > format_emapper_output_1.sh

cat eggnog_mapping/*mapping.1 > allspecies_1.mapping
wc -l allspecies_1.mapping # 2881617 genes mapped

# the swissprot tree species codes don't match the uniprot codes for the following organisms:
# ST=SPETR, UP=ICTTR, TAXID=43179
# ST=CANFA, UP=CANLF, TAXID=9615
# ST=SULSO, UP=SACS2, TAXID=273057

# Also just found this one, doesn't match the UniProt controlled vocabulary of species
# but the proteome for 13642 still denotes 'POLPA' even though it matches 'HETPA' in spec list
# ST=POLPA, UP=HETPA, TAXID=13642

# let's edit the species tree to make it match the new UP IDs..
vim species_tree/speciestree_fmt.nwk
:%s/SPETR/ICTTR/
:%s/CANFA/CANLF/
:%s/SULSO/SACS2/

#### compute dollo parsimony:

wget http://www.iro.umontreal.ca/~csuros/gene_content/Count.tgz
tar -xvzf Count.tgz
rm Count.tgz

# path to Java executable: /project/rmcox/programs/Count/Count.jar

# from within /project/rmcox/LECA/qfo_tree:
# IMPORTANT: remember to (1) start running XMing, (2) set X11 to froward, then (3) start the application
java -Xmx1024M -jar /project/rmcox/programs/Count/Count.jar  # might need to give it more memory

# vim command for formatting the .nwk species tree:
:%s/__\w\+_\w\+__\d\+//g

#######################################################################################
# Internal CFMS data (January 2021)
#######################################################################################

#### create directory structure
awk -F',' '{print tolower($5)}' speciesinfo_01252021.csv | tail -n +2 > speciescodes.txt

while read id; do echo "mkdir /project/rmcox/LECA/ms/cfms/raw/${id}"; done < codes.txt > /project/rmcox/LECA/records/make_ms_dirs.sh
while read id; do echo "mkdir /project/rmcox/LECA/ms/cfms/processed/${id}"; done < codes.txt > /project/rmcox/LECA/records/make_ms_dirs.sh
bash /project/rmcox/LECA/records/make_ms_dirs.sh

#### msblender dirs:

# copy paste paths from MS master spreadsheet into vim
# for mzXML:
:%s!^!mkdir -p !
:w! make_mzxml_dirs.sh

# for DB, output, working:
:%s/mzXML/working/
:w! make_working_dirs.sh
:%s/working/output/
:w! make_output_dirs.sh
:%s/output/DB/
:w! make_db_dirs.sh

# execute
cat make_mzxml_dirs.sh | parallel -j8
cat make_output_dirs.sh | parallel -j8
cat make_working_dirs.sh | parallel -j8
cat make_db_dirs.sh | parallel -j8

# symlink data

# problem children 2/22/2021:
../ms/cfms/processed/caeel/beads_iex_1/mzXML
../ms/cfms/processed/caeel/beads_iex_9/mzXML
../ms/cfms/processed/camnx/bng_10/mzXML
../ms/cfms/processed/camnx/bng_1/mzXML
../ms/cfms/processed/camnx/bng_2/mzXML
../ms/cfms/processed/camnx/bng_3/mzXML
../ms/cfms/processed/camnx/bng_4/mzXML
../ms/cfms/processed/camnx/bng_5/mzXML
../ms/cfms/processed/camnx/bng_6/mzXML
../ms/cfms/processed/camnx/bng_7/mzXML
../ms/cfms/processed/camnx/bng_8/mzXML
../ms/cfms/processed/camnx/bng_9/mzXML
../ms/cfms/processed/euggr/sec_2/mzXML
../ms/cfms/processed/human/iex_1/mzXML
../ms/cfms/processed/human/iex_22/mzXML
../ms/cfms/processed/human/iex_23/mzXML
../ms/cfms/processed/xenla/sec_3/mzXML

# problem children 3/1/2021
/project/rmcox/LECA/ms/cfms/processed/mouse/iex_1/mzXML  # fixed
/project/rmcox/LECA/ms/cfms/processed/xenla/sec_3/mzXML  # fixed
/project/rmcox/LECA/ms/cfms/processed/caeel/beads_iex_1/mzXML  # fixed
/project/rmcox/LECA/ms/cfms/processed/human/iex_24/mzXML  # fixed
/project/rmcox/LECA/ms/cfms/processed/human/iex_16/mzXML  # fixed ... I think
/project/rmcox/LECA/ms/cfms/processed/human/iex_1/mzXML  # fixed
/project/rmcox/LECA/ms/cfms/processed/human/iex_23/mzXML  # fixed
/project/rmcox/LECA/ms/cfms/processed/euggr/sec_2/mzXML  # this is fine; haven't run yet

# to count the actualy number of fractions
# inside /project/rmcox/LECA/ms/cfms/processed:
ls -x */*/mzXML/* > fraction_count.txt

#### look for missing data
find ./ -type f -name "Pf.*\.raw" -print # look for file
find /MS/processed/ -iname 'WAN100329_OT2_HS3NE_HCW*mzXML' -type f
find /MS/processed/ -iname 'wgIEF*mzXML' -type f
find /MS/processed/ -iname 'MES_IEX_mixedbed*mzXML' -type f
find /MS/ -iname 'SEC_sputum_buffer1_SEC_02192016*mzXML' -type f

find /MS/ -iname '2012_06_12_SHIMADZUUHPLC*mzXML' -type f
find /MS/ -iname 'PH100503_HS3NE_TCS*' -type f
find /MS/ -iname 'PH081204_HS3NE_AF*' -type f
find /MS/ -iname 'WAN100329_HS3NE_HCW*raw' -type f

# search for a lot of missing data
while read exp; do echo "find /MS/ -iname '${exp}*mzXML' -type f >> loc.txt"; done < missing_exps.txt > find_data.sh
cat find_data.sh | parallel -j8

du -a | cut -d/ -f2 | sort | uniq -c | sort -nr # count number of files in each subdirectory (i.e., # of fractions)

# extract compressed mzXMLs
while read file; do echo "tar -xvf ${file}"; done < tarlist.txt > extract_tar.sh
while read file; do echo "tar -zxvf ${file}"; done < targzlist.txt >> extract_tar.sh
cat extract_tar.sh | parallel -j18

#######################################################################################
# External CFMS data (Jan-Feb 2021)
#######################################################################################

# Hollier et al (Choudhary lab)
# Blue native PAGE fractionations for 3 species of Plasmodium (870 total fractions)
# Publication: https://doi.org/10.1016/j.celrep.2019.07.019
# CEDAR: https://www3.cmbi.umcn.nl/cedar/browse/experiments/CRX20
# ProteomeXchange: https://www.ebi.ac.uk/pride/archive/projects/PXD001220
curl -l  > sample_list.txt
while read sample; do echo "wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/08/PXD009039/${sample}"; done < sample_list.txt > get_data.sh
cat get_data.sh | parallel -j32

# Gerovac et al (Vogel lab)
# Glycerol gradient fractionation of P. aeruginosa (22 total fractions)
# Publication: https://doi.org/10.1128/mBio.03454-20
# ProteomeXchange: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD021359
curl -l ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/01/PXD021359/ > sample_list.txt
while read sample; do echo "wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/01/PXD021359/${sample}"; done < sample_list.txt > get_data.sh
cat get_data.sh | parallel -j10

# Crozier et al (Ferguson lab)
# SEC and IEX fractionation of T. brucei
# Publication: https://doi.org/10.1074/mcp.O117.068122
# ProteomeXchange: ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD005968
curl -l ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD005968/ > sample_list.txt
while read sample; do echo "wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/04/PXD005968/${sample}"; done < sample_list.txt > get_data.sh
cat get_data.sh | parallel -j10

# Berger et al (Brandt lab)
# Blue native PAGE fractionations for Ca. Methanoperedens
# Publication: https://doi.org/10.1016/j.bbabio.2020.148308
# ProteomeXchange: None; Data provided via SURFfilesender (link expired 1 week after sent)
# Raw data transferred to:
/MS/submit/Brandt_BNs_Methanoperedens/
# Processed to mzXML:
/MS/processed/Brandt_BNs_Methanoperedens/

# Shatsky et al (Biggin lab)
# SEC and HIC fractionation of D. vulgaris (5273 fractions)
# Publication: https://doi.org/10.1074/mcp.M115.057117
# ProteomeXchange: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD003392
curl -l ftp://massive.ucsd.edu/MSV000079440/raw/ > sample_list.txt
# these data were in a weird format not easily combined with our CFMS dataset

#######################################################################################
# Downloading, creating, and mapping proteomes
#######################################################################################

# rotifer did not have a reference proteome;
# transcriptome that has been convering to CDS was downloaded:
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/GI/NZ/GINZ01/GINZ01.1.fsa_nt.gz
gunzip GINZ01.1.fsa_nt.gz
grep -c '^>' GINZ01.1.fsa_nt # 21118 CDS

# without pfam and blastp
TransDecoder.Predict -t GINZ01.1.fsa_nt
grep -c '^>' GINZ01.1.fsa_nt.transdecoder.pep # 20895 proteins

# forced inclusion of pfam and blastp significant hits
blastp -query  GINZ01.1.fsa_nt.transdecoder_dir/longest_orfs.pep -db /project/rmcox/databases/uniprot/uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 32 > rotifer.blastp.outfmt6
hmmscan --cpu 32 --domtblout rotifer.pfam.domtblout /project/rmcox/databases/pfam/Pfam-A.hmm GINZ01.1.fsa_nt.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t GINZ01.1.fsa_nt --retain_pfam_hits rotifer.pfam.domtblout --retain_blastp_hits rotifer.blastp.outfmt6
grep -c '^>' GINZ01.1.fsa_nt.transdecoder.pep # 21062 proteins
cp GINZ01.1.fsa_nt.transdecoder.pep ../brart.fasta

# all other proteomes were downloaded as indicated on master spreadsheet:
/project/rmcox/LECA/ms/meta/speciesinfo_02112021.csv

# map to eggnog groups
while read code; do echo "nice -5 python /project/rmcox/programs/eggnog-mapper/emapper.py -i ${code}.fasta --output /project/rmcox/LECA/ms/og_proteomes/nog_mapping/${code}.euNOG.diamond -d euNOG --override --resume --scratch_dir /project/rmcox/LECA/ms/og_proteomes/nog_scratch --no_file_comments --keep_mapping_files -m diamond"; done < codes.txt > map_leca_orgs.sh
cat map_leca_orgs.sh | parallel -j30

# format emapper output for all files
while read code; do echo "Rscript /project/rmcox/LECA/scripts/format_emapper_output.R -f nog_mapping/${code}.euNOG.diamond.emapper.annotations -o nog_mapping/${code}.euNOG.diamond.mapping.1 -s diamond -l 1 -v 2 -p ENOG50"; done < codes.txt > format_emapper_output_1.sh # root

while read code; do echo "Rscript /project/rmcox/LECA/scripts/format_emapper_output.R -f nog_mapping/${code}.euNOG.diamond.emapper.annotations -o nog_mapping/${code}.euNOG.diamond.mapping.2759 -s diamond -l 2759 -v 2 -p ENOG50"; done < codes.txt > format_emapper_output_2759.sh # eukarya

# generate orthocollapsed proteomes for each species
while read code; do echo "python /project/rmcox/LECA/scripts/concat_ortho_proteins.py -m og_proteomes/nog_mapping/${code}.euNOG.diamond.mapping.1 -f proteomes/${code}.fasta -o /project/rmcox/LECA/ms/og_proteomes/${code}.1.collapsed.fasta --limit_length"; done < prok_codes.txt > concat_proks_1.sh

while read code; do echo "python /project/rmcox/LECA/scripts/concat_ortho_proteins.py -m og_proteomes/nog_mapping/${code}.euNOG.diamond.mapping.1 -f proteomes/${code}.fasta -o /project/rmcox/LECA/ms/og_proteomes/${code}.2759.collapsed.fasta --limit_length"; done < euk_codes.txt > concat_euks_2759.sh # CONTAINS ERROR

# count proteins mapped and OGs
for x in *mapping.1; do echo "${x}"; awk -F'\t' '{print $1}' ${x} | tail -n +2 | sort -u | wc -l; done >> num_proteins_mapped_1.txt
for x in *mapping.1; do echo "${x}"; awk -F'\t' '{print $2}' ${x} | tail -n +2 | sort -u | wc -l; done >> num_ogs_mapped_1.txt
for x in *mapping.2759; do echo "${x}"; awk -F'\t' '{print $1}' ${x} | tail -n +2 | sort -u | wc -l; done >> num_proteins_mapped_2759.txt
for x in *mapping.2759; do echo "${x}"; awk -F'\t' '{print $2}' ${x} | tail -n +2 | sort -u | wc -l; done >> num_ogs_mapped_2759.txt


#######################################################################################
# Running all experiments through MSblender
#######################################################################################

# check that all mzXML files have been successfully symlinked
find /project/rmcox/LECA/ms/cfms/processed/ -type d -empty | grep 'mzXML' # only second euglena fractionation missing

# run MSblender on hopper
vim explist.txt  # names of folders containing DB, mzXML, working and output directories

# put contam in all dirs
vim make_db_dirs.sh
:%s/mkdir -p/cp contam.fasta/
:w! cp_contam_fasta.sh
cat cp_contam_fasta.sh | parallel -j8

# concat dbs and contam files
while read code; do echo "cat ${code}*collapsed.fasta ../../records/contam.fasta > contam_combined/${code}.collapsed_contam.combined.fasta"; done < ../codes.txt > combine_contams.sh

# run msblender
bash /project/rmcox/programs/run_msblender/scripts/create_commands.sh explist.txt /project/rmcox/programs/MSblender/runMS2.sh /project/rmcox/programs/SearchGUI-3.3.5/
find /project/rmcox/LECA/ms/cfms/processed/ -iname '*COMMANDS.txt' -type f > msblender_cmds.txt

while read cmd; do sed -i -e 's/^/nice -5 /' $cmd; done < msblender_cmds.txt

while read cmd; do cat $cmd | parallel -j32; done < msblender_cmds.txt  # hmm this didn't really work

while read cmd; do cat $cmd; done < msblender_cmds.txt & | parallel -j8

nohup bash -c 'for cmd in `cat msblender_cmds.txt`; do cat $cmd | parallel -j64; done' > msblender_cmds.out &

# check if all experiments finished
find /project/rmcox/LECA/ms/cfms/processed/ -type d -empty | grep 'output' | wc -l

python /stor/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/msblender2elution.py --prot_count_files /project/rmcox/programs/run_msblender/halobacterium/output/*.group --output_filename halobacterium.prot_count_mFDRpsm001.unique.elut --fraction_name_from_filename --parse_uniprot_id --remove_zero_unique

python /project/cmcwhite/github/protein_complex_maps/protein_complex_maps/preprocessing_util/elutionTidy.py --input_elution /stor/work/Marcotte/project/rmcox/programs/run_msblender/halobacterium/halobacterium.prot_count_mFDRpsm001.unique.elut --outfile halobacterium.unique.protcount.tidy --firstcol ProteinID --valuecol SpectralCounts

#######################################################################################
# Running all experiments through MSblender: TACC Special Edition
#######################################################################################

##### set up directory structure

# species dirs
while read id; do echo "mkdir /scratch/05819/rmcox/cfms/raw/${id}"; done < codes.txt > /scratch/05819/rmcox/cfms/records/make_ms_dirs.sh
while read id; do echo "mkdir /scratch/05819/rmcox/cfms/processed/${id}"; done < codes.txt >> /scratch/05819/rmcox/cfms/records/make_ms_dirs.sh
bash /project/rmcox/LECA/records/make_ms_dirs.sh

# copy paste paths from MS master spreadsheet into vim
# for mzXML:
:%s/\/project\/rmcox\/LECA\/ms/\/scratch\/05819\/rmcox
:%s!^!mkdir -p !
:w! make_mzxml_dirs.sh

# for DB, output, working:
:%s/mzXML/working/
:w! make_working_dirs.sh
:%s/working/output/
:w! make_output_dirs.sh
:%s/output/DB/
:w! make_db_dirs.sh

# make em all
cat make_ms_dirs.sh make_mzxml_dirs.sh make_output_dirs.sh make_working_dirs.sh make_db_dirs.sh > make_all_dirs.sh
bash make_all_dirs.sh

# make contam combined fastas
for x in *fasta; do echo "cat ${x} contam.fasta > contam_combined/${x%.fasta}_contam.combined.fasta"; done > contam_combine.sh
bash contam_combine.sh

# symlink data to DB dirs
# step 1, copy symlink cmds from MS master sheet
# step 2, fix file paths and write
:%s/\/project\/rmcox\/LECA\/ms\/og_proteomes/\/scratch\/05819\/rmcox\/cfms\/proteomes
:%s/\/project\/rmcox\/LECA\/ms/\/scratch\/05819\/rmcox
:w! symlink_dbs.sh

# symlink data to mzXML dirs
# step 1, copy mzxml symlink cmds from MS master sheet
# step 2, fix file paths and write
:%s/\/MS\/processed/\/scratch\/04393\/rctftacc\/rmcox\/MS\/processed
:%s/\/MS\.processed/\/scratch\/04393\/rctftacc\/rmcox\/MS\/processed
:%s/\/project\/rmcox\/LECA\/ms/\/scratch\/05819\/rmcox
:w! symlink_mzxmls.sh

# NOTE: make explist.txt file by copying file paths from last column of MS master sheet
# generate commands and put them all in one file
bash /work/05819/rmcox/lonestar/run_msblender/scripts/create_commands.sh exp_list.txt /work/projects/BioITeam/ls5/gitdir/MSBlender/MSblender_consistent/runMS2.sh /work/05819/rmcox/lonestar/SearchGUI-3.3.5
find /scratch/05819/rmcox/cfms/processed/ -iname '*COMMANDS.txt' -type f > msblender_cmds.txt

# Anna's new "create commands" on TACC
# ensures working dir is always printed
echo "bash $runscript $mzxml_file $proteome $working_dir $output_dir $searchgui_dir > $working_dir/$pfx.msblender.log 2>&1" >> $EXP_COMMANDS

bash /work/projects/BioITeam/ls5/gitdir/MSBlender/run_msblender/scripts/create_commands.sh exp_list.txt /work/projects/BioITeam/ls5/gitdir/MSBlender/MSblender_consistent/runMS2.sh /work/05819/rmcox/lonestar/SearchGUI-3.3.5
find /scratch/05819/rmcox/cfms/processed/ -iname '*COMMANDS.txt' -type f > msblender_cmds.txt

##### set up commands for parallel processing

# all fractions
for exp in `cat msblender_cmds.txt`; do echo "`cat $exp`"; done >> all_ms_exps.sh

# how many?
wc -l all_ms_exps.sh # 12692

# split into sets of 40
split -l 40 ../all_ms_exps.sh --numeric-suffixes cfms --additional-suffix=.cmds

# also, Claire gave me a template file for MPI
# can I loop through, edit the template for each batch, and save it as a new file?
for batch in ms_exp*; do echo "sed 's/job1/${batch}/g; s/commands\.txt/\/scratch\/05819\/rmcox\/cfms\/${batch}/g' ls5_template.slurm > ${batch}.slurm"; done > make_slurm_files.sh

# other things:
# how to keep execution loop going when disconnected from LS5? ------> tmux or screen
# how to track completions/failures?

# MSblender dir:
/work/projects/BioITeam/ls5/gitdir/MSBlender/MSblender_consistent/runMS2.sh

# try Anna's launcher script:
/work/projects/BioITeam/common/script/launcher_maker.py -t 2 -n ms_exp_batch9227 -w 8 -a A-cm10 -q dev -e rachaelcox@utexas.edu
sbatch ms_exp_batch9227.slurm

for exp in cfms*cmds; do echo "/work/projects/BioITeam/common/script/launcher_maker.py -t 2 -n ${exp} -w 8 -a A-cm10 -q dev -e rachaelcox@utexas.edu"; done > make_slurm_files.sh

# ANNA EXAMPLE DIRECTORY:/project/abattenh/rmcox_msblender/

# \\\ TO DO LIST
# 1. update run_MS2 path to Anna's edited script in BioITeam ----- DONE
# 2. fix all proteome symlinks ----- DONE
# 3. create commands with Anna's edited script, or manually add on a msblender.log redirect to each experiment ----- DONEls

# problem with Anna's create commmands script
# this is what it does now:
bash /work/projects/BioITeam/ls5/gitdir/MSBlender/MSblender_consistent/runMS2.sh /scratch/05819/rmcox/cfms/processed/arath/wwc_1/mzXML/Athaliana_sproutsWWC_43a_02192017.mzXML /scratch/05819/rmcox/cfms/processed/arath/wwc_1/DB/arath.2759.collapsed_contam.combined.fasta
/scratch/05819/rmcox/cfms/processed/arath/wwc_1/DB/arath.collapsed_contam.combined.fasta /scratch/05819/rmcox/cfms/processed/arath/wwc_1/working /scratch/05819/rmcox/cfms/processed/arath/wwc_1/output /work/05819/rmcox/lonestar/SearchGUI-3.3.5

# need it to do this:
bash /work/projects/BioITeam/ls5/gitdir/MSBlender/MSblender_consistent/runMS2.sh /scratch/05819/rmcox/cfms/processed/arath/wwc_1/mzXML/Athaliana_sproutsWWC_43a_02192017.mzXML /scratch/05819/rmcox/cfms/processed/arath/wwc_1/DB/arath.2759.collapsed_contam.combined.fasta /scratch/05819/rmcox/cfms/processed/arath/wwc_1/working /scratch/05819/rmcox/cfms/processed/arath/wwc_1/output /work/05819/rmcox/lonestar/SearchGUI-3.3.5

# ah the problem is not with Anna's script. it's that I have multiple proteomes in the DB directories
# clean it up
find /scratch/05819/rmcox/cfms/processed/ -iname '*contam.combined.fasta' -print -exec rm '{}' ';'

# relink
bash /scratch/05819/rmcox/cfms/records/symlink_dbs.sh

# regenerate commands files
bash /work/projects/BioITeam/ls5/gitdir/MSBlender/run_msblender/scripts/create_commands.sh exp_list.txt /work/projects/BioITeam/ls5/gitdir/MSBlender/MSblender_consistent/runMS2.sh /work/05819/rmcox/lonestar/SearchGUI-3.3.5
find /scratch/05819/rmcox/cfms/processed/ -iname '*COMMANDS.txt' -type f > msblender_cmds.txt
touch all_ms_exps.sh | for exp in `cat msblender_cmds.txt`; do echo "`cat $exp`"; done >> all_ms_exps.sh
# make sure it's right this time
head all_ms_exps.sh
# preview:
bash /work/projects/BioITeam/ls5/gitdir/MSBlender/MSblender_consistent/runMS2.sh /scratch/05819/rmcox/cfms/processed/orysj/iex_1/mzXML/RiceL_IEX_01-1a_20150611.mzXML /scratch/05819/rmcox/cfms/processed/orysj/iex_1/DB/orysj.2759.collapsed_contam.combined.fasta /scratch/05819/rmcox/cfms/processed/orysj/iex_1/working /scratch/05819/rmcox/cfms/processed/orysj/iex_1/output /work/05819/rmcox/lonestar/SearchGUI-3.3.5 > /scratch/05819/rmcox/cfms/processed/orysj/iex_1/working/RiceL_IEX_01-1a_20150611.msblender.log 2>&1 # looks good
# /scratch/05819/rmcox/cfms$ cd batched_cmds/
split -l 40 ../all_ms_exps.sh --numeric-suffixes cfms --additional-suffix=.cmds
for exp in cfms*cmds; do echo "/work/projects/BioITeam/common/script/launcher_maker.py -t 2 -n ${exp} -w 4 -a A-cm10 -q dev -e rachaelcox@utexas.edu"; done > make_slurm_files.sh
bash make_slurm_files.sh

# test one
sbatch cfms9227.slurm 

# doesn't work. let's troubleshoot in an idev session
idev -p development -m 60 -A A-cm10

# the problem is here:
$LAUNCHER_DIR/paramrun

# trying to test in idev...
export LAUNCHER_PPN=40
export LAUNCHER_NHOSTS=5
export LAUNCHER_RMI_HOSTFILE=`mktemp -t launcher.cfms00.hostlist.XXXXXXXX`
scontrol show hostname $SLURM_NODELIST > $LAUNCHER_RMI_HOSTFILE

# quick reference for killing this monster:
kill `ps -aux | grep 'rmcox' | awk '{print $2}'`

# some commands are failing:
for i in `cat failed_cmds.txt`; do printf ">>>>>>>>FILE:${i}<<<<<<<<\n`tail ${i}`\n\n"; done > failed_reasons.txt
for i in `cat failed_cmds.txt`; do echo "`tail ${i}`"; done > failed_reasons.txt

# make sbatch loop
mkdir test_loop
cd test_loop
split -l 20 ../all_ms_exps.sh --numeric-suffixes cfms --additional-suffix=.cmds
rm cfms[1-9]* # leave only 10 files to test with
for exp in cfms*cmds; do echo "/work/projects/BioITeam/common/script/launcher_maker.py -t 1 -n ${exp} -w 8 -a A-cm10 -q dev -e rachaelcox@utexas.edu"; done > make_slurm_files.sh
bash make_slurm_files.sh
ls *slurm > all_slurm_files.txt

vim serial_submit.sh
#!/bin/bash
runlist=$1
userid=$2
finished_dir=$3
while read slrm; do
        echo -e "\nSubmitting ${slrm} at $(date)."
        sbatch $slrm
        running=$(showq -u ${userid} | grep 'Running\|Waiting')
        until [ -z "$running" ]; do
                echo "Waiting on previous run(s) to finish ..."
                echo "Currently running or in queue: ${running}"
                sleep 60s
                running=$(showq -u ${userid} | grep 'Running\|Waiting')
        done
        files=${slrm%.*}
        echo "Moving ${files} to ${finished_dir}."
        mv $files* $finished_dir
        echo "No running or waiting jobs detected."
done < $runlist
echo "All jobs submitted."

bash serial_submit.sh all_slurm_files.txt rmcox finished

# divide up cmd files for me and Anna
find . -maxdepth 1 -type f -iname "*cmds" | tail -316 | xargs cp -t ../anna_fractions
find . -maxdepth 1 -type f -iname "*cmds" | head -319 | xargs cp -t ../rachael_fractions

# make my slurm files
for exp in cfms*cmds; do echo "/work/projects/BioITeam/common/script/launcher_maker.py -t 1 -n ${exp} -w 4 -a A-cm10 -q dev -e rachaelcox@utexas.edu"; done > make_slurm_files.sh
bash make_slurm_files.sh
ls *slurm > all_slurm_files.txt
bash ../scripts/serial_submit.sh all_slurm_files.txt rmcox finished
bash ../scripts/serial_submit.sh all_slurm_files.txt rmcox finished > serialsubmit.log 2>&1

# code for organizing logs by error code
rmcox/cfms/records

# anna has edited the serial submission script in
rmcox/cfms/anna_fractions

# failed to submit:
grep -A5 'Too many simultaneous jobs in queue.' serialsubmit.log | grep 'Moving cfms*' # about 1700 have completed, 9 have failed to submit as of 8:30AM 4/20/21 (24 hours of run time)

# make list of files tp rerun
grep -A5 'Too many simultaneous jobs in queue.' serialsubmit_reruns.log | grep 'Moving cfms*' > rerun_slurms.txt

# move them out of finished
while read files; do mv finished/${files}* .; done < rerun_slurms.txt
ls *slurm > all_slurm_reruns.txt
bash ../scripts/serial_submit.sh all_slurm_reruns.txt rmcox finished > serialsubmit_reruns2.log 2>&1

# Anna's LS5 scratch dir:
/scratch/01063/abattenh/cfms

# how many runs had problems?
ls5:/scratch/05819/rmcox/cfms/rachael_fractions/logs$ find . -type f | cut -d/ -f2 | sort | uniq -c

      1 comet_xml  # just wgSEC_26 w/ syntax error parsing XML, tried rerunning, no luck; claire doesn't have this fraction in her final matrix
   5224 good  # successes
     98 msb_in  # these are all 'division by zero' or 'list index out of range' errors (make-msblender-in.py)
    490 msgf_plus  # one of these is wgSEC_14, which is in claire's final matrix
    550 no_prots  # no proteins detected
     10 tandemk  # these are all errors with input spectra (low quality?)

# Anna has modified runMS2 such that it will proceed if 1 of the search aglorithms fails and successfully rerun most of the fractions with errors
# Her notes on this are here (including some **really nice** bash scripting examples for handling complicated directory structures and generating summary info):
/project/abattenh/rmcox_msblender/msblender_local_reruns.sh

# Run status of all 12,691 fractions is here:
/project/rmcox/LECA/ms/cfms/meta/all_reruns.info.txt

# Run statuses are summarized here:
/project/rmcox/LECA/ms/cfms/meta/all_stats_by_result.txt
# --- good w/ 3 search alg: 10,550
# --- good w/ 2 search algs: 1,012
# --- no proteins detected: 1,081 (almost all species have at least 1)
# --- no viable spectra: 26 (19 human + 4 arath + 2 nemve + 1 wheat)
# --- parse xml error: 22 (21 plaf7 + 1 wheat)

#######################################################################################
# Finding all MSblender results
#######################################################################################

find processed/ -type f -iname *.group > elut_files/groupfile_paths.txt  # 11,801 cfms fractions, but only 11,562 "good fractions"

# edit so we can glob for experiments
vim groupfile_paths.txt
:%s/output.*/output\/\*\.group/  # hmm. this caught 10,716 fractions, with 1,085 in weird tmp subdirectories

# get the problem files into a .txt, all problem files are in /working/ subdirectories:
grep '/output/' groupfile_paths.txt | wc -l # 10716 fractions
grep -v '/output/' groupfile_paths.txt | wc -l # 1085 fractions
grep -v '/output/' groupfile_paths.txt > problem_fractions.txt

# edit to get problem fraction names
vim problem_fractions.txt
:%s/.*working\///
:%s/\.temp.*//
:w problem_fraction_names.txt

# after further investigation almost all .groups in tmp dirs are the fractions with no prots detected, except for 4:
#   /project/rmcox/LECA/ms/cfms/processed/pseae/sec_4/, Ps_sputum_SEC_buffer2_membrane_tritonX_02102016_7
#   /project/rmcox/LECA/ms/cfms/processed/orysj/iex_2/, Rice_leaf_IEX_fraction_42a_20150518
#   /project/rmcox/LECA/ms/cfms/processed/orysj/iex_2/, Rice_leaf_IEX_fraction_45a_20150518
#   /project/rmcox/LECA/ms/cfms/processed/orysj/iex_2/, Rice_leaf_IEX_fraction_59a_20150521

# !!!
# 11562 - 10720 = 842
# where are the missing 842 "good" fractions??
# !!!

# make a file of good fraction IDs to compare to summary list
grep '/output/' groupfile_paths.txt > good_fractions.txt

# initial investigation of CAMNX missing frac "20190816_Methanop_6-12_S1_SB_ACO29" shows that it did indeed run Comet and TandemK, but errored out on MSGF+. No .group file was written, according to this:
/project/rmcox/LECA/ms/cfms/processed/camnx/bng_9/working/20190816_Methanop_6-12_S1_SB_ACO29.msblender.log

# !!! MISSING !!! 
# 323 CAMNX fractions, 511 TRYB2 fractions, 17 PSEAE fractions are all labeled "good_2" but have no .group output
# !!!		   !!!

# 96 .groups were written in CAMNX output dirs, all are "good_2", and all are labeled "rerun.2"
# all 853 missing fractions are labeled "rerun.1"
# something related to rerun.1? check:
/project/rmcox/LECA/ms/cfms/rerun/rerun.1/logs/good_2/

# after looking into the fraction logs in /project/rmcox/LECA/ms/cfms/rerun/rerun.1/logs, I can see that for the missing fraction "20190816_Methanop_6-12_S1_SB_ACO29" that the output directory was originally:
/stor/work/Marcotte/project/rmcox/LECA/ms/cfms/reprocess/camnx/bng_9/output/

# the log indicates that a .group file was written for this fraction
# as far as I can tell, this dir no longer exists. did Anna move the output from here?
# found a reprocess directory here:
/project/rmcox/LECA/ms/cfms/tmp/reprocess/

# can we find .groups here?
find . -type f -iname *.group > potential_missing_fracs.txt  # 253 total
find camnx/ -type f -iname *.group > potential_missing_camnx.txt  # 96 camnx, but these are all accounted for

# let's just search all of hopper for these missing files
# made a file called missing_msb_output.txt in R and wrote it to /project/rmcox/LECA/ms
while read exp; do echo "find /project/rmcox/ -iname '${exp}*group' -type f >> locate_output.txt"; done < missing_msb_output.txt > find_data.sh
touch locate_output.txt
cat find_data.sh | parallel -j16  # found 30

# no dice. but I did find this in Anna's notes:
rsync -avrP --exclude=tacc.logs \
  abattenh@hfogstor01:/project/rmcox/LECA/ms/cfms/rerun/ \
  /project/rmcox/LECA/ms/cfms/rerun/rerun.1/hopefog/

# so rerun.1 is on hopefog
ssh rmcox@hfogcomp01.ccbb.utexas.edu
cdw
cd /project/rmcox/LECA/ms/cfms/
find . -type f -iname *.group  # well tryb2 fracs are here
find . -type f -iname *.group | grep 'tryb2' | wc -l  # 458 tryb2
find . -type f -iname *.group | grep 'camnx' | wc -l  # 0 camnx
find . -type f -iname *.group | grep 'pseae' | wc -l  # 0 pseae

# ok, Anna has rerun the missing fractions
# let's make sure we can find them all
find processed/ -type f -iname *.group > elut_files/groupfile_paths_060721.txt

# edit so we can glob for experiments
vim groupfile_paths.txt
:%s/output.*/output\/\*\.group/  # 11,569 "good" fractions

#######################################################################################
# Formatting all MSblender results
#######################################################################################

# template to make .eluts
python /stor/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/msblender2elution.py --prot_count_files /project/rmcox/programs/run_msblender/halobacterium/output/*.group --output_filename elut_files/halobacterium.prot_count_mFDRpsm001.unique.elut --fraction_name_from_filename --parse_uniprot_id --remove_zero_unique

# template to make .tidy.eluts
python /project/cmcwhite/github/protein_complex_maps/protein_complex_maps/preprocessing_util/elutionTidy.py --input_elution /stor/work/Marcotte/project/rmcox/programs/run_msblender/halobacterium/halobacterium.prot_count_mFDRpsm001.unique.elut --outfile halobacterium.unique.protcount.tidy --firstcol ProteinID --valuecol SpectralCounts

# template to count unique protein IDs
awk -F',' '{print $2}' halobacterium.unique.protcount.tidy | sort -u | wc -l  # 2337 unique protein IDs


# make commands for all experiments
zsh
code=(`awk -F'/' '{print $2}' /project/rmcox/LECA/ms/cfms/explist.txt`) \
exp=(`awk -F'/' '{print $3}' /project/rmcox/LECA/ms/cfms/explist.txt`)
for i j (${code:^exp}) echo "i: $i, j: $j"  # I think it works

for i j (${code:^exp}) echo "python /stor/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/msblender2elution.py --prot_count_files /project/rmcox/LECA/ms/cfms/processed/${i}/${j}/output/*.group --output_filename /project/rmcox/LECA/ms/cfms/elut_files/${i}.${j}.prot_count_mFDRpsm001.unique.elut --fraction_name_from_filename --parse_uniprot_id --remove_zero_unique" >> group_to_elut_cmds.sh # nevermind this does not work, the array order ends up being semi-random

# finally just ended up doing it with concatenate in the google doc
# example to make sure that's actually working:
python /stor/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/msblender2elution.py --prot_count_files /project/rmcox/LECA/ms/cfms/processed/arath/iex_1/output/*.group --output_filename /project/rmcox/LECA/ms/cfms/elut_files/arath.iex_1.prot_count_mFDRpsm001.unique.elut --fraction_name_from_filename --parse_uniprot_id --remove_zero_unique

# ok, think we're good to go
cat group_to_elut_cmds.sh | parallel -j16 2>&1 | tee group_to_elut.log

# see if any didn't work
grep 'IOError' group_to_elut.log

# these 2 didn't work
IOError: File /project/rmcox/LECA/ms/cfms/processed/camnx/bng_10/output/*.group does not exist # no prots for whole experiment
IOError: File /project/rmcox/LECA/ms/cfms/processed/human/iex_17/output/*.group does not exist # no mzXMLs... symlink didn't work

# i think it's fine, we have human data out the wazoo but keep human/iex_17 in mind the next time i do a batch of reruns

# let's make it tidy so i can actually do something with this

# template
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THIS NO LONGER EXISTS
python /project/cmcwhite/github/protein_complex_maps/protein_complex_maps/preprocessing_util/elutionTidy.py --input_elution /stor/work/Marcotte/project/rmcox/programs/run_msblender/halobacterium/${code}.${exp}.prot_count_mFDRpsm001.unique.elut --outfile ${code}.${exp}.unique.protcount.tidy --firstcol Orthogroup --valuecol SpectralCounts
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THIS DOES NOT WORK
python /project/rmcox/github_repos/protein_complex_maps/protein_complex_maps/util --input_elutions /stor/work/Marcotte/project/rmcox/programs/run_msblender/halobacterium/${code}.${exp}.prot_count_mFDRpsm001.unique.elut --outfile ${code}.${exp}.unique.protcount.tidy --sep '\t'
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

# made my own script here for processing to tidy format:
/project/rmcox/LECA/scripts/elut2tidy.R

# test
Rscript /project/rmcox/LECA/scripts/elut2tidy.R --elut_file arath.iex_1.prot_count_mFDRpsm001.unique.elut --out_file arath.iex_1.prot_count_mFDRpsm001.unique.elut.tidy # it works

# convert all .elut files to tidy format
rm elut_to_tidy_cmds.sh
rm elut_to_tidy.log
rm *unique.tidy
while IFS= read line; do
  line=`echo $line | perl -pe '~s|\r?\n||' `
  if [[ "$line" == "" ]]; then continue; fi
  code=`echo "$line" | awk '{print $1}'`
  exp=`echo "$line" | awk '{print $2}'`
  echo "$code $exp"
  rootDir="/stor/work/Marcotte/project/rmcox"
  scriptDir="LECA/scripts"
  elutionDir="LECA/ms/cfms2/elut_files"
  infile="${code}.${exp}.prot_count_mFDRpsm001.unique.filtdollo.norm.elut"
  outfile="${code}.${exp}.prot_count_mFDRpsm001.unique.filtdollo.norm.tidy"
  cmd="Rscript $rootDir/$scriptDir/elut2tidy.R"
  cmd="$cmd --elut_file $rootDir/$elutionDir/$infile"
  cmd="$cmd --out_file $rootDir/$elutionDir/$outfile"
  cmd="$cmd 2>&1 | tee elutnorm_to_tidy.log"
  echo "$cmd" | tee -a elutnorm_to_tidy_cmds.sh
done < <(tail -n +2 explist.txt)

cat elutnorm_to_tidy_cmds.sh | parallel -j16

# cluster results
ls *unique.elut > all_eluts.txt

# let's just try 1
python /stor/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/elutionCluster.py \
--input_elutions /project/rmcox/LECA/ms/cfms/elut_files/tetts.sec_3.prot_count_mFDRpsm001.unique.elut \
--outfile tetts.sec_3.unique.protcount.hclust.xlsx --hclust_metric pearson \
--hclust_method average --threshold 10 \
--path_to_annotations /project/rmcox/LECA/ms/cfms/annotations/eunog_refs.human_genes.tsv # this works

# try with all -- does not work
python /stor/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/elutionCluster.py \
--input_elutions all_eluts.csv \
--outfile ../analysis/leca.unique.protcount.hclust.xlsx --hclust_metric pearson \
--hclust_method average --threshold 10 \
--path_to_annotations /project/rmcox/LECA/ms/cfms/annotations/eunog_refs.human_genes.tsv

# this works with tail 2, it's pulling 1 xenla and 1 yeast elut, maybe the problem is the prefix
python /stor/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/elutionCluster.py \
--input_elutions `cat all_eluts.txt | tail -2` \
--outfile ../analysis/leca.unique.protcount.hclust.xlsx --hclust_metric pearson \
--hclust_method average --threshold 10 \
--path_to_annotations /project/rmcox/LECA/ms/cfms/annotations/eunog_refs.human_genes.tsv


grep(!"halsa|camnx|pseae|ecoli|staa8", elut_files, value = TRUE)

# clustering
hclust_elut(concat_elut_file, outfile)

hclust_elut("/project/rmcox/LECA/ms/cfms/results/leca_15exp_test_elut.filt20p.raw.txt", 
	"/project/rmcox/LECA/ms/cfms/results/leca_15exp_test_elut.filt20p.hclust.txt") #test
# start at 11:50am, finished around 1:15pm

hclust_elut("/project/rmcox/LECA/ms/cfms/results/leca_70exp_test_elut.filt10p.raw.txt", 
	"/project/rmcox/LECA/ms/cfms/results/leca_70exp_test_elut.filt10p.hclust.txt") # started around 1:30pm 6/13/21

hclust_elut("/project/rmcox/LECA/ms/cfms/results/leca_euk_elut.filt10p.raw.txt", 
	"/project/rmcox/LECA/ms/cfms/results/leca_euk_elut.filt10p.hclust.txt") # started around 1:30pm 6/13/21


# need to filter before clustering
# at least 100 PSMs observed across all species

#######################################################################################
# Rerunning all eukaryotic MSblender experiments
#######################################################################################

# well it turns out that the search proteomes are not euk-mapped but in fact are root-mapped due to this error:
while read code; do echo "python /project/rmcox/LECA/scripts/concat_ortho_proteins.py -m og_proteomes/nog_mapping/${code}.euNOG.diamond.mapping.1 -f proteomes/${code}.fasta -o /project/rmcox/LECA/ms/og_proteomes/${code}.2759.collapsed.fasta --limit_length"; done < euk_codes.txt > concat_euks_2759.sh

# this is the correct version
while read code; do echo "python /project/rmcox/LECA/scripts/concat_ortho_proteins.py -m og_proteomes/nog_mapping/${code}.euNOG.diamond.mapping.2759 -f proteomes/${code}.fasta -o /project/rmcox/LECA/ms/og_proteomes/${code}.2759.collapsed.fasta --limit_length"; done < euk_codes.txt > concat_euks_2759.sh

cat concat_euks_2759.sh | parallel -j16


##### move rootNOG MSblender results to new location
cp -r processed/ processed_root/
cp -r elut_files/ elut_files_root/

##### reconstruct MSblender directory tree for all euk experiments w/ euNOG proteomes

bash /work/projects/BioITeam/ls5/gitdir/MSBlender/run_msblender/scripts/create_commands.sh exp_list.txt /work/projects/BioITeam/ls5/gitdir/MSBlender/MSblender_consistent/runMS2.sh /work/05819/rmcox/lonestar/SearchGUI-3.3.5
find /scratch/05819/rmcox/cfms/processed/ -iname '*COMMANDS.txt' -type f > msblender_cmds.txt

##### set up commands for parallel processing

# all fractions
for exp in `cat msblender_cmds.txt`; do echo "`cat $exp`"; done >> all_ms_exps.sh

# how many?
wc -l all_ms_exps.sh # 12692

# split into sets of 40
split -l 40 ../all_ms_exps.sh --numeric-suffixes cfms --additional-suffix=.cmds

##### rsync reconfigured tree to Hopefog to divide workload in 1/2
rsync -avrWL processed/ \
 rmcox@hfogcomp01.ccbb.utexas.edu:/project/rmcox/LECA/ms/cfms/processed/

rsync -avrWL processed/ \
 rmcox@hfogcomp02.ccbb.utexas.edu:/project/rmcox/LECA/ms/cfms/processed/

## JOB LIMITS
# REMEMBER TO ADJUST OUTPUT DIRS (for each computer)
 # hopper = 4
 # hfog1 = 3
 # hfog2 = 4
 # hfog3 = 4

# check for CCT and 20S proteasome counts
cd /project/rmcox/LECA/ms/cfms2/qc

grep -Ff cct_ogs.txt ../process_leca/*/*/output/*group | awk -F' ' '{print $1}' > tmp.txt
sed 's/\//\t/g; s/:/\t/' tmp.txt | cut -f 2- > cct_obs_allspec.txt
rm tmp.txt

grep -Ff proteasome_ogs.txt ../process_leca/*/*/output/*group | awk -F' ' '{print $1}' > tmp.txt
sed 's/\//\t/g; s/:/\t/' tmp.txt | cut -f 2- > proteasome_obs_allspec.txt
rm tmp.txt

#######################################################################################
# Mapping rootNOGs to euNOGs (from Jaime)

#######################################################################################
wget http://eggnog5.embl.de/download/eggnog_5.0/raw_data/groups_parents_enog.tar.gz

# check degree of rootNOG to eukNOG redundancy using LGL
nohup /stor/work/Marcotte/project/rmcox/programs/LGL/bin/lgl.pl -c conf_eunog_map &> nohup_conf_file_eunogs.out&

# keep getting an error; LGL does not like things that are circular, e.g.
level_2759 level_1
KOG3118 KOG3117
KOG3117 KOG3118

# had to remove all instances of a euNOG in the rootNOG column, now it works

#######################################################################################
# Clustering elution matrix
#######################################################################################

# hierarchically cluster filtered matrix (sum of PSMs >= 300)
python /project/rmcox/LECA/scripts/elutionCluster.py --input_elutions leca_euks_elut.filt150p.raw.txt --outfile leca_euks_elut.filt150p.hclust.csv --hclust_metric pearson --hclust_method average

sort -k3 -n leca_euks_elut.filt150p.hclust.csv > leca_euks_elut.filt150p.hclust.sorted.csv
cut --complement -d',' -f2,3 leca_euks_elut.filt150p.hclust.sorted.csv > leca_euks_elut.filt150p.hclust.nosums.csv

# hierarchically cluster filtered matrix (filtered for LECA genes & PSMs >= 150)
python /project/rmcox/LECA/scripts/elutionCluster.py --input_elutions leca_euks_elut.filtdollo.filt150p.raw.txt --outfile leca_euks_elut.filtdollo.filt150p.hclust.csv --hclust_metric pearson --hclust_method average
sort -k3 -n leca_euks_elut.filtdollo.filt150p.hclust.csv > leca_euks_elut.filtdollo.filt150p.hclust.sorted.csv

python /project/rmcox/LECA/scripts/elutionCluster.py --input_elutions leca_euks_elut.filtdollo.raw.txt --outfile leca_euks_elut.filtdollo.hclust.csv --hclust_metric pearson --hclust_method average

#######################################################################################
# Gold standard PPIs from Complex Portal
#######################################################################################
curl -l ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/ > gold_ids.txt
while read sample; do echo "wget http://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/${sample}"; done < gold_ids.txt > get_data.sh
cat get_data.sh | parallel -j32

# retrieve proteomes from uniprot
while read id; do echo "perl /project/rmcox/LECA/scripts/retrieve_proteomes.pl ${id}"; done < gold_ids.txt > retrieve_proteomes.sh

# map to eggnog groups
while read id; do echo "nice -5 python /project/rmcox/programs/eggnog-mapper/emapper.py -i ${id}.fasta --output /project/rmcox/LECA/ms/gold_stds/proteomes/nog_mapping/${id}.euNOG.diamond --override --resume --scratch_dir /project/rmcox/LECA/ms/gold_stds/proteomes/nog_scratch/ --no_file_comments --keep_mapping_files -m diamond"; done < ../gold_ids.txt > map_gold_stds.sh
cat map_gold_stds.sh | parallel -j30

# format emapper output for all files
while read code; do echo "Rscript /project/rmcox/LECA/scripts/format_emapper_output.R -f nog_mapping/${code}.euNOG.diamond.emapper.annotations -o nog_mapping/${code}.euNOG.diamond.mapping.2759 -s diamond -l 2759 -v 2 -p ENOG50"; done < ../gold_ids.txt > format_emapper_output_2759.sh # eunog

# format protein IDs into file for input to Claire's mapping script
Rscript /project/rmcox/LECA/scripts/format_comportal.R -d /project/rmcox/LECA/ms/gold_stds/comportal -p '*tsv' --write_annotations TRUE

# run Claire's complex mapping script
python /project/cmcwhite/data/gold_standards/corum_complexes/translate_corum2.py --complexes /project/rmcox/LECA/ms/gold_stds/comportal/9606.gold.cmplx.txt --mapping /project/rmcox/LECA/ms/gold_stds/ac_mappings/9606.euNOG.diamond.mapping.2759 --outfile 9606.euNOG.gold.cmplx.txt # test # nope didn't work

# need to re-format emapper output such that: (1) no header, (2) space-separated
while read code; do echo "Rscript /project/rmcox/LECA/scripts/format_uniprot_ids.R -f /project/rmcox/LECA/ms/gold_stds/proteomes/nog_mapping/${code}.euNOG.diamond.mapping.2759"; done < ../gold_ids.txt > format_uniprot_ids.sh
cat format_uniprot_ids.sh | parallel -j12

while read code; do echo "tail -n +2 /project/rmcox/LECA/ms/gold_stds/proteomes/nog_mapping/${code}.euNOG.diamond.mapping.2759 | sed 's/\t/ /g' > ac_mappings/${code}.euNOG.diamond.mapping.2759.fmt"; done < ../gold_ids.txt > format_ac_mapping.sh

# now let's try Claire's script
python /project/cmcwhite/data/gold_standards/corum_complexes/translate_corum2.py --complexes /project/rmcox/LECA/ms/gold_stds/ppis/9606.gold.cmplx.txt --mapping /project/rmcox/LECA/ms/gold_stds/ac_mappings/9606.euNOG.diamond.mapping.2759.fmt --outfile 9606.euNOG.gold.cmplx.txt # yes

while read code; do echo "python /project/cmcwhite/data/gold_standards/corum_complexes/translate_corum2.py --complexes /project/rmcox/LECA/ms/gold_stds/ppis/${code}.gold.cmplx.txt --mapping /project/rmcox/LECA/ms/gold_stds/ac_mappings/${code}.euNOG.diamond.mapping.2759.fmt --outfile ${code}.euNOG.gold.cmplx.txt"; done < ../gold_ids.txt > map_gold_ppis.sh

# some of these didn't work, potentially because: (1) prokaryotic and/or (2) file contains only 1 gold std complex
find -size 0
./9615.euNOG.gold.cmplx.txt # Canis lupus familiaris, 1 complex
./562.euNOG.gold.cmplx.txt # Escherichia coli, prokaryote
./7788.euNOG.gold.cmplx.txt # Torpedo marmorata, 1 complex
./7787.euNOG.gold.cmplx.txt # Tetronarce californica, 1 complex
./83333.euNOG.gold.cmplx.txt # Escherichia coli (strain K12), prokaryote
./9940.euNOG.gold.cmplx.txt # Ovis aries, 1 complex 
./6523.euNOG.gold.cmplx.txt # Lymnaea stagnalis, 1 complex

# get unique list of all gold standard euNOG identifiers
# make sure euNOG complex file isn't empty, then append it to file all.euNOG.cmplx.txt
rm all.euNOG.gold.cmplx.txt
touch all.euNOG.gold.cmplx.txt
for file in *[[:digit:]]*.euNOG.gold.cmplx.txt; do
	if [ -s $file ]; then cat ${file} >> all.euNOG.gold.cmplx.txt; fi
done

sed 's/ /\n/g' all.euNOG.gold.cmplx.txt | sort | uniq | wc -l # 2599 unique gold standard orthogroups

sed 's/ /\n/g' all.euNOG.gold.cmplx.txt | sort -u > all.euNOG.gold_prots.unique.txt

# let's add a few more databases:
# CORUM (mammals) + TAIR (arabidopsis) + Claire's 6 manually curated complexes (arabidopsis)
# previously acquired ARATH data in /project/rmcox/LECA/ms/gold_stds/:
./proteomes/3702.fasta
./nog_mapping/3702.euNOG.diamond.mapping.2759
./ac_mappings/3702.euNOG.diamond.mapping.2759.fmt
./comportal/3702.annotations.gold.cmplx.txt

# make eggNOG file for CORUM complexes
# 10090 = mouse, 10116 = rat, 9606 = human, 9913 = bovine, 9823 = pig
touch corum.euNOG.diamond.mapping.2759.fmt
while read x; do cat ${x}.euNOG.diamond.mapping.2759.fmt >> corum.euNOG.diamond.mapping.2759.fmt; done < ../corum_codes.txt

# run claire's script on the cmplx files
python /project/cmcwhite/data/gold_standards/corum_complexes/translate_corum2.py --complexes /project/rmcox/LECA/ms/gold_stds/corum/corum.cmplx.txt --mapping /project/rmcox/LECA/ms/gold_stds/ac_mappings/corum.euNOG.diamond.mapping.2759.fmt --outfile /project/rmcox/LECA/ms/gold_stds/corum/corum.euNOG.gold.cmplx.txt

python /project/cmcwhite/data/gold_standards/corum_complexes/translate_corum2.py --complexes /project/rmcox/LECA/ms/gold_stds/tair/tair.cmplx.txt --mapping /project/rmcox/LECA/ms/gold_stds/ac_mappings/3702.euNOG.diamond.mapping.2759.fmt --outfile /project/rmcox/LECA/ms/gold_stds/tair/tair.euNOG.gold.cmplx.txt

python /project/cmcwhite/data/gold_standards/corum_complexes/translate_corum2.py --complexes /project/rmcox/LECA/ms/gold_stds/claire/claire.cmplx.txt --mapping /project/rmcox/LECA/ms/gold_stds/ac_mappings/3702.euNOG.diamond.mapping.2759.fmt --outfile /project/rmcox/LECA/ms/gold_stds/claire/claire.euNOG.gold.cmplx.txt

# merge based on jaccard coefficient
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

# my python script isn't working
# something fundamentally wrong with the recursion
# first fix weird formatting issues, e.g. quotes in random places
sed -i 's/"//g' all.gold.cmplx.txt

# let's try kevin's script before we go any further
python /project/rmcox/github_repos/protein_complex_maps/protein_complex_maps/preprocessing_util/complexes/complex_merge.py --help

python /project/rmcox/github_repos/protein_complex_maps/protein_complex_maps/preprocessing_util/complexes/complex_merge.py --cluster_filename all.gold.cmplx.txt --output_filename all.gold.cmplx.merged.humap.txt --merge_threshold 0.6 --complex_size_threshold 30
# yup this works
# it's also built into the cfms script

#######################################################################################
# Supervised learning for protein interactions (cfmsflow)
#######################################################################################

nextflow main.nf -params-file example_params/example_wholepipeline.json  # this works now
NXF_VER=20.04.1.5335 nextflow main.nf –params-file example_params/example_wholepipeline.json
/home/cmcwhite/programs/nextflow main.nf -params-file example_params/example_wholepipeline.json

# make feature matrix
../nextflow main.nf -params-file generate_params.json --entrypoint 1 --exitpoint 2

# train model
nextflow main.nf -params-file generate_params.json --entrypoint 3 --exitpoint 5

# gold standards:
/project/cmcwhite/data/gold_standards/corum_complexes/*

# do normalized elution profiles
ls -d $PWD/elut_files/*filtdollo.norm.elut > elutlist_filtdollo_norm.txt
tmux new -s featmat_filtdollo_norm
cd /project/rmcox/programs/cfmsflow/
cp parameters_featmat_filtdollo.json parameters_featmat_filtdollo_norm.json # edited to refer to correct paths
nextflow main.nf -params-file parameters_featmat_filtdollo_norm.json
# process is failing at cfmsiner_corr:alph_process (human.iex_15.prot_count_mFDRpsm001.unique.filtdollo.norm.elut)
# memory error?
# 2:05PM 8/31/21 -- [ 38%] 225 of 596

# lecaNOGs only, raw counts
ls -d $PWD/elut_files/*filtdollo.elut > elutlist_filtdollo.txt
cd /project/rmcox/programs/cfmsflow2/
nextflow main.nf -params-file parameters_featmat_filtdollo.json
tmux new -s featmat_filtdollo
nextflow main.nf -params-file parameters_featmat_filtdollo.json -resume
# on hopefog
../nextflow main.nf -params-file parameters_featmat_filtdollo.json -resume
/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/elutlist_filtdollo.txt

# gold standard proteins only, raw counts
ls -d $PWD/elut_files/*unique.gold.raw.elut > elutlist_goldstd.txt
cd /project/rmcox/programs/cfmsflow3/
nextflow main.nf -params-file generate_params.json --entrypoint 1 --exitpoint 2
mv parameters.json parameters_featmat_goldstd.json
tmux new -s featmat_goldstd
nextflow main.nf -params-file parameters_featmat_goldstd.json -resume
# on hopefog
../nextflow main.nf -params-file parameters_featmat_goldstd.json -resume
/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/elutlist_filtdollo.txt

# gold standard proteins only, normalized counts
Rscript /project/rmcox/LECA/scripts/normalize_elut.R --elut_dir /project/rmcox/LECA/ms/cfms2/elut_files --elut_pattern '*unique.gold.raw.elut'
ls -d $PWD/elut_files/*unique.gold.norm.elut > elutlist_goldstd_norm.txt
../nextflow main.nf -params-file parameters_featmat_goldstd_norm.json

# featmat locations
./cfmsflow_gn/work/77/ea29b66059d2d74a31af961534fbd5/  # for normalized
./cfmsflow_gr/work/2e/80f3c6c0a7cdad08a9e05362f79fc1/  # for raw data

# this is taking way too long... need to move it to hopefog
rmcox@hfogcomp01:/stor/work/Marcotte/project/rmcox/LECA/ms$ rsync -r rmcox@hopper.icmb.utexas.edu:/project/rmcox/LECA/ms/cfms2 .
rmcox@hfogcomp01:/stor/work/Marcotte/project/rmcox/programs$ rsync -r rmcox@hopper.icmb.utexas.edu:/project/rmcox/programs/cfmsflow3 .

# ok, this is going pretty fast, but there is some problem with the output being written to the specified directory
# I wonder if it has something to do with the file paths on Hopefog??
# maybe I can extract the files manually with the find command...
find /stor/work/Marcotte/project/rmcox/programs/cfmsflow_gn/ -iname '*gold*feat' -type f > feat_files.txt  # this works
find /stor/work/Marcotte/project/rmcox/programs/cfmsflow_gr/ -iname '*gold*feat' -type f > feat_files.txt 

while read file; do cp $file /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow/by_exp_goldstd_norm/; done < feat_files.txt
while read file; do cp $file /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow/by_exp_goldstd_raw/; done < feat_files.txt

# copy computed featmat to correct directory
cp /stor/work/Marcotte/project/rmcox/programs/cfmsflow_gn/work/77/ea29b66059d2d74a31af961534fbd5/featmat /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow/by_exp_goldstd_norm/
cp /stor/work/Marcotte/project/rmcox/programs/cfmsflow_gr/work/2e/80f3c6c0a7cdad08a9e05362f79fc1/featmat /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow/by_exp_goldstd_raw/

# copy computed features and feature matrix to hopper
rmcox@marccomp01:/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2$ rsync -avP rmcox@hfogcomp01.ccbb.utexas.edu:/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow/by_exp_goldstd_raw .
rmcox@marccomp01:/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2$ rsync -avP rmcox@hfogcomp01.ccbb.utexas.edu:/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow/by_exp_goldstd_norm .

# from the gold standard feature matrix, use this code to select features:
python /project/cmcwhite/data/protein_complexes_plant/scripts/feature_selections.py --help

# example
python /project/cmcwhite/data/protein_complexes_plant/scripts/feature_selections.py --feature_names ../feature_matrix/allplants_all_features.csv --plainfeatmat ../feature_matrix/allplants_feature_matrix_missing1.featmat.labeled1.scaled --output_file allplants_feature_matrix_missing1.featmat.labeled1.scaled.featselect

# so it looks like that "feature_names" file is just a comma-separated list of all the feature file names
# to make the --feature_names file:
ls by_exp_goldstd_raw/*gold*raw*feat | sed 's/by_exp.*\///' | sed ':a;N;$!ba;s/\n/,/g' > alleuks_goldstd_raw_features.csv
ls by_exp_goldstd_norm/*gold*norm*feat | sed 's/by_exp.*\///' | sed ':a;N;$!ba;s/\n/,/g' > alleuks_goldstd_norm_features.csv

ls *norm*feat | sed ':a;N;$!ba;s/\n/,/g' > /project/dbarth/feat_selection/all_archaea_features.csv

python /project/cmcwhite/data/protein_complexes_plant/scripts/feature_selections.py --feature_names alleuks_goldstd_norm_features.csv --plainfeatmat by_exp_goldstd_norm/featmat --output_file feature_selection/gold.norm.featselect

python /project/cmcwhite/data/protein_complexes_plant/scripts/feature_selections.py --feature_names /project/dbarth/feat_selection/all_archaea_features.csv --plainfeatmat /project/dbarth/software/archaea_cfmsflow_featmat1/output/featmat --output_file /project/dbarth/feat_selection/archaea.featselect

# daryl's stuff
python /project/cmcwhite/data/protein_complexes_plant/scripts/feature_selections.py --feature_names /project/dbarth/feat_selection/all_archaea_features.csv --plainfeatmat /project/dbarth/feat_selection/featmat.norm.labeled.csv --output_file feature_selection/gold.norm.featselect

# got an error about labels
# noticed that the label orders in the featmat file do not match the order of the .csv file
# also, it looks like Claire's input as extra columns for ID1, ID2, and label (true or false ppi)
# can I make this from cfmsflow, or should I write my own script?
../nextflow main.nf -params-file generate_params.json --entrypoint 3 --exitpoint 3
../nextflow main.nf -params-file parameters_featsel_goldstd_norm.json

# i ended up writing my own script
# on hopefog:
python /stor/work/Marcotte/project/rmcox/LECA/scripts/feature_selections.py --feature_names alleuks_goldstd_norm_features.csv --plainfeatmat featmat.gold.norm.labeled.csv --output_file feature_selection/gold.norm.101821.featselect  # this worked
python /stor/work/Marcotte/project/rmcox/LECA/scripts/feature_selections.py --feature_names alleuks_goldstd_raw_features.csv --plainfeatmat featmat.gold.raw.labeled.csv --output_file feature_selection/gold.raw.101821.featselect

# merged gold std assemblies w/ kevin's script
python /stor/work/Marcotte/project/rmcox/LECA/scripts/feature_selections.py --feature_names alleuks_goldstd_norm_features.csv --plainfeatmat featmat.gold.norm.labeled.csv --output_file feature_selection/gold.norm.102021.featselect  # this worked
python /stor/work/Marcotte/project/rmcox/LECA/scripts/feature_selections.py --feature_names alleuks_goldstd_raw_features.csv --plainfeatmat featmat.gold.raw.labeled.csv --output_file feature_selection/gold.raw.102021.featselect

# right now looks like the gradient boosting classifier might be best
# lets try building a model with what we've got
# this is on hopefog
../nextflow main.nf -params-file parameters3to5.json  # this doesn't work, will have to run each step separately?

../nextflow main.nf -params-file generate_params.json --entrypoint 3 --exitpoint 3
mv parameters.json parameters3.json
../nextflow main.nf -params-file parameters3.json

../nextflow main.nf -params-file generate_params.json --entrypoint 4 --exitpoint 4
mv parameters.json parameters4.json
../nextflow main.nf -params-file parameters4.json

../nextflow main.nf -params-file generate_params.json --entrypoint 5 --exitpoint 5
mv parameters.json parameters5.json
../nextflow main.nf -params-file parameters5.json


#######################################################################################
# Adding in Kevin's huMAP feature matrix
#######################################################################################
# this works
python3 /project/rmcox/LECA/scripts/replace_featmat_ids.py --orthomap /project/rmcox/LECA/ms/og_proteomes/nog_mapping/human.euNOG.diamond.mapping.2759.fmt --featmat /project/rmcox/LECA/ms/apms/humap2_ppis_ACC_20200821.pairsWprob.test

# now a bigger featmat
python3 /project/rmcox/LECA/scripts/replace_featmat_ids.py --orthomap /project/rmcox/LECA/ms/og_proteomes/nog_mapping/human.euNOG.diamond.mapping.2759.fmt --featmat /project/rmcox/LECA/ms/apms/humap2_ppis.upids.featmat.test

# now the full thing
python3 /project/rmcox/LECA/scripts/replace_featmat_ids.py --orthomap /project/rmcox/LECA/ms/og_proteomes/nog_mapping/human.euNOG.diamond.mapping.2759.fmt --featmat /project/rmcox/LECA/ms/apms/humap2_ppis.upids.featmat

# features that we actually want:
ext_Dm_guru,ext_Hs_malo,entropy_orig9k,zscore_orig9k,nwdscore_orig9k,plate_zscore_orig9k,uPeps_orig9k,ratio_orig9k,total_psms_orig9k,ratioTotalPSMs_orig9k,UtoTratio_orig9k,neg_ln_pval,pair_count,prey.bait.correlation,valid.values,log10.prey.bait.ratio,log10.prey.bait.expression.ratio,hein_neg_ln_pval,hein_pair_count,ave_apsm,nwdscore_bioplex2,zscore_bioplex2,plate_zscore_bioplex2,entropy_bioplex2,uPeps_bioplex2,ratio_bioplex2,total_psms_bioplex2,ratioTotalPSMs_bioplex2,UtoTratio_bioplex2,neg_ln_pval_bioplex2_Z4,pair_count_bioplex2_Z4,neg_ln_pval_bioplex2_Z2,pair_count_bioplex2_Z2,AvgSpec,AvgP,MaxP,Fold_Change,BFDR,AvgSpec_nonciliated_bioid,AvgP_nonciliated_bioid,MaxP_nonciliated_bioid,Fold_Change_nonciliated_bioid,BFDR_nonciliated_bioid,neg_ln_pval_cilium_hygeo,pair_count_cilium_hygeo,neg_ln_pval_cilium_hygeo_avgspec2,pair_count_cilium_hygeo_avgspec2,neg_ln_pval_cilium_hygeo_avgspec4,pair_count_cilium_hygeo_avgspec4,SAij,Sij,Sji,Mij,neg_ln_pval_boldt_apms_hygeo,pair_count_boldt_apms_hygeo,neg_ln_pval_boldt_apms_hygeo_gt4,pair_count_boldt_apms_hygeo_gt4,neg_ln_pval_treiber_hygeo_gt4,pair_count_treiber_hygeo_gt4,AvgSpec_youn,AvgP_youn,MaxP_youn,FoldChange,BFDR_youn,neg_ln_pval_youn_hygeo,pair_count_youn_hygeo,neg_ln_pval_youn_hygeo_gt2,pair_count_youn_hygeo_gt2,neg_ln_pval_youn_hygeo_gt4,pair_count_youn_hygeo_gt4,neg_ln_pval_treiber_hygeo_gt2,pair_count_treiber_hygeo_gt2

# replace all empty values with 0
sed -i 's/^,/0,/; :a;s/,,/,0,/g;ta' humap_ext_feats.csv

# get min and max for each feature
python3 ../../scripts/find_feat_minmax.py -f humap_ext_feats.csv

#######################################################################################
# Generating clade-specific featmats
#######################################################################################
touch sort_cmds.sh
while read x; do echo "cp elut_files/${x}*elut elut_animals/"; done < animal_codes.txt >> sort_cmds.sh
while read x; do echo "cp elut_files/${x}*elut elut_plants/"; done < plant_codes.txt >> sort_cmds.sh
while read x; do echo "cp elut_files/${x}*elut elut_tsar/"; done < tsar_codes.txt >> sort_cmds.sh
while read x; do echo "cp elut_files/${x}*elut elut_excavate/"; done < excavate_codes.txt >> sort_cmds.sh
cat sort_cmds.sh | parallel -j10

# let's just do these all with raw counts for now

# --- animals --- 
ls -d $PWD/elut_animals/animals_concat.raw.150p.gold.elut > elutlist_animals.txt
cd /stor/work/Marcotte/project/rmcox/programs/cfmsflow_animals/
../nextflow main.nf -params-file generate_params.json --entrypoint 1 --exitpoint 2
mv parameters.json parameters_animals_12.json
tmux new -s animals_feats
../nextflow main.nf -params-file parameters_animals_12.json -resume

cd /stor/work/Marcotte/project/rmcox/programs/cfmsflow_animals/output/
mv featmat featmat_animals
cp featmat_animals ../../../LECA/ms/cfms2/cfmsflow_120721/featmat/

# --- plants --- 
ls -d $PWD/elut_plants/plants_concat.raw.150p.elut > elutlist_plants.txt
cd /stor/work/Marcotte/project/rmcox/programs/cfmsflow_plants/
cp ../cfmsflow_animals/parameters_animals_12.json parameters_plants_12.json
sed -i 's/animals/plants/g' parameters_plants_12.json
tmux new -s plants_feats
../nextflow main.nf -params-file parameters_plants_12.json -resume

cd /stor/work/Marcotte/project/rmcox/programs/cfmsflow_plants/output/
mv featmat featmat_plants
cp featmat_plants ../../../LECA/ms/cfms2/cfmsflow_120721/featmat/

# --- tsar --- 
ls -d $PWD/elut_tsar/tsar_concat.raw.150p.elut > elutlist_tsar.txt
cd /stor/work/Marcotte/project/rmcox/programs/cfmsflow_tsar/
cp ../cfmsflow_animals/parameters_animals_12.json parameters_tsar_12.json
sed -i 's/animals/tsar/g' parameters_tsar_12.json
tmux new -s tsar_feats
../nextflow main.nf -params-file parameters_tsar_12.json -resume

cd /stor/work/Marcotte/project/rmcox/programs/cfmsflow_tsar/output/
mv featmat featmat_tsar
cp featmat_tsar ../../../LECA/ms/cfms2/cfmsflow_120721/featmat/

# --- excavate --- 
ls -d $PWD/elut_excavate/excavate_concat.raw.150p.elut > elutlist_excavate.txt
cd /stor/work/Marcotte/project/rmcox/programs/cfmsflow_excavate/
cp ../cfmsflow_animals/parameters_animals_12.json parameters_excavate_12.json
sed -i 's/animals/excavate/g' parameters_excavate_12.json
tmux new -s excavate_feats
../nextflow main.nf -params-file parameters_excavate_12.json -resume

cd /stor/work/Marcotte/project/rmcox/programs/cfmsflow_excavate/output/
mv featmat featmat_excavate
cp featmat_excavate ../../../LECA/ms/cfms2/cfmsflow_120721/featmat/

/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/results/leca_euks_elut.filtdollo.150p.raw.csv
leca_euks_elut.filtdollo.filt150p.raw.csv


#######################################################################################
# Filter for pearson's R
#######################################################################################

pr_files <- dir(path = "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_firstpass/by_exp_filtdollo_raw", pattern = "*pearsonR.feat$", full.names = TRUE)

pr_tidy <- pr_files %>%
  map(read_csv) %>%
  reduce(rbind)

pr_tidy[is.na(pr_tidy)] <- 0

print("Writing tidy correlations ...")
system.time(
  data.table::fwrite(pr_tidy, 
                     file.path(workdir, 
                               "cfms2/cfmsflow_120721/pearR_feats_allprots.csv"), 
                     sep = ",")
)
print("Done!")

python3 /stor/work/Marcotte/project/rmcox/LECA/scripts/filter_coefficients2.py -f pearR_feats_allprots.csv -t 0.3

for x in `ls *pearsonR.feat`; do echo "python3 /stor/work/Marcotte/project/rmcox/LECA/scripts/filter_coefficients2.py -f ${x} -t 0.3"; done > filter_coefficients.sh
cat filter_coefficients.sh | parallel -j16

# 96% reduction in size
cat *pearsonR.feat | wc -l
823289295
cat *corr_filt.feat | wc -l
50234264

# get all species, exps, & ids in 1 file
awk -F',' '{print FILENAME (NF?",":"") $0}' *corr_filt.feat > ids_exps_species.raw.150p.p3.csv
wc -l ids_exps_species.raw.150p.p3.csv
50234264

# filter for ids observed in 2 clades & retain max correlation
# testing
shuf -n 10000 ids_exps_species.raw.150p.p3.csv > ids_exps_species.raw.150p.p3.test.csv
python3 ../../../../scripts/filter_clades.py -f ids_exps_species.raw.150p.p3.test.csv -t 2

# full thing
python3 ../../../../scripts/filter_clades.py -f ids_exps_species.raw.150p.p3.csv -t 2

rmcox@hfogcomp01:/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/sql$ bash do_load_featmats.sh &> load_featmats.log
#######################################################################################
# Second pass w/ new filters/features
#######################################################################################
# in rmcox@hfogcomp01:/stor/work/Marcotte/project/rmcox/programs/cfmsflow_120721
../nextflow main.nf -params-file parameters3.json
../nextflow main.nf -params-file parameters4.json
../nextflow main.nf -params-file parameters5.json

# here are the top feats from RFE
# 69 total features
rmcox@hfogcomp01:/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_120721/model_training/topfeats_rfe.txt

# need to organize all potential features into one folder
# from top level cfmsflow github dirs:
find ./ -type f -name "*\.feat" -print
find cfmsflow_*/output -type f -name "*\.feat" -exec cp {} /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_120721/allfeats/ \;

while read x; do echo "cp ../cfmsflow_firstpass/by_exp_filtdollo_raw/${x} topfeats_rfe/"; done < model_training/topfeats_rfe.txt > get_cfms_feats.sh
cat get_cfms_feats.sh | parallel -j16 # only 6/68 from individual cfms, 5 are from arath

grep -c 'concat' topfeats_rfe.txt # 16 feats from concatenated clade matrices

# 22 total cfms features
# 47 total apms features

# need to regnerate elut iles filtered at 150 psms
# ran parse elut in R
# now setting up to run through cfmsflow
ls -d $PWD/elut_files/*filtdollo.raw.150p.elut > elutlist_150p_raw.txt
../nextflow main.nf -params-file generate_params.json --entrypoint 1 --exitpoint 2
mv parameters.json parameters12.json
../nextflow main.nf -params-file parameters12.json


#######################################################################################
# Third pass w/ all protein-protein PPIs w/ >= 0.3 pearson R and observed in >2 clades 
#######################################################################################
# in rmcox@hfogcomp01:/stor/work/Marcotte/project/rmcox/programs/cfmsflow_012022
../nextflow main.nf -params-file parameters3.json
../nextflow main.nf -params-file parameters4.json
../nextflow main.nf -params-file parameters5.json


# need to remove ribosomal proteins from gold standard file
grep '60S' leca_eunog_annots.020722.tsv | awk -F'\t' '{print $1}'

"_goldstandard_complexes": "/stor/work/Marcotte/project/rmcox/LECA/ms/gold_stds/all.gold.cmplx.noRibos.txt"
"postrain": "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.train_ppis.ordered"
"negtrain": "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/negtrain",
"postest": "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.test_ppis.ordered",
"negtest": "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/negtest"
"annotation_file": "/stor/work/Marcotte/project/rmcox/LECA/ms/annotations/leca_eunog_annots.020722.tsv"

#######################################################################################
# Debugging weird gold standard behavior
#######################################################################################

# SCORED OUTPUT:
- ID1	ID2	set
- KOG0477	KOG0480	train
- KOG0359	KOG0362	test
- KOG0359	KOG0364	test
- KOG0362	KOG0364	predicted

# input gold standards:
all.gold.cmplx.noRibos.txt: * KOG0363 "KOG0362" KOG0361 KOG0360 "KOG0364" KOG0357 KOG0358 "KOG0359"
all.gold.cmplx.noRibos.txt: * KOG3313 KOG0363 "KOG0362" KOG0361 KOG0360 KOG1760 "KOG0364" KOG3048 KOG0357 KOG3501 KOG0358 "KOG0359"
all.gold.cmplx.noRibos.txt: KOG0363 "KOG0362" KOG0361 KOG0357 KOG0358 "KOG0359"
all.gold.cmplx.noRibos.txt: KOG0363 KOG0360 "KOG0364" KOG0357 KOG0358 "KOG0359"
all.gold.cmplx.noRibos.txt: * KOG0363 "KOG0362" "KOG0364" KOG0357 KOG0358 "KOG0359"

# gold standards after pipeling performs merging:
grep 'KOG0359' *filt.t*ordered
goldstandard_filt.test_ppis.ordered: KOG0359 KOG0362
goldstandard_filt.test_ppis.ordered: KOG0359 KOG0364
goldstandard_filt.train_ppis.ordered: KOG0359 KOG3048
goldstandard_filt.train_ppis.ordered: KOG0359 KOG3501
goldstandard_filt.train_ppis.ordered: KOG0359 KOG1760
goldstandard_filt.train_ppis.ordered: KOG0359 KOG3313

grep 'KOG0362' *filt.t*ordered
goldstandard_filt.test_ppis.ordered: KOG0359 KOG0362
goldstandard_filt.test_ppis.ordered: KOG0357 KOG0362
goldstandard_filt.test_ppis.ordered: KOG0360 KOG0362
goldstandard_filt.test_ppis.ordered: ENOG502QPS5 KOG0362
goldstandard_filt.test_ppis.ordered: ENOG502QUYD KOG0362
goldstandard_filt.train_ppis.ordered: KOG0362 KOG1760
goldstandard_filt.train_ppis.ordered: KOG0362 KOG3501
goldstandard_filt.train_ppis.ordered: KOG0362 KOG3048
goldstandard_filt.train_ppis.ordered: KOG0362 KOG3313

grep 'KOG0364' *filt.t*ordered
goldstandard_filt.test_ppis.ordered: KOG0363 KOG0364
goldstandard_filt.test_ppis.ordered: ENOG502QUYD KOG0364
goldstandard_filt.test_ppis.ordered: ENOG502QPS5 KOG0364
goldstandard_filt.test_ppis.ordered: KOG0360 KOG0364
goldstandard_filt.test_ppis.ordered: KOG0359 KOG0364
goldstandard_filt.train_ppis.ordered: KOG0364 KOG3048
goldstandard_filt.train_ppis.ordered: KOG0364 KOG1760
goldstandard_filt.train_ppis.ordered: KOG0364 KOG3313
goldstandard_filt.train_ppis.ordered: KOG0364 KOG3501
goldstandard_filt.train_ppis.ordered: KOG0361 KOG0364

#######################################################################################
# Important miscellaneous notes
#######################################################################################

# 154 euk experiments, 29 prok experiments
# prokaryotic codes
proks <- c("halsa", "camnx", "pseae", "ecoli", "staa8")

# number of abundance measurements
wc -l */*/output/*group | tail -n 1 >> num_meas.txt  # 13236665 total

# syncing LECA project dir from hopper to hopefog
# the ms subdirectory still pretty massive
# anna got pretty upset about this, better not try it again
rsync -avrWL rmcox@hopper.icmb.utexas.edu:/project/rmcox/LECA/ms/ ms/

### // TO DO LIST //
# - ML model tuning:
# -- weight protein complexes in gold standard set
# -- add normalized concatenated features
# - make circle plot
# - make phylogenetic PPI tree
# - add volcano plot to diffprot
# - LECA interactome visualization
# -- make static interactive HTML

######################################################################################
# AlphaFold template
######################################################################################
rm jobfile.cmds
touch jobfile.cmds
while read jobname; do echo "module load alphafold/2.1.0-ctr; run_alphafold.sh --fasta_paths=/work/05819/rmcox/maverick2/alphafold/fastas/monomer/${jobname}.fasta --output_dir=/work/05819/rmcox/maverick2/alphafold/output --model_preset=monomer --flagfile=/work/projects/tacc/bio/alphafold/test/flags/full_dbs.ff"; done < fastas/monomer/monomer_files.txt >> jobfile.cmds


######################################################################################
# Generating elution profile figures by clade
######################################################################################
workdir <- "LECA/ms/cfms2/"
dir_animals <- "elut_animals"
dir_plants <- "elut_plants"
dir_tsar <- "elut_tsar"
dir_excavate <- "elut_excavate" 

# input_dir = paste0(workdir, dir_animals)
# identifier = str_extract(input_dir, "(?<=_)(\\w+)")

# cat_xtract_group <- function(data_path, identifier){
  
#   print("Locating data ...")
#   files <- dir(data_path, pattern = "*.filtdollo.norm.elut")  
  
#   print("Joining data while retaining filename information ...")
#   data <- data_frame(filename = files) %>% 
#     mutate(file_contents = map(filename,          
#                                ~read_csv(file.path(data_path, .)))) %>%
#     unnest(cols = c(file_contents))
  
#   print("Formatting species & experiment info ...")
#   data_fmt <- data %>%
#     separate(filename, into = c("species", "experiment"), sep="\\.") %>%
#     unique() %>% 
#     mutate(clade == identifier)
  
#   print("Done!")
  
# }

# sort tidy clade files
touch sort_tidy_elut_cmds.sh
while read x; do echo "cp elut_files/${x}*tidy elut_animals/"; done < meta/animal_codes.txt >> sort_tidy_elut_cmds.sh
while read x; do echo "cp elut_files/${x}*tidy elut_plants/"; done < meta/plant_codes.txt >> sort_tidy_elut_cmds.sh
while read x; do echo "cp elut_files/${x}*tidy elut_tsar/"; done < meta/tsar_codes.txt >> sort_tidy_elut_cmds.sh
while read x; do echo "cp elut_files/${x}*tidy elut_excavate/"; done < meta/excavate_codes.txt >> sort_tidy_elut_cmds.sh
cat sort_tidy_elut_cmds.sh | parallel -j16

# rscript for generating sparkline plots
Rscript /stor/work/Marcotte/project/rmcox/LECA/scripts/plot_sparklines.R -i Coatmer-2.csv  # template

# all input data
ls *csv > cmplx_list.txt
mv cmplx_list.txt input_data/

# -----> full sparklines
rm cmds/generate_sparklines_full.sh
while read cmplx_file; do echo "Rscript /stor/work/Marcotte/project/rmcox/LECA/scripts/plot_sparklines.R -i /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/sparklines/${cmplx_file} --best_exp FALSE --sample_clades FALSE 2>&1 | tee /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/sparklines/logs/${cmplx_file%.*}.full.log"; done < input_data/cmplx_list.txt > cmds/generate_sparklines_full.sh
cat cmds/generate_sparklines_full.sh | parallel -j8

# -----> best sampled experiment
rm cmds/generate_sparklines_bestexp.sh
while read cmplx_file; do echo "Rscript /stor/work/Marcotte/project/rmcox/LECA/scripts/plot_sparklines.R -i /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/sparklines/${cmplx_file} --best_exp TRUE --sample_clades FALSE 2>&1 | tee /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/sparklines/logs/${cmplx_file%.*}.bestexp.log"; done < input_data/cmplx_list.txt > cmds/generate_sparklines_bestexp.sh
cat cmds/generate_sparklines_bestexp.sh | parallel -j8

# -----> sampling of clades
rm cmds/generate_sparklines_bestexp_sampleclades.sh
while read cmplx_file; do echo "Rscript /stor/work/Marcotte/project/rmcox/LECA/scripts/plot_sparklines.R -i /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/sparklines/${cmplx_file} --best_exp TRUE --sample_clades TRUE 2>&1 | tee /stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/sparklines/logs/${cmplx_file%.*}.best_exp.sample_clades.log"; done < input_data/cmplx_list.txt > cmds/generate_sparklines_bestexp_sampleclades.sh
cat cmds/generate_sparklines_bestexp_sampleclades.sh | parallel -j8

# -----> run em all
cat cmds/* > cmds/make_all_sparklines.sh
cat cmds/make_all_sparklines.sh | parallel -j24

# -----> reset for a new run
rm *png
rm *pdf
rm logs/*log

######################################################################################
# Generating draft protein complex annotations
######################################################################################

# get 1 file with all eggNOG mapping
# in rmcox@hfogcomp01:/stor/work/Marcotte/project/rmcox/LECA/ms/og_proteomes/nog_mapping:
cat *euNOG.diamond.mapping.2759 >all.euNOG.diamond.mapping.2759


# testing multimapping
# KOG0079 = RAB35/IFT27
# tetrahymena IDs:
grep 'KOG0079' tetts.euNOG.diamond.mapping.2759
I7MM63  KOG0079
Q24D93  KOG0079

grep -A5 'I7MM63' tetts.fasta | grep -A5 'Q24D93' tetts.fasta
>tr|I7MM63|I7MM63_TETTS Ras family small GTPase OS=Tetrahymena thermophila (strain SB210) OX=312017 GN=TTHERM_00298510 PE=4 SV=2
MYTTNGQTFLKDYNMTQAAEVSSKIIPFPEKDKDVELFFFDISGQDCYQSITSELLGQAD
LLVLVYDCTSQESYNKLQTWYDKVKQKNSKSLQGVLISTKNDLSGARQVDPQQARDLAKK
LNLQFFEVSSARNTNIDLPFQNLAERLV
>tr|Q24D93|Q24D93_TETTS Small guanosine triphosphatase family Ras family protein OS=Tetrahymena thermophila (strain SB210) OX=312017 GN=TTHERM_01041920 PE=4 SV=1
MNGGASDGIAKVKILTLGESGVGKSSLLLRYKDDKFAGNFVTTLGVEYKQKQIQIENIPL
TVQVWDTAGQERFKTITPNYYRSVDGALLVFDISELETFDKIEYWIQDLQKDADLTKVIM
VLVGNKSDLEEKRKVDFQTAQKKAQEYQIEYIETSAKTNSNVDQVFEMIVRKVIDKKGGI
EKFKAQYQEQIAQNLNLKKANTSNGASGSNKNKQNNNNGGCC

# human IDs:
grep 'KOG0079' human.euNOG.diamond.mapping.2759
sp|Q15286|RAB35_HUMAN   KOG0079
sp|Q9BW83|IFT27_HUMAN   KOG0079

grep -A5 'sp|Q15286|RAB35_HUMAN' human.fasta | grep -A5 'sp|Q9BW83|IFT27_HUMAN' human.fasta
>sp|Q15286|RAB35_HUMAN Ras-related protein Rab-35 OS=Homo sapiens (Human) OX=9606 GN=RAB35 PE=1 SV=1
MARDYDHLFKLLIIGDSGVGKSSLLLRFADNTFSGSYITTIGVDFKIRTVEINGEKVKLQ
IWDTAGQERFRTITSTYYRGTHGVIVVYDVTSAESFVNVKRWLHEINQNCDDVCRILVGN
KNDDPERKVVETEDAYKFAGQMGIQLFETSAKENVNVEEMFNCITELVLRAKKDNLAKQQ
QQQQNDVVKLTKNSKRKKRCC
>sp|Q9BW83|IFT27_HUMAN Intraflagellar transport protein 27 homolog OS=Homo sapiens (Human) OX=9606 GN=IFT27 PE=1 SV=1
MVKLAAKCILAGDPAVGKTALAQIFRSDGAHFQKSYTLTTGMDLVVKTVPVPDTGDSVEL
FIFDSAGKELFSEMLDKLWESPNVLCLVYDVTNEESFNNCSKWLEKARSQAPGISLPGVL
VGNKTDLAGRRAVDSAEARAWALGQGLECFETSVKEMENFEAPFHCLAKQFHQLYREKVE
VFRALA

# global alignment output:
(
(
sp|Q15286|RAB35_HUMAN:0.31592
,
tr|Q24D93|Q24D93_TETTS:0.31592
):0.068025
,
(
sp|Q9BW83|IFT27_HUMAN:0.334459
,
tr|I7MM63|I7MM63_TETTS:0.334459
):0.049486
)
;

# so Q24D93 = RAB35, I7MM63 = IFT27 ... probably.


# get summed pep counts
python2 /stor/work/Marcotte/project/rmcox/LECA/scripts/msblender2elution.py --prot_count_files *.1  --output_filename tetts_iex1.pepcount --fraction_name_from_filename --msblender_format --spectral_count_type TotalCount --pepcount


####
# path to example slurm for colabfold
####

/work/08481/epenning/shared/launchers/sample_cf_launcher.slurm


###
# xlink mapping code
###

https://github.com/yaviddang20/cilia_xlink_new/tree/master/results


### localization_ml

# get quickgo annotations using their API
curl -X GET --header 'Accept:text/tsv' 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?taxonId=3702%2C9606%2C559292%2C312017&aspect=cellular_component&qualifier=located_in' -o quickgo_annots_all.tsv  # this is missing GO term column, and for some reason capped out at 10,000 ...

curl -X GET --header 'Accept:application/json' 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?proteome=%20gcrpCan&includeFields=goName&includeFields=taxonName&includeFields=name&taxonId=3702%2C9606%2C559292%2C312017&aspect=cellular_component&qualifier=located_in'

curl -X GET --header 'Accept:text/tsv' 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?proteome=%20gcrpCan&includeFields=goName&includeFields=taxonName&includeFields=name&taxonId=3702%2C9606%2C559292%2C312017&aspect=cellular_component&qualifier=located_in' # error: "messages":["Annotation proteome_unsorted requires 'geneProductType_unsorted=protein' to be set."

# ok, trying that:
curl -X GET --header 'Accept:text/tsv' 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?proteome=%20gcrpCan&geneProductType_unsorted=protein&includeFields=goName&includeFields=taxonName&includeFields=name&taxonId=3702%2C9606%2C559292%2C312017&aspect=cellular_component&qualifier=located_in' # same error

# another attempt with &geneProductType=protein and &downloadLimit=700000
curl -X GET --header 'Accept:text/tsv' 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?proteome=%20gcrpCan&includeFields=goName&includeFields=taxonName&includeFields=name&taxonId=3702%2C9606%2C559292%2C312017&aspect=cellular_component&qualifier=located_in&geneProductType=protein&downloadLimit=700000'

# David's output:
ls /stor/work/Marcotte/project/dy4652/local_ml_go_sepClassifierWeights


# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Regenerate heat map w/ all data
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# ::::::::::::::::::::::::::::::::::::::::::
# Concatenate all elution profiles
# ::::::::::::::::::::::::::::::::::::::::::

# had to reorder the T. brucei fractions, so need to make new concatenated elution matrix
# also rewrote join script in python while I was at it

# one for each clade
python3 scripts/join_eluts.py --input_dir ppi_ml/results/elutions/elut_animals/ --ext mFDRpsm001.unique.elut --outfile ppi_ml/results/elutions/elut_ordered/amorphea.unique.ordered.elut

python3 scripts/join_eluts.py --input_dir ppi_ml/results/elutions/elut_plants/ --ext mFDRpsm001.unique.elut --outfile ppi_ml/results/elutions/elut_ordered/viridiplantae.unique.ordered.elut

python3 scripts/join_eluts.py --input_dir ppi_ml/results/elutions/elut_tsar/ --ext mFDRpsm001.unique.elut --outfile ppi_ml/results/elutions/elut_ordered/tsar.unique.ordered.elut

python3 scripts/join_eluts.py --input_dir ppi_ml/results/elutions/elut_excavate/ --ext mFDRpsm001.unique.elut --outfile ppi_ml/results/elutions/elut_ordered/excavate.ordered.unique.elut

# this takes about an hour for the full set
python3 scripts/join_eluts.py --input_dir ppi_ml/results/elutions/elut_ordered/ --ext unique.elut.ordered --outfile ppi_ml/results/elutions/elut_ordered/leca.unique.ordered.elut

# ::::::::::::::::::::::::::::::::::::::::::
# Cluster concatenated elution profile
# ::::::::::::::::::::::::::::::::::::::::::

# want to sweep clustering metrics to see which one looks best

# test:
python3 scripts/cluster.py --infile ppi_ml/results/elutions/pkl/leca.unique.filtdollo.150p.norm.ordered.elut.pkl --distance_metric euclidean --cluster_method average --output_dir ppi_ml/results/elutions/clustered/ --pickle_output --plot

# generate parallelizable sweeps:
while read metric; do echo "python3 scripts/cluster.py --infile ppi_ml/results/elutions/pkl/leca.unique.filtdollo.150p.norm.ordered.elut.pkl --distance_metric ${metric} --cluster_method average --output_dir ppi_ml/results/elutions/clustered/test_params --plot"; done < ppi_ml/annotations/lists/distance_metrics.txt > ppi_ml/records/sweep_metrics_norm_avg_cmds.sh 
cat ppi_ml/records/sweep_metrics_norm_avg_cmds.sh | parallel -j24

for m in ward complete average single; do echo "python3 scripts/cluster.py --infile ppi_ml/results/elutions/pkl/leca.unique.filtdollo.150p.norm.ordered.elut.pkl --distance_metric euclidean --cluster_method ${m} --output_dir ppi_ml/results/elutions/clustered/test_params/ --plot"; done> ppi_ml/records/sweep_methods_norm_avg_cmds.sh

# best clustering:
# - braycurtis
# - correlation**
# - jaccard 

# generate heat maps for each step of filtering
python3 scripts/cluster.py --infile ppi_ml/results/elutions/pkl/leca.unique.filtdollo.0p.norm.elut.ordered.pkl --distance_metric correlation --cluster_method average --output_dir ppi_ml/results/elut_clustering/ --pickle_output --plot

python3 scripts/cluster.py --infile ppi_ml/results/elutions/pkl/leca.unique.150p.norm.ordered.elut.pkl --distance_metric correlation --cluster_method average --output_dir ppi_ml/results/elut_clustering/ --pickle_output --plot

python3 scripts/cluster.py --infile ppi_ml/results/elutions/pkl/leca.unique.filtdollo.150p.norm.ordered.elut.pkl --distance_metric correlation --cluster_method average --output_dir ppi_ml/results/elut_clustering/ --pickle_output --plot

# -------------------------
# -- generate sparklines --
# -------------------------
cd /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/sparklines
rm /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/sparklines/*
rm /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/sparklines/logs/*
# full sparklines
rm cmds/generate_sparklines_full.sh
for x in `ls /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/`; do echo "Rscript /stor/work/Marcotte/project/rmcox/leca/scripts/plot_sparklines.R \
-i /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/${x} \
-o /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/sparklines/${x%.*} \
--best_exp FALSE --sample_clades FALSE 2>&1 | tee /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/sparklines/logs/${x%.*}.full.log"; done > cmds/generate_sparklines_full.sh
cat cmds/generate_sparklines_full.sh | parallel -j8

# best sampled experiment "--best_exp TRUE --sample_clades FALSE", "${cmplx_file%.*}.bestexp.log"
rm cmds/generate_sparklines_bestexp.sh
for x in `ls /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/`; do echo "Rscript /stor/work/Marcotte/project/rmcox/leca/scripts/plot_sparklines.R \
-i /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/${x} \
-o /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/sparklines/${x%.*}.best_exp \
--best_exp TRUE --sample_clades FALSE 2>&1 | tee /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/sparklines/logs/${x%.*}.best_exp.log"; done > cmds/generate_sparklines_bestexp.sh
cat cmds/generate_sparklines_bestexp.sh | parallel -j8

# sampling of clades "--best_exp TRUE --sample_clades TRUE", "${cmplx_file%.*}.best_exp.sample_clades.log"
rm cmds/generate_sparklines_bestexp_sampleclades.sh
for x in `ls /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/`; do echo "Rscript /stor/work/Marcotte/project/rmcox/leca/scripts/plot_sparklines.R \
-i /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/${x} \
-o /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/sparklines/${x%.*}.best_exp.sample_clades \
--best_exp TRUE --sample_clades TRUE 2>&1 | tee /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/sparklines/logs/${x%.*}.best_exp.sample_clades.log"; done > cmds/generate_sparklines_bestexp_sampleclades.sh
cat cmds/generate_sparklines_bestexp_sampleclades.sh | parallel -j8

# ------------------------
# -- generate dot plots --
# ------------------------
cd /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/dotplots
rm /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/dotplots/*
rm /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/dotplots/logs/*
rm cmds/generate_dotplots.sh
for x in `ls /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/`; do echo "Rscript /stor/work/Marcotte/project/rmcox/leca/scripts/plot_phylo_dots.R -i /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/${x} -o /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/dotplots/${x%.*}_phyloplot 2>&1 | tee /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/dotplots/logs/${x%.*}.phylodots.log"; done > cmds/generate_dotplots.sh
cat cmds/generate_dotplots.sh | parallel -j8

# ----------------------------
# -- generate network plots --
# ----------------------------
cd /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/networks
rm /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/networks/*
rm /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/networks/logs/*
rm cmds/generate_networks.sh
for x in `ls /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/`; do echo "Rscript /stor/work/Marcotte/project/rmcox/leca/scripts/plot_ppi_networks.R \
--cmplx /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/${x} \
--scores /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/results_lj.noribo/annotatedscores_top50k.fmt.csv \
--annotations /stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/leca_euNOGs_plot_annots.csv \
--outfile /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/networks/${x%.*}_ppigraph \
--scan_layouts 2>&1 | tee /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/networks/logs/${x%.*}.ppigraph.log"; done > cmds/generate_networks.sh
cat cmds/generate_dotplots.sh | parallel -j8


# ::::::::::::::::::::::::::::::::::::::::::
# Get pep counts for each experiment
# ::::::::::::::::::::::::::::::::::::::::::

# template
python2 /stor/work/Marcotte/project/rmcox/leca/scripts/msblender2elution.py --prot_count_files *1 --output_filename arath_iex1.pepcount --fraction_name_from_filename --msblender_format --spectral_count_type TotalCount --pepcount

# get exp info
awk -F',' '{print $3}' leca_master_gsheet.csv | sed 's/\/project\/rmcox\/LECA\/ms\/cfms\/processed\///' | sed 's/\/mzXML//' | sed 's/\//\t/' | tail -n +2  > ms_exp_list.txt

# line by line info
while IFS=$'\t' read -r -a x; do
    echo "species: ${x[0]} | exp: ${x[1]}"
done < ms_exp_list.txt

# pep file path
for x in `ls /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/msblender/leca_level/*/*/output/*pep_count*1`; do echo $x; done

# combine
while IFS=$'\t' read -r -a x; do
    echo "python2 /stor/work/Marcotte/project/rmcox/leca/scripts/msblender2elution.py \
    --prot_count_files /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/msblender/leca_level/${x[0]}/${x[1]}/output/*pep_count*1 \
    --output_filename /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/pepcounts/${x[0]}_${x[1]}.pepcount \
    --fraction_name_from_filename \
    --msblender_format \
    --spectral_count_type TotalCount --pepcount"
done < /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/meta/ms_exp_list.txt > records/get_pepcounts.sh

cat records/get_pepcounts.sh | parallel -j16

# make sure output is good (in other words, check if file is empty)
for file in *.pepcount; do if [ ! -s $file ]; then echo $file; fi; done

# each .pepcount: 1 col for each fraction, 1 row for each peptide; 2nd col is total counts
# need to convert to 1 col for each experiment (total count), 1 row for each peptide
# probs do this in python: /stor/work/Marcotte/project/rmcox/leca/notebooks/gather_pep_data.ipynb


# need to join peptide <-> ID assignments
# e.g., from /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/msblender/leca_level/arath/iex_1/output/AT_indark_IEX_fraction_94a_20150503.prot_list

# example:
AT_indark_IEX_fraction_94a_20150503.00090.3.VLHTLLNRSGKQLR      KOG4341
AT_indark_IEX_fraction_94a_20150503.00097.7.KAPAERPNAKELLKHRFIKNARKSPKLLERIRERPK        KOG0201(pre=K,post=Y)

# is it found in peptotals files?
grep 'KAPAERPNAKELLKHRFIKNARKSPKLLERIRERPK' arath_peptotals.csv # nope
grep 'APAERPNAKELLKHRFIKNARKSPKLLERIRERP' arath_peptotals.csv # nope
grep 'VLHTLLNRSGKQLR' arath_peptotals.csv # nope

# back to the drawing board ...?
# potentially use this script I wrote forever ago...
# /stor/work/Marcotte/project/rmcox/leca/scripts/assign_peptides.py 

# new script for getting pep totals for each OG assignment:
python3 /stor/work/Marcotte/project/rmcox/leca/scripts/get_pep_assignments.py --data_dir /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/msblender/leca_level/arath/wwc_1/output/ --outfile /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/pep_assign/arath_wwc_1.pep_assign

# generate commands for each experiment
while IFS=$'\t' read -r -a x; do
    echo "python3 /stor/work/Marcotte/project/rmcox/leca/scripts/get_pep_assignments.py \
    --data_dir /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/msblender/leca_level/${x[0]}/${x[1]}/output/ \
    --outfile /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/pep_assign/${x[0]}_${x[1]}.pep_assign"
done < /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/meta/ms_exp_list.txt > records/get_pep_assignments.sh

cat records/get_pep_assignments.sh | parallel -j8 2>&1 | tee logs/get_pep_assignments.log

# did totals for each species manually at the end of this notebook: /stor/work/Marcotte/project/rmcox/leca/notebooks/gather_prot_lists.ipynb

# now make back assignment files:
python3 /stor/work/Marcotte/project/rmcox/leca/scripts/back_assign_peptides.py --peptides /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/pep_assign_totals/dicdi_pep_assign_totals.csv --fasta /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/proteomes/dicdi.fasta --mapping /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/og_proteomes/nog_mapping/dicdi.euNOG.diamond.mapping.2759 --outfile /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/pep_assign_posthoc/dicdi_ogs.back_assign_peps

# line by line info
while IFS=$'\t' read -r -a x; do
    echo "species: ${x[0]} | exp: ${x[1]}"
done < ms_exp_list.txt

# make command file
while IFS=$'\t' read -r -a x; do
    echo "python3 /stor/work/Marcotte/project/rmcox/leca/scripts/back_assign_peptides.py \
    --peptides /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/pep_assign_totals/${x[0]}_pep_assign_totals.csv \
    --fasta /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/proteomes/${x[0]}.fasta \
    --mapping /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/og_proteomes/nog_mapping/${x[0]}.euNOG.diamond.mapping.2759 \
    --outfile /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/cfms/pep_assign_posthoc/${x[0]}_ogs.back_assign_peps > logs/${x[0]}.back_assign_peps.log"
done < /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/meta/ms_exp_list.txt > records/back_assign_peps.sh

# run for each species (this takes a very, very, very long time)
uniq records/back_assign_peps.sh | cat | parallel -j16

# next step: choose sequence based on (1) 1:1 complex members, (2) OG subunit evidence, (3) subunit sequence length
# this file is potentially helpful: /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/meta/avg_famsize_per_og

# ::::::::::::::::::::::::::::::::::::::::::
# Get CFMS and APMS score for each PPI
# ::::::::::::::::::::::::::::::::::::::::::

# need to filter this file for paired scores:
/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/featmat/featmat_humap

# this file is massive; turned it into a pickled dictionary:
python3 scripts/make_apms_dict.py # -> /project/rmcox/leca/ppi_ml/data/apms/apms_dict.pkl

# performing this data munging in this notebook: /stor/work/Marcotte/project/rmcox/leca/notebooks/get_apms_scores.ipynb

# top feats (/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/model_training_lj.noribo/tpot_fitted_model.p.featureimportances):
0,0.0331,excavate_concat.raw.150p.spearmanR.feat,rfe
1,0.0276,plants_concat.raw.150p.spearmanR.feat,rfe
2,0.0258,plants_concat.raw.150p.pearsonR.feat,rfe
3,0.0243,tsar_concat.raw.150p.pearsonR.feat,rfe
4,0.0209,animals_concat.raw.150p.pearsonR.feat,rfe
5,0.0184,tsar_concat.raw.150p.spearmanR.feat,rfe
6,0.016,pair_count_cilium_hygeo,rfe  # col 26; WMM gupta(?)
7,0.0138,neg_ln_pval_youn_hygeo_gt2,rfe  # col 43; ; WMM youn
8,0.0137,neg_ln_pval_youn_hygeo,rfe  # col 41; WMM youn
9,0.0136,ave_apsm,rfe  # col 15; ?
10,0.013,neg_ln_pval_treiber_hygeo_gt2,rfe  # col 47; WMM youn

# min max file: /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/apms/humap/humap_ext_feats.minmax.csv


# -------------------------------
# -- generate apms table plots --
# -------------------------------
cd /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/apms_scores
rm /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/apms_scores/*
rm /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/apms_scores/logs/*
rm cmds/generate_apms_scores.sh
for x in `ls /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/`; do echo "Rscript /stor/work/Marcotte/project/rmcox/leca/scripts/plot_scores.R -i /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/${x} -o /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/apms_scores/${x%.*}_scores 2>&1 | tee /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/apms_scores/logs/${x%.*}.scores.log"; done > cmds/generate_apms_scores.sh
cat cmds/generate_apms_scores.sh | parallel -j16

# test set
shuf -n 30 cmds/generate_apms_scores.sh > cmds/test.sh
cat cmds/test.sh | parallel -j10

# hmmmm something weird is happening with the scores
# there are a large number of missing pairs
# eg: PSMA4 + PSMA5 / P25789 + P28066 / KOG0178 + KOG0176

# get all ids in the featmat pre-og assignment
tail -n +2 humap2_feature_matrix_20200820.featmat | awk -F',' '{print $1}' | uniq > all_ids1.txt
tail -n +2 humap2_feature_matrix_20200820.featmat | awk -F',' '{print $2}' | uniq > all_ids2.txt
cat all_ids1.txt all_ids2.txt | uniq > all_ids.txt

# let's get some numbers...
wc -l all_ids.txt # 17,579,788 total ids, but a lot of them are just like '1, 2, 3, 4, 5...'
grep 'sp\|tr' all_ids.txt | wc -l # 112,098 uniprot IDs

grep 'P25789\|P28066' all_ids.txt  # but neither PSMA4 or PSMA5 are in there ...

# how many ids ended up in the actual featmat?
tail -n +2 featmat_humap | awk -F',' '{print $1}' featmat_humap > fmat_ids.txt
tail -n +2 fmat_ids.txt | awk -F' ' '{print $1}' | uniq > fmat_ids1.txt
tail -n +2 fmat_ids.txt | awk -F' ' '{print $2}' | uniq > fmat_ids2.txt
cat fmat_ids1.txt fmat_ids2.txt | sort -u > all_fmat_ids.txt
wc -l all_fmat_ids.txt # 9886

# how many unique ogs?
grep 'KOG\|ENOG' all_fmat_ids.txt | sort -u > all_fmat_ogs.txt
wc -l all_fmat_ogs.txt # 6771 out of 20505 total possible human OGs

# these units are definitely in these humap files:
grep 'P25789\|P28066' humap2_ppis_ACC_20200821.pairsWprob
grep 'KOG0178' humap2_ppis_ACC_20200821.eggnogIDs.pairsWprob
grep 'KOG0176' humap2_ppis_ACC_20200821.eggnogIDs.pairsWprob

# individually, the ogs are in this file:
grep 'KOG0178\|KOG0176' humap_ext_feats_final.csv
	# but the pair + score are not
	grep 'KOG0178,KOG0176' humap_ext_feats_final.csv # nope
	grep 'KOG0176,KOG0178' humap_ext_feats_final.csv # nope
	# num pairs in this file:
	wc -l humap_ext_feats_final.csv  # 17,516,684


# what about this file?
wc -l orig9k_bioplex2_hygeo_bioid_hygeo_boldt_apms_hygeo_treiber_hygeo_wgt2_youn_hygeo_trimCols_groupbyMean_ACCs.featmat # 17,516,684

# these units are in this file
grep 'P25789\|P28066' orig9k_bioplex2_hygeo_bioid_hygeo_boldt_apms_hygeo_treiber_hygeo_wgt2_youn_hygeo_trimCols_groupbyMean_ACCs.featmat

# wait let's just check PSMA4
grep 'P25789' orig9k_bioplex2_hygeo_bioid_hygeo_boldt_apms_hygeo_treiber_hygeo_wgt2_youn_hygeo_trimCols_groupbyMean_ACCs.featmat  # oh PSMA4 straight up is not in here

# but the pair + score are not
grep 'P25789,P28066' orig9k_bioplex2_hygeo_bioid_hygeo_boldt_apms_hygeo_treiber_hygeo_wgt2_youn_hygeo_trimCols_groupbyMean_ACCs.featmat  # nope
grep 'P28066,P25789' orig9k_bioplex2_hygeo_bioid_hygeo_boldt_apms_hygeo_treiber_hygeo_wgt2_youn_hygeo_trimCols_groupbyMean_ACCs.featmat  # nope

# summary:
# kevin's humap2 probability scores file has all the relevant proteins
# the actual featmat files are missing UniProt IDs for some number of proteins (e.g., PSMA4/P25789)
# there's a possibility PSMA4 is in the featmat with a different ID
# if that's true... how do I find it?
# in a previous email chain with Kevin, we talked about ID mapping issues... he directed me to this notebook directory: https://github.com/marcottelab/protein_complex_maps/tree/humap2/protein_complex_maps/notebooks
# and even if I fix this problem... should I regenerate my model?
# (╯°□°)╯︵ ┻━┻

# --------------------------
# -- fix apms ID mapping --
# --------------------------

# entrez/ensembl --> uniprot
head /stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/uniprot/EntrezToUniprot_clean.csv

# i'm going to need to build in a check for non-unique pairs...
# eg, a pair might be unique pre-ID replacement, but is now duplicate post-replacement
# new featmat with uniprot IDs replaced:
head /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/apms/humap/orig9k_featmat.upids.csv
grep 'P25789' orig9k_featmat.upids.csv # wtf, still nothing
grep 'P25789' humap2_20200820_id1s_replaced.csv # yes
grep 'KOG0178,KOG0176' humap2_featmat_20200820.euNOGs.csv # yes

# ok, so we'll probably move forward with the featmat downloaded from the website
head /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/apms/humap/humap2_featmat_20200820.euNOGs.target_scores.group_mean.csv

# wait. after regenerating the apms score tables, it looks like there is still missing data
grep 'KOG0178' results/cmplx_files/* # 20S
grep 'KOG0178' data/apms/clustered_og_scores.csv # ok it is there. just need to regenerate the scores again

# how much did we improve?
# old apms fmat:
grep 'KOG\|ENOG' all_fmat_ids.txt | sort -u > all_fmat_ogs.txt
wc -l all_fmat_ogs.txt # 6771 out of 20505 total possible human OGs
# new apms fmat:
tail -n +2 featmat_humap | awk -F' ' '{print $1}' | uniq > fmat_ids1_new.txt
tail -n +2 featmat_humap | awk -F' ' '{print $2}' | uniq > fmat_ids2_new.txt
sed -i 's/,.*//g' fmat_ids2_new.txt
cat fmat_ids1_new.txt fmat_ids2_new.txt | sort -u > all_fmat_ids_new.txt
wc -l all_fmat_ids_new.txt # 7598

# so added 827 orthogroups (+12%);

# the real question is how many protein pairs were added
# and how many of them were gold standard

cd /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/apms_scores
rm /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/apms_scores/*
rm /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/apms_scores/logs/*
rm cmds/generate_apms_scores.sh
for x in `ls /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/`; do echo "Rscript /stor/work/Marcotte/project/rmcox/leca/scripts/plot_scores.R -i /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cmplx_files/${x} -o /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/apms_scores/${x%.*}_scores 2>&1 | tee /stor/work/Marcotte/project/rmcox/leca/ppi_ml/figures/apms_scores/logs/${x%.*}.scores.log"; done > cmds/generate_apms_scores.sh
cat cmds/generate_apms_scores.sh | parallel -j8

# ****TO DO: FIGURE THIS OUT:
# the real question is how many protein pairs were added
# and how many of them were gold standard

# --------------------------
# -- regenerate PPI model --
# --    (fourth pass)     --
# --------------------------

# let's review all our features:
* "featmat_allexps_p3c2" = 4X correlation metrics for each experiment with protein pairs with pearsons R greater than or equal to 0.3 (based on raw counts/our left-most df)
* "featmat_allconcat" = 4X correlation metrics for all raw elution profiles concatenated (150PSMs)
* "featmat_animals" = 4X correlation metrics for amorphea raw elution profiles concatenated (150PSMs)
* "featmat_excavate" = 4X correlation metrics for excavate raw elution profiles concatenated (150PSMs)
* "featmat_tsar" = 4X correlation metrics for tsar raw elution profiles concatenated (150PSMs)
* "featmat_plants" = 4X correlation metrics for archaeplastida raw elution profiles concatenated (150PSMs)
* "featmat_humap" = 47 APMS metrics from humap2 where IDs were converted to orthogroups and the metric became the group mean

# other input files
"_goldstandard_complexes": "/stor/work/Marcotte/project/rmcox/LECA/ms/gold_stds/all.gold.cmplx.noRibos.txt"
"postrain": "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.train_ppis.ordered"
"negtrain": "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/negtrain",
"postest": "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.test_ppis.ordered",
"negtest": "/stor/work/Marcotte/project/rmcox/LECA/ms/cfms2/cfmsflow_012022/model_training_lj.noribo/negtest"
"annotation_file": "/stor/work/Marcotte/project/rmcox/LECA/ms/annotations/leca_eunog_annots.020722.tsv"

# low key might re-do this without cfmsflow..
# getting kevin's scripts:
git clone https://github.com/marcottelab/protein_complex_maps.git

# complex merging:
python /stor/work/Marcotte/project/rmcox/programs/protein_complex_maps/protein_complex_maps/preprocessing_util/complexes/complex_merge.py --cluster_filename /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/gold_stds/all.gold.cmplx.noRibos.txt --output_filename /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/gold_stds/all.gold.cmplx.noRibos.merged.txt --complex_size_threshold 30 --merge_threshold 0.6

# how many did we merge?
wc -l *noRibos*
	# 1499 all.gold.cmplx.noRibos.merged.txt
	# 2508 all.gold.cmplx.noRibos.txt

# mmmm
# the feature extraction scripts actually do look kinda annoying
# step 1/2: just need to make an elut path file for target elution matrix & edit it into the params file
# i think i can do this in a bash script
# note: do not add poisson reps to normalized data

for exp in animals excavate tsar plants all150p; do
	new_name=cfmsflow_${exp}_norm
	git clone git@github.com:marcottelab/cfmsflow.git
	mv cfmsflow ${new_name}
done

for exp in animals excavate tsar plants all150p; do
	cd cfmsflow_${exp}_norm
	../nextflow main.nf -params-file generate_params.json --entrypoint 1 --exitpoint 2
	file=${exp}_norm_eluts_path.txt
	sed -i '/"input_elution_pattern": "",/d' parameters.json
	sed -i "s/\\\"input_elution_file\\\": \\\"\\\",/\\\"input_elution_file\\\": \\\"\/stor\/work\/Marcotte\/project\/rmcox\/programs\/data4cfmsflow\/elutions\/elut_paths\/${file}\\\",/" parameters.json
	cd ../
done

for exp in animals excavate tsar plants all150p; do
	echo "${exp}:"
	cat cfmsflow_${exp}_norm/parameters.json
	echo ""
done

for exp in animals excavate tsar plants all150p; do
	echo "cd cfmsflow_${exp}_norm && ../nextflow main.nf -params-file parameters.json"
done > get_normalized_feats.cmds
cat get_normalized_feats.cmds | parallel -j5

# permissions error again
# try symlinking the data into the same dir:
ln -s /stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/elutions/ .

# regenerate elut path:
for exp in animals excavate tsar plants all150p; do
	cd cfmsflow_${exp}_norm
	../nextflow main.nf -params-file generate_params.json --entrypoint 1 --exitpoint 2
	file=${exp}_norm_eluts_path.txt
	sed -i '/"input_elution_pattern": "",/d' parameters.json
	sed -i "s/\\\"input_elution_file\\\": \\\"\\\",/\\\"input_elution_file\\\": \\\"\/stor\/work\/Marcotte\/project\/rmcox\/programs\/data4cfmsflow\/elutions\/elut_paths\/${file}\\\",/" parameters.json
	cd ../
done

# still didn't work
# (╯°□°)╯︵ ┻━┻

# calculate pearsonR
python /indocker_repos/protein_complex_maps/protein_complex_maps/features/ExtractFeatures/canned_scripts/extract_features.py --format csv --normalize row_max -f pearsonR -o excavate.unique.norm.ordered.pearsonR.feat excavate.unique.norm.ordered.elut -r poisson_noise -i 0
# calculate braycurtis
python /indocker_repos/protein_complex_maps/protein_complex_maps/features/ExtractFeatures/canned_scripts/extract_features.py --format csv --normalize row_max -f braycurtis -o excavate.unique.norm.ordered.braycurtis.feat excavate.unique.norm.ordered.elut

# from protein_complex_maps:
python protein_complex_maps/protein_complex_maps/features/ExtractFeatures/canned_scripts/extract_features.py
# this error --> ImportError: No module named protein_complex_maps.features.ExtractFeatures.Features
# (╯°□°)╯︵ ┻━┻

# ok it works when i move it up to the right level
# claire must have messed with this when she made the nextflow pipeline
# cool cool cool
python /stor/work/Marcotte/project/rmcox/programs/protein_complex_maps/extract_features.py --help

# template
# feats: pearsonR, spearmanR, euclidean, braycurtis, sum_difference, spearmanR_weighted
python /stor/work/Marcotte/project/rmcox/programs/protein_complex_maps/extract_features.py ${fmat} --format csv --feature ${feat} --as_pickle

# test; this works
python /stor/work/Marcotte/project/rmcox/programs/protein_complex_maps/extract_features.py ../elutions/elut_ordered/leca.unique.filtdollo.150p.norm.ordered.elut --format csv --feature pearsonR --outfile all.norm.150p.pearsonR --as_pickle

for x in pearsonR spearmanR euclidean braycurtis covariance spearmanR_weighted; do
	echo "python /stor/work/Marcotte/project/rmcox/programs/protein_complex_maps/extract_features.py ../elutions/elut_ordered/leca.unique.filtdollo.150p.norm.ordered.elut --format csv --feature ${x} --outfile all.norm.150p.${x}.feat"
done > cmds/extract_feats_all.sh
cat cmds/extract_feats_all.sh  | parallel -j3

# for the rest of the clades
for x in amorphea excavate tsar viridiplantae; do for y in pearsonR spearmanR euclidean braycurtis covariance spearmanR_weighted; do echo $x $y; done; done > cmds/explist_clades.txt

# had to edit extract_features.py to put in a manually supplied outfile argument
while IFS= read line; do
	line=`echo $line | perl -pe '~s|\r?\n||' `
	code=`echo "$line" | awk '{print $1}'`
	feat=`echo "$line" | awk '{print $2}'`
	rootDir="/stor/work/Marcotte/project/rmcox"
	scriptDir="programs/protein_complex_maps"
	elutionDir="leca/ppi_ml/data/elutions/elut_ordered"
	featDir="leca/ppi_ml/data/calc_feats"
	infile="${code}.unique.filtdollo.norm.ordered.elut"
	outfile="${code}.filtdollo.norm.150p.${feat}.feat"
	cmd="python $rootDir/$scriptDir/extract_features.py"
	cmd="$cmd $rootDir/$elutionDir/$infile --format csv"
	cmd="$cmd --outfile $rootDir/$featDir/$outfile"
	cmd="$cmd --feature $feat"
	echo $cmd
done > cmds/extract_feats_clades.cmds < cmds/explist_clades.txt
cat cmds/extract_feats_clades.cmds | parallel -j8

# keep this in mind, since i'm doing this manually:
# (from cfmsflow modules)
	# if (corr == "euclidean" || corr == "braycurtis")
	#python ${params.protein_complex_maps_dir}/protein_complex_maps/features/normalize_features.py --input_feature_matrix $feat --output_filename ${feat}.rescaled --features $corr --min 0 --sep , --inverse
# eg eucl & brays require normalization (was done automatically on raw feats)

# this is theoreically the next step?
python ${params.protein_complex_maps_dir}/protein_complex_maps/features/alphabetize_pairs_chunks.py --feature_pairs $feat --outfile ${feat}.ordered --sep $sep --chunksize 1000000

# but the feats look unique to me? either way, i'll join them together use frozen sets so i guess it doesn't matter

# this is how kevin did it:
python ${params.protein_complex_maps_dir}/protein_complex_maps/features/build_feature_matrix.py --input_pairs_files $features --store_interval 10 --output_file featmat --sep ','

# followed by:
python ${params.protein_complex_maps_dir}/protein_complex_maps/features/add_label_chunks.py --input_feature_matrix $featmat --input_positives $positives --input_negatives $negatives --sep , --ppi_sep ' ' --id_column ID --output_file featmat_labeled --fillna 0 --id_sep ' ' --chunksize 100000

# kevin's labeling script
protein_complex_maps/preprocessing_util/complexes/split_complexes.py

### how i'm actually going to do it:

# pickle feats
for x in `ls /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/features/`; do echo "python3 /stor/work/Marcotte/project/rmcox/leca/scripts/make_pkl.py -i ${x}; done" > pickle_feat_files.sh

cat pickle_feat_files.sh | parallel -j16

# build fmat
python3 /stor/work/Marcotte/project/rmcox/leca/scripts/build_featmat.py --directory /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/features/ --pickle --left_join_file featmat_allexps_p3c2.pkl --outfile_name /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/featmats/featmat

# label fmat
python3 /stor/work/Marcotte/project/rmcox/leca/scripts/label_featmat.py --featmat /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/featmats/featmat.pkl --gold_std_file /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/gold_stds/all.gold.cmplx.noRibos.merged.txt --outfile_name /stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/featmats/featmat_labeled --seed 13 --shuffle_feats

# running TPOT:
python /indocker_repos/run_TPOT/train_TPOT.py --training_data featmat_labeled1 --outfile pipeline.py --template Selector-Classifier --selector_subset $SELECTORS_FORMATTED --classifier_subset $CLASSIFIERS_FORMATTED --id_cols 0 --n_jobs 20 --generations 10 --population_size 20 --labels -1 1 --temp_dir auto --groupcol traincomplexgroups --max_features_to_select 50

# ::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::

# go ahead and run my featmat through cfmsflow
# template to copy a working cfmsflow directory:
rsync -av --progress sourcefolder /destinationfolder --exclude thefoldertoexclude --exclude anotherfoldertoexclud

rsync -av --progress /stor/work/Marcotte/project/rmcox/programs/cfmsflow_010222 /stor/work/Marcotte/project/rmcox/programs/cfmsflow_012623 --exclude output --exclude work

# ::::::::::::::::::::::::::::::::::::::::::
# cfmsflow - step 3
# ::::::::::::::::::::::::::::::::::::::::::

# featmat path:
/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/featmats/featmat_final

# gold std path:
/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/gold_stds/all.gold.cmplx.noRibos.txt

# run step 3:
../nextflow main.nf -params-file parameters3.json

# ::::::::::::::::::::::::::::::::::::::::::
# cfmsflow - step 4
# ::::::::::::::::::::::::::::::::::::::::::

# existing featmat:
/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/featmats/featmat_final

# labeled featmat:
/stor/work/Marcotte/project/rmcox/programs/cfmsflow_012623/output/featmat_labeled

# train/test paths:
/stor/work/Marcotte/project/rmcox/programs/cfmsflow_012623/output/goldstandard_filt.neg_test_ppis.ordered

/stor/work/Marcotte/project/rmcox/programs/cfmsflow_012623/output/goldstandard_filt.neg_train_ppis.ordered

/stor/work/Marcotte/project/rmcox/programs/cfmsflow_012623/output/goldstandard_filt.test_ppis.ordered

/stor/work/Marcotte/project/rmcox/programs/cfmsflow_012623/output/goldstandard_filt.train_ppis.ordered

# CV groups:
/stor/work/Marcotte/project/rmcox/programs/cfmsflow_012623/output/traincomplexgroups

# run step 4:
../nextflow main.nf -params-file parameters4.json

# ::::::::::::::::::::::::::::::::::::::::::
# cfmsflow - step 5
# ::::::::::::::::::::::::::::::::::::::::::

# scored interactions:
/stor/work/Marcotte/project/rmcox/programs/cfmsflow_012623/output/scored_interactions

# annotation file:
/stor/work/Marcotte/project/rmcox/leca/ppi_ml/annotations/leca_euNOGs_human-arath_annotated.tsv

# elution file:
/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/elutions/leca_euks_elut.filtdollo.filt150p.raw.noheaders.csv

# run step 5:
../nextflow main.nf -params-file parameters5.json

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# rerun with test/train split from first run
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# ::::::::::::::::::::::::::::::::::::::::::
# cfmsflow - step 4
# ::::::::::::::::::::::::::::::::::::::::::

# train/test paths:
/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.train_ppis.ordered

/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.neg_train_ppis.ordered

/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.test_ppis.ordered

/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.train_ppis.ordered

# groups:
/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/model_training_lj.noribo/traincomplexgroups

# run step 4:
../nextflow main.nf -params-file parameters4_oldTTsplit.json

#
"postrain": "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.train_ppis.ordered",
"negtrain": "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.neg_train_ppis.ordered",
"postest": "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.test_ppis.ordered",
"negtest": "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/cfmsflow_012022/model_training_lj.noribo/goldstandard_filt.neg_test_ppis.ordered",

ls /stor/work/Marcotte/project/rmcox/programs/cfmsflow_012623/featmat_labeled