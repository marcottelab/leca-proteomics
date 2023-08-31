# clean list of commands used for LECA processing & analyses
# major sections are found by searching for  '# ||||' within this file
# minor sections are found by searching for '# ::::' within this file

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# LECA GENOME CONSTRUCTION
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# these are the commands run to compute dollo parsimony on QoF/SwissProt reference species orthogroups [May 2020 - June 2020]

# ::::::::::::::::::::::::::::::::::::::::::
# Download data
# ::::::::::::::::::::::::::::::::::::::::::

# 78 reference organisms (RELEASE 22-Apr-2020; RETRIEVED 4-May-2020)
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

readlink -f *fasta >> all_proteome_paths.txt

# ::::::::::::::::::::::::::::::::::::::::::
# Map eggNOG groups
# ::::::::::::::::::::::::::::::::::::::::::

# create commands
cat all_proteome_paths.txt | { while read path; do proteome=`basename $path .fasta`; echo "nice -5 python /project/rmcox/programs/eggnog-mapper/emapper.py -i ${path} --output /project/rmcox/LECA/qfo_tree/eggnog_mapping/${proteome}.euNOG.diamond -d euNOG --override --resume --scratch_dir /project/rmcox/LECA/qfo_tree/eggnog_mapping/scratch --no_file_comments --keep_mapping_files -m diamond"; done } > map_qfo_organisms.sh

# execute eggnog mapper
cat map_qfo_organisms.sh | parallel -j32

# format eggnog output for LECA-level taxonomic scope (2759)
cat all_proteome_paths.txt | { while read path; do proteome=`basename $path .fasta`; echo "Rscript ../scripts/format_emapper_output.R -f eggnog_mapping/${proteome}.euNOG.diamond.emapper.annotations -o eggnog_mapping/${proteome}.euNOG.diamond.mapping.2759 -s diamond -l 2759 -v 2 -p ENOG50"; done } > format_emapper_output_2759.sh

# format eggnog output for LUCA-level taxnomic scope (1)
cat all_proteome_paths.txt | { while read path; do proteome=`basename $path .fasta`; echo "Rscript ../scripts/format_emapper_output.R -f eggnog_mapping/${proteome}.euNOG.diamond.emapper.annotations -o eggnog_mapping/${proteome}.euNOG.diamond.mapping.1 -s diamond -l 1 -v 2 -p ENOG50"; done } > format_emapper_output_1.sh

cat records/format_emapper_output_* | parallel -j32

cat eggnog_mapping/*mapping.1 > allspecies_1.mapping
wc -l allspecies_1.mapping # 842063 genes mapped

cat eggnog_mapping/*mapping.2759 > allspecies_2759.mapping
wc -l allspecies_2759.mapping # 745713 genes mapped

# format uniprot IDs for future purposes:
for x in *mapping.2759; do Rscript /project/rmcox/LECA/scripts/format_uniprot_ids.R -f $x; done

# ::::::::::::::::::::::::::::::::::::::::::
# Evaluate Dollo parsimony
# ::::::::::::::::::::::::::::::::::::::::::

wget http://www.iro.umontreal.ca/~csuros/gene_content/Count.tgz
tar -xvzf Count.tgz
rm Count.tgz

# path to Java executable: /project/rmcox/programs/Count/Count.jar

# from within /project/rmcox/LECA/qfo_tree:
# IMPORTANT: remember to (1) start running XMing, (2) set X11 to froward, then (3) start the application
java -Xmx1024M -jar /project/rmcox/programs/Count/Count.jar  # might need to give it more memory

# vim command for formatting the .nwk species tree:
:%s/__\w\+_\w\+__\d\+//g


# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# LECA CFMS PROCESSING
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# these are the commands run to process internal and external co-fractionation mass spectrometry (CFMS) data through MSblender [January 2021 - May 2021]

# ::::::::::::::::::::::::::::::::::::::::::
# Download external data
# ::::::::::::::::::::::::::::::::::::::::::

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

# ::::::::::::::::::::::::::::::::::::::::::
# Get proteomes for each species
# ::::::::::::::::::::::::::::::::::::::::::

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

# ::::::::::::::::::::::::::::::::::::::::::
# Map proteomes to eggNOG groups
# ::::::::::::::::::::::::::::::::::::::::::

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

# ::::::::::::::::::::::::::::::::::::::::::
# Create directories for each species
# ::::::::::::::::::::::::::::::::::::::::::

awk -F',' '{print tolower($5)}' speciesinfo_02112021.csv | tail -n +2 > speciescodes.txt

while read id; do echo "mkdir /project/rmcox/LECA/ms/cfms/raw/${id}"; done < codes.txt > /project/rmcox/LECA/records/make_ms_dirs.sh
while read id; do echo "mkdir /project/rmcox/LECA/ms/cfms/processed/${id}"; done < codes.txt > /project/rmcox/LECA/records/make_ms_dirs.sh
bash /project/rmcox/LECA/records/make_ms_dirs.sh

# ::::::::::::::::::::::::::::::::::::::::::
# Create MSblender directories
# ::::::::::::::::::::::::::::::::::::::::::

# master MS spreadsheet: https://docs.google.com/spreadsheets/d/1IR60UjWbeUV42kCOP8jS-JiQUr8Y8fefBBUjsUPIH9o/edit?usp=sharing

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

# symlink data using master spreadsheet
# check that all mzXML files have been successfully symlinked
find /project/rmcox/LECA/ms/cfms/processed/ -type d -empty | grep 'mzXML' # only second euglena fractionation missing

# put contam in all dirs
vim make_db_dirs.sh
:%s/mkdir -p/cp contam.fasta/
:w! cp_contam_fasta.sh
cat cp_contam_fasta.sh | parallel -j8

# concat dbs and contam files
while read code; do echo "cat ${code}*collapsed.fasta ../../records/contam.fasta > contam_combined/${code}.collapsed_contam.combined.fasta"; done < ../codes.txt > combine_contams.sh

# ::::::::::::::::::::::::::::::::::::::::::
# Run MSblender (TACC)
# ::::::::::::::::::::::::::::::::::::::::::

# first attempt at processing ~13k fractions on Hopper crashed the server
# --> worked with Anna Battenhouse to create a pipeline on TACC
# --> used same commands/proteomes described above
# --> however, paths described in this section will refer to my project directory on TACC instead of Hopper

# MSblender dir:
/work/projects/BioITeam/ls5/gitdir/MSBlender/MSblender_consistent/runMS2.sh

# Anna's new "create commands" on TACC
# ensures working dir is always printed
echo "bash $runscript $mzxml_file $proteome $working_dir $output_dir $searchgui_dir > $working_dir/$pfx.msblender.log 2>&1" >> $EXP_COMMANDS

# generate commands:
bash /work/projects/BioITeam/ls5/gitdir/MSBlender/run_msblender/scripts/create_commands.sh exp_list.txt /work/projects/BioITeam/ls5/gitdir/MSBlender/MSblender_consistent/runMS2.sh /work/05819/rmcox/lonestar/SearchGUI-3.3.5
find /scratch/05819/rmcox/cfms/processed/ -iname '*COMMANDS.txt' -type f > msblender_cmds.txt

# set up commands for parallel processing
# all fractions
for exp in `cat msblender_cmds.txt`; do echo "`cat $exp`"; done >> all_ms_exps.sh

# how many?
wc -l all_ms_exps.sh # 12692

# split into sets of 40
split -l 40 ../all_ms_exps.sh --numeric-suffixes cfms --additional-suffix=.cmds

# put cmds into sbatch files
for exp in cfms*cmds; do echo "/work/projects/BioITeam/common/script/launcher_maker.py -t 1 -n ${exp} -w 8 -a A-cm10 -q dev -e rachaelcox@utexas.edu"; done > make_slurm_files.sh
bash make_slurm_files.sh
ls *slurm > all_slurm_files.txt

# script for serial slurm submission:

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

# run msblender on all fractions
bash serial_submit.sh all_slurm_files.txt rmcox finished > serialsubmit.log 2>&1

# ::::::::::::::::::::::::::::::::::::::::::
# QC MSblender output
# ::::::::::::::::::::::::::::::::::::::::::

# to count the actualy number of fractions
# inside /project/rmcox/LECA/ms/cfms/processed:
ls -x */*/mzXML/* > fraction_count.txt # 12692

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

# run status of all 12,691 fractions is here:
/project/rmcox/LECA/ms/cfms/meta/all_reruns.info.txt

# run statuses are summarized here:
/project/rmcox/LECA/ms/cfms/meta/all_stats_by_result.txt
# --- good w/ 3 search alg: 10,550
# --- good w/ 2 search algs: 1,012
# --- no proteins detected: 1,081 (almost all species have at least 1)
# --- no viable spectra: 26 (19 human + 4 arath + 2 nemve + 1 wheat)
# --- parse xml error: 22 (21 plaf7 + 1 wheat)

# locating all results:
find processed/ -type f -iname *.group > elut_files/groupfile_paths_060721.txt

# edit so we can glob for experiments
vim groupfile_paths.txt
:%s/output.*/output\/\*\.group/  # 11,569 "good" fractions

# counting unique protein IDs
awk -F',' '{print $2}' halobacterium.unique.protcount.tidy | sort -u | wc -l  # 2337 unique protein IDs

# ::::::::::::::::::::::::::::::::::::::::::
# Format MSblender output
# ::::::::::::::::::::::::::::::::::::::::::

# template to make .eluts
python /stor/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/preprocessing_util/msblender2elution.py --prot_count_files /project/rmcox/programs/run_msblender/halobacterium/output/*.group --output_filename elut_files/halobacterium.prot_count_mFDRpsm001.unique.elut --fraction_name_from_filename --parse_uniprot_id --remove_zero_unique

# ended up making the group_to_elut commands in the master google doc (I was not very good at bash scripting at this point in my life)
# run these commands:
cat group_to_elut_cmds.sh | parallel -j16 2>&1 | tee group_to_elut.log

# see if any didn't work
grep 'IOError' group_to_elut.log

# made my own script here for processing to tidy format:
/project/rmcox/LECA/scripts/elut2tidy.R

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

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# LECA CFMS ANALYSIS
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# ::::::::::::::::::::::::::::::::::::::::::
# Concatenate elution profiles
# ::::::::::::::::::::::::::::::::::::::::::

# function for concatenating elut data
concat_elut <- function(elut_list, psm_thres, outfile_prefix){
  
  # function for reading in data
  read_elut <- function(elut_file){
    
    elut_df <- read_delim(elut_file, delim = ",")
    
    return(elut_df)
  }
  
  results <- mapply(read_elut, 
                    elut_file = elut_list)
  results[1] <- NULL
  
  # takes a long time
  elut_concat <- results %>%
    reduce(full_join, by = 'X1')
  elut_concat[is.na(elut_concat)] <- 0
  
  path <- dirname(normalizePath(elut_list[1]))
  
  # these take a long time to write
  write_tsv(elut_concat, paste0(path, "/", outfile_prefix, ".raw.txt"))
}

# ::::::::::::::::::::::::::::::::::::::::::
# Hierarchical clustering of elution matrix
# ::::::::::::::::::::::::::::::::::::::::::

# for preliminary visualization;
# clustered concatenated elution profiles using morpheus:
# https://software.broadinstitute.org/morpheus/

# want to sweep clustering metrics to see which one looks best

# test my own:
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
python3 scripts/cluster.py --infile ppi_ml/results/elutions/pkl/leca.unique.150p.norm.ordered.elut.pkl --distance_metric correlation --cluster_method average --output_dir ppi_ml/results/elutions/clustered/ --pickle_output --plot

python3 scripts/cluster.py --infile ppi_ml/results/elutions/pkl/leca.unique.filtdollo.150p.norm.ordered.elut.pkl --distance_metric correlation --cluster_method average --output_dir ppi_ml/results/elutions/clustered/ --pickle_output --plot

# ::::::::::::::::::::::::::::::::::::::::::
# Gold standard PPIs (Complex Portal)
# ::::::::::::::::::::::::::::::::::::::::::

# download data
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

# need to re-format emapper output such that: (1) no header, (2) space-separated
while read code; do echo "Rscript /project/rmcox/LECA/scripts/format_uniprot_ids.R -f /project/rmcox/LECA/ms/gold_stds/proteomes/nog_mapping/${code}.euNOG.diamond.mapping.2759"; done < ../gold_ids.txt > format_uniprot_ids.sh
cat format_uniprot_ids.sh | parallel -j12

while read code; do echo "tail -n +2 /project/rmcox/LECA/ms/gold_stds/proteomes/nog_mapping/${code}.euNOG.diamond.mapping.2759 | sed 's/\t/ /g' > ac_mappings/${code}.euNOG.diamond.mapping.2759.fmt"; done < ../gold_ids.txt > format_ac_mapping.sh

# run Claire's complex mapping script
python /project/cmcwhite/data/gold_standards/corum_complexes/translate_corum2.py --complexes /project/rmcox/LECA/ms/gold_stds/ppis/9606.gold.cmplx.txt --mapping /project/rmcox/LECA/ms/gold_stds/ac_mappings/9606.euNOG.diamond.mapping.2759.fmt --outfile 9606.euNOG.gold.cmplx.txt

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

# ::::::::::::::::::::::::::::::::::::::::::
# Gold standard PPIs (CORUM)
# ::::::::::::::::::::::::::::::::::::::::::

# make eggNOG file for CORUM complexes
# 10090 = mouse, 10116 = rat, 9606 = human, 9913 = bovine, 9823 = pig
touch corum.euNOG.diamond.mapping.2759.fmt
while read x; do cat ${x}.euNOG.diamond.mapping.2759.fmt >> corum.euNOG.diamond.mapping.2759.fmt; done < ../corum_codes.txt

# run complex mapping script
python /project/cmcwhite/data/gold_standards/corum_complexes/translate_corum2.py --complexes /project/rmcox/LECA/ms/gold_stds/corum/corum.cmplx.txt --mapping /project/rmcox/LECA/ms/gold_stds/ac_mappings/corum.euNOG.diamond.mapping.2759.fmt --outfile /project/rmcox/LECA/ms/gold_stds/corum/corum.euNOG.gold.cmplx.txt

# ::::::::::::::::::::::::::::::::::::::::::
# Gold standard PPIs (manual ARATH data)
# ::::::::::::::::::::::::::::::::::::::::::

# previously acquired ARATH data in /project/rmcox/LECA/ms/gold_stds/:
./proteomes/3702.fasta
./nog_mapping/3702.euNOG.diamond.mapping.2759
./ac_mappings/3702.euNOG.diamond.mapping.2759.fmt
./comportal/3702.annotations.gold.cmplx.txt

# run complex mapping script
python /project/cmcwhite/data/gold_standards/corum_complexes/translate_corum2.py --complexes /project/rmcox/LECA/ms/gold_stds/claire/claire.cmplx.txt --mapping /project/rmcox/LECA/ms/gold_stds/ac_mappings/3702.euNOG.diamond.mapping.2759.fmt --outfile /project/rmcox/LECA/ms/gold_stds/claire/claire.euNOG.gold.cmplx.txt

# ::::::::::::::::::::::::::::::::::::::::::
# Gold standard PPIs (merging similar complexes)
# ::::::::::::::::::::::::::::::::::::::::::

# kevin's script
python /project/rmcox/github_repos/protein_complex_maps/protein_complex_maps/preprocessing_util/complexes/complex_merge.py --help

python /project/rmcox/github_repos/protein_complex_maps/protein_complex_maps/preprocessing_util/complexes/complex_merge.py --cluster_filename all.gold.cmplx.txt --output_filename all.gold.cmplx.merged.humap.txt --merge_threshold 0.6 --complex_size_threshold 30

# :::::::::::::::::::::::::::::::::::::::::::::::::::
# rachael's cirque du machine learning
# :::::::::::::::::::::::::::::::::::::::::::::::::::

# sync to latest version:
git pull
conda env update --file environment.yml --prune

# set directory vars to make cmds more readable
# also makes it easier for:
#       - moving the ML pipeline scripts to a new repo or different clone
#       - changing project directory without editing a bunch of cmds
script_dir="/stor/work/Marcotte/project/rmcox/leca/scripts"
proj_dir="/stor/work/Marcotte/project/rmcox/leca/ppi_ml"

# if setting up a new project directory:
cd $proj_dir
mkdir {data,results,logs,cmds}
mkdir data/{feature_files,feature_subsets,feature_matrix}
mkdir results/{tpot,feature_selection,ppi_predict,walktrap}
tree -d * # confirm directory structure

# ---------------------------------------------------
# build feature matrix
# ---------------------------------------------------

# pickle feats (run time ~5min)
[ ! -e cmds/pickle_feat_files.sh ] || rm cmds/pickle_feat_files.sh; touch cmds/pickle_feat_files.sh
for x in `ls --ignore=*.pkl ${proj_dir}/data/features/`; do echo "python3 ${script_dir}/make_pkl.py -i ${proj_dir}/data/features/${x}"; done >> cmds/pickle_feat_files.sh
cat cmds/pickle_feat_files.sh | parallel -j24

# combine subsets (run time ~45min)
for dir in all chicken dolphin mouse pig rabbit; do echo \
"python3 ${script_dir}/build_featmat.py \
--directory ${proj_dir}/data/feature_files/${dir}/ \
--pickle --outfile_name ${proj_dir}/data/feature_subsets/${dir}_features"; done > cmds/join_subsets.sh
cat cmds/join_subsets.sh | parallel -j6 2>&1 | tee ${proj_dir}/logs/join_subsets.log

# combine all (run time ~20min)
python3 ${script_dir}/build_featmat.py --directory ${proj_dir}/data/feature_subsets/ \
--pickle --outfile_name ${proj_dir}/data/feature_matrix/feature_matrix 2>&1 | tee ${proj_dir}/logs/build_featmat.log

# ---------------------------------------------------
# label feature matrix
# ---------------------------------------------------

# label positive/negative PPIs (run time ~2 hours)
# seed makes negative PPI labels reproducible
python3 ${script_dir}/label_featmat.py \
--featmat ${proj_dir}/data/featmats/featmat.pkl \
--gold_std_file ${proj_dir}/data/gold_stds/all.gold.cmplx.noRibos.merged.txt \
--outfile_name ${proj_dir}/data/featmats/featmat_labeled --seed 13 2>&1 | tee ${proj_dir}/logs/label_featmat.log

# ---------------------------------------------------
# get optimized ML pipeline
# ---------------------------------------------------

# tpot "generations" and "pop_size" set here to default cfmsflow params
# seed makes group splits, train/test splits, & model reproducible
# extremely slow (will probably run overnight)
# to speed things up reduce num_splits, generations, & pop_size
python3 ${script_dir}/run_tpot.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--outdir ${proj_dir}/results/tpot/gss/ \
--group_split_method GroupShuffleSplit --num_splits 5 --train_size 0.75 \
--generations 10 --pop_size 20 --seed 17 2>&1 | tee ${proj_dir}/logs/run_tpot_gss.log

python3 ${script_dir}/run_tpot.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--outdir ${proj_dir}/results/tpot/gkfold/ \
--group_split_method GroupKFold --num_splits 5 --train_size 0.75 \
--generations 10 --pop_size 20 --seed 17 2>&1 | tee ${proj_dir}/logs/run_tpot_gkfold.log

python3 ${script_dir}/run_tpot.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--outdir ${proj_dir}/results/tpot/sgkfold/ \
--group_split_method StratifiedGroupKFold --num_splits 5 --train_size 0.75 \
--generations 10 --pop_size 20 --seed 13 2>&1 | tee ${proj_dir}/logs/run_tpot_sgkfold.log

# get all test scores
grep --with-filename 'Test set score' ../../logs/run_tpot*log > all_test_scores.txt

# get all pipelines
sed -n '/Average/,/Fix random/p' gss/tpot_pipeline_1

[ ! -e all_pipelines.txt ] || rm all_pipelines.txt; touch all_pipelines.txt
for dir in gkfold gss sgkfold; do 
        echo "-------------------------------------------------------"
        echo "-------------------- ${dir} ---------------------------" 
        echo "-------------------------------------------------------"
        sed -n '/Average/,/Fix random/p' ${dir}/tpot_pipeline* |
        sed '/# Fix random state for all the steps in exported pipeline/d' |
        sed '/# Fix random state in exported estimator/d' |
        sed 's/# Average/\n# Average/'
done >> all_pipelines.txt

# best models:
gkfold/tpot_model_3  # Linear SVC (3 steps)
gkfold/tpot_model_5  # SGD (3 steps)
gss/tpot_model_1     # KNN (3 steps)
gss/tpot_model_4     # ExtraTrees (1 step)

# ---------------------------------------------------
# run recursive feature elimination
# ---------------------------------------------------

# choose best tpot models as input
# note: KNN does not provide feature selection logic

# seed makes group splits, train/test splits, & model reproducible
# highly variable run time (will take 0.25-6+ hours depending on # of features)
# to speed things up increase remove_per_step value

# [[ LinearSVC ]]
python3 ${script_dir}/select_features.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_3.pkl \
--group_split_method GroupKFold \
--outdir ${proj_dir}/results/feature_selection/ \
--threads 12 --remove_per_step 3 --num_splits 5 --seed 17 2>&1 | tee ${proj_dir}/logs/select_features_SVC.log

# [[ SGDClassifier ]]
python3 ${script_dir}/select_features.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_5.pkl \
--group_split_method GroupKFold \
--outdir ${proj_dir}/results/feature_selection/ \
--threads 12 --remove_per_step 3 --num_splits 5 --seed 17 2>&1 | tee ${proj_dir}/logs/select_features_SGD.log

# [[ ExtraTrees ]]
python3 ${script_dir}/select_features.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--model ${proj_dir}/results/tpot/gss/tpot_model_4.pkl \
--group_split_method GroupKFold \
--outdir ${proj_dir}/results/feature_selection/ \
--threads 12 --remove_per_step 3 --num_splits 5 --seed 17 2>&1 | tee ${proj_dir}/logs/select_features_ExtraTrees.log

# ---------------------------------------------------
# predict PPIs
# ---------------------------------------------------

# seed makes group splits, train/test splits, & model reproducible
# fairly quick (run time ~5-15min)

script_dir="/stor/work/Marcotte/project/rmcox/leca/scripts"
proj_dir="/stor/work/Marcotte/project/rmcox/leca/ppi_ml"

# with seed=777 & train_size=0.5, actual train/test split is 75.6/24.4

# [[ LinearSVC ]]
# generate predictions using all features
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_3.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/all/ \
--fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_LinearSVC_all_feats.log

# generate predictions across a sweep of features
for x in 5 10 25 50 100 250; do echo "python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_3.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/${x} \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_LinearSVC_RFECV_all.csv \
--num_feats ${x} --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_LinearSVC_${x}_feats.log"; done > cmds/sweep_feats_LinearSVC.sh
cat cmds/sweep_feats_LinearSVC.sh | parallel -j5


# [[ SGDClassifier ]]
# generate predictions using all features
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_5.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/all/ \
--fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_SGDClassifier_all_feats.log

# generate predictions across a sweep of features
for x in 5 10 25 50 100 250; do echo "python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_5.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/${x} \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_SGDClassifier_RFECV_all.csv \
--num_feats ${x} --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_SGDClassifier_${x}_feats.log"; done > cmds/sweep_feats_SGDClassifier.sh
cat cmds/sweep_feats_SGDClassifier.sh | parallel -j5


# [[ ExtraTrees ]]
# generate predictions using all features
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gss/tpot_model_4.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/all/ \
--fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_ExtraTreesClassifier_all_feats.log

# generate predictions across a sweep of features
for x in 5 10 25 50 100 250; do echo "python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gss/tpot_model_4.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/${x} \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_ExtraTreesClassifier_RFECV_all.csv \
--num_feats ${x} --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_ExtraTreesClassifier_${x}_feats.log"; done > cmds/sweep_feats_ExtraTreesClassifier.sh
cat cmds/sweep_feats_ExtraTreesClassifier.sh | parallel -j5

# ---------------------------------------------------
# run cross-validation
# ---------------------------------------------------

# [[ ExtraTrees ]]
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gss/tpot_model_4.pkl \
--outdir ${proj_dir}/results/cross_val \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_ExtraTreesClassifier_RFECV_all.csv \
--group_split_method GroupKFold --num_splits 5 \
--num_feats 25 --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/crossval_ExtraTreesClassifier_25_feats.log

# [[ SGDClassifier ]]
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_5.pkl \
--outdir ${proj_dir}/results/cross_val \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_SGDClassifier_RFECV_all.csv \
--group_split_method GroupKFold --num_splits 5 \
--num_feats 25 --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/crossval_SGDClassifier_25_feats.log

# [[ LinearSVC ]]
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_3.pkl \
--outdir ${proj_dir}/results/cross_val \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_LinearSVC_RFECV_all.csv \
--group_split_method GroupKFold --num_splits 5 \
--num_feats 100 --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/crossval_LinearSVC_100_feats.log

# ---------------------------------------------------
# cluster PPIs
# ---------------------------------------------------

# best models:
# ExtraTrees (100 feats) -> ppi_predict/feature_sweep/100/scored_interactions_fdr10_ExtraTreesClassifier.csv
# SGD (100 feats) -> ppi_predict/feature_sweep/100/scored_interactions_fdr10_SGDClassifier.csv

# seed makes results reproducible
# very fast (run time ~1min)

# [[ LinearSVC ]]
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/detect_communities.py --scores ${proj_dir}/results/ppi_predict/feature_sweep/100/scored_interactions_fdr10_LinearSVC.csv --steps ${x} --outfile ${proj_dir}/results/walktrap/step_sweep/LinearSVC_100feats_fdr10_${x}steps --annotations ${proj_dir}/annotations/leca_eunog_annots_complete.030721.csv --seed 777 2>&1 | tee ${proj_dir}/logs/walktrap_LinearSVC_100feats_fdr10_${x}steps.log"; done > cmds/sweep_walktrap_steps_LinearSVC_100feats_fdr10.sh
cat cmds/sweep_walktrap_steps_LinearSVC_100feats_fdr10.sh | parallel -j5

# calculate precision-recall per cut
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/calc_cluster_pr.py --cluster_file ${proj_dir}/results/walktrap/step_sweep/LinearSVC_100feats_fdr10_${x}steps.csv --positive_ppi_dict ${proj_dir}/data/featmats/positive_ppi_dict.pkl --negative_ppi_dict ${proj_dir}/data/featmats/negative_ppi_dict.pkl --outfile ${proj_dir}/results/walktrap/step_sweep/LinearSVC_100feats_fdr10_${x}steps_precision_recall.csv 2>&1 | tee ${proj_dir}/logs/walktrap_LinearSVC_100feats_fdr10_${x}steps_pr.log"; done > cmds/walktrap_pr_LinearSVC_100feats_fdr10.sh
cat cmds/walktrap_pr_LinearSVC_100feats_fdr10.sh | parallel -j5

# [[ ExtraTrees ]]
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/detect_communities.py --scores ${proj_dir}/results/ppi_predict/feature_sweep/25/scored_interactions_fdr10_ExtraTreesClassifier.csv --steps ${x} --outfile ${proj_dir}/results/walktrap/step_sweep/ExtraTreesClassifier_25feats_fdr10_${x}steps --annotations ${proj_dir}/annotations/leca_eunog_annots_complete.030721.csv --seed 777 2>&1 | tee ${proj_dir}/logs/walktrap_ExtraTreesClassifier_25feats_fdr10_${x}steps.log"; done > cmds/sweep_walktrap_steps_ExtraTreesClassifier_25feats_fdr10.sh
cat cmds/sweep_walktrap_steps_ExtraTreesClassifier_25feats_fdr10.sh | parallel -j5

# calculate precision-recall per cut
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/calc_cluster_pr.py --cluster_file ${proj_dir}/results/walktrap/step_sweep/ExtraTreesClassifier_25feats_fdr10_${x}steps.csv --positive_ppi_dict ${proj_dir}/data/featmats/positive_ppi_dict.pkl --negative_ppi_dict ${proj_dir}/data/featmats/negative_ppi_dict.pkl --outfile ${proj_dir}/results/walktrap/step_sweep/ExtraTreesClassifier_25feats_fdr10_${x}steps_precision_recall.csv 2>&1 | tee ${proj_dir}/logs/walktrap_ExtraTreesClassifier_25feats_fdr10_${x}steps_pr.log"; done > cmds/walktrap_pr_ExtraTreesClassifier_25feats_fdr10.sh
cat cmds/walktrap_pr_ExtraTreesClassifier_25feats_fdr10.sh | parallel -j5

# [[ SGDClassifier ]]
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/detect_communities.py --scores ${proj_dir}/results/ppi_predict/feature_sweep/25/scored_interactions_fdr10_SGDClassifier.csv --steps ${x} --outfile ${proj_dir}/results/walktrap/step_sweep/SGDClassifier_25feats_fdr10_${x}steps --annotations ${proj_dir}/annotations/leca_eunog_annots_complete.030721.csv --seed 777 2>&1 | tee ${proj_dir}/logs/walktrap_SGDClassifier_25feats_fdr10_${x}steps.log"; done > cmds/sweep_walktrap_steps_SGDClassifier_25feats_fdr10.sh
cat cmds/sweep_walktrap_steps_SGDClassifier_25feats_fdr10.sh | parallel -j5


# calculate precision-recall per cut
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/calc_cluster_pr.py --cluster_file ${proj_dir}/results/walktrap/step_sweep/SGDClassifier_25feats_fdr10_${x}steps.csv  --positive_ppi_dict ${proj_dir}/data/featmats/positive_ppi_dict.pkl --negative_ppi_dict ${proj_dir}/data/featmats/negative_ppi_dict.pkl --outfile ${proj_dir}/results/walktrap/step_sweep/SGDClassifier_25feats_fdr10_${x}steps_precision_recall.csv 2>&1 | tee ${proj_dir}/logs/walktrap_SGDClassifier_25feats_fdr10_${x}steps_pr.log"; done > cmds/walktrap_pr_SGDClassifier_25feats_fdr10.sh
cat cmds/walktrap_pr_SGDClassifier_25feats_fdr10.sh | parallel -j5

# ---------------------------------------------------
# sanity check the pipeline
# ---------------------------------------------------

# pretty quick (run time ~5min)
# will make a nice table in the log file
python3 ${script_dir}/sanity_checks.py \
--featmat_file ${proj_dir}/data/feature_matrix/feature_matrix_labeled.pkl \
--results_file ${proj_dir}/results/ppi_predict/all_feats/scored_interactions_all_ExtraTreesClassifier.csv 2>&1 | tee /project/rmcox/misc/vy.pipeline/logs/sanity_checks.log

# can also check group merge stats
python3 ${script_dir}/label_featmat.py \
--featmat ${proj_dir}/data/feature_matrix/feature_matrix.pkl \
--gold_std_file /project/vyqtdang/MS_brain/nextflow_brain/cfmsflow/gold_std/verNOG/CORUM4.1_core_merged.txt \
--outfile_name ${proj_dir}/data/feature_matrix/feature_matrix_labeled_allgroups \
--keep_cmplx_overlap --seed 13 2>&1 | tee ${proj_dir}/logs/label_featmat_allgroups.log

python3 ${script_dir}/sanity_checks.py \
--featmat_file ${proj_dir}/data/feature_matrix/feature_matrix_labeled_allgroups_traintest.pkl \
--results_file ${proj_dir}/results/ppi_predict/all_feats/scored_interactions_all_ExtraTreesClassifier.csv 2>&1 | tee ${proj_dir}/logs/sanity_checks_groups.log


# ---------------------------------------------------
# final model
# ---------------------------------------------------


# ---------------------------------------------------
# calculate enrichment for external ppis
# ---------------------------------------------------

script_dir="/stor/work/Marcotte/project/rmcox/leca/scripts"
proj_dir="/stor/work/Marcotte/project/rmcox/leca/ppi_ml"

# template
python3 ${script_dir}/calc_enrichment.py \
--int_file ${proj_dir}/results/ppi_predict/feature_sweep/100/scored_interactions_all_LinearSVC.csv \
--ext_file ${proj_dir}/data/validation/y2h_fmt.csv \
--int_filter_file ${proj_dir}/annotations/lists/leca_nogs.txt \
--bin_width 10000 \
--outfile_name ${proj_dir}/results/ext_val/y2h_ext_enrichment_10k_bw.csv

# create commands
for exp in xlms y2h cx; do
        echo "python3 ${script_dir}/calc_enrichment.py \
--int_file ${proj_dir}/results/ppi_predict/feature_sweep/100/scored_interactions_all_LinearSVC.csv \
--ext_file ${proj_dir}/data/validation/${exp}_fmt.csv \
--int_filter_file ${proj_dir}/annotations/lists/leca_nogs.txt \
--bin_width 10000 \
--outfile_name ${proj_dir}/results/ext_val/${exp}_ext_enrichment_10k_bw.csv"
done > cmds/calc_ext_enrich_10k.sh

for exp in xlms y2h cx; do
        echo "python3 ${script_dir}/calc_enrichment.py \
--int_file ${proj_dir}/results/ppi_predict/feature_sweep/100/scored_interactions_all_LinearSVC.csv \
--ext_file ${proj_dir}/data/validation/${exp}_fmt.csv \
--int_filter_file ${proj_dir}/annotations/lists/leca_nogs.txt \
--bin_width 1000 \
--outfile_name ${proj_dir}/results/ext_val/${exp}_ext_enrichment_1k_bw.csv"
done > cmds/calc_ext_enrich_1k.sh

for exp in xlms y2h cx; do
        echo "python3 ${script_dir}/calc_enrichment.py \
--int_file ${proj_dir}/results/ppi_predict/feature_sweep/100/scored_interactions_all_LinearSVC.csv \
--ext_file ${proj_dir}/data/validation/${exp}_fmt.csv \
--int_filter_file ${proj_dir}/annotations/lists/leca_nogs.txt \
--bin_width 45000 \
--outfile_name ${proj_dir}/results/ext_val/${exp}_ext_enrichment_45k_bw.csv"
done > cmds/calc_ext_enrich_45k.sh


# ::::::::::::::::::::::::::::::::::::::::::
# Crosslinking analysis
# ::::::::::::::::::::::::::::::::::::::::::

# map custom annotated fasta to eggnog groups
/stor/work/Marcotte/project/rmcox/programs/eggnog-mapper-2.0.5/emapper.py -i xlink/data/fastas/tett_marcotte_annotations_fixed.fasta --output xlink/data/tetts_custom.diamond --no_file_comments -m diamond --cpu 16

# format results
python3 ../scripts/format_emapper_output.py -i data/eggnog/tetts_custom.diamond.emapper.annotations -o data/eggnog/tetts_custom.diamond.mapping.2759 -l 2759


# ::::::::::::::::::::::::::::::::::::::::::
# Disease network propagation
# ::::::::::::::::::::::::::::::::::::::::::

proj_dir="/stor/work/Marcotte/project/rmcox/leca/human_disease/network_propagation"
python3 ../scripts/network_propagation.py --ppi_network ${proj_dir}/data/leca_map.csv --disease_network ${proj_dir}/data/disease_network.tsv --annotations ${proj_dir}/data/leca_id_map.tsv --outfile_name ${proj_dir}/results/leca_diseases_leave1out_auroc.tsv