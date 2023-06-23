# ::::::::::::::::::::::::::::::::::::::::::
# eggNOG mapper for orthology assigment
# ::::::::::::::::::::::::::::::::::::::::::

## -----------------------
## step 0: install emapper
## -----------------------
wget https://github.com/eggnogdb/eggnog-mapper/archive/refs/tags/2.0.5.tar.gz
tar -xvxf 2.0.5.tar.gz
cd eggnog-mapper-2.0.5

# download eggNOG taxa and diamond databases
# answer y to all questions (or use -y to do that automatically)
python3 download_eggnog_data.py --data_dir data

# download hmmer databases (takes a long time):
python3 download_eggnog_data.py --data_dir data -H -d 2759 # leca
python3 download_eggnog_data.py --data_dir data -H -d 40674 # mammalia
python3 download_eggnog_data.py --data_dir data -H -d 4751 # fungi

# !!! important !!!
# -> when asked to enter a "non-empty name" for the hmmer db, let the script default to the tax ID number or enter the tax ID number yourself; anything else will result in problems
# -> hmmer mode specifically requires python 3.10.8

## -----------------------
## step 1: run emapper
## -----------------------

# human (diamond)
python3 eggnog-mapper-2.0.5/emapper.py -i test/human.fasta --output test/human.diamond --no_file_comments -m diamond --cpu 16
# human (hmmer)
python3 eggnog-mapper-2.0.5/emapper.py -i test/human.fasta --output test/human.hmmer --no_file_comments -m hmmer -d 40674 --cpu 16

# yeast (diamond)
python3 eggnog-mapper-2.0.5/emapper.py -i test/yeast.fasta --output test/yeast.diamond --no_file_comments -m diamond --cpu 16
# yeast (hmmer)
python3 eggnog-mapper-2.0.5/emapper.py -i test/yeast.fasta --output test/yeast.hmmer --no_file_comments -m hmmer -d 4751 --cpu 16


## -----------------------
## step 2: format results
## -----------------------

# get scripts from git
git clone https://github.com/marcottelab/leca-proteomics.git
# this has a billion things in it (sorry ¯\_(ツ)_/¯) but for eggnog stuff you'll only need 2 scripts: format_emapper_output.py and bin_orthogroups.py

# both diamond and hmmer mode will output a ".annotations" file which is used as input to the format_emapper_output script below

python3 leca-proteomics/scripts/format_emapper_output.py -i test/human.diamond.emapper.annotations -o test/human.diamond.mapping.40674 -l 40674 # extract mammalia-level assignments

python3 leca-proteomics/scripts/format_emapper_output.py -i test/yeast.hmmer.emapper.annotations -o test/yeast.hmmer.mapping.40674 -l 4751 # extract fungi-level assignments

# this script also has a --format_uniprot_id arg, which you can use to return a formatted ID column that e.g. converts "sp|A2P2R3|YM084_YEAST" to "A2P2R3"; it's useful if you're doing annotation table joins

## -----------------------
## step 3: bin proteome
## -----------------------

# the format_emapper_script outputs a ".mapping" file, which is used as input to the bin_orthogroups script below
# ensure the ProteinID column in the .mapping files matches the ID format of the original input FASTA file you used in step 1

python3 leca-proteomics/scripts/bin_orthogroups.py -m test/human.diamond.mapping.40674 -f test/human.fasta -o test/human.diamond.40674.collapsed.fasta --limit_length

python3 leca-proteomics/scripts/bin_orthogroups.py -m test/yeast.hmmer.mapping.4751 -f test/yeast.fasta -o test/yeast.hmmer.4751.collapsed.fasta --limit_length

# --limit_length caps binned FASTA entry length at 100,000aa;
# if a collapsed orthogroup exceeds this, it will be broken back out into it's individual original FASTA entries
# good example: xenopus olfactory protein family
# it was added because proteome discoverer breaks on FASTA entries larger than this
# even if i'm not using PD, i always keep this enabled

