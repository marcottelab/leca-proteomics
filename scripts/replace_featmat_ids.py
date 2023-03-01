import csv
import argparse
import time
from pathlib import Path

def main():
    print("Reading in protein ID <-> eggNOG ID file ...")
    id_dict = {}
    with open(args.orthomap, 'r') as f:
        for mapping in csv.reader(f, delimiter="\t"):
            proteinID = mapping[0]
            eggnogID = mapping[1]
            id_dict[proteinID] = eggnogID

    print("Converting protein IDs to eggNOG IDs ...")
    with open (args.featmat, "r") as f_in:
        with open(writefile, "w") as f_out:
            for line in f_in:
                id1, id2, feats = line.split(",", 2)
                id1_new = (id_dict.get(id1) if id1 in id_dict else id1)
                id2_new = (id_dict.get(id2) if id2 in id_dict else id2)
                new_line = "{},{},{}".format(id1_new, id2_new, feats)
                f_out.write(new_line)

    print("Done!")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Replaces a set of ID1-ID2 accessions in a feature matrix with their corresponding eggNOG1-eggNOG2 accessions.")
    parser.add_argument("-m", "--orthomap", action="store", required=True,
                                        help="Tab-separated file mapping protein accessions to eggNOG groups where \
                                        the columns are labeled 'ID' for the eggNOG groups and 'ProteinID' for the protein accession. \
                                        (This is the output from the format_emapper script).")
    parser.add_argument("-f", "--featmat", action="store", required=True,
                                        help="Comma-separated feature file containing ID1-ID2 accessions as the first two columns")
    parser.add_argument("-o", "--out", action="store", required=False,
                                        help="Name for outfile (optional)")

    args = parser.parse_args()
    if args.out == None:
        filepath = Path(args.featmat)
        filename = str(filepath.with_suffix(""))
        writefile = filename + ".eggnogIDs" + filepath.suffixes[-1]
    else:
        writefile = args.out
    start_time = time.time()
    main()
    print("Runtime = {} seconds.".format(time.time() - start_time))





