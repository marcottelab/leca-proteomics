import csv
import argparse
import time
from time import sleep
from tqdm import tqdm
from pathlib import Path

def main():
    print("Filtering tidy feature table for correlations >= {} ...".format(threshold))
    with open(corr_tbl, 'r') as f_in:
        with open(writefile, "w") as f_out:
            freader = csv.reader(f_in, delimiter=",")
            next(freader)
            for row in tqdm(freader):
                #sleep(0.00001)
                proteins = row[0]
                info = row[1]
                species, exp = info.split(".", 1)
                corr = abs(float(row[2]))
                if corr >= threshold:
                    new_line = "{},{},{}\n".format(proteins, species, corr)
                    f_out.write(new_line)
    print("Done!")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filters a tidy feature table based on a given threshold.")
    parser.add_argument("-f", "--file", action="store", required=True,
                                        help="Comma-separated feature file containing columns \
                                        'ID', 'exp', 'correlation'.")
    parser.add_argument("-t", "--threshold", action="store", required=True, help="Threshold to filter at.")
    parser.add_argument("-o", "--out", action="store", required=False,
                                        help="Name for outfile (optional).")

    args = parser.parse_args()
    if args.out == None:
        filepath = Path(args.file)
        filename = str(filepath.with_suffix(""))
        writefile = filename + ".corr_filt" + filepath.suffixes[-1]
    else:
        writefile = args.out
    corr_tbl = args.file
    threshold = float(args.threshold)
    start_time = time.time()
    main()
    print("Runtime = {} seconds.".format(time.time() - start_time))
