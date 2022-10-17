import csv
import argparse
import time
from time import sleep
from tqdm import tqdm
from pathlib import Path

def main():

	clade_list = ["animals", "plants", "tsar", "excavate"]
	animals = ["brart", "caeel", "dicdi", "drome", "human", "mouse", \
	"nemve", "pig", "strpu", "xenla", "yeast"]
	plants = ["arath", "braol", "cansa", "cerri", "cheqi", "chlre", \
	"cocnu", "maize", "orysj", "selml", "sollc", "soybn", "wheat"]
	tsar = ["phatc", "plaba", "plaf7", "plafo", "plakh", "tetts"]
	excavate = ["euggr", "tryb2"]

	print("Formating input for clade information ...")
	clade_dict = dict()
	#nested_defaults = {'animals': 0, 'animals_max': 0, 'plants': 0, 'plants_max': 0, 'tsar': 0, 'tsar_max': 0, 'excavate': 0, 'excavate_max': 0}
	with open(corr_tbl, 'r') as f_in:
		freader = csv.reader(f_in, delimiter=",")
		next(freader)
		for row in tqdm(freader):
			info = row[0]
			species, exp = info.split(".", 1)
			ids = row[1]+' '+row[2]
			corr = abs(float(row[3]))

			if ids not in clade_dict.keys():
				clade_dict[ids] = {'animals': 0, 'animals_max': 0, 'plants': 0, 'plants_max': 0, 'tsar': 0, 'tsar_max': 0, 'excavate': 0, 'excavate_max': 0}

			if species in animals:
				clade_dict[ids]['animals'] = 1
				if corr > clade_dict[ids]['animals_max']:
					clade_dict[ids]['animals_max'] = corr
			elif species in plants:
				clade_dict[ids]['plants'] = 1
				if corr > clade_dict[ids]['plants_max']:
					clade_dict[ids]['plants_max'] = corr
			elif species in tsar:
				clade_dict[ids]['tsar'] = 1
				if corr > clade_dict[ids]['tsar_max']:
					clade_dict[ids]['tsar_max'] = corr
			elif species in excavate:
				clade_dict[ids]['excavate'] = 1
				if corr > clade_dict[ids]['excavate_max']:
					clade_dict[ids]['excavate_max'] = corr

			# eventually need to convert the code above to a nice function
			# for clade in clade_list:
			# 	if species in clade:
			# 		clade_dict[ids][clade] = 1
			# 		if corr > clade_dict
	
	print("Filtering for ID1-ID2 present in >= {} clades ...".format(threshold))
	final_dict = {}
	for ids in clade_dict:
		clade_sum = clade_dict[ids]['animals'] + clade_dict[ids]['plants'] + clade_dict[ids]['tsar'] + clade_dict[ids]['excavate']
		if clade_sum >= threshold:
			final_dict[ids] = clade_dict[ids]

	print("Writing results to file {} ...".format(writefile))
	fields = ['ID', 'animals', 'animals_max', 'plants', 'plants_max', 'tsar', 'tsar_max', 'excavate', 'excavate_max']
	with open(writefile, "w") as f_out:
		writer = csv.DictWriter(f_out, fields)
		writer.writeheader()
		for key in final_dict:
			final_dict[key]['ID'] = key
			writer.writerow(final_dict[key])

	print("Done!")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filters a table for PPIs observed in at least 2 clades.")
    parser.add_argument("-f", "--file", action="store", required=True,
                                        help="Comma-separated feature file containing 4 columns (1st column = experiment info, 2nd col = ID1, 3rd col = ID2, 4th col = Pearson's R). File should not have headers.")
    parser.add_argument("-t", "--threshold", action="store", required=True, help="Number of clades to require.")
    parser.add_argument("-o", "--out", action="store", required=False,
                                        help="Name for outfile (optional).")

    args = parser.parse_args()
    if args.out == None:
        filepath = Path(args.file)
        filename = str(filepath.with_suffix(""))
        writefile = filename + ".clade_filt" + filepath.suffixes[-1]
    else:
        writefile = args.out
    corr_tbl = args.file
    threshold = float(args.threshold)
    start_time = time.time()
    main()
    print("Runtime = {} seconds.".format(time.time() - start_time))