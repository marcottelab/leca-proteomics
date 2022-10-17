import csv
import argparse
import time
from pathlib import Path

def main():
	print("Finding the minimum and maximum for each feature ...")
	with open(matrix, 'r') as f_in:
		reader = csv.reader(f_in,delimiter=',')
		feats = list(next(reader))
		print("Feature names:", feats)
		nfeats = len(next(reader))
		print("Number of features:", nfeats)
		maxlist = [float('-inf') for x in range(nfeats)]
		print("Empty 'max' array:\n", maxlist)
		minlist = [float('inf') for x in range(nfeats)]
		print("Empty 'min' array:\n", maxlist)
		next(reader, None)
		for row in reader:
			i = 0
			for value in row:
				if float(value) < minlist[i]:
					minlist[i] = float(value)
				if float(value) > maxlist[i]:
					maxlist[i] = float(value)
				i += 1
		print("Length of max values:", len(maxlist))
		print("Populated 'max' array", maxlist)
		print("Length of min values:", len(minlist))
		print("Populated 'min' array", minlist)

	print("Writing results to {} ...".format(writefile))
	with open(writefile, "w") as f_out:
		writer = csv.writer(f_out)
		writer.writerow(['feature', 'min', 'max'])
		values = zip(feats, minlist, maxlist)
		writer.writerows(values)

	print("Done!")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Finds min/max of each feature in a feature matrix.")
    parser.add_argument("-f", "--file", action="store", required=True,
                                        help="Comma-separated feature file containing.")
    parser.add_argument("-o", "--out", action="store", required=False,
                                        help="Name for outfile (optional).")

    args = parser.parse_args()
    if args.out == None:
        filepath = Path(args.file)
        filename = str(filepath.with_suffix(""))
        writefile = filename + ".minmax" + filepath.suffixes[-1]
    else:
        writefile = args.out
    matrix = args.file
    start_time = time.time()
    main()
    print("Runtime = {} seconds.".format(time.time() - start_time))

