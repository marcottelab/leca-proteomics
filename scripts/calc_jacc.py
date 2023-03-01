import csv
import argparse
import numpy as np

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute Jaccard similarity for gold standard protein complexes, merging where J > 0.6. Outputs merged complex file.")
    parser.add_argument("-i", "--input_file", action="store", required=True,
                                        help="File ('.txt') containing new line-separated protein complexes, with 1 protein complex per line where the subunits are space separated")
    parser.add_argument("-o", "--out", action="store", required=False,
                                        help="Name for outfile of merged complexes (optional)")
    args = parser.parse_args()

    input_file = args.input_file

class Counter(object) :
    def __init__(self, fun) :
        self._fun = fun
        self.counter=0
    def __call__(self,*args, **kwargs) :
        self.counter += 1
        print("# of iterations = {}".format(self.counter))
        return self._fun(*args, **kwargs)

def jaccard(a, b):
	intersection = len(list(set(a).intersection(set(b))))
	union = len(a) + len(b) - intersection
	jacc = float(intersection)/union
	return(jacc)

def merge(a, b):
	cmplx_1 = set(a)
	cmplx_2 = set(b)
	unique_mems = list(cmplx_2 - cmplx_1)
	merged_cmplx = a + unique_mems
	return(merged_cmplx)

def set_outname(out):
	if args.out == None:
	    new_ext = '.merged.txt'
	    writefile = new_ext.join(input_file.rsplit('.txt', 1))
	    print("Outfile name: {}".format(writefile))
	else:
	    writefile = args.out
	    print("Outfile name: {}".format(writefile))
	return(writefile)

def read_cmplx(input_file):
	cmplx_file = open(input_file).read().splitlines()

	cmplx_list = []
	for string in cmplx_file:
		cmplx = string.strip('"').split(' ')
		cmplx_list.append(cmplx)
	return(cmplx_list)

def eval_cmplx(cmplx_list):
	to_remove = []
	merged_list = []
	merged_count = 0
	for i in range(len(cmplx_list)):
		for j in range(i + 1, len(cmplx_list)):
			c1 = cmplx_list[i]
			c2 = cmplx_list[j]
			jacc = jaccard(c1, c2)
			#print('\nComplex1:{}\nComplex2:{}\nJaccard:{}'.format(c1, c2, jacc))
			if jacc > 0.6:
				to_remove.append(c1)
				to_remove.append(c2)
				merged_cmplx = merge(c1, c2)
				if merged_cmplx not in merged_list:
					#print('The merged complex is:{}\n'.format(merged_cmplx))
					merged_list.append(merged_cmplx)
					merged_count = merged_count + 1
			else:
				continue

	final_list = [x for x in cmplx_list if x not in to_remove]
	#print("Complexes not merged: {}".format(final_list))

	for cmplx in merged_list:
		final_list.append(cmplx)
	#print("Final list: {}".format(final_list))
	#print("# in final list: {}".format(len(final_list)))

	final_unique = [list(x) for x in set(tuple(x) for x in final_list)]
	#final_unique = np.unique(np.array(final_list))
	#print("Final list, unique: {}".format(final_unique))
	#print("# in final list, unique: {}".format(len(final_unique)))

	#print('-------------------')
	print('TOTAL # MERGED = {}'.format(merged_count))
	print('-------------------')
	
	# recursively merge complexes until they converge
	if merged_count > 0:
		eval_cmplx(final_unique)
	else:
		print("Complex similarity has converged; writing results to {}.".format(writefile))
	
	return(final_unique)

eval_cmplx = Counter(eval_cmplx)

writefile = set_outname(args.out)
cmplx_list = read_cmplx(args.input_file)
print("Merging complexes with Jaccard coefficient > 0.6...")
merged = eval_cmplx(cmplx_list)

with open(writefile, "w") as outfile:
	for cmplx in merged:
		cmplx_fmt = " ".join(map(str, cmplx))
		outfile.write('{}\n'.format(cmplx_fmt))

# if __name__ == '__main__':
#     main()