import pickle
import csv
from tqdm import tqdm

apms_fmat = "/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/apms/featmat_humap"
apms_dict = dict()
with open(apms_fmat, "r") as f:
    reader = csv.reader(f, delimiter=",")
    headers = next(reader, None)
    for line in tqdm(reader, total=12321600):
        pair = frozenset(line[0].split(' '))
        scores = []
        scores = [line[1:]]
        print(scores)
        if pair not in apms_dict.keys():
            apms_dict.update({pair: scores})

with open('/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/apms/apms_scores_dict.pkl', 'wb') as handle:
    pickle.dump(apms_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
