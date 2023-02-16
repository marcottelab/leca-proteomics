import pickle
import numpy as np
import datetime as dt
from tpot import TPOTClassifier
from sklearn.model_selection import GroupShuffleSplit

# future argparse vars, probably
tpot_dir = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/tpot/'
timestamp = dt.datetime.now().strftime('%Y-%m-%d')
outfile = tpot_dir+'tpot_pipeline'
fmatfile = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/featmats/featmat_labeled_traintest.pkl'
train_size = 0.7
num_splits = 5
seed = 13

# define auto ML objects
gss = GroupShuffleSplit(n_splits = num_splits, train_size=train_size, random_state=seed)
pipeline_opt = TPOTClassifier()
pipeline_opt = TPOTClassifier(generations=5, population_size=20, cv=5, random_state=13, verbosity=2)

# read in data
print(f'Reading in {fmatfile} ...')
with open(fmatfile, 'rb') as handle:
    fmat = pickle.load(handle)
    
# define cols
label_cols = ['ID', 'label', 'super_group']
data_cols = [c for c in fmat.columns.values.tolist() if c not in label_cols]

# make data, target, and group arrays
print(f'Formatting arrays ...')
X = fmat[data_cols].to_numpy()
y = fmat[label_cols[1]].to_numpy()
groups = fmat[label_cols[2]].to_numpy()

# get gss splits for each iteration
print(f'Getting test/train splits (# splits = {num_splits}) ...')
for i, (test_idx, train_idx) in enumerate(gss.split(X, y, groups)):
    
    X_train = X[train_idx]
    y_train = y[train_idx]
    X_test = X[test_idx]
    y_test = y[test_idx]
    
    print(f'Running TPOT for split #{i+1}...')
    print(f"--> # train PPIs = {len(X[train_idx])}")
    print(f"--> # test PPIs = {len(X[test_idx])}")
    
    pipeline_opt.fit(X_train, y_train) 
    print(pipeline_opt.score(X_test, y_test))
    pipeline_opt.export(outfile+'_'+str(i+1))

# for train_idx, test_idx in gss.split(X, y, groups):
#     X_train = X[train_idx]
#     y_train = y[train_idx]
#     X_test = X[test_idx]
#     y_test = y[test_idx]
    
# get optimized pipeline
# print(f'Running TPOT ...')
# pipeline_opt.fit(X_train, y_train)
# print(pipeline_opt.score(X_test, y_test))
# pipeline_opt.export(outfile)

