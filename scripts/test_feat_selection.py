import pandas  as pd
import numpy as np
import pickle
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import RFECV
from sklearn.metrics import PrecisionRecallDisplay
import matplotlib.pyplot as plt
from functools import reduce

# future argparse vars, probably
fmatfile = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/featmats/featmat_labeled_traintest.pkl'
out_dir = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/results/feature_selection/'
seed = 13
train_size = 0.7
num_splits = 5
min_features_to_select = 1

# load data
print(f'Reading in {fmatfile} ...')
with open(fmatfile, 'rb') as handle:
    fmat = pickle.load(handle)

# define cols
label_cols = ['ID', 'label', 'super_group']
data_cols = [c for c in fmat.columns.values.tolist() if c not in label_cols]

# convert cols to arrays
print(f'Formatting arrays ...')
X = fmat[data_cols].to_numpy()
y = fmat[label_cols[1]].to_numpy()
groups = fmat[label_cols[2]].to_numpy()

# group split strategy
gss = GroupShuffleSplit(n_splits = num_splits, train_size=train_size, random_state=seed)

# define model to train on
clf = ExtraTreesClassifier(n_estimators=100, random_state=seed)
print(f'Selected model: {clf}')
print(f"Total # positive labels: {len(fmat[fmat['label']==1])}")
print(f"Total # negative labels: {len(fmat[fmat['label']==-1])}")

# feature selector (rfe w/ cross-validation)
cv = StratifiedKFold(5)
rfecv = RFECV(
    estimator=clf,
    step=7,
    cv=cv,
    scoring="accuracy",
    min_features_to_select=min_features_to_select,
    n_jobs=5,
)

# get gss splits for each iteration
print(f'----- Running recursive feature elimination for {num_splits} GSS splits -----')
df_list = []
# NOTE: i am getting some weird inconsistent train set sizes with GSS
# need to look into this further
for i, (test_idx, train_idx) in enumerate(gss.split(X, y, groups)):
    
    # define test/train splits
    X_train = X[train_idx]
    y_train = y[train_idx]
    X_test = X[test_idx]
    y_test = y[test_idx]
    
    # define output vars
    suffix = 'gssfold'+str(i+1)
    outname = f'sel_feats_xtrees_{suffix}'
    
    # run feature selection
    print(f'Executing RFE for GSS split #{i+1} ...')
    print(f"Total # train PPIs = {len(X[train_idx])}")
    print(f"Total # test PPIs = {len(X[test_idx])}")
    
    rfecv.fit(X_train, y_train)
    print(f"Optimal number of features: {rfecv.n_features_}")
    
    # plot results
    print(f'Generating feature selection evaluation plots ...')
    
    # mean test accuracy vs number of feats selected
    n_scores = len(rfecv.cv_results_["mean_test_score"])
    plt.figure()
    plt.xlabel(f"Number of features selected\n(optimal={rfecv.n_features_} features)")
    plt.ylabel("Mean CV test accuracy")
    plt.errorbar(
        range(min_features_to_select, n_scores + min_features_to_select),
        rfecv.cv_results_["mean_test_score"],
        yerr=rfecv.cv_results_["std_test_score"],
    )
    plt.title(f"Recursive feature elimination\nwith correlated features (fold #{i+1})")
    print(f'Saving plot to {out_dir+outname}_cv-test_nfeats-vs-acc.png ...')
    plt.savefig(f'{out_dir+outname}_cv-test_nfeats-vs-acc.png', dpi=300, transparent=True)
    
    # PR curve
    from sklearn.metrics import PrecisionRecallDisplay
    PrecisionRecallDisplay.from_estimator(rfecv, X_test, y_test)
    plt.savefig(f'{out_dir+outname}_holdout-test_prcurve.png', dpi=300, transparent=True)
    print(f'Saving plot to {out_dir+outname}_holdout-test_prcurve.png ...')
    
    # get result table
    print(f'Extracting selected features & scores ...')
    results = pd.DataFrame({'feature':data_cols, 'rank':rfecv.ranking_, 'support':rfecv.support_})
    sel_feats = results[results['support'] == True]
    sel_feats_scored = sel_feats.head(rfecv.n_features_)
    sel_feats_scored.drop(['support'], axis=1, inplace=True)
    sel_feats_scored.drop(['rank'], axis=1, inplace=True)
    sel_feats_scored['mdi'] = rfecv.estimator_.feature_importances_
    sel_feats_scored = sel_feats_scored.sort_values('mdi')
    sel_feats_scored['fold'] = int(i+1)
    df_list.append(sel_feats_scored)
    
    # write output
    print(f'Writing results to {out_dir+outname}.csv ...')
    sel_feats_scored.to_csv(f'{out_dir+outname}.csv', index=False)

print(f'Getting selected feature representation & summary stats across all folds ...')
all_res = pd.concat(df_list)
gb = all_res.groupby(['feature'])
counts = gb.size().to_frame(name='counts')
agg_res = (counts
           .join(gb.agg({'fold': lambda x: ', '.join(set(x.astype(str).dropna()))}).rename(columns={'fold': 'gss_folds'}))
           .join(gb.agg({'mdi': 'mean'}).rename(columns={'mdi': 'mean_mdi'}))
           .join(gb.agg({'mdi': 'min'}).rename(columns={'mdi': 'min_mdi'}))
           .join(gb.agg({'mdi': 'max'}).rename(columns={'mdi': 'max_mdi'}))
           .sort_values(['counts', 'mean_mdi'], ascending=[False, False])
           .reset_index()
          )
feat_intxn = agg_res[agg_res['counts'] == num_splits]
print(f'# common features: {len(feat_intxn)}')
print(f'Writing results to {out_dir} ...')
agg_res.to_csv(f'{out_dir}sel_feats_xtrees_allres.csv', index=False)
feat_intxn.to_csv(f'{out_dir}sel_feats_xtrees_intxn.csv', index=False)