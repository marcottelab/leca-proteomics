import pandas
import pickle
import numpy

leftfile = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/features/featmat_allexps_p3c2.pkl'
fmatfile = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/data/featmats/featmat.pkl'
outfile = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml/debug/debug_ljoin_ids.txt'

def make_fset(x):
    if len(set(x.split(' '))) < 2:
        print(f"WARNING: Correlation metrics for '{x}' (self-self PPI) detected; make sure you mean for this to be included ...")
        x1 = x.split(' ')[0]
        fset = frozenset({x1})
        return(fset)
    else:
        x1 = x.split(' ')[0]
        x2 = x.split(' ')[1]
        fset = frozenset({x1,x2})
        return(fset)

def difference(list1, list2):
    return list(set(list1) ^ set(list2))

print(f'Loading {fmatfile} ...')
with open(fmatfile, 'rb') as handle:
    fmat = pickle.load(handle)
print(f'Loaded feature matrix shape: {fmat.shape}')
fmat_ids = [make_fset(i) for i in fmat['ID']]
counts = fmat.groupby(['ID']).size().sort_values(ascending=False)
print(counts)

# ok so there are duplicate rows
# idk why
# so theoretically we can just drop them ..?
fmat_uniq = fmat.drop_duplicates()
print(f'De-duplicated feature matrix shape: {fmat_uniq.shape}')


# print(f'Loading {leftfile} ...')
# with open(leftfile, 'rb') as handle:
#     left_df = pickle.load(handle)
# print(f'Left join seed data frame shape: {left_df.shape}')
# left_ids = [make_fset(i) for i in left_df['ID']]
# del left_df

# print(f'Finding differences in the ID column ...')
# diff = difference(fmat_ids, left_ids)

# print(f'Writing results to {outfile} ...')
# with open(outfile, 'w') as f:
#     f.write('\n'.join(map(str, diff)))
    
print('Done!')
