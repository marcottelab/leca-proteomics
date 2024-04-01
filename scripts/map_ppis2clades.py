''' 
lazy script for parallelizing the creation of a ppi:phyla dictionary 
'''

import sys
import pandas as pd
import pickle
import re

''' functions '''
def make_fset(x, drop=True):
    if len(set(x.split(' '))) < 2:
        print(f"WARNING: Features for '{x}-{x}' (self-self PPI) detected ...")
        if drop == False:
            x1 = x.split(' ')[0]
            fset = frozenset({x1,x1})
            return(fset)
        else:
            return(None)
    else:
        x1 = x.split(' ')[0]
        x2 = x.split(' ')[1]
        fset = frozenset({x1,x2})
        return(fset)
    
def make_clade_dict(clade_file):
    spec_df = pd.read_csv(clade_file)
    spec_df['code'] = [i.lower() for i in spec_df['code']]
    spec_df['clade'] = ['Excavata' if i == 'Excavate' else i for i in spec_df['clade']]
    clade_dict = dict(zip(spec_df['code'], spec_df['clade']))
    return(clade_dict)

def make_phyla_dict(fmat, clade_dict):
    # map ppis to phyla
    phyla_dict = dict()
    for i in range(len((fmat))):
        ppi_id = fmat['fs'][i]
        # get non-zero features
        feats = list(fmat[fmat['fs']==ppi_id].dropna(axis=1,how='all'))[1:]
        species_set = set()
        clade_set = set()
        # get species names from features
        for f in feats:
            code = f.split('.', 1)[0]
            if code != 'fs':
                species_set.add(code)
                clade_set.add(clade_dict[code])
        clades_sorted = sorted(list(clade_set))
        if 'Archaeplastida' in clade_set:
            clades_sorted.remove('Archaeplastida')
            clades_sorted.append('Archaeplastida')
        if ppi_id not in phyla_dict.keys():
            phyla_dict[ppi_id] = {'n_phyla':{}, 'phyla':{}}
            phyla_dict[ppi_id]['n_phyla'] = len(clade_set)
            phyla_dict[ppi_id]['phyla'] = ', '.join(clades_sorted)
    return(phyla_dict)



''' args '''
chunk_file = sys.argv[1]
out_dict = sys.argv[2]


''' main '''
fmat = pd.read_pickle(chunk_file)
proj_dir = '/stor/work/Marcotte/project/rmcox/leca/ppi_ml'
clade_dict = make_clade_dict(f'{proj_dir}/data/meta/speciesinfo_clades.csv')
phyla_dict = make_phyla_dict(fmat, clade_dict)
with open(out_dict, 'wb') as file:
    pickle.dump(phyla_dict, file)