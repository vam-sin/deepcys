# libraries
import pandas as pd 
from sklearn.utils import resample
from sklearn.preprocessing import LabelEncoder
import numpy as np

ds = pd.read_excel('../data/new_ds_data.xlsx')

ds.columns = ['no', 'pdb', 'res', 'ch', 'pka', 'bf', 'rhpy', 'mod']

encoder = LabelEncoder()
encoder.fit(ds['mod'])
ds['mod'] = encoder.transform(ds['mod'])
print(encoder.inverse_transform([0,1,2,3]))

counts = ds['mod'].value_counts()
counts = list(counts)
number_of_samples = min(counts)
print("Number of Samples ", number_of_samples)
''' Biased highly towards Thioether
thioether         3244
Disulphide        2600
metal-binding     2342
Sulphenylation    2315

3    3244
0    2600
2    2342
1    2315

Name: mod, dtype: int64
[0,1,2,3]: ['Disulphide' 'Sulphenylation' 'metal-binding' 'thioether']
y:
thio
dis
mb
sul
'''

ds_majority1 = ds[ds['mod'] == 3]
ds_majority2 = ds[ds['mod'] == 0]
ds_majority3 = ds[ds['mod'] == 2]
ds_majority4 = ds[ds['mod'] == 1]

# Undersample

ds_majority1 = resample(ds_majority1, 
                                 replace=False,    # sample without replacement
                                 n_samples=number_of_samples,     # to match minority class
                                 random_state=42) # reproducible results

ds_majority2 = resample(ds_majority2, 
                                 replace=False,    # sample without replacement
                                 n_samples=number_of_samples,     # to match minority class
                                 random_state=42) # reproducible results

ds_majority3 = resample(ds_majority3, 
                                 replace=False,    # sample without replacement
                                 n_samples=number_of_samples,     # to match minority class
                                 random_state=42) # reproducible results

ds_majority4 = resample(ds_majority4, 
                                 replace=False,    # sample without replacement
                                 n_samples=number_of_samples,     # to match minority class
                                 random_state=42) # reproducible results

ds_downsampled = pd.concat([ds_majority1, ds_majority2, ds_majority3, ds_majority4])

# print("ds_downsampled")
# print(ds_downsampled['mod'].value_counts())

pdb = []
for i in ds_downsampled['pdb']:
	pdb.append(i)
pdb = pd.DataFrame(pdb)

res = []
for i in ds_downsampled['res']:
	res.append(i)
res = pd.DataFrame(res)

ch = []
for i in ds_downsampled['ch']:
	ch.append(i)
ch = pd.DataFrame(ch)

pka = []
for i in ds_downsampled['pka']:
	pka.append(i)
pka = pd.DataFrame(pka)

bf = []
for i in ds_downsampled['bf']:
	bf.append(i)
bf = pd.DataFrame(bf)

rhpy = []
for i in ds_downsampled['rhpy']:
	rhpy.append(i)
rhpy = pd.DataFrame(rhpy)

mod = []
for i in ds_downsampled['mod']:
	mod.append(i)
mod = pd.DataFrame(mod)

ds = pd.concat([pdb, res, ch, pka, bf, rhpy, mod], axis=1)
ds.columns = ['pdb', 'resid', 'chain', 'pka', 'bf', 'rhpy', 'mod']
print(ds)

# ds.to_excel('../data/balanced_dataset.xlsx')