# Protein Class
import pandas as pd
import pickle
from sklearn.preprocessing import LabelEncoder

# dataset import and processing
ds = pd.read_csv('../../dataset/dataset.csv')
pdb = list(ds.iloc[:,1].values)
new_pdb = []
for i in pdb:
	i = i.replace('.pdb', '')
	new_pdb.append(i)

pdb = new_pdb
res = ds.iloc[:,2]
chain = ds.iloc[:,3]

func = []
for i in range(len(pdb)):
	try:
		path_ = '../../../../pdb/' + str(pdb[i]).lower() + '.pdb'
		f = open(path_, "r")
	except:
		path_ = '../../../../pdb/' + str(pdb[i]).upper() + '.pdb'
		f = open(path_, "r")

	for x in f:
		x = x.replace("HEADER    ", "")
		x = x.split(' ')
		ind_func = ''
		for j in range(5):
			ind_func += x[j]
		print(ind_func)
		func.append(ind_func)
		break

filename = 'feature/NF3.pickle'
outfile = open(filename,'wb')
pickle.dump(func ,outfile)
outfile.close()

le = LabelEncoder()
func = le.fit_transform(func)
print(le.classes_, len(le.classes_))

filename = 'feature/NF3_le.pickle'
outfile = open(filename,'wb')
pickle.dump(func ,outfile)
outfile.close()
