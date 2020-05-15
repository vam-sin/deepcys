# Protein Class
import pandas as pd
import pickle
from sklearn.preprocessing import LabelEncoder

ds = pd.read_excel('../../data/dataset.xlsx')
pdb = ds.iloc[:,1].values

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
			# ind_func += ' '
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

filename = 'new_thio/NF3.pickle'
outfile = open(filename,'wb')
pickle.dump(func ,outfile)
outfile.close()
