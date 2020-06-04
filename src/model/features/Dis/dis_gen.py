import pandas as pd 
import pickle 
ds = pd.read_excel("../../data/balanced_dataset.xlsx")
print(ds)
pka = list(ds.iloc[:,4])

dis = []
for i in pka:
	if i == 99.9:
		dis.append(1)
	else:
		dis.append(0)

# print(dis)

# print(list(pka))
filename = 'feature/dis.pickle'
outfile = open(filename,'wb')
pickle.dump(dis,outfile)
outfile.close()