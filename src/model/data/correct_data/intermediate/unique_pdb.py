import pandas as pd 

ds = pd.read_csv('get_three_features.csv')
pdb = list(ds.iloc[:,1:2].values)

new_pdb = []

for i in pdb:
	j = i[0]
	j = j.replace('.pdb', '')
	# print(j)
	new_pdb.append(j)

uniq = []
for j in new_pdb:
	if j not in uniq:
		uniq.append(j)

string = ''

for i in uniq:
	string += str(i)
	string += ', '

print(string)

# 12, 360 unique pdb structures in total.