# Checks if the mentioned disulphides have SSBOND tags

import pandas as pd 

ds = pd.read_excel('../data/new_ds_data.xlsx')
ds.columns = ['no', 'pdb', 'res', 'ch', 'pka', 'bf', 'rhpy', 'mod']

ds = ds.loc[ds['mod'] == 'disulphide']

pdb = list(ds['pdb'].values)
res = list(ds['res'].values)
ch = list(ds['ch'].values)
print(ds)

y = ds['mod']

dis = []

for i in range(1, len(pdb)):
	print(i+1, " Out of ", len(pdb))
	try:
		filename_pdb = '../../../pdb/' + pdb[i].lower() + '.pdb'
		f = open(filename_pdb, "r")
		for line in f:
			line = line.split(' ')
			new_line = []
			for j in line:
				if j != '':
					new_line.append(j)
			if new_line[0] == 'SSBOND':
				# 2,3,4 | 5,6,7 (CYS, CH, RES)
				if (new_line[2] == 'CYS' and new_line[3] == ch[i] and new_line[4] == str(res[i])) or (new_line[5] == 'CYS' and new_line[6] == ch[i] and new_line[7] == str(res[i])): 
					add = []
					add.append(pdb[i])
					add.append(res[i])
					add.append(ch[i])
					dis.append(add)

	except:
		pass

	try:
		filename_pdb = '../../../pdb/' + pdb[i].upper() + '.pdb'
		f = open(filename_pdb, "r")
		for line in f:
			line = line.split(' ')
			new_line = []
			for j in line:
				if j != '':
					new_line.append(j)
			if new_line[0] == 'SSBOND':
				# 2,3,4 | 5,6,7 (CYS, CH, RES)
				if (new_line[2] == 'CYS' and new_line[3] == ch[i] and new_line[4] == str(res[i])) or (new_line[5] == 'CYS' and new_line[6] == ch[i] and new_line[7] == str(res[i])): 
					add = []
					add.append(pdb[i])
					add.append(res[i])
					add.append(ch[i])
					dis.append(add)

	except:
		pass

# print(dis, len(dis))

print("Number of entries with SSBOND: ", len(dis))
print("Number of Disulphide Entires: ", len(pdb))