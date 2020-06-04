import pandas as pd 

ds = pd.read_excel('../data/new_ds_data.xlsx')
ds.columns = ['no', 'pdb', 'res', 'ch', 'pka', 'bf', 'rhpy', 'mod']

ds = ds.loc[ds['mod'] == 'meta-binding']

pdb = list(ds['pdb'].values)
res = list(ds['res'].values)
ch = list(ds['ch'].values)
mod = list(ds['mod'].values)

# print(pdb)

y = ds['mod']

met = []

for i in range(1, len(pdb)):
	print(i)
	try:
		filename_pdb = '../../../pdb/' + pdb[i].lower() + '.pdb'
		f = open(filename_pdb, "r")
		for line in f:
			line = line.split(' ')
			new_line = []
			for x in line:
				if x != '':
					new_line.append(x)
			if new_line[0] == 'LINK':
				# LINK        ZN    ZN A 375                 SG  CYS A 110     1555   1555  2.24  
				if (new_line[5] == 'SG' and new_line[6] == 'CYS' and new_line[7] == ch[i] and str(res[i]) == new_line[8] and ('ZN' in new_line[2] or 'FE' in new_line[2] or 'CD' in new_line[2] or 'CU' in new_line[2])) or (new_line[1] == 'SG' and new_line[2] == 'CYS' and new_line[3] == ch[i] and str(res[i]) == new_line[4] and ('ZN' in new_line[5] or 'FE' in new_line[5] or 'CD' in new_line[5] or 'CU' in new_line[5])):
					print(new_line)
					add = str(pdb[i]) +', ' + str(res[i]) + ', ' + str(ch[i] + ", " + str(mod[i])) 
					met.append(add)
	except:
		pass

	try:
		filename_pdb = '../../../pdb/' + pdb[i].upper() + '.pdb'
		f = open(filename_pdb, "r")
		for line in f:
			line = line.split(' ')
			new_line = []
			for x in line:
				if x != '':
					new_line.append(x)
			if new_line[0] == 'LINK':
				# LINK        ZN    ZN A 375                 SG  CYS A 110     1555   1555  2.24  
				if (new_line[5] == 'SG' and new_line[6] == 'CYS' and new_line[7] == ch[i] and str(res[i]) == new_line[8] and ('ZN' in new_line[2] or 'FE' in new_line[2] or 'CD' in new_line[2] or 'CU' in new_line[2])) or (new_line[1] == 'SG' and new_line[2] == 'CYS' and new_line[3] == ch[i] and str(res[i]) == new_line[4] and ('ZN' in new_line[5] or 'FE' in new_line[5] or 'CD' in new_line[5] or 'CU' in new_line[5])):
					print(new_line)
					add = str(pdb[i]) +', ' + str(res[i]) + ', ' + str(ch[i] + ", " + str(mod[i])) 
					met.append(add)
	except:
		pass

print("Number of entries with Match: ", len(dis))
print("Number of Metal-binding Entires: ", len(pdb))