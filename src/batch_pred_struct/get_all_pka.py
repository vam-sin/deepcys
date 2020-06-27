import propka.lib, propka.molecular_container
import os
from pathlib import Path
import pandas as pd 
import pickle

def get_pka(pdb, residid, chain, mod):
	residid = residid[0]
	chain = chain[0]
	mod = mod[0]
	print(pdb, residid, chain, mod)
	if mod == 'disulphide':
		print("99.9")
		return 99.9
	else:
		try:
			# PROJECT_PATH = os.path.dirname(__file__) + "/"
			# print(PROJECT_PATH)
			filename = 'pKa_Data/' + pdb + '.pka'
			print(filename)

			# my_molecule = propka.molecular_container.Molecular_container(pdbfile, options)
			# my_molecule.calculate_pka()
			# my_molecule.write_pka()
			# print(filename)
			# print(cur_path)
			# print(options, pdbfile)
			pka_val_file = open(filename, "r")
			search_results = []
			for x in pka_val_file:
				search_text1 = 'CYS'
				search_text2 = str(residid)
				search_text3 = str(chain)
				if (search_text1 in x) and (search_text2 in x) and (search_text3 in x):
					search_results.append(x)
			req_text = search_results[len(search_results)-1] 
			req_text = req_text.split(' ')
			values_req_text = []
			for j in req_text:
				if j != "":
					values_req_text.append(j)
			print(float(values_req_text[3]))
			pka_val = float(values_req_text[3])
			pka_val_file.close()
			print(pka_val)

			return pka_val
		except:
			print("Propka error")

			return 0.0

if __name__ == '__main__':
	ds = pd.read_csv('get_three_features.csv')
	pdb = list(ds.iloc[:,1:2].values)
	res = list(ds.iloc[:,2:3].values)
	chain = list(ds.iloc[:,3:4].values)
	mod = list(ds.iloc[:,4:5].values)
	pka = []

	# new_pdb = []

	# for i in pdb:
	# 	j = i[0]
	# 	j = j.replace('.pdb', '')
	# 	# print(j)
	# 	new_pdb.append(j)

	# uniq = []
	# for j in new_pdb:
	# 	if j not in uniq:
	# 		uniq.append(j)

	# print(uniq)

	for i in range(len(pdb)):
		print("######################")
		print(i+1, " Out of ", len(pdb))
		print("######################")
		# try:
		pka.append(get_pka(pdb[i][0].replace('.pdb',''), res[i], chain[i], mod[i]))
		# except:
		# 	pass

filename = 'pka.pickle'
outfile = open(filename,'wb')
pickle.dump(pka,outfile)
outfile.close()
