# Secondary Structure Folds

# Libraries
import pandas as pd
import numpy as np
import wget
import re
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.DSSP import DSSP
from Bio import SeqIO
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import pickle
from sklearn.preprocessing import StandardScaler, LabelEncoder

# Tasks
# Separate window sizes (3, 5, 7, 9, 11, 13)
window = 5

# dataset import and preprocessing
ds = pd.read_excel('../../data/balanced_dataset.xlsx')
pdb = ds.iloc[:,1]
print(pdb)
res = ds.iloc[:,2]
chain = ds.iloc[:,3]

# Structures
# H,G,I: 1 
# T: 2 (T, S)
# S: 3
# B: 4
# E: 5
# - 6
# Exception: 7

ssf_list = []
p = PDBParser()
last_file = '../../../../pdb/1b2l.pdb'
last_dssp = dssp_dict_from_pdb_file(last_file)
for i in range(len(pdb)):
	pdb_id = str(pdb[i])
	print(pdb_id)
	try:
		file = '../../../../pdb/' + pdb_id.lower() +'.pdb'
		if file == last_file:
			dssp = last_dssp
		else:
			last_file = file 
			dssp = dssp_dict_from_pdb_file(file)
			last_dssp = dssp
	except:
		file = '../../../../pdb/' + pdb_id.upper() +'.pdb'
		if file == last_file:
			dssp = last_dssp
		else:
			last_file = file 
			dssp = dssp_dict_from_pdb_file(file)
			last_dssp = dssp
	dssp = dssp[0]
	# print(dssp)
	ssf = []
	start = res[i] - window
	end = res[i] + window
	structure = ''
	for k, v in dssp:
		chain = k
		break
	for j in range(start-1, end):
		try:
			structure = dssp[chain, (' ', j, ' ')][1]
			if structure == 'H' or structure == 'G' or structure == 'I':
				ssf.append(1)
			elif structure == 'T' or structure == 'S':
				ssf.append(2)
			elif structure == 'B':
				ssf.append(3)
			elif structure == 'E':
				ssf.append(4)
			else:
				ssf.append(5)
		except:
			ssf.append(6)
	print(ssf, i)
	ssf_list.append(ssf)

# # Pickle
filename = 'feature/NF1_' + str(window) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(ssf_list,outfile)
outfile.close()
