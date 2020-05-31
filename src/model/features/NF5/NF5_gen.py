# Primary Sequence

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
window = 13

# dataset import and preprocessing
ds = pd.read_excel('../../data/dataset.xlsx')
pdb = ds.iloc[:,1]
# print(pdb)
res = ds.iloc[:,2]
chain = ds.iloc[:,3]

prim = []
for i in range(len(pdb)):
	pdb_id = str(pdb[i])
	string = []
	print(pdb_id, i)
	list_ind = 0
	try:
		file = '../../../../fasta/' + pdb_id.upper() +'.fasta.txt'
	except:
		file = '../../../../fasta/' + pdb_id.lower() +'.fasta.txt'

	record = list(SeqIO.parse(file, "fasta"))
	for j in range(len(record)):
		if type(chain[i]) == str:	
			ind = record[j].id[5]
			if ind == chain[i]:
				list_ind = j
				break
		else:
			list_ind = chain[i] - 1

	seq = record[list_ind].seq
	start = res[i] - window
	end = res[i] + window
	# print(seq[res[i]-1])
	for j in range(start-1, end):
		try:
			if seq[j] == 'A':
				string.append(1)
			elif seq[j] == 'R':
				string.append(2)
			elif seq[j] == 'N':
				string.append(3)
			elif seq[j] == 'D':
				string.append(4)
			elif seq[j] == 'C':
				string.append(5)
			elif seq[j] == 'G':
				string.append(6)
			elif seq[j] == 'Q':
				string.append(7)
			elif seq[j] == 'E':
				string.append(8)
			elif seq[j] == 'H':
				string.append(9)
			elif seq[j] == 'I':
				string.append(10)
			elif seq[j] == 'L':
				string.append(11)
			elif seq[j] == 'K':
				string.append(12)
			elif seq[j] == 'M':
				string.append(13)
			elif seq[j] == 'F':
				string.append(14)
			elif seq[j] == 'P':
				string.append(15)
			elif seq[j] == 'S':
				string.append(16)
			elif seq[j] == 'T':
				string.append(17)
			elif seq[j] == 'W':
				string.append(18)
			elif seq[j] == 'Y':
				string.append(19)
			elif seq[j] == 'V':
				string.append(20)
		except:
			string.append(21)
	print(string)
	prim.append(string)

prim = pd.DataFrame(prim)
prim = prim.replace(np.nan,0)
print(prim.isnull().sum().sum())

filename = 'feature/NF5_' + str(window) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(prim,outfile)
outfile.close()


