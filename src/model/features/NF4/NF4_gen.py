# Motifs

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
window = 3

# dataset import and preprocessing
ds = pd.read_excel('../../data/balanced_dataset.xlsx')
pdb = ds.iloc[:,1]
# print(pdb)
res = ds.iloc[:,2]
chain = ds.iloc[:,3]

# Identify Motifs: CC, CXC, CX4C, CX3C, CX2C, CX2CX5C, CX2CX2C, CX2CX2CX3C
motif = []
for i in range(len(pdb)):
	sing_motif = np.zeros(8)
	pdb_id = str(pdb[i])
	string = ''
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
	for j in range(start-1, end):
		try:
			string += seq[j]
		except:
			string += '-'
	
	# # Motif check: CC, CXC, CX4C, CX3C, CX2C, CX2CX5C, CX2CX2C, CX2CX2CX3C
	if len(re.findall(r"CC", string)) > 0:
		sing_motif[0] = 1.0
		print("0")
	if len(re.findall(r"C[A-Z]C", string)) > 0:
		sing_motif[1] = 1.0
		print("1")
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string)) > 0:
		sing_motif[2] = 1.0
		print("2")
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string)) > 0:
		sing_motif[3] = 1.0
		print("3")
	if len(re.findall(r"C[A-Z][A-Z]C", string)) > 0:
		sing_motif[4] = 1.0
		print("4")
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string)) > 0:
		sing_motif[5] = 1.0
		print("5")
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string)) > 0:
		sing_motif[6] = 1.0
		print("6")
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string)) > 0:
		sing_motif[7] = 1.0
		print("7")

	print(sing_motif)
	motif.append(sing_motif)

filename = 'feature/NF4_' + str(window) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif,outfile)
outfile.close()
