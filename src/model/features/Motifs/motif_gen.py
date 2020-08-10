# Motifs
# Takes Sequence from the PDB Structure now.
# Libraries
import pandas as pd
import numpy as np
import wget
import os
import re
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.DSSP import DSSP
from Bio import SeqIO
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import pickle
from sklearn.preprocessing import StandardScaler, LabelEncoder
from seq_extract import get_sequence

# Tasks
# Separate window sizes (3, 5, 7, 9, 11, 13)
window = 13

# read structure from file
ds = pd.read_csv('../../dataset/dataset.csv')
pdb = list(ds.iloc[:,1].values)
new_pdb = []
for i in pdb:
	i = i.replace('.pdb', '')
	new_pdb.append(i)

pdb = new_pdb
res = ds.iloc[:,2]
chain = ds.iloc[:,3]

# Identify Motifs: CC, CXC, CX4C, CX3C, CX2C, CX2CX5C, CX2CX2C, CX2CX2CX3C
motif13 = []
motif11 = []
motif9 = []
motif7 = []
motif5 = []
motif3 = []

for i in range(len(pdb)):
	print(i, len(motif13))
	sing_motif13 = np.zeros(8)
	sing_motif11 = np.zeros(8)
	sing_motif9 = np.zeros(8)
	sing_motif7 = np.zeros(8)
	sing_motif5 = np.zeros(8)
	sing_motif3 = np.zeros(8)
	
	pdb_id = str(pdb[i])
	string = ''
	print(pdb_id, i, len(pdb))
	list_ind = 0
	file1 = '../../../../pdb/' + pdb_id.upper() +'.pdb'
	file2 = '../../../../pdb/' + pdb_id.lower() +'.pdb'
	if os.path.isfile(file1) == True:
		string13, string11, string9, string7, string5, string3 = get_sequence(file1, res[i], chain[i], window)
	else:
		string13, string11, string9, string7, string5, string3 = get_sequence(file2, res[i], chain[i], window)
	
	# Motif check: CC, CXC, CX4C, CX3C, CX2C, CX2CX5C, CX2CX2C, CX2CX2CX3C
	if len(re.findall(r"CC", string13)) > 0:
		sing_motif13[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string13)) > 0:
		sing_motif13[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string13)) > 0:
		sing_motif13[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string13)) > 0:
		sing_motif13[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string13)) > 0:
		sing_motif13[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string13)) > 0:
		sing_motif13[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string13)) > 0:
		sing_motif13[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string13)) > 0:
		sing_motif13[7] = 1.0

	if len(re.findall(r"CC", string11)) > 0:
		sing_motif11[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string11)) > 0:
		sing_motif11[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string11)) > 0:
		sing_motif11[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string11)) > 0:
		sing_motif11[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string11)) > 0:
		sing_motif11[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string11)) > 0:
		sing_motif11[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string11)) > 0:
		sing_motif11[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string11)) > 0:
		sing_motif11[7] = 1.0


	if len(re.findall(r"CC", string9)) > 0:
		sing_motif9[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string9)) > 0:
		sing_motif9[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string9)) > 0:
		sing_motif9[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string9)) > 0:
		sing_motif9[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string9)) > 0:
		sing_motif9[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string9)) > 0:
		sing_motif9[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string9)) > 0:
		sing_motif9[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string9)) > 0:
		sing_motif9[7] = 1.0


	if len(re.findall(r"CC", string7)) > 0:
		sing_motif7[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string7)) > 0:
		sing_motif7[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string7)) > 0:
		sing_motif7[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string7)) > 0:
		sing_motif7[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string7)) > 0:
		sing_motif7[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string7)) > 0:
		sing_motif7[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string7)) > 0:
		sing_motif7[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string7)) > 0:
		sing_motif7[7] = 1.0

	if len(re.findall(r"CC", string5)) > 0:
		sing_motif5[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string5)) > 0:
		sing_motif5[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string5)) > 0:
		sing_motif5[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string5)) > 0:
		sing_motif5[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string5)) > 0:
		sing_motif5[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string5)) > 0:
		sing_motif5[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string5)) > 0:
		sing_motif5[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string5)) > 0:
		sing_motif5[7] = 1.0


	if len(re.findall(r"CC", string3)) > 0:
		sing_motif3[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string3)) > 0:
		sing_motif3[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string3)) > 0:
		sing_motif3[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string3)) > 0:
		sing_motif3[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string3)) > 0:
		sing_motif3[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string3)) > 0:
		sing_motif3[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string3)) > 0:
		sing_motif3[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string3)) > 0:
		sing_motif3[7] = 1.0

	print(sing_motif13)
	print(sing_motif11)
	print(sing_motif9)
	print(sing_motif7)
	print(sing_motif5)
	print(sing_motif3)
	motif13.append(sing_motif13)
	motif11.append(sing_motif11)
	motif9.append(sing_motif9)
	motif7.append(sing_motif7)
	motif5.append(sing_motif5)
	motif3.append(sing_motif3)

filename = 'feature/NF4_' + str(13) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif13,outfile)
outfile.close()

filename = 'feature/NF4_' + str(11) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif11,outfile)
outfile.close()

filename = 'feature/NF4_' + str(9) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif9,outfile)
outfile.close()

filename = 'feature/NF4_' + str(7) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif7,outfile)
outfile.close()

filename = 'feature/NF4_' + str(5) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif5,outfile)
outfile.close()

filename = 'feature/NF4_' + str(3) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif3,outfile)
outfile.close()
