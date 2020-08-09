# Libraries
import requests 
import numpy as np
import pickle
import os
from feature_gen import get_nf1, get_nf2, get_nf3, get_nf4
from pka import get_pka
from Bf_rhpy import get_bf_rhpy
import os.path

# Take in PDB ID and residue ID. (Ex: 1b2l, 137, A, Output: Disulfide)

def get_features(pdb, res, chain):
	# Parameters:
	res = int(res)

	# Get FASTA and PDB.
	PROJECT_PATH = os.path.dirname(__file__) + "/"
	print(PROJECT_PATH)
	print("\nSteps.")
	filename_pdb = 'PDB_Data/' + pdb.replace(' ', '') + '.pdb'
	if os.path.isfile(filename_pdb) == False:
		url = 'https://files.rcsb.org/download/' + pdb.upper() + '.pdb'
		r = requests.get(url)
		f = open(filename_pdb, 'wb')
		f.write(r.content)
		f.close()
		print("Obtained PDB. ", res)

	# BF_rHpy
	BF, rHpy = get_bf_rhpy(pdb, res, chain)
	print("Calculated BF and rHpy: " + str(BF) + ", " + str(rHpy))


	# Secondary Structure Folds
	nf1_7 = get_nf1(pdb, res, chain, 7)

	print("Calculated NF1.")

	# Amino Acid Signatures in Interaction Shells
	nf2_8, nf2_7, nf2_6, nf2_5 = get_nf2(pdb, res, chain)

	print("Calculated NF2.")

	# Enzyme Class
	nf3 = get_nf3(pdb)
	print(nf3)
	print("Calculated NF3")
	
	# Motifs
	nf4_3 = get_nf4(pdb, res, chain, 3)
	nf4_5 = get_nf4(pdb, res, chain, 5)
	nf4_7 = get_nf4(pdb, res, chain, 7)
	nf4_9 = get_nf4(pdb, res, chain, 9)
	nf4_11 = get_nf4(pdb, res, chain, 11)
	nf4_13 = get_nf4(pdb, res, chain, 13)

	print("Calculated NF4")

	# # Compile X
	X = []
	X.append(BF)
	X.append(rHpy)

	for i in nf1_7:
		X.append(i)

	for i in nf2_5:
		X.append(i)
	for i in nf2_6:
		X.append(i)
	for i in nf2_7:
		X.append(i)
	for i in nf2_8:
		X.append(i)

	for i in nf3:
		X.append(i)

	for i in nf4_3:
		X.append(i)
	for i in nf4_5:
		X.append(i)
	for i in nf4_7:
		X.append(i)
	for i in nf4_9:
		X.append(i)
	for i in nf4_11:
		X.append(i)
	for i in nf4_13:
		X.append(i)


	X = np.asarray(X)
	print(X.shape)

	return X



