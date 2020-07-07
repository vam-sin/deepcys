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
		open(filename_pdb, 'wb').write(r.content)
		print("Obtained PDB.")

	# # pKa 
	pKa = get_pka(pdb, res, chain)
	print("pKa Calculation Done: " + str(pKa))

	# # BF_RHPY
	BF, rHpy = get_bf_rhpy(pdb, res, chain)
	print("Calculated BF and rHpy: " + str(BF) + ", " + str(rHpy))

	# # NF1
	nf1_9 = get_nf1(pdb, res, chain, 9)

	print("Calculated NF1.")

	# # NF2
	nf2_8, nf2_7, nf2_6, nf2_5 = get_nf2(pdb, res, chain)

	print("Calculated NF2.")

	# # NF3
	nf3 = get_nf3(pdb)
	print(nf3)
	print("Calculated NF3")
	
	# # NF4
	nf4_3 = get_nf4(pdb, res, chain, 3)
	nf4_5 = get_nf4(pdb, res, chain, 5)
	nf4_7 = get_nf4(pdb, res, chain, 7)
	nf4_9 = get_nf4(pdb, res, chain, 9)
	nf4_11 = get_nf4(pdb, res, chain, 11)
	nf4_13 = get_nf4(pdb, res, chain, 13)

	print("Calculated NF4")

	# # Compile X
	X = []
	X.append(pKa)
	X.append(BF)
	X.append(rHpy)

	for i in nf1_9:
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


	# print(len(nf1_13), len(nf2_8), len(nf2_7), len(nf2_6), len(nf2_5), len(nf3), len(nf4_13), len(nf4_11), len(nf4_9), len(nf4_7), len(nf4_5), len(nf4_3), len(nf5_13), len(X))
	X = np.asarray(X)
	print(X.shape)

	return X



