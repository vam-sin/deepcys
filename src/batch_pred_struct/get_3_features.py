# Libraries
import requests 
import numpy as np
import pickle
import os
from pka_website import get_pka
from Bf_rhpy_website import get_bf_rhpy

# Take in PDB ID and residue ID. (Ex: 1b2l, 137, A, Output: Disulfide)

def get_features(pdb, res, chain, mod):
	# Parameters:
	res = int(res)

	# Get FASTA and PDB.
	PROJECT_PATH = os.path.dirname(__file__) + "/"
	print(PROJECT_PATH)
	print("\nSteps.")
	# url = 'https://files.rcsb.org/download/' + pdb.upper() + '.pdb'
	# r = requests.get(url)
	filename_pdb = PROJECT_PATH + '/PDB_Data/' + pdb + '.pdb'
	# open(filename_pdb, 'wb').write(r.content)
	print("Obtained PDB.")

	# url = 'https://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=fastachain&compression=NO&structureId=' + pdb.lower() + '&chainId=' + chain
	# r = requests.get(url)
	# filename_fasta = PROJECT_PATH + '/PDB_fasta/' + pdb + '.fasta'
	# open(filename_fasta, 'wb').write(r.content)
	# print("Obtained fasta.")

	# pKa 
	if mod == 'disulphide':
		pka = 99.9
	else:
		pKa = get_pka(pdb, res, chain)
	print("pKa Calculation Done: " + str(pKa))

	# BF_RHPY
	BF, rHpy = get_bf_rhpy(pdb, res, chain)
	print("Calculated BF and rHpy: " + str(BF) + ", " + str(rHpy))

	# Compile X
	X = []
	X.append(pKa)
	X.append(BF)
	X.append(rHpy)

	# print(len(nf1_13), len(nf2_8), len(nf2_7), len(nf2_6), len(nf2_5), len(nf3), len(nf4_13), len(nf4_11), len(nf4_9), len(nf4_7), len(nf4_5), len(nf4_3), len(nf5_13), len(X))
	X = np.asarray(X)

	return X



