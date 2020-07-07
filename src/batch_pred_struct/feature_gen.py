import os
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import numpy as np
from Bio.PDB import *
from Bio.PDB import PDBParser
from Bio import SeqIO
import pickle
from sklearn.preprocessing import LabelEncoder
from seq_extract import get_sequence
import re

def get_nf1(pdb, res, chain, nf1_window):
	PROJECT_PATH = os.path.dirname(__file__) + "/"
	filename_pdb = PROJECT_PATH + '/PDB_Data/' + pdb + '.pdb'
	dssp = dssp_dict_from_pdb_file(filename_pdb)
	dssp = dssp[0]
	nf1 = []
	start = res - nf1_window
	end = res + nf1_window
	structure = ''
	for k, v in dssp:
		chain = k
		break
	for j in range(start-1, end):
		try:
			structure = dssp[chain, (' ', j, ' ')][1]
			if structure == 'H' or structure == 'G' or structure == 'I':
				nf1.append(1)
			elif structure == 'T' or structure == 'S':
				nf1.append(2)
			elif structure == 'B':
				nf1.append(3)
			elif structure == 'E':
				nf1.append(4)
			else:
				nf1.append(5)
		except:
			nf1.append(6)

	print("NF1_" + str(nf1_window) + ": " + str(nf1))

	return nf1

def get_nf2(pdb, res, chain):
	nf2_8_single = np.zeros(20, dtype = int)
	nf2_7_single = np.zeros(20, dtype = int)
	nf2_6_single = np.zeros(20, dtype = int)
	nf2_5_single = np.zeros(20, dtype = int) 
	PROJECT_PATH = os.path.dirname(__file__) + "/"
	filename_pdb = PROJECT_PATH + '/PDB_Data/' + pdb + '.pdb'
	parser = PDBParser()
	structure = parser.get_structure('PHA-L', filename_pdb)
	model = structure[0]
	try:
		# Iterate for all chains
		for chain in model:
			residue1 = chain[res] 
			for residue2 in chain:
				if residue1 != residue2:
					try:
						distance = residue1['CA'] - residue2['CA']
					except KeyError:
						continue
					if distance < 5: 
						if residue2.get_resname() == 'ALA':
							nf2_5_single[0] += 1
						elif residue2.get_resname() == 'ARG':
							nf2_5_single[1] += 1
						elif residue2.get_resname() == 'ASN':
							nf2_5_single[2] += 1
						elif residue2.get_resname() == 'ASP':
							nf2_5_single[3] += 1
						elif residue2.get_resname() == 'CYS':
							nf2_5_single[4] += 1
						elif residue2.get_resname() == 'GLU':
							nf2_5_single[5] += 1
						elif residue2.get_resname() == 'GLN':
							nf2_5_single[6] += 1
						elif residue2.get_resname() == 'GLY':
							nf2_5_single[7] += 1
						elif residue2.get_resname() == 'HIS':
							nf2_5_single[8] += 1
						elif residue2.get_resname() == 'ILE':
							nf2_5_single[9] += 1
						elif residue2.get_resname() == 'LEU':
							nf2_5_single[10] += 1
						elif residue2.get_resname() == 'LYS':
							nf2_5_single[11] += 1
						elif residue2.get_resname() == 'MET':
							nf2_5_single[12] += 1
						elif residue2.get_resname() == 'PHE':
							nf2_5_single[13] += 1
						elif residue2.get_resname() == 'PRO':
							nf2_5_single[14] += 1
						elif residue2.get_resname() == 'SER':
							nf2_5_single[15] += 1
						elif residue2.get_resname() == 'THR':
							nf2_5_single[16] += 1
						elif residue2.get_resname() == 'TRP':
							nf2_5_single[17] += 1
						elif residue2.get_resname() == 'TYR':
							nf2_5_single[18] += 1
						elif residue2.get_resname() == 'VAL':
							nf2_5_single[19] += 1
					if distance < 6:
						if residue2.get_resname() == 'ALA':
							nf2_6_single[0] += 1
						elif residue2.get_resname() == 'ARG':
							nf2_6_single[1] += 1
						elif residue2.get_resname() == 'ASN':
							nf2_6_single[2] += 1
						elif residue2.get_resname() == 'ASP':
							nf2_6_single[3] += 1
						elif residue2.get_resname() == 'CYS':
							nf2_6_single[4] += 1
						elif residue2.get_resname() == 'GLU':
							nf2_6_single[5] += 1
						elif residue2.get_resname() == 'GLN':
							nf2_6_single[6] += 1
						elif residue2.get_resname() == 'GLY':
							nf2_6_single[7] += 1
						elif residue2.get_resname() == 'HIS':
							nf2_6_single[8] += 1
						elif residue2.get_resname() == 'ILE':
							nf2_6_single[9] += 1
						elif residue2.get_resname() == 'LEU':
							nf2_6_single[10] += 1
						elif residue2.get_resname() == 'LYS':
							nf2_6_single[11] += 1
						elif residue2.get_resname() == 'MET':
							nf2_6_single[12] += 1
						elif residue2.get_resname() == 'PHE':
							nf2_6_single[13] += 1
						elif residue2.get_resname() == 'PRO':
							nf2_6_single[14] += 1
						elif residue2.get_resname() == 'SER':
							nf2_6_single[15] += 1
						elif residue2.get_resname() == 'THR':
							nf2_6_single[16] += 1
						elif residue2.get_resname() == 'TRP':
							nf2_6_single[17] += 1
						elif residue2.get_resname() == 'TYR':
							nf2_6_single[18] += 1
						elif residue2.get_resname() == 'VAL':
							nf2_6_single[19] += 1
					if distance < 7:
						if residue2.get_resname() == 'ALA':
							nf2_7_single[0] += 1
						elif residue2.get_resname() == 'ARG':
							nf2_7_single[1] += 1
						elif residue2.get_resname() == 'ASN':
							nf2_7_single[2] += 1
						elif residue2.get_resname() == 'ASP':
							nf2_7_single[3] += 1
						elif residue2.get_resname() == 'CYS':
							nf2_7_single[4] += 1
						elif residue2.get_resname() == 'GLU':
							nf2_7_single[5] += 1
						elif residue2.get_resname() == 'GLN':
							nf2_7_single[6] += 1
						elif residue2.get_resname() == 'GLY':
							nf2_7_single[7] += 1
						elif residue2.get_resname() == 'HIS':
							nf2_7_single[8] += 1
						elif residue2.get_resname() == 'ILE':
							nf2_7_single[9] += 1
						elif residue2.get_resname() == 'LEU':
							nf2_7_single[10] += 1
						elif residue2.get_resname() == 'LYS':
							nf2_7_single[11] += 1
						elif residue2.get_resname() == 'MET':
							nf2_7_single[12] += 1
						elif residue2.get_resname() == 'PHE':
							nf2_7_single[13] += 1
						elif residue2.get_resname() == 'PRO':
							nf2_7_single[14] += 1
						elif residue2.get_resname() == 'SER':
							nf2_7_single[15] += 1
						elif residue2.get_resname() == 'THR':
							nf2_7_single[16] += 1
						elif residue2.get_resname() == 'TRP':
							nf2_7_single[17] += 1
						elif residue2.get_resname() == 'TYR':
							nf2_7_single[18] += 1
						elif residue2.get_resname() == 'VAL':
							nf2_7_single[19] += 1
					if distance < 8:
						if residue2.get_resname() == 'ALA':
							nf2_8_single[0] += 1
						elif residue2.get_resname() == 'ARG':
							nf2_8_single[1] += 1
						elif residue2.get_resname() == 'ASN':
							nf2_8_single[2] += 1
						elif residue2.get_resname() == 'ASP':
							nf2_8_single[3] += 1
						elif residue2.get_resname() == 'CYS':
							nf2_8_single[4] += 1
						elif residue2.get_resname() == 'GLU':
							nf2_8_single[5] += 1
						elif residue2.get_resname() == 'GLN':
							nf2_8_single[6] += 1
						elif residue2.get_resname() == 'GLY':
							nf2_8_single[7] += 1
						elif residue2.get_resname() == 'HIS':
							nf2_8_single[8] += 1
						elif residue2.get_resname() == 'ILE':
							nf2_8_single[9] += 1
						elif residue2.get_resname() == 'LEU':
							nf2_8_single[10] += 1
						elif residue2.get_resname() == 'LYS':
							nf2_8_single[11] += 1
						elif residue2.get_resname() == 'MET':
							nf2_8_single[12] += 1
						elif residue2.get_resname() == 'PHE':
							nf2_8_single[13] += 1
						elif residue2.get_resname() == 'PRO':
							nf2_8_single[14] += 1
						elif residue2.get_resname() == 'SER':
							nf2_8_single[15] += 1
						elif residue2.get_resname() == 'THR':
							nf2_8_single[16] += 1
						elif residue2.get_resname() == 'TRP':
							nf2_8_single[17] += 1
						elif residue2.get_resname() == 'TYR':
							nf2_8_single[18] += 1
						elif residue2.get_resname() == 'VAL':
							nf2_8_single[19] += 1

		print(nf2_5_single)
		print(nf2_6_single)
		print(nf2_7_single)
		print(nf2_8_single)
	except:
		print("NF2 Production Failed")
		print(nf2_5_single)
		print(nf2_6_single)
		print(nf2_7_single)
		print(nf2_8_single)

	return nf2_8_single, nf2_7_single, nf2_6_single, nf2_5_single

def get_nf3(pdb):
	try:
		PROJECT_PATH = os.path.dirname(__file__) + "/"
		path_ = PROJECT_PATH + '/PDB_Data/' + pdb + '.pdb'
		infile_path = PROJECT_PATH + '/NF3.pickle'
		infile = open(infile_path,'rb')
		func = pickle.load(infile)
		infile.close()

		le = LabelEncoder()
		func = le.fit(func)
		f = open(path_, "r")

		for x in f:
			x = x.replace("HEADER    ", "")
			x = x.split(' ')
			ind_func = ''
			for j in range(5):
				ind_func += x[j]
				# ind_func += ' '
			print(ind_func)
			nf3 = le.transform([ind_func])

			break
		
		print("NF3: ", nf3)

		return nf3 

	except:
		print("NF3: ", 999)

		return [999]

def get_nf4(pdb, res, chain, window):
	sing_motif = np.zeros(8)
	try:
		file = 'PDB_Data/' + pdb.upper() +'.pdb'
		string = get_sequence(file, res, chain, window)
	except:
		file = 'PDB_Data/' + pdb.lower() +'.pdb'
		string = get_sequence(file, res, chain, window)
	print(string)
	
	if len(re.findall(r"CC", string)) > 0:
		sing_motif[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string)) > 0:
		sing_motif[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string)) > 0:
		sing_motif[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string)) > 0:
		sing_motif[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string)) > 0:
		sing_motif[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string)) > 0:
		sing_motif[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string)) > 0:
		sing_motif[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string)) > 0:
		sing_motif[7] = 1.0

	print("NF4_" + str(window) + ": " + str(sing_motif))
	
	return sing_motif

