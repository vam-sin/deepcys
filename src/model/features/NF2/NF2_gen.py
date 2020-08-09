# Proteins in the interaction shell

from Bio.PDB import PDBParser
import pandas as pd
import numpy as np
import pickle

# create parser
parser = PDBParser()

# parameters
radius = 6

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

# Iterate for the protein
inshell_proteins_5 = []
inshell_proteins_6 = []
inshell_proteins_7 = []
inshell_proteins_8 = []
for i in range(len(pdb)):
	print(i, len(pdb))
	inshell_proteins_single_5 = np.zeros(20, dtype = int) 
	inshell_proteins_single_6 = np.zeros(20, dtype = int) 
	inshell_proteins_single_7 = np.zeros(20, dtype = int) 
	inshell_proteins_single_8 = np.zeros(20, dtype = int) 
	try:
		try:
			path_ = '../../../../pdb/' + str(pdb[i]).lower() + '.pdb'
			structure = parser.get_structure('PHA-L', path_)
		except:
			path_ = '../../../../pdb/' + str(pdb[i]).upper() + '.pdb'
			structure = parser.get_structure('PHA-L', path_)
		model = structure[0]
		try:
			for chain in model:
				residue1 = chain[resid[i]] 
				for residue2 in chain:
					if residue1 != residue2:
						try:
							distance = residue1['CA'] - residue2['CA']
						except KeyError:
							continue
						if distance < 5: 
							if residue2.get_resname() == 'ALA':
								inshell_proteins_single_5[0] += 1
							elif residue2.get_resname() == 'ARG':
								inshell_proteins_single_5[1] += 1
							elif residue2.get_resname() == 'ASN':
								inshell_proteins_single_5[2] += 1
							elif residue2.get_resname() == 'ASP':
								inshell_proteins_single_5[3] += 1
							elif residue2.get_resname() == 'CYS':
								inshell_proteins_single_5[4] += 1
							elif residue2.get_resname() == 'GLU':
								inshell_proteins_single_5[5] += 1
							elif residue2.get_resname() == 'GLN':
								inshell_proteins_single_5[6] += 1
							elif residue2.get_resname() == 'GLY':
								inshell_proteins_single_5[7] += 1
							elif residue2.get_resname() == 'HIS':
								inshell_proteins_single_5[8] += 1
							elif residue2.get_resname() == 'ILE':
								inshell_proteins_single_5[9] += 1
							elif residue2.get_resname() == 'LEU':
								inshell_proteins_single_5[10] += 1
							elif residue2.get_resname() == 'LYS':
								inshell_proteins_single_5[11] += 1
							elif residue2.get_resname() == 'MET':
								inshell_proteins_single_5[12] += 1
							elif residue2.get_resname() == 'PHE':
								inshell_proteins_single_5[13] += 1
							elif residue2.get_resname() == 'PRO':
								inshell_proteins_single_5[14] += 1
							elif residue2.get_resname() == 'SER':
								inshell_proteins_single_5[15] += 1
							elif residue2.get_resname() == 'THR':
								inshell_proteins_single_5[16] += 1
							elif residue2.get_resname() == 'TRP':
								inshell_proteins_single_5[17] += 1
							elif residue2.get_resname() == 'TYR':
								inshell_proteins_single_5[18] += 1
							elif residue2.get_resname() == 'VAL':
								inshell_proteins_single_5[19] += 1
						if distance < 6:
							if residue2.get_resname() == 'ALA':
								inshell_proteins_single_6[0] += 1
							elif residue2.get_resname() == 'ARG':
								inshell_proteins_single_6[1] += 1
							elif residue2.get_resname() == 'ASN':
								inshell_proteins_single_6[2] += 1
							elif residue2.get_resname() == 'ASP':
								inshell_proteins_single_6[3] += 1
							elif residue2.get_resname() == 'CYS':
								inshell_proteins_single_6[4] += 1
							elif residue2.get_resname() == 'GLU':
								inshell_proteins_single_6[5] += 1
							elif residue2.get_resname() == 'GLN':
								inshell_proteins_single_6[6] += 1
							elif residue2.get_resname() == 'GLY':
								inshell_proteins_single_6[7] += 1
							elif residue2.get_resname() == 'HIS':
								inshell_proteins_single_6[8] += 1
							elif residue2.get_resname() == 'ILE':
								inshell_proteins_single_6[9] += 1
							elif residue2.get_resname() == 'LEU':
								inshell_proteins_single_6[10] += 1
							elif residue2.get_resname() == 'LYS':
								inshell_proteins_single_6[11] += 1
							elif residue2.get_resname() == 'MET':
								inshell_proteins_single_6[12] += 1
							elif residue2.get_resname() == 'PHE':
								inshell_proteins_single_6[13] += 1
							elif residue2.get_resname() == 'PRO':
								inshell_proteins_single_6[14] += 1
							elif residue2.get_resname() == 'SER':
								inshell_proteins_single_6[15] += 1
							elif residue2.get_resname() == 'THR':
								inshell_proteins_single_6[16] += 1
							elif residue2.get_resname() == 'TRP':
								inshell_proteins_single_6[17] += 1
							elif residue2.get_resname() == 'TYR':
								inshell_proteins_single_6[18] += 1
							elif residue2.get_resname() == 'VAL':
								inshell_proteins_single_6[19] += 1
						if distance < 7:
							if residue2.get_resname() == 'ALA':
								inshell_proteins_single_7[0] += 1
							elif residue2.get_resname() == 'ARG':
								inshell_proteins_single_7[1] += 1
							elif residue2.get_resname() == 'ASN':
								inshell_proteins_single_7[2] += 1
							elif residue2.get_resname() == 'ASP':
								inshell_proteins_single_7[3] += 1
							elif residue2.get_resname() == 'CYS':
								inshell_proteins_single_7[4] += 1
							elif residue2.get_resname() == 'GLU':
								inshell_proteins_single_7[5] += 1
							elif residue2.get_resname() == 'GLN':
								inshell_proteins_single_7[6] += 1
							elif residue2.get_resname() == 'GLY':
								inshell_proteins_single_7[7] += 1
							elif residue2.get_resname() == 'HIS':
								inshell_proteins_single_7[8] += 1
							elif residue2.get_resname() == 'ILE':
								inshell_proteins_single_7[9] += 1
							elif residue2.get_resname() == 'LEU':
								inshell_proteins_single_7[10] += 1
							elif residue2.get_resname() == 'LYS':
								inshell_proteins_single_7[11] += 1
							elif residue2.get_resname() == 'MET':
								inshell_proteins_single_7[12] += 1
							elif residue2.get_resname() == 'PHE':
								inshell_proteins_single_7[13] += 1
							elif residue2.get_resname() == 'PRO':
								inshell_proteins_single_7[14] += 1
							elif residue2.get_resname() == 'SER':
								inshell_proteins_single_7[15] += 1
							elif residue2.get_resname() == 'THR':
								inshell_proteins_single_7[16] += 1
							elif residue2.get_resname() == 'TRP':
								inshell_proteins_single_7[17] += 1
							elif residue2.get_resname() == 'TYR':
								inshell_proteins_single_7[18] += 1
							elif residue2.get_resname() == 'VAL':
								inshell_proteins_single_7[19] += 1
						if distance < 8:
							if residue2.get_resname() == 'ALA':
								inshell_proteins_single_8[0] += 1
							elif residue2.get_resname() == 'ARG':
								inshell_proteins_single_8[1] += 1
							elif residue2.get_resname() == 'ASN':
								inshell_proteins_single_8[2] += 1
							elif residue2.get_resname() == 'ASP':
								inshell_proteins_single_8[3] += 1
							elif residue2.get_resname() == 'CYS':
								inshell_proteins_single_8[4] += 1
							elif residue2.get_resname() == 'GLU':
								inshell_proteins_single_8[5] += 1
							elif residue2.get_resname() == 'GLN':
								inshell_proteins_single_8[6] += 1
							elif residue2.get_resname() == 'GLY':
								inshell_proteins_single_8[7] += 1
							elif residue2.get_resname() == 'HIS':
								inshell_proteins_single_8[8] += 1
							elif residue2.get_resname() == 'ILE':
								inshell_proteins_single_8[9] += 1
							elif residue2.get_resname() == 'LEU':
								inshell_proteins_single_8[10] += 1
							elif residue2.get_resname() == 'LYS':
								inshell_proteins_single_8[11] += 1
							elif residue2.get_resname() == 'MET':
								inshell_proteins_single_8[12] += 1
							elif residue2.get_resname() == 'PHE':
								inshell_proteins_single_8[13] += 1
							elif residue2.get_resname() == 'PRO':
								inshell_proteins_single_8[14] += 1
							elif residue2.get_resname() == 'SER':
								inshell_proteins_single_8[15] += 1
							elif residue2.get_resname() == 'THR':
								inshell_proteins_single_8[16] += 1
							elif residue2.get_resname() == 'TRP':
								inshell_proteins_single_8[17] += 1
							elif residue2.get_resname() == 'TYR':
								inshell_proteins_single_8[18] += 1
							elif residue2.get_resname() == 'VAL':
								inshell_proteins_single_8[19] += 1
								
			print(inshell_proteins_single_5)
			print(inshell_proteins_single_6)
			print(inshell_proteins_single_7)
			print(inshell_proteins_single_8)
		except:
			print("Error")
	except:
		print("Error2")
		pass
	inshell_proteins_5.append(inshell_proteins_single_5)
	inshell_proteins_6.append(inshell_proteins_single_6)
	inshell_proteins_7.append(inshell_proteins_single_7)
	inshell_proteins_8.append(inshell_proteins_single_8)

filename = 'NF2_5.pickle'
outfile = open(filename,'wb')
pickle.dump(inshell_proteins_5 ,outfile)
outfile.close()

filename = 'NF2_6.pickle'
outfile = open(filename,'wb')
pickle.dump(inshell_proteins_6 ,outfile)
outfile.close()

filename = 'NF2_7.pickle'
outfile = open(filename,'wb')
pickle.dump(inshell_proteins_7 ,outfile)
outfile.close()

filename = 'NF2_8.pickle'
outfile = open(filename,'wb')
pickle.dump(inshell_proteins_8 ,outfile)
outfile.close()