import propka.lib, propka.molecular_container
import os
from pathlib import Path
import os.path
import pandas as pd 
import numpy as np 
import random
import pickle

def preprocess_pdb(pdb):
	pdbfile = '../../pdb/' + pdb.replace(' ', '') + '.pdb'
	f = open(pdbfile, 'r')
	g = open('PDB_Data/' + pdb.replace(' ', '') + '_1.pdb', 'w')
	for line in f:
		line = line.replace('HETATM', 'ATOM  ')
		line = line.replace('CSO', 'CYS')
		line = line.replace('OCY', 'CYS')
		line = line.replace('5CA', 'CYS')
		line = line.replace('OCS', 'CYS')
		line = line.replace('CSO', 'CYS')
		line = line.replace('PRS', 'CYS')
		line = line.replace('BCS', 'CYS')
		line = line.replace('CMH', 'CYS')
		line = line.replace('YCM', 'CYS')
		line = line.replace('DCY', 'CYS')

		g.write(line)

def get_pka(pdb, residid, chain):
	pdb = pdb.lower()
	print(pdb, residid, chain)
	pdb = pdb.replace('.pdb', '')
	filename = pdb + '_1.pka'
	cur_path = os.path.dirname(__file__)
	if os.path.isfile(filename) == False:
		PROJECT_PATH = os.path.dirname(__file__) + "/"
		options, _ = propka.lib.loadOptions()
		preprocess_pdb(pdb)
		pdbfile = 'PDB_Data/' + pdb.replace(' ', '') + '_1.pdb'
		filename = pdb + '_1.pka'
		print(options, pdbfile)
		try:
			my_molecule = propka.molecular_container.Molecular_container(pdbfile, options)
			my_molecule.calculate_pka()
			my_molecule.write_pka()
		except:
			print("Random Return")
			return random.uniform(1, 14)
	print(filename)
	print(cur_path)
	
	pka_val_file = open(filename, "r")
	search_results = []
	for x in pka_val_file:
		search_text1 = 'CYS'
		search_text2 = str(residid)
		search_text3 = str(chain)
		if (search_text1 in x) and (search_text2 in x) and (search_text3 in x):
			search_results.append(x)
	try:
		req_text = search_results[len(search_results)-1] 
		req_text = req_text.split(' ')
		values_req_text = []
		for j in req_text:
			if j != "":
				values_req_text.append(j)
		print(float(values_req_text[3]))
		pka_val = float(values_req_text[3])
		pka_val_file.close()
	except:
		print("Random Return")
		return random.uniform(1, 14)
	try:
		cmd = 'rm ' + str(pdb).lower() + '_1.propka_input'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm ' + str(pdb).upper() + '_1.propka_input'
		os.system(cmd)
	except:
		pass

	print("Actual Answer")

	return pka_val



