import os
import random

def get_bf_rhpy(pdb, res, chain):
	PROJECT_PATH = os.path.dirname(__file__) + "/"
	pdbfile = PROJECT_PATH + 'PDB_Data/' + pdb + '.pdb'

	list_file = os.path.join(PROJECT_PATH, 'menv_server', 'list')
	num_file = os.path.join(PROJECT_PATH, 'menv_server', 'num')

	f = open(list_file, "w+")
	f.write(pdbfile)

	f = open(num_file, "w+")
	f.write(str(res))

	# Run script
	script_file = os.path.join(PROJECT_PATH, 'menv_server', 'script3')
	print(script_file)
	cmd = 'csh ' + script_file + ' ' + list_file + ' ' + num_file 
	print(cmd)
	os.system(cmd)

	# Extract from result file
	try:
		result_file = os.path.join(PROJECT_PATH, 'PDB_Data', pdb + '-new.menv')
		f = open(result_file, "r")
		for x in f:
			x = x.split(',')
			residue_num = x[1].split(' ')
			residue = int(residue_num[len(residue_num)-3])
			atom = x[3]
			ch = residue_num[len(residue_num)-6]
			bf = float(x[4])
			rhpy = float(x[7])

			if (residue == res) and (ch == chain) and ((atom == 'TO') or (atom == 'TD')):
				print("New Menv")
				return bf, rhpy

	except:
		print("New Menv Failed")

	try:
		print("Try Running Renumber")
		renumber_script = os.path.join(PROJECT_PATH, 'menv_server', 'newrenumber.sh')
		psf_file = PROJECT_PATH + 'PDB_Data/' + pdb + '-psf.menv'
		cmd = 'csh ' + renumber_script + ' ' + pdbfile + ' ' + num_file 
		print(cmd)
		os.system(cmd)
		result_file = os.path.join(PROJECT_PATH, 'PDB_Data', pdb + '-new.menv')
		f = open(result_file, "r")
		for x in f:
			x = x.split(',')
			residue_num = x[1].split(' ')
			residue = int(residue_num[len(residue_num)-3])
			atom = x[3]
			ch = residue_num[len(residue_num)-6]
			bf = float(x[4])
			rhpy = float(x[7])

			if (residue == res) and (ch == chain) and ((atom == 'TO') or (atom == 'TD')):
				print("New Menv")
				return bf, rhpy
	except:
		print("Renumber failed")

	try:
		result_file = os.path.join(PROJECT_PATH, 'PDB_Data', pdb + '-psf.menv')
		f = open(result_file, "r")
		for x in f:
			x = x.split(',')
			residue = int(x[1])
			aa = x[2]
			atom = x[4]
			bf = float(x[5])
			rhpy = float(x[8])
			if (aa == '  CYS') and ((atom == '  TO') or (atom == '  TD')):
				print("PSF Menv")
				return bf, rhpy
	except:
		print("No files, Sending random")
		return random.uniform(0, 1), random.uniform(0, 1)

	print("send random")
	return random.uniform(0, 1), random.uniform(0, 1)