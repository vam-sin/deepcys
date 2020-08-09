import os
import random
import os.path

def remove_files(pdb):
	try:
		cmd = 'rm combined.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm commaremoved.menv'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm Crdresults.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm finalresult.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm finalappended.menv'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm finalcut.menv'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm finalfinal.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm finalresultwithnum.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm fort.71'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm fort.72'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm fort.88'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm menvfilewithnum.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm Menvresults.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm outfile.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm output.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm pythonresults.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm requiredoutput.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm revisedoutput.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm totalresidues.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm uniqueMenvresidues.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm uniqueresidues.txt'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm temp'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).upper() + '-new.menv'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).lower() + '-new.menv'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).upper() + '-psf.crd'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).lower() + '-psf.crd'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).upper() + '-psf.menv'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).lower() + '-psf.menv'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).upper() + '-psf.out'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).lower() + '-psf.out'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).upper() + '-psf.psf'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).lower() + '-psf.psf'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).upper() + '-psf.pdb'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm PDB_Data/' + str(pdb).lower() + '-psf.pdb'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm result.pdb'
		os.system(cmd)
	except:
		pass
	try:
		cmd = 'rm uniqueresult.pdb'
		os.system(cmd)
	except:
		pass


def get_bf_rhpy(pdb, res, chain):
	try:
		print(res, chain)
		PROJECT_PATH = os.path.dirname(__file__) + "/"
		pdbfile = PROJECT_PATH + 'PDB_Data/' + pdb + '.pdb'

		list_file = os.path.join(PROJECT_PATH, 'menv_server', 'list')
		num_file = os.path.join(PROJECT_PATH, 'menv_server', 'num')

		f = open(list_file, "w+")
		f.write(pdbfile)

		f = open(num_file, "w+")
		f.write(str(res))

		result_file = os.path.join(PROJECT_PATH, 'PDB_Data', pdb + '-new.menv')
		if os.path.isfile(result_file) == False:
			# Run script
			script_file = os.path.join(PROJECT_PATH, 'menv_server', 'script3')
			print(script_file)
			cmd = 'csh ' + script_file + ' ' + list_file + ' ' + num_file 
			print(cmd)
			os.system(cmd)

			# Extract from result file
			try:
				result_file = os.path.join(PROJECT_PATH, 'PDB_Data', pdb + '-new.menv')
				print("First time New Menv")
				f = open(result_file, "r")
				for x in f:
					x = x.split(',')
					residue_num = x[1].split(' ')
					residue = int(residue_num[len(residue_num)-3])
					atom = x[3]
					atom = atom.replace(' ', '')
					ch = residue_num[5].replace(' ', '')
					bf = float(x[4])
					rhpy = float(x[7])

					if (int(residue) == int(res)) and (ch == chain.replace('\n', '')):
						print("Half Match, Check Atom")
						if atom == 'TO' or atom == 'TD':
							print("New Menv")
							f.close()
							remove_files(pdb)
							return bf, rhpy

			except:
				print("New Menv Failed")

			try:
				print("Try Running Renumber")
				renumber_script = os.path.join(PROJECT_PATH, 'menv_server', 'newrenumber.sh')
				psf_file = PROJECT_PATH + 'PDB_Data/' + pdb + '-psf.menv'
				cmd = 'sh ' + renumber_script + ' ' + pdbfile + ' ' + psf_file
				print(cmd)
				os.system(cmd)
				result_file = os.path.join(PROJECT_PATH, 'PDB_Data', pdb + '-new.menv')
				f = open(result_file, "r")
				for x in f:
					x = x.split(',')
					residue_num = x[1].split(' ')
					residue = int(residue_num[len(residue_num)-3])
					aa = ''
					atom = x[3]
					atom = atom.replace(' ', '')
					ch = residue_num[5]
					bf = float(x[4])
					rhpy = float(x[7])

					if (residue == res) and (ch == chain.replace('\n', '')) and ((atom == 'TO') or (atom == 'TD')):
						print("Renumbered New Menv")
						remove_files(pdb)
						return bf, rhpy
			except:
				print("Renumber failed")

			try:
				print("Basic psf because new menv renumber also failed")
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
						remove_files(pdb)
						return bf, rhpy
			except:
				print("No files, Sending random")
				remove_files(pdb)
				return random.uniform(0, 1), random.uniform(0, 1)

			print("send random")
			remove_files(pdb)
			return random.uniform(0, 1), random.uniform(0, 1)

		else:
			print("First time New Menv")
			f = open(result_file, "r")
			for x in f:
				x = x.split(',')
				residue_num = x[1].split(' ')
				residue = int(residue_num[len(residue_num)-3])
				atom = x[3]
				atom = atom.replace(' ', '')
				ch = residue_num[5].replace(' ', '')
				bf = float(x[4])
				rhpy = float(x[7])

				if (int(residue) == int(res)) and (ch == chain.replace('\n', '')):
					print("Half Match, Check Atom")
					if atom == 'TO' or atom == 'TD':
						print("New Menv")
						f.close()
						# remove_files(pdb)
						return bf, rhpy


		print("Random Return")
		remove_files(pdb)
		return random.uniform(0, 1), random.uniform(0, 1)

	except:
		print("Random Return")
		remove_files(pdb)
		return random.uniform(0, 1), random.uniform(0, 1)