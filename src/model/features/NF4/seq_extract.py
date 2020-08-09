# Get required sequence from the pdb structure

def three_to_one(aa):
	d = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER':'*',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

	for i in d.keys():
		if i in aa:
			return d[i]

def get_sequence(pdb, res, chain, window):
	sequence = ''
	start = res - window
	instart = start
	end = res + window
	mod_types = []
	if start <= 0:
		for i in range(abs(start)+1):
			sequence += 'Z'
		start = 1
	for i in range(start, end+1):
		try:
			pdb_file = open(pdb, "r")
			for line in pdb_file:
				line = line.split(' ')
				new_line = []
				for seg in line:
					if seg != '':
						new_line.append(seg)

				if new_line[0] == 'ATOM' or new_line[0] == 'HETATM':
					try:
						try:
							res_num = int(new_line[4].replace(chain, ''))
							ch = new_line[4][0]
							aa = new_line[3]
						except:
							res_num = int(new_line[5])
							aa = new_line[3]
							ch = new_line[4]
					except:
						res_num = 0
						aa = 'A'
						ch = 'A'

					if int(i) == res_num and chain == ch:
						try:
							try:
								sequence += three_to_one(aa)
							except:
								if 'CSO' in aa:
									sequence += 'C'
								elif 'CYS' in aa:
									sequence += 'C'
								elif 'OCY' in aa:
									sequence += 'C'
								elif '5CA' in aa:
									sequence += 'C'
								elif 'OCS' in aa:
									sequence += 'C'
								elif 'PRS' in aa:
									sequence += 'C'
								elif 'BCS' in aa:
									sequence += 'C'
								elif 'CMH' in aa:
									sequence += 'C'
								elif 'YCM' in aa:
									sequence += 'C'
								elif 'DCY' in aa:
									sequence += 'C'
								else:
									for mod in mod_types:
										if mod in aa:
											sequence += 'C'
											print("MOD RES SLAYING")
											break
							break
						except:
							pass

				elif new_line[0] == 'MODRES':
					if new_line[5] == 'CYS':
						mod_types.append(new_line[2]) 
			if len(sequence) != int(i) - instart + 1:
				sequence += 'Z'						 
		except:
			sequence += 'Z'

	total_len = window*2 + 1
	if len(sequence) < total_len:
		left = total_len - len(sequence)
		for j in range(left):
			sequence += 'Z'

	return sequence, sequence[2:25], sequence[4:23], sequence[6:21], sequence[8:19], sequence[10: 17]
