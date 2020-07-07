# Get required sequence from the pdb structure

def three_to_one(aa):
	d = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER':'*',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

	return d[aa]

def get_sequence(pdb, res, chain, window):

	pdb_file = open(pdb, "r")
	sequence = ''
	start = res - window
	end = res + window
	if start <= 0:
		print("Start error")
		for i in range(abs(start)+1):
			sequence += 'Z'
		start = 1

	for i in range(start, end+1):
	# try:
		for line in pdb_file:
			line = line.split(' ')
			new_line = []
			for seg in line:
				if seg != '':
					new_line.append(seg)

			if new_line[0] == 'ATOM':
				try:
					res_num = int(new_line[5])
					aa = new_line[3]
					ch = new_line[4]
					if i == res_num and chain == ch:
						# print(i)
						sequence += three_to_one(aa)
						break 
				except:
					pass
				
					# except:
					# 	pass
		# except:
		# 	sequence.append('X')

	total_len = window*2 + 1
	if len(sequence) < total_len:
		print("End Error")
		left = total_len - len(sequence)
		for j in range(left):
			sequence += 'Z'

	return sequence

if __name__ == '__main__':
	window = 13
	sequence = get_sequence('2H4W.pdb', 285, 'A', window)
	print(sequence, len(sequence), sequence[window])