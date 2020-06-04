# Structure to Sequence residue number variations

import pandas as pd 
from Bio import SeqIO

ds = pd.read_excel('../data/dataset.xlsx')
# print(ds)
ds.columns = ['no', 'pdb', 'res', 'ch', 'pka', 'bf', 'rhpy', 'mod']

pdb = list(ds['pdb'].values)
# print(pdb)
res = list(ds['res'].values)
# print(res)
ch = list(ds['ch'].values)
# print(ch)
mod = list(ds['mod'].values)
changed = 0
for i in range(1, len(pdb)):
	filename_fasta = '../../../fasta/' + pdb[i].upper() + '.fasta.txt'
	record = list(SeqIO.parse(filename_fasta, "fasta"))
	list_ind = 0
	for j in range(len(record)):
		if type(ch[i]) == str:	
			ind = record[j].id[5]
			if ind == ch[i]:
				list_ind = j
				break
		else:
			list_ind = int(ch[i]) - 1

	seq = record[list_ind].seq

	try:
		nex = 0
	# print(seq[res[i]-1])
		if seq[res[i]-1] != 'C':
			print(pdb[i], res[i], ch[i])
			filename_pdb = '../../../pdb/' + pdb[i].lower() + '.pdb'
			pdb_file = open(filename_pdb, "r")
			ind = 0
			ssqind = 0
			for line in pdb_file:
				if nex == 1:
					line = line.split(' ')
					switch = int(line[10] + line[11] + line[12])
					print(switch, seq[res[i]-switch], seq[res[i]-switch-5:res[i]-switch+5])
					if seq[res[i] - switch] == 'C':
						res[i] = res[i] - switch + 1
						changed += 1
					break
				if 'M RES C SSSEQI' in line and nex == 0:
					ssqind = ind
					nex = 1
					print("YAAS " + str(ind))
				ind += 1

	except:
		# print("Cant Help")
		pass

print("Total Fixes: " + str(changed))

# print(res)
res = pd.DataFrame(res)

ds['res'] = res
ds = ds.iloc[:,1:]
# print(ds)
# ds.to_excel('../data/dataset.xlsx')
