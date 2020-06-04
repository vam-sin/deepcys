import pandas as pd 

from Bio import SeqIO

ds = pd.read_excel('../data/dataset.xlsx')
print(ds)
ds.columns = ['no', 'pdb', 'res', 'ch', 'pka', 'bf', 'rhpy', 'mod']

pdb = list(ds['pdb'].values)
# print(pdb)
res = list(ds['res'].values)
# print(res)
ch = list(ds['ch'].values)
# print(ch)
mod = list(ds['mod'].values)
changed = 0
tot = 0
g = open('residue_out_of_range.txt', "w")
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
		if res[i] >= len(seq):
			print(pdb[i], res[i], ch[i])
			g.write(pdb[i])
			g.write(', ')
			g.write(str(res[i]))
			g.write(', ')
			g.write(ch[i])
			g.write("\n")
			tot+=1

	except:
		# print("Cant Help")
		pass

print("Total Fixes: " + str(changed))
print("Total: " + str(tot))

# print(res)
# res = pd.DataFrame(res)

# ds['res'] = res
# ds = ds.iloc[:,1:]
# print(ds)
# ds.to_excel('../data/dataset.xlsx')

# four