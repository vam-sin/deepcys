import pandas as pd 
from Bio import SeqIO

ds = pd.read_excel('../data/dataset.xlsx')
ds.columns = ['no', 'pdb', 'res', 'ch', 'pka', 'bf', 'rhpy', 'mod']

pdb = list(ds['pdb'].values)
# print(pdb)
res = list(ds['res'].values)
# print(res)
ch = list(ds['ch'].values)
# print(ch)
mod = list(ds['mod'].values)
# print(mod)

not_a_cys_pdb = []
not_a_cys_res = []
not_a_cys_ch = []
not_a_cys_mod = []
residue_out_of_range_pdb = []
residue_out_of_range_res = []
residue_out_of_range_ch = []
residue_out_of_range_mod = []

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
		# print(seq[res[i]-1])
		if seq[res[i]-1] != 'C':
			not_a_cys_pdb.append(pdb[i])
			not_a_cys_res.append(res[i])
			not_a_cys_ch.append(ch[i])
			not_a_cys_mod.append(mod[i])
	except:
		pass
	try:
		if res[i] > len(seq):
			residue_out_of_range_pdb.append(pdb[i])
			residue_out_of_range_res.append(res[i])
			residue_out_of_range_ch.append(ch[i])
			residue_out_of_range_mod.append(mod[i])
	except:
		pass

print(len(not_a_cys_pdb), len(residue_out_of_range_pdb))
# print(not_a_cys_mod)
f_cys = open("Not_a_Cysteine.txt", "w")

for i in range(len(not_a_cys_pdb)):
	# print(not_a_cys_pdb[i], not_a_cys_res[i], not_a_cys_ch[i])
#	print(i)
	f_cys.write(not_a_cys_pdb[i])
	f_cys.write(", ")
	f_cys.write(str(not_a_cys_res[i]))
	f_cys.write(", ")
	f_cys.write(str(not_a_cys_ch[i]))
	f_cys.write(", ")
	f_cys.write(str(not_a_cys_mod[i]))
	f_cys.write("\n")

f_cys = open("residue_out_of_range.txt", "w")

for i in range(len(residue_out_of_range_pdb)):
	# print(not_a_cys_pdb[i], not_a_cys_res[i], not_a_cys_ch[i])
	f_cys.write(residue_out_of_range_pdb[i])
	f_cys.write(", ")
	f_cys.write(str(residue_out_of_range_res[i]))
	f_cys.write(", ")
	f_cys.write(str(residue_out_of_range_ch[i]))
	f_cys.write(", ")
	f_cys.write(str(not_a_cys_mod[i]))
	f_cys.write("\n")