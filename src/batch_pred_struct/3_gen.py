from get_3_features import get_features
import pandas as pd 

# pdb = []
# res = []
# chain = []

# f = open("data.txt", "r")
# for x in f:
# 	x = x.split(',')
# 	print(x[0], x[1], x[2])
# 	pdbid = x[0]
# 	pdbid = pdbid.replace(" ", "")
# 	ch = x[2]
# 	ch = ch.replace(" ", "")
# 	pdb.append(pdbid)
# 	res.append(int(x[1]))
# 	chain.append(ch)

ds = pd.read_csv('get_three_features.csv')
pdb = list(ds.iloc[:,1:2].values)
res = list(ds.iloc[:,2:3].values)
chain = list(ds.iloc[:,3:4].values)
mod = list(ds.iloc[:,4:5].values)
print(mod)

new_pdb = []

for i in pdb:
	j = i[0]
	j = j.replace('.pdb', '')
	# print(j)
	new_pdb.append(j)

features_df = []

for i in range(len(pdb)):
	print("########################")
	print(i, " Out of ", len(pdb))
	print("########################")
	features = get_features(new_pdb[i], res[i][0], chain[i][0], mod[i][0])
	features_df.append(features)

# Save as xlsx
features_df = pd.DataFrame(features_df)
print(features_df)
features_df.to_csv("Batch_Pred_Features.csv")

# 3 per minute
