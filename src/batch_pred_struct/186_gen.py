from get_186_features import get_features
import pandas as pd 

pdb = []
res = []
chain = []
f = open("data.txt", "r")
for x in f:
	x = x.split(',')
	print(x[0], x[1], x[2])
	pdbid = x[0]
	pdbid = pdbid.replace(" ", "")
	ch = x[2]
	ch = ch.replace(" ", "")
	pdb.append(pdbid)
	res.append(int(x[1]))
	chain.append(ch)

features_df = []

for i in range(len(pdb)):
	features = get_features(pdb[i], res[i], chain[i])
	features_df.append(features)

# Save as xlsx
features_df = pd.DataFrame(features_df)
print(features_df)
features_df.to_csv("Batch_Pred_Features.csv")
