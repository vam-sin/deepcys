from app_backend import predict_class

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

g = open("results.txt", "w")

g.write("PDB, Residue Number, Chain Letter, Cysteine Modification \n\n")

for i in range(len(pdb)):
	classes, prediction, pred_max = predict_class(pdb[i], res[i], chain[i])
	val = str(classes[pred_max])
	pred_acc = prediction[pred_max]
	ch = chain[i].replace("\n", "")
	g.write(pdb[i] + ", " + str(res[i]) + ", " + ch + ": " + val + " , Prediction Accuracy: " + str(pred_acc))
	g.write("\n")

g.close()
