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

for i in range(len(pdb)):
	classes, prediction, pred_max = predict_class(pdb[i], res[i], chain[i])
	val = str(classes[pred_max])
	print(val)
	g.write(val + "\n")

g.close()