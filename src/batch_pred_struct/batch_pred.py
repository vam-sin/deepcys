from tensorflow.keras.models import load_model
import pandas as pd
import numpy as np 

# Uncomment for GPU
from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession


config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

ds = pd.read_csv("Batch_Pred_Features.csv")

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

ds = ds.iloc[:,1:]
print(ds)
model_path = 'ann.h5'
classifier = load_model(model_path)
ds = np.expand_dims(ds.values,axis=2)
prediction = classifier.predict(ds)

g = open("results.txt", "w")
m = open("mod.txt", "w")
g.write("Batch Prediction Results \n\n")
g.write("PDB-ID | Residue Number | Chain | Modification \n\n")

for i in range(len(prediction)):
	pred = list(prediction[i])
	for j in range(len(pred)):
		pred[j] *= 100 
		pred[j] = float("{:.2f}".format(pred[j]))
	# Print Results
	classes = ['Disulphide', 'Sulphenylation', 'Metal-Binding', 'Thioether']
	pred_max = pred.index(max(pred))
	print(i, classes[pred_max])
	
	g.write(pdb[i])
	g.write(', ')
	g.write(str(res[i]))
	g.write(', ')
	g.write(ch[i])
	g.write(', ')
	g.write(str(classes[pred_max]))
	g.write("\n")

	m.write(str(classes[pred_max]))
	m.write("\n")

g.close()
m.close()
