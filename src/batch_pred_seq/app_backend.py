# Libraries
import requests 
import numpy as np
from tensorflow.keras.models import load_model
import pickle
import os
from feature_gen import get_nf1, get_nf2, get_nf3, get_nf4, get_nf5
from pka_website import get_pka
from Bf_rhpy_website import get_bf_rhpy
import keras.backend.tensorflow_backend as tb
tb._SYMBOLIC_SCOPE.value = True

# Take in PDB ID and residue ID. (Ex: 1b2l, 137, A, Output: Disulfide)

def predict_class(pdb, res, chain):
	# Parameters:
	res = int(res)

	# Get FASTA and PDB.
	PROJECT_PATH = os.path.dirname(__file__) + "/"
	print(PROJECT_PATH)
	print("\nSteps.")

	url = 'https://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=fastachain&compression=NO&structureId=' + pdb.lower() + '&chainId=' + chain
	r = requests.get(url)
	filename_fasta = PROJECT_PATH + '/PDB_fasta/' + pdb + '.fasta'
	open(filename_fasta, 'wb').write(r.content)
	print("Obtained fasta.")

	# NF5 (Works)
	nf5_13 = get_nf5(pdb, res, chain, 13)

	# Compile X
	X = []

	for i in nf5_13:
		X.append(i)

	X = np.asarray(X)
	dim1 = len(X)
	dim2 = 13*2 + 1

	X = X.reshape(1, dim2, 1)
	X = X.astype(np.float32)

	# # Load Model and Predict
	model_path = PROJECT_PATH + '/gru.h5'
	classifier = load_model(model_path)
	prediction = classifier.predict(X)
	prediction = prediction[0]
	prediction = list(prediction)
	for i in range(len(prediction)):
		prediction[i] *= 100 
		prediction[i] = float("{:.2f}".format(prediction[i]))
	# Print Results
	classes = ['Disulphide', 'Sulphenylation', 'Metal-Binding', 'Thioether']
	print("\nProbaility Results:-")
	print("Disulphide: " + str(prediction[0]) + '%')
	print("Selphenylation: " + str(prediction[1]) + '%')
	print("Metal-Binding: " + str(prediction[2]) + '%')
	print("Thioether: " + str(prediction[3]) + '%')
	print("\n")
	pred_max = prediction.index(max(prediction))
	print("Highest Probability to be " + classes[pred_max]) 

	return classes, prediction, pred_max

if __name__ == '__main__':
	pdb = input("Enter PDB ID: ").upper()
	res = int(input("Enter Residue ID: "))
	chain = input("Enter Chain Alphabet: ").upper()

	predict_class(pdb, res, chain)


