# Libraries
import pandas as pd
import numpy as np
import keras
from keras.models import Sequential
from tensorflow.keras.models import load_model
import keras.backend as K
from keras import optimizers
from keras.layers import Dense, Dropout, BatchNormalization, Conv1D, Flatten, MaxPooling1D, Reshape, GRU, SpatialDropout1D
from keras.layers.embeddings import Embedding
from keras.preprocessing import sequence
from sklearn.metrics import confusion_matrix, accuracy_score
from keras.utils import to_categorical
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.utils import shuffle
from keras.utils import np_utils
from sklearn.model_selection import train_test_split, KFold, cross_val_score, StratifiedKFold
import pickle

# Notes
# GRU > LSTM

# Dataset
ds = pd.read_excel('../data/dataset.xlsx')
# X = ds.iloc[:,3:7]
y = ds.iloc[:,7]

# Seed
seed = 1337
np.random.seed(1337)

window = 13
filename = '../features/NF5/feature/NF5_13.pickle'
infile = open(filename,'rb')
pos = pickle.load(infile)
infile.close()
pos = np.asarray(pos)
X = pd.DataFrame(pos)

dim1 = len(pos)
dim2 = window*2 + 1

# Neural Network model
def NN():
	model = Sequential()
	model.add(GRU(100)) # Try with Spatial Dropout
	model.add(Dropout(0.5))
	model.add(Dropout(0.5))
	model.add(Dense(4, activation='softmax'))
	model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

	return model

kfold = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
cvscores = []
for train, test in kfold.split(X, y):

	classifier = NN()

	# mcp_save = keras.callbacks.callbacks.ModelCheckpoint('gru.h5', save_best_only=True, monitor='val_accuracy', verbose = 1)
	# reduce_lr = keras.callbacks.callbacks.ReduceLROnPlateau(monitor='accuracy', factor=0.1, patience=10, verbose=1, mode='auto', min_delta=0.0001, cooldown=0, min_lr=0)
	# callbacks_list = [mcp_save, reduce_lr]

	X = pd.DataFrame(X)

	# 	# evaluate the model

	encoder = LabelEncoder()
	encoder.fit(y)
	y = encoder.transform(y)
	y = pd.DataFrame(y)
	test_ind = []
	for k in test:
		test_ind.append(k)
	train_ind = []
	for k in test:
		train_ind.append(k)
	# print(test_ind)
	X_train = np.asarray(X.iloc[train])
	X_train = X_train.reshape(len(X_train), dim2, 1)
	y_train = np_utils.to_categorical(y.iloc[train])
	X_test = np.asarray(X.iloc[test])
	X_test = X_test.reshape(len(X_test), dim2, 1)
	y_test = np_utils.to_categorical(y.iloc[test])

	mcp_save = keras.callbacks.callbacks.ModelCheckpoint('gru.h5', save_best_only=True, monitor='val_accuracy', verbose=1)
	reduce_lr = keras.callbacks.callbacks.ReduceLROnPlateau(monitor='val_accuracy', factor=0.1, patience=10, verbose=1, mode='auto', min_delta=0.0001, cooldown=0, min_lr=0)
	callbacks_list = [reduce_lr, mcp_save]

	classifier.fit(X_train, y_train, epochs=100, batch_size=16, verbose=1, callbacks = callbacks_list)

	# classifier = load_model('gru.h5')
	# X_test = X_test.astype(np.float32)
	scores = classifier.evaluate(X_test, y_test, verbose=1)
	print("%s: %.2f%%" % (classifier.metrics_names[1], scores[1]*100))
	cvscores.append(scores[1] * 100)

print("%.2f%% (+/- %.2f%%)" % (np.mean(cvscores), np.std(cvscores))) 


