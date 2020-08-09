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
from sklearn.metrics import confusion_matrix, accuracy_score, f1_score
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.utils import shuffle
from keras.utils import np_utils
from sklearn.model_selection import train_test_split, KFold, cross_val_score, StratifiedKFold
import pickle
from sklearn.utils import class_weight

# Notes
# GRU > LSTM
# GPU
import tensorflow as tf
from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

tf.keras.backend.clear_session()



config = ConfigProto()
config.gpu_options.allow_growth = True
gpu_options = tf.compat.v1.GPUOptions(per_process_gpu_memory_fraction=0.333)

sess = tf.compat.v1.Session(config=tf.compat.v1.ConfigProto(gpu_options=gpu_options))

LIMIT = 3 * 1024
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        tf.config.experimental.set_virtual_device_configuration(
            gpus[0],
            [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=LIMIT)])
        logical_gpus = tf.config.experimental.list_logical_devices('GPU')
        print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
    except RuntimeError as e:
        # Virtual devices must be set before GPUs have been initialized
        print(e)

# Dataset
ds = pd.read_csv('../data/correct_data/dataset.csv')
y = ds.iloc[:,7]
y_w = y

# Seed
seed = 42
np.random.seed(seed)

window = 11
filename = '../features/NF5/feature/NF5_11.pickle'
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

classifier = NN()

encoder = LabelEncoder()
encoder.fit(y)
y = encoder.transform(y)
y = np_utils.to_categorical(y)

# Make Test Data
metal_X = X.iloc[18759:18859,:]
X = X.drop(X.index[18759:18859])
y_ = pd.DataFrame(y)
metal_y = y_.iloc[18759:18859]
y_ = y_.drop(y_.index[18759:18859])
y_w = y_w.drop(y_w.index[18759:18859])

sulph_X = X.iloc[18961:19061,:]
sulph_y = y_.iloc[18961:19061]
X = X.drop(X.index[18961:19061])
y_ = y_.drop(y_.index[18961:19061])
y_w = y_w.drop(y_w.index[18961:19061])

dis_X = X.iloc[55170:55270,:]
dis_y = y_.iloc[55170:55270]
X = X.drop(X.index[55170:55270])
y_ = y_.drop(y_.index[55170:55270])
y_w = y_w.drop(y_w.index[55170:55270])

thio_X = X.iloc[107260:107360,:]
thio_y = y_.iloc[107260:107360]
X = X.drop(X.index[107260:107360])
y_ = y_.drop(y_.index[107260:107360])
y_w = y_w.drop(y_w.index[107260:107360])

X_test = pd.concat([metal_X, sulph_X, dis_X, thio_X], axis = 0)
y_test = np.asarray(pd.concat([metal_y, sulph_y, dis_y, thio_y], axis = 0))
y_train = np.asarray(y_) 
X_train = X

y_w = list(y_w.values)
encoder = LabelEncoder()
encoder.fit(y_w)
y_w = encoder.transform(y_w)
class_weights = class_weight.compute_class_weight('balanced',
                                                 np.unique(y_w),
                                                 y_w)
print(class_weights, np.unique(y_w))

X_train, y_train = shuffle(X_train, y_train, random_state = 42)

y_ = pd.DataFrame(y_)
X_train = np.expand_dims(X_train,axis=2)
X_test = np.expand_dims(X_test,axis=2)

mcp_save = keras.callbacks.callbacks.ModelCheckpoint('gru.h5', save_best_only=True, monitor='val_accuracy', verbose=1)
reduce_lr = keras.callbacks.callbacks.ReduceLROnPlateau(monitor='val_accuracy', factor=0.1, patience=10, verbose=1, mode='auto', min_delta=0.0001, cooldown=0, min_lr=0)
callbacks_list = [reduce_lr, mcp_save]

weights = {0: 0.31614373, 1: 1.43080227, 2: 46.60362694, 3: 8.58253817}
history = classifier.fit(X_train, y_train, epochs=100, batch_size=256, verbose=1, validation_data = (X_test, y_test), callbacks = callbacks_list, class_weight = weights)

classifier = load_model('gru.h5')
scores = classifier.evaluate(X_test, y_test, verbose=1)
print("Loss: " + str(scores[0]) + ", Accuracy: " + str(scores[1]))
print(scores)

y_pred = classifier.predict(X_test)

# Metrics
print("Confusion Matrix")
matrix = confusion_matrix(y_test.argmax(axis=1), y_pred.argmax(axis=1))
print(matrix)

print("F1 Score")
print(f1_score(y_test.argmax(axis=1), y_pred.argmax(axis=1), average = 'weighted'))

'''
Results:
WS, Epochs, BS, VA
3, 50, 256, loss: 0.4601 - accuracy: 0.7677 - val_loss: 0.6939 - val_accuracy: 0.7825
5, 50, 256, loss: 0.3065 - accuracy: 0.8245 - val_loss: 0.8791 - val_accuracy: 0.8050
7, 50, 256, loss: 0.4637 - accuracy: 0.7285 - val_loss: 0.7555 - val_accuracy: 0.7375
9, 50, 256, loss: 0.3675 - accuracy: 0.7923 - val_loss: 0.7294 - val_accuracy: 0.7800
11, 50, 256, 
13, 50, 256,
'''


