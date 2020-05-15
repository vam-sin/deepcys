# Libraries
import pandas as pd
import numpy as np
import keras
from keras.models import Model
from tensorflow.keras.models import load_model
import keras.backend as K
from keras import optimizers
from keras.layers import Dense, Dropout, BatchNormalization, Conv1D, Flatten, Input, GaussianNoise, LeakyReLU, Add
from keras.utils import to_categorical, np_utils
from keras.regularizers import l2
from sklearn.model_selection import train_test_split, KFold, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder, normalize
from sklearn.utils import shuffle
from sklearn.metrics import confusion_matrix, accuracy_score
import matplotlib.pyplot as plt
import pickle

# Tasks

# dataset import and preprocessing
ds = pd.read_excel('../data/dataset.xlsx')
X1 = ds.iloc[:,4:6] # ('pKa', 'BF')
X1 = normalize(X1, norm = 'l2')
X1 = pd.DataFrame(X1)
X2 = ds.iloc[:,6:7] # rHpy
X = pd.concat([X1, X2], axis=1, sort=False)
X = X.fillna(X.mean())
y = ds.iloc[:,7]

# Seed
seed = 1337
np.random.seed(1337)

# Features
# Secondary Structure Folds (NF1)
infile = open('../features/NF1/feature/NF1_13.pickle','rb')
nf1_13 = pickle.load(infile)
infile.close()

# # Amino Acid Signatures in the Interaction Shells (NF2)
infile = open('../features/NF1/feature/NF2_8.pickle','rb')
nf2_8 = pickle.load(infile)
infile.close()
infile = open('../features/NF1/feature/NF2_7.pickle','rb')
nf2_7 = pickle.load(infile)
infile.close()
infile = open('../features/NF1/feature/NF2_6.pickle','rb')
nf2_6 = pickle.load(infile)
infile.close()
infile = open('../features/NF1/feature/NF2_5.pickle','rb')
nf2_5 = pickle.load(infile)
infile.close()

# # # Protein Class features (NF3)
infile = open('../features/NF1/feature/NF3.pickle','rb')
nf3 = pickle.load(infile)
infile.close()

# # # # Motif (NF4)
infile = open('../features/NF1/feature/NF4_13.pickle','rb')
nf4_13 = pickle.load(infile)
infile.close()

infile = open('../features/NF1/feature/NF4_11.pickle','rb')
nf4_11 = pickle.load(infile)
infile.close()

infile = open('../features/NF1/feature/NF4_9.pickle','rb')
nf4_9 = pickle.load(infile)
infile.close()

infile = open('../features/NF1/feature/NF4_7.pickle','rb')
nf4_7 = pickle.load(infile)
infile.close()

infile = open('../features/NF1/feature/NF4_5.pickle','rb')
nf4_5 = pickle.load(infile)
infile.close()

infile = open('../features/NF1/feature/NF4_3.pickle','rb')
nf4_3 = pickle.load(infile)
infile.close()

# # # # Protein Primary Sequence features (NF5)
infile = open('../features/NF1/feature/NF5_13.pickle','rb')
nf5_13 = pickle.load(infile)
infile.close()

# Feature Selection
nf1_13 = pd.DataFrame(nf1_13)
nf2_8 = pd.DataFrame(nf2_8)
nf2_7 = pd.DataFrame(nf2_7)
nf2_6 = pd.DataFrame(nf2_6)
nf2_5 = pd.DataFrame(nf2_5)
nf3 = pd.DataFrame(nf3)
nf4_3 = pd.DataFrame(nf4_3)
nf4_5 = pd.DataFrame(nf4_5)
nf4_7 = pd.DataFrame(nf4_7)
nf4_9 = pd.DataFrame(nf4_9)
nf4_11 = pd.DataFrame(nf4_11)
nf4_13 = pd.DataFrame(nf4_13)
nf5_13 = pd.DataFrame(nf5_13)

X = pd.concat([X, nf1_13, nf2_5, nf2_6, nf2_7, nf2_8, nf3, nf4_13, nf4_11, nf4_9, nf4_7, nf4_5, nf4_3, nf5_13], axis=1, sort=False)

# convert integers to dummy variables (i.e. one hot encoded)
encoder = LabelEncoder()
encoder.fit(y)
encoded_Y = encoder.transform(y)
y = np_utils.to_categorical(encoded_Y)
print("The classes in y are: " + str(encoder.classes_))

X, y = shuffle(X, y, random_state = 42)

X_train, X_test, y_train, y_test = train_test_split(X.values, y, test_size=0.1, random_state = 42)

# Conv1D Layers
X_train = np.expand_dims(X_train,axis=2)
X_test = np.expand_dims(X_test,axis=2)

NN (Skip Connections) Model
input_ = Input(shape = (len(X.columns),1,))
# x = GaussianNoise(0.5)(input_)
x = Conv1D(128, (3), padding = 'same')(input_)
x = LeakyReLU(alpha = 0.05)(x)
x = Conv1D(64, (3), padding = 'same')(x)
x = GaussianNoise(0.5)(x)
x = LeakyReLU(alpha = 0.05)(x)
x = Conv1D(32, (3), padding = 'same')(x)
x = LeakyReLU(alpha = 0.05)(x)

# Skip Connection #1
input_c = Conv1D(32, (3), padding = 'same')(input_)
skc = Add()([input_c, x])
skc_f = Flatten()(skc)

x = Dense(512, kernel_initializer = 'glorot_uniform')(skc_f)
# x = GaussianNoise(0.5)(x) # For Higher Overfitting
x = LeakyReLU(alpha = 0.05)(x)
x = BatchNormalization()(x)

# Block 2
z = Dense(256, kernel_initializer = 'glorot_uniform')(x)
z = LeakyReLU(alpha = 0.05)(z)
z = Dropout(0.2)(z)
z = BatchNormalization()(z)

# Block 3
x = Dense(128, kernel_initializer = 'glorot_uniform')(z)
x = LeakyReLU(alpha = 0.05)(x)
x = Dropout(0.2)(x)
x = BatchNormalization()(x) 

# Block 4
# Skip connection #2
z = Dense(128, kernel_initializer = 'glorot_uniform')(skc_f)
sk1 = Add()([x, z])
z = Dense(64, kernel_initializer = 'glorot_uniform')(sk1)
z = LeakyReLU(alpha = 0.05)(z)
z = Dropout(0.2)(z)
z = BatchNormalization()(z)

# Block 5
x = Dense(32, kernel_initializer = 'glorot_uniform')(z)
x = LeakyReLU(alpha = 0.05)(x)
x = Dropout(0.2)(x)
x = BatchNormalization()(x) 

# Block 6
z = Dense(16, kernel_initializer = 'glorot_uniform')(x)
z = LeakyReLU(alpha = 0.05)(z)
z = Dropout(0.2)(z)
z = BatchNormalization()(z)

# Block 7
# Skip Connection #3
x = Dense(16, kernel_initializer = 'glorot_uniform')(skc_f)
sk2 = Add()([z, x])
x = Dense(8, kernel_initializer = 'glorot_uniform')(sk2)
x = LeakyReLU(alpha = 0.05)(x)
x = BatchNormalization()(x) 
x = Dropout(0.5)(x)
x = Dense(4, activation = 'softmax')(x)

model = Model(input_, x)
model.compile(optimizer = 'adam', loss = 'categorical_crossentropy', metrics=['accuracy'])

# callbacks
mcp_save = keras.callbacks.callbacks.ModelCheckpoint('ann.h5', save_best_only=True, monitor='val_accuracy', verbose=1)
reduce_lr = keras.callbacks.callbacks.ReduceLROnPlateau(monitor='val_accuracy', factor=0.1, patience=10, verbose=1, mode='auto', min_delta=0.0001, cooldown=0, min_lr=0)
callbacks_list = [reduce_lr, mcp_save]

# Training
history = model.fit(X_train, y_train, batch_size = 256, epochs = 100, validation_data = (X_test, y_test), shuffle = False, callbacks = callbacks_list)

# Testing
model = load_model('ann.h5')
eval = model.evaluate(x = X, y = y)
print("Loss: " + str(eval[0]) + ", Accuracy: " + str(eval[1]))

# # Plot History
# # summarize history for accuracy
plt.plot(history.history['accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()
# summarize history for loss
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()