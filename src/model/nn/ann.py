# Libraries
import pandas as pd
import numpy as np
import keras
import tensorflow as tf
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
from sklearn.metrics import confusion_matrix, accuracy_score, f1_score
import matplotlib.pyplot as plt
import pickle
from keras import regularizers
from keras import backend as K
from sklearn.utils import class_weight

# GPU
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

# dataset import and preprocessing
ds = pd.read_csv("../dataset/dataset.csv")

X1 = ds.iloc[:,5:6] # (BF)
X1 = pd.DataFrame(X1)
X2 = ds.iloc[:,6:7] # rHpy
X = pd.concat([X1, X2], axis=1)
y = ds.iloc[:,7]
y_w = y

# Seed
seed = 1337
np.random.seed(1337)

# Features
# Secondary Structure Folds (NF1)
infile = open('../features/NF1/feature/NF1_7.pickle','rb')
nf1_9 = pickle.load(infile)
infile.close()

# # Amino Acid Signatures in the Interaction Shells (NF2)
infile = open('../features/NF2/feature/NF2_8.pickle','rb')
nf2_8 = pickle.load(infile)
infile.close()
infile = open('../features/NF2/feature/NF2_7.pickle','rb')
nf2_7 = pickle.load(infile)
infile.close()
infile = open('../features/NF2/feature/NF2_6.pickle','rb')
nf2_6 = pickle.load(infile)
infile.close()
infile = open('../features/NF2/feature/NF2_5.pickle','rb')
nf2_5 = pickle.load(infile)
infile.close()

# # # Protein Class features (NF3)
infile = open('../features/NF3/feature/NF3_le.pickle','rb')
nf3 = pickle.load(infile)
infile.close()

# # # # Motif (NF4)
infile = open('../features/NF4/feature/NF4_13.pickle','rb')
nf4_13 = pickle.load(infile)
infile.close()

infile = open('../features/NF4/feature/NF4_11.pickle','rb')
nf4_11 = pickle.load(infile)
infile.close()

infile = open('../features/NF4/feature/NF4_9.pickle','rb')
nf4_9 = pickle.load(infile)
infile.close()

infile = open('../features/NF4/feature/NF4_7.pickle','rb')
nf4_7 = pickle.load(infile)
infile.close()

infile = open('../features/NF4/feature/NF4_5.pickle','rb')
nf4_5 = pickle.load(infile)
infile.close()

infile = open('../features/NF4/feature/NF4_3.pickle','rb')
nf4_3 = pickle.load(infile)
infile.close()

# Feature Selection
nf1_9 = pd.DataFrame(nf1_9)
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

X = pd.concat([X, nf1_9, nf2_5, nf2_6, nf2_7, nf2_8, nf3, nf4_13, nf4_11, nf4_9, nf4_7, nf4_5, nf4_3], axis=1, sort=False)

# convert integers to dummy variables (i.e. one hot encoded)
encoder = LabelEncoder()
encoder.fit(y)
encoded_Y = encoder.transform(y)
counts_df = pd.DataFrame(encoded_Y)
y = np_utils.to_categorical(encoded_Y)
print("The classes in y are: " + str(encoder.classes_))

# Make Test Data
one = datapoints.iloc[18759:18859,:]
metal_X = X.iloc[18759:18859,:]
X = X.drop(X.index[18759:18859])
y_ = pd.DataFrame(y)
metal_y = y_.iloc[18759:18859]
y_ = y_.drop(y_.index[18759:18859])
y_w = y_w.drop(y_w.index[18759:18859])

two = datapoints.iloc[18961:19061,:]
sulph_X = X.iloc[18961:19061,:]
sulph_y = y_.iloc[18961:19061]
X = X.drop(X.index[18961:19061])
y_ = y_.drop(y_.index[18961:19061])
y_w = y_w.drop(y_w.index[18961:19061])

three = datapoints.iloc[55170:55270,:]
dis_X = X.iloc[55170:55270,:]
dis_y = y_.iloc[55170:55270]
X = X.drop(X.index[55170:55270])
y_ = y_.drop(y_.index[55170:55270])
y_w = y_w.drop(y_w.index[55170:55270])

four = datapoints.iloc[107260:107360,:]
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

# NN (Skip Connections) Model
input_ = Input(shape = (len(X.columns),1,))
x = Conv1D(128, (3), padding = 'same', kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(input_)
x = LeakyReLU(alpha = 0.05)(x)
x = Conv1D(64, (3), padding = 'same', kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)
x = LeakyReLU(alpha = 0.05)(x)
x = Conv1D(32, (3), padding = 'same', kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)
x = LeakyReLU(alpha = 0.05)(x)

# Skip Connection #1
input_c = Conv1D(32, (3), padding = 'same', kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(input_)
skc = Add()([input_c, x])
skc_f = Flatten()(skc)

x = Dense(512, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(skc_f)
x = LeakyReLU(alpha = 0.05)(x)
x = BatchNormalization()(x)

# Block 2
z = Dense(256, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)
z = LeakyReLU(alpha = 0.05)(z)
# z = Dropout(0.1)(z)
z = BatchNormalization()(z)

# Block 3
x = Dense(128, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(z)
x = LeakyReLU(alpha = 0.05)(x)
# x = Dropout(0.1)(x)
x = BatchNormalization()(x) 

# Block 4
# Skip connection #2
z = Dense(128, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(skc_f)
sk1 = Add()([x, z])
z = Dense(64, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(sk1)
z = LeakyReLU(alpha = 0.05)(z)
# z = Dropout(0.1)(z)
z = BatchNormalization()(z)

# Block 5
x = Dense(32, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(z)
x = LeakyReLU(alpha = 0.05)(x)
# x = Dropout(0.1)(x)
x = BatchNormalization()(x) 

# Block 6
z = Dense(16, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)
z = LeakyReLU(alpha = 0.05)(z)
# z = Dropout(0.1)(z)
z = BatchNormalization()(z)

# Block 7
# Skip Connection #3
x = Dense(16, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(skc_f)
sk2 = Add()([z, x])
x = Dense(8, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(sk2)
x = LeakyReLU(alpha = 0.05)(x)
x = BatchNormalization()(x) 
x = Dropout(0.5)(x)
x = Dense(4, activation = 'softmax', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)

model = Model(input_, x)

opt = keras.optimizers.Adam(learning_rate = 1e-4)
model.compile(optimizer = opt, loss = 'categorical_crossentropy', metrics=['accuracy'])

# callbacks
mcp_save = keras.callbacks.callbacks.ModelCheckpoint('ann.h5', save_best_only=True, monitor='val_accuracy', verbose=1)
reduce_lr = keras.callbacks.callbacks.ReduceLROnPlateau(monitor='val_accuracy', factor=0.1, patience=5, verbose=1, mode='auto', min_delta=0.0001, cooldown=0, min_lr=0)
callbacks_list = [reduce_lr, mcp_save]

# Training
weights = {0: 0.9, 1: 3.5, 2: 48.5, 3: 19.5}
history = model.fit(X_train, y_train, batch_size = 256, epochs = 50, validation_data = (X_test, y_test), shuffle = False, callbacks = callbacks_list, class_weight = weights)

# Testing
model = load_model('ann.h5')
eval_ = model.evaluate(x = X_test, y = y_test)
print("Loss: " + str(eval_[0]) + ", Accuracy: " + str(eval_[1]))
print(eval_)

y_pred = model.predict(X_test)

# Metrics
print("Confusion Matrix")
matrix = confusion_matrix(y_test.argmax(axis=1), y_pred.argmax(axis=1))
print(matrix)

print("F1 Score")
print(f1_score(y_test.argmax(axis=1), y_pred.argmax(axis=1), average = 'weighted'))

# Plot History
# summarize history for accuracy
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


# Benchmark
# infile = open('benchmark/dis2to25.pickle','rb')
# features_X = pickle.load(infile)
# infile.close()

# print(features_X.shape)

# features_X = features_X[:,1:151]
# features_X = np.delete(features_X, 2, 1)
# features_X = np.delete(features_X, 2, 1)
# features_X = np.delete(features_X, 17, 1)
# features_X = np.delete(features_X, 17, 1)
# drop_cols = [2,3,19,20]

# print(features_X.shape)

# features_X = np.expand_dims(features_X, axis=2)

# classes_dict = {"0": "disulphide", "1": "metal-binding", "2": "sulphenylation", "3": "thioether"}

# y_pred = model.predict(features_X)

# y_pred = np.asarray(y_pred)

# for i in y_pred:
#     pred_num = np.argmax(i)

# g = open("benchmark/results_dis.txt", "w")
# # g.write("##### DeepCys Batch Prediction Results #####\n\n")
# for i in range(len(y_pred)):
#     text = classes_dict[str(np.argmax(y_pred[i]))]
#     g.write(str(text))
#     g.write("\n")

# g.close()

# infile = open('benchmark/thioether2to25.pickle','rb')
# features_X = pickle.load(infile)
# infile.close()
# features_X = features_X[:,1:151]
# features_X = np.delete(features_X, 2, 1)
# features_X = np.delete(features_X, 2, 1)
# features_X = np.delete(features_X, 17, 1)
# features_X = np.delete(features_X, 17, 1)
# features_X = np.expand_dims(features_X, axis=2)

# classes_dict = {"0": "disulphide", "1": "metal-binding", "2": "sulphenylation", "3": "thioether"}

# y_pred = model.predict(features_X)

# y_pred = np.asarray(y_pred)

# for i in y_pred:
#     pred_num = np.argmax(i)

# g = open("benchmark/results_thio.txt", "w")
# # g.write("##### DeepCys Batch Prediction Results #####\n\n")
# for i in range(len(y_pred)):
#     text = classes_dict[str(np.argmax(y_pred[i]))]
#     g.write(str(text))
#     g.write("\n")

# g.close()

# infile = open('benchmark/CSO2to25.pickle','rb')
# features_X = pickle.load(infile)
# infile.close()
# features_X = features_X[:,1:151]
# features_X = np.delete(features_X, 2, 1)
# features_X = np.delete(features_X, 2, 1)
# features_X = np.delete(features_X, 17, 1)
# features_X = np.delete(features_X, 17, 1)
# features_X = np.expand_dims(features_X, axis=2)

# classes_dict = {"0": "disulphide", "1": "metal-binding", "2": "sulphenylation", "3": "thioether"}

# y_pred = model.predict(features_X)

# y_pred = np.asarray(y_pred)

# for i in y_pred:
#     pred_num = np.argmax(i)

# g = open("benchmark/results_cso.txt", "w")
# # g.write("##### DeepCys Batch Prediction Results #####\n\n")
# for i in range(len(y_pred)):
#     text = classes_dict[str(np.argmax(y_pred[i]))]
#     g.write(str(text))
#     g.write("\n")

# g.close()

# infile = open('benchmark/metalbinding2to25.pickle','rb')
# features_X = pickle.load(infile)
# infile.close()
# features_X = features_X[:,1:151]
# features_X = np.delete(features_X, 2, 1)
# features_X = np.delete(features_X, 2, 1)
# features_X = np.delete(features_X, 17, 1)
# features_X = np.delete(features_X, 17, 1)
# features_X = np.expand_dims(features_X, axis=2)

# classes_dict = {"0": "disulphide", "1": "metal-binding", "2": "sulphenylation", "3": "thioether"}

# y_pred = model.predict(features_X)

# y_pred = np.asarray(y_pred)

# for i in y_pred:
#     pred_num = np.argmax(i)

# g = open("benchmark/results_mb.txt", "w")
# # g.write("##### DeepCys Batch Prediction Results #####\n\n")
# for i in range(len(y_pred)):
#     text = classes_dict[str(np.argmax(y_pred[i]))]
#     g.write(str(text))
#     g.write("\n")

# g.close()

# infile = open('benchmark/DUF.pickle','rb')
# features_X = pickle.load(infile)
# infile.close()
# features_X = features_X[:,1:151]
# features_X = np.delete(features_X, 2, 1)
# features_X = np.delete(features_X, 2, 1)
# features_X = np.delete(features_X, 17, 1)
# features_X = np.delete(features_X, 17, 1)
# features_X = np.expand_dims(features_X, axis=2)

# classes_dict = {"0": "disulphide", "1": "metal-binding", "2": "sulphenylation", "3": "thioether"}

# y_pred = model.predict(features_X)

# y_pred = np.asarray(y_pred)

# for i in y_pred:
#     pred_num = np.argmax(i)

# g = open("results_duf.txt", "w")
# g.write("##### DeepCys Batch Prediction Results #####\n\n")
# for i in range(len(y_pred)):
#     text = classes_dict[str(np.argmax(y_pred[i]))]
#     g.write(str(text))
#     g.write("\n")

# g.close()

# # Benchmark Accuracy values
# f = open("benchmark/results_dis.txt","r")
# dis_b = []
# count = 0
# for line in f:
#     l = line.replace('\n','')
#     dis_b.append(l)
#     if l == 'disulphide':
#         count += 1

# print("Disulphide Benchmark: ", str(count/len(dis_b)))

# f = open("benchmark/results_mb.txt","r")
# mb_b = []
# count = 0
# for line in f:
#     l = line.replace('\n','')
#     mb_b.append(l)
#     if l == 'metal-binding':
#         count += 1

# print("Metal-Binding Benchmark: ", str(count/len(mb_b)))

# f = open("benchmark/results_thio.txt","r")
# thio_b = []
# count = 0
# for line in f:
#     l = line.replace('\n','')
#     thio_b.append(l)
#     if l == 'thioether':
#         count += 1

# print("Thioether Benchmark: ", str(count/len(thio_b)))

# f = open("benchmark/results_cso.txt","r")
# sul_b = []
# count = 0
# for line in f:
#     l = line.replace('\n','')
#     sul_b.append(l)
#     if l == 'sulphenylation':
#         count += 1

# f = open("benchmark/results_sulph.txt","r")
# for line in f:
#     l = line.replace('\n','')
#     sul_b.append(l)
#     if l == 'sulphenylation':
#         count += 1

# print("Sulphenylation Benchmark: ", str(count/len(sul_b)))

