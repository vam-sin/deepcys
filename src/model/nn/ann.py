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

# Tasks

# dataset import and preprocessing
ds = pd.read_excel("../data/balanced_dataset.xlsx")
X1 = ds.iloc[:,4:6] # ('pKa', 'BF')
# X1 = normalize(X1, norm = 'l2')
X1 = pd.DataFrame(X1)
X2 = ds.iloc[:,6:7] # rHpy
X = pd.concat([X1, X2], axis=1)
# print(X)
# X = X.fillna(X.mean())
y = ds.iloc[:,7]
# print(X, y)
# Seed
seed = 1337
np.random.seed(1337)

# Probable error with met
# met = []
# for i in range(len(y)):
#     if i < 6945 and i > 4452:
#         met.append(1)
#     else:
#         met.append(0)

# met = pd.DataFrame(met)

# Features
# Dis
infile = open('../features/Dis/feature/dis.pickle','rb')
dis = pickle.load(infile)
infile.close()

# Secondary Structure Folds (NF1)
infile = open('../features/NF1/feature/NF1_13.pickle','rb')
nf1_13 = pickle.load(infile)
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
infile = open('../features/NF3/feature/NF3.pickle','rb')
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

# # # # Protein Primary Sequence features (NF5)
infile = open('../features/NF5/feature/NF5_13.pickle','rb')
nf5_13 = pickle.load(infile)
infile.close()

# Feature Selection
dis = pd.DataFrame(dis)
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

# Met error
# X = pd.concat([X, met, nf1_13, nf2_5, nf2_6, nf2_7, nf2_8, nf3, nf4_13, nf4_11, nf4_9, nf4_7, nf4_5, nf4_3, nf5_13], axis=1, sort=False)
# print(X)

# convert integers to dummy variables (i.e. one hot encoded)
encoder = LabelEncoder()
encoder.fit(y)
encoded_Y = encoder.transform(y)
y = np_utils.to_categorical(encoded_Y)
print("The classes in y are: " + str(encoder.classes_))
# print(list(encoder.inverse_transform(encoded_Y)))
# print(X, y)
X, y = shuffle(X, y, random_state = 42)

X_train, X_test, y_train, y_test = train_test_split(X.values, y, test_size=0.1, random_state = 42)

# print(y_test)
# Conv1D Layers
X_train = np.expand_dims(X_train,axis=2)
X_test = np.expand_dims(X_test,axis=2)

# NN (Skip Connections) Model
input_ = Input(shape = (len(X.columns),1,))
# x = GaussianNoise(0.5)(input_)
x = Conv1D(128, (3), padding = 'same', kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(input_)
x = LeakyReLU(alpha = 0.05)(x)
x = Conv1D(64, (3), padding = 'same', kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)
x = GaussianNoise(0.5)(x)
x = LeakyReLU(alpha = 0.05)(x)
x = Conv1D(32, (3), padding = 'same', kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)
x = LeakyReLU(alpha = 0.05)(x)

# Skip Connection #1
input_c = Conv1D(32, (3), padding = 'same', kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(input_)
skc = Add()([input_c, x])
skc_f = Flatten()(skc)

x = Dense(512, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(skc_f)
# x = GaussianNoise(0.5)(x) # For Higher Overfitting
x = LeakyReLU(alpha = 0.05)(x)
x = BatchNormalization()(x)

# Block 2
z = Dense(256, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)
z = LeakyReLU(alpha = 0.05)(z)
z = Dropout(0.5)(z)
z = BatchNormalization()(z)

# Block 3
x = Dense(128, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(z)
x = LeakyReLU(alpha = 0.05)(x)
x = Dropout(0.5)(x)
x = BatchNormalization()(x) 

# Block 4
# Skip connection #2
z = Dense(128, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(skc_f)
sk1 = Add()([x, z])
z = Dense(64, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(sk1)
z = LeakyReLU(alpha = 0.05)(z)
z = Dropout(0.5)(z)
z = BatchNormalization()(z)

# Block 5
x = Dense(32, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(z)
x = LeakyReLU(alpha = 0.05)(x)
x = Dropout(0.5)(x)
x = BatchNormalization()(x) 

# Block 6
z = Dense(16, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)
z = LeakyReLU(alpha = 0.05)(z)
z = Dropout(0.5)(z)
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
model.compile(optimizer = 'adam', loss = 'categorical_crossentropy', metrics=['accuracy'])

# callbacks
mcp_save = keras.callbacks.callbacks.ModelCheckpoint('ann.h5', save_best_only=True, monitor='val_accuracy', verbose=1)
reduce_lr = keras.callbacks.callbacks.ReduceLROnPlateau(monitor='val_accuracy', factor=0.1, patience=10, verbose=1, mode='auto', min_delta=0.0001, cooldown=0, min_lr=0)
callbacks_list = [reduce_lr, mcp_save]

# Training
history = model.fit(X_train, y_train, batch_size = 32, epochs = 150, validation_data = (X_test, y_test), shuffle = False, callbacks = callbacks_list)

# Testing
model = load_model('ann.h5')
eval = model.evaluate(x = X_test, y = y_test)
print("Loss: " + str(eval[0]) + ", Accuracy: " + str(eval[1]))

y_pred = model.predict(X_test)
count = 0
for i in range(len(y_pred)):
    if np.where(y_pred[i] == max(y_pred[i]))[0][0] == np.where(y_test[i] == (max(y_test[i])))[0][0]:
        count = count + 1
print("Accuracy: ", count/len(y_pred))

# Metrics
print("Confusion Matrix")
matrix = confusion_matrix(y_test.argmax(axis=1), y_pred.argmax(axis=1))
print(matrix)

print("F1 Score")
print(f1_score(y_test.argmax(axis=1), y_pred.argmax(axis=1), average = 'weighted'))



# # Plot History
# # # summarize history for accuracy
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

''' 
[[152   0  58  33]
 [  0 212   2   0]
 [ 63   0 127  25]
 [ 19   0  23 204]]

loss: 0.7292 - accuracy: 0.6842 - val_loss: 0.6676 - val_accuracy: 0.7571
[0,1,2,3]: ['Disulphide' 'Sulphenylation' 'metal-binding' 'thioether']
'''