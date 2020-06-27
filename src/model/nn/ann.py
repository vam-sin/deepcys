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

# Tasks
def recall_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall

def precision_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision

def f1_m(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))

# dataset import and preprocessing
ds = pd.read_csv("../data/correct_data/dataset.csv")
X1 = ds.iloc[:,4:6] # ('pKa', 'BF')
# print(X1)
# X1 = normalize(X1, norm = 'l2')
X1 = pd.DataFrame(X1)
X2 = ds.iloc[:,6:7] # rHpy
X = pd.concat([X1, X2], axis=1)
# print(X)
# X = X.fillna(X.mean())
y = ds.iloc[:,7]
y_w = y
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

# Make Test Data
metal_X = X.iloc[18759:18859,:]
X = X.drop(X.index[18759:18859])
y_ = pd.DataFrame(y)
metal_y = y_.iloc[18759:18859]
y_ = y_.drop(y_.index[18759:18859])
y_w = y_w.drop(y_w.index[18759:18859])
# print(metal_y)

sulph_X = X.iloc[18961:19061,:]
sulph_y = y_.iloc[18961:19061]
X = X.drop(X.index[18961:19061])
y_ = y_.drop(y_.index[18961:19061])
y_w = y_w.drop(y_w.index[18961:19061])
# print(sulph_y)

dis_X = X.iloc[55170:55270,:]
dis_y = y_.iloc[55170:55270]
X = X.drop(X.index[55170:55270])
y_ = y_.drop(y_.index[55170:55270])
y_w = y_w.drop(y_w.index[55170:55270])
# print(dis_y)

thio_X = X.iloc[107260:107360,:]
thio_y = y_.iloc[107260:107360]
X = X.drop(X.index[107260:107360])
y_ = y_.drop(y_.index[107260:107360])
y_w = y_w.drop(y_w.index[107260:107360])
# print(thio_y)
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

# print(X_test, y_test)
# print(list(encoder.inverse_transform(encoded_Y)))
# print(X, y)
X_train, y_train = shuffle(X_train, y_train, random_state = 42)

# X_train, X_test, y_train, y_test = train_test_split(X.values, y, test_size=0.1, random_state = 42)

# print(y_test)
# Conv1D Layers
y_ = pd.DataFrame(y_)
X_train = np.expand_dims(X_train,axis=2)
X_test = np.expand_dims(X_test,axis=2)

# NN (Skip Connections) Model
input_ = Input(shape = (len(X.columns),1,))
# x = GaussianNoise(0.5)(input_)
x = Conv1D(128, (3), padding = 'same', kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(input_)
x = LeakyReLU(alpha = 0.05)(x)
x = Conv1D(64, (3), padding = 'same', kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)
# x = GaussianNoise(0.5)(x)
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
z = Dropout(0.7)(z)
z = BatchNormalization()(z)

# Block 3
x = Dense(128, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(z)
x = LeakyReLU(alpha = 0.05)(x)
x = Dropout(0.7)(x)
x = BatchNormalization()(x) 

# Block 4
# Skip connection #2
z = Dense(128, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(skc_f)
sk1 = Add()([x, z])
z = Dense(64, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(sk1)
z = LeakyReLU(alpha = 0.05)(z)
z = Dropout(0.7)(z)
z = BatchNormalization()(z)

# Block 5
x = Dense(32, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(z)
x = LeakyReLU(alpha = 0.05)(x)
x = Dropout(0.7)(x)
x = BatchNormalization()(x) 

# Block 6
z = Dense(16, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)
z = LeakyReLU(alpha = 0.05)(z)
z = Dropout(0.7)(z)
z = BatchNormalization()(z)

# Block 7
# Skip Connection #3
x = Dense(16, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(skc_f)
sk2 = Add()([z, x])
x = Dense(8, kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(sk2)
x = LeakyReLU(alpha = 0.05)(x)
x = BatchNormalization()(x) 
x = Dropout(0.7)(x)
x = Dense(4, activation = 'softmax', kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-4), bias_regularizer=regularizers.l2(1e-4), activity_regularizer=regularizers.l2(1e-5))(x)

model = Model(input_, x)
model.compile(optimizer = 'adam', loss = 'categorical_crossentropy', metrics=['accuracy'])

# callbacks
mcp_save = keras.callbacks.callbacks.ModelCheckpoint('ann.h5', save_best_only=True, monitor='val_accuracy', verbose=1)
reduce_lr = keras.callbacks.callbacks.ReduceLROnPlateau(monitor='val_accuracy', factor=0.1, patience=2, verbose=1, mode='auto', min_delta=0.0001, cooldown=0, min_lr=0)
callbacks_list = [reduce_lr, mcp_save]

# Training
# weights = [ 0.31614373  1.43080227 46.60362694  8.58253817]
weights = {0: 0.31614373, 1: 1.43080227, 2: 46.60362694, 3: 8.58253817}
history = model.fit(X_train, y_train, batch_size = 512, epochs = 20, validation_data = (X_test, y_test), shuffle = False, callbacks = callbacks_list, class_weight = weights)

# Testing
model = load_model('ann.h5')
eval_ = model.evaluate(x = X_test, y = y_test)
print("Loss: " + str(eval_[0]) + ", Accuracy: " + str(eval_[1]))
print(eval_)

y_pred = model.predict(X_test)
# count = 0
# for i in range(len(y_pred)):
#     if np.where(y_pred[i] == max(y_pred[i]))[0][0] == np.where(y_test[i] == (max(y_test[i])))[0][0]:
#         count = count + 1
# print("Accuracy: ", count/len(y_pred))

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
######## Tasks ######$$$$
1. Split into specific train and test. (Done)
Take 100 of each. 
2. Convert to weighted crossentropy (Done)
If you are talking about the regular case, where your network produces only one output, 
then your assumption is correct. In order to force your algorithm to treat every instance of class 1 as 50 instances of class 0 you have to.
"treat every instance of class 1 as 50 instances of class 0" means that in your loss function you assign higher value to these instances. 
Hence, the loss becomes a weighted average, where the weight of each sample is specified by class_weight and its corresponding class. 
3. Choose a different performance metric to maximize - Done. Fixed Class Imbalance issue.
Precision, Recall, F1 Score
Undersampling Majority
Oversampling Majority
SMOTE
####### Results #########

loss: 1.6696 - accuracy: 0.8350 - val_loss: 1.2270 - val_accuracy: 0.8625
Confusion Matrix
[[100   0   0   0]
 [  0  86   1  13]
 [  0   1  99   0]
 [  0  30  10  60]]
F1 Score
0.857281372366213

Grid Search:
Batch Size: 256
Epochs: 20

####### Parameters #######

[0,1,2,3]: ['disulphide' 'metal-binding' 'sulphenylation' 'thioether']: (weights) [0.31694402, 1.42852999, 39.88733432, 8.34879778]: (frequency)

'''