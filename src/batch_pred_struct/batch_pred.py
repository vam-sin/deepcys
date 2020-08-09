import pandas as pd 
import numpy as np 
from get_151_features import get_features
from tensorflow.keras.models import load_model
from sklearn.preprocessing import StandardScaler, LabelEncoder, normalize
from sklearn.metrics import confusion_matrix, accuracy_score, f1_score
from keras.utils import to_categorical, np_utils
from tqdm import tqdm
import pickle

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

# # read from data.txt
pdb = []
res = []
chain = []
f = open("data.txt", "r")
for line in f:
    line_split = line.split(',')
    new_line_split = []
    for el in line_split:
        if el != ' ' and el != '':
            new_line_split.append(el)
    line_split = new_line_split
    pdb.append(line_split[0])
    res.append(int(line_split[1].replace(' ', '')))
    chain.append(line_split[2].replace("\n", '').replace(' ', ''))

# feature gen
features_X = []

for i in tqdm(range(len(pdb))):
	print("####################")
	print(i+1, " Out of ", len(pdb))
	print("####################")
	features_X.append(get_features(pdb[i].replace(' ',''), res[i], chain[i]))

features_X = np.asarray(features_X)

features_X = np.expand_dims(features_X, axis=2)

classes_dict = {"0": "disulphide", "1": "metal-binding", "2": "sulphenylation", "3": "thioether"}

# load model and predictions
model = load_model('ann_90.h5')

y_pred = model.predict(features_X)

y_pred = np.asarray(y_pred)

for i in y_pred:
    pred_num = np.argmax(i)
    print(classes_dict[str(pred_num)])

g = open("results.txt", "w")
g.write("##### DeepCys Batch Prediction Results #####\n\n")
for i in range(len(y_pred)):
    text = classes_dict[str(np.argmax(y_pred[i]))]
    print(str(text))
    g.write(str(text))
    g.write("\n")

g.close()

