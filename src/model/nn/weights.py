import pandas as pd 
from sklearn.utils import class_weight
import numpy as np
from sklearn.preprocessing import StandardScaler, LabelEncoder, normalize

ds = pd.read_csv("../data/correct_data/dataset.csv")

y = ds.iloc[:,7]
y = list(y.values)
encoder = LabelEncoder()
encoder.fit(y)
y = encoder.transform(y)
class_weights = class_weight.compute_class_weight('balanced',
                                                 np.unique(y),
                                                 y)
print(class_weights, np.unique(y))

'''
['disulphide' 'metal-binding' 'sulphenylation' 'thioether'] : [ 0.31694402  1.42852999 39.88733432  8.34879778]
[0,1,2,3]: ['Disulphide' 'Sulphenylation' 'metal-binding' 'thioether']
weights = [0.31694402, 39.88733432, 1.42852999, 8.34879778]
'''