# Libraries
import numpy as np 
import pandas as pd 

# Load all the data

# metal-binding
metal2 = pd.read_csv('metal-2+-1.5-2.0.csv')
metal2 = metal2.iloc[:,2:5]

y = []
for i in range(len(metal2)):
	y.append('metal-binding')

y = pd.DataFrame(y)
metal2 = pd.concat([metal2, y], axis=1)

metal2.columns = ['pdb', 'res', 'chain', 'mod']

# print(metal2)

metal3 = pd.read_csv('metal-3+-1.5-2.0.csv')
metal3 = metal3.iloc[:,2:5]

y = []
for i in range(len(metal3)):
	y.append('metal-binding')

y = pd.DataFrame(y)
metal3 = pd.concat([metal3, y], axis=1)

metal3.columns = ['pdb', 'res', 'chain', 'mod']

# print(metal3)

metal = pd.concat([metal2, metal3], axis=0)

metal = metal.reset_index()
metal = metal.iloc[:,1:]
# print(metal)

# thioether (thio has all three already calculated)
thio = pd.read_excel('Thioether_1.5-2.0.xlsx')
# print(thio)

# sulph
sulph = pd.read_csv('suphenylation1.5-2.0.csv')
sulph = sulph.iloc[:,1:4]

y = []
for i in range(len(sulph)):
	y.append('sulphenylation')

y = pd.DataFrame(y)
sulph = pd.concat([sulph, y], axis=1)

sulph.columns = ['pdb', 'res', 'chain', 'mod']
# print(sulph)

# disulphide
dis = pd.read_csv('disulphide1.5-2.0.csv')
dis = dis.iloc[:,1:4]
dis = dis.drop_duplicates()
dis = dis.reset_index()
dis = dis.iloc[:,1:]

y = []
for i in range(len(dis)):
	y.append('disulphide')

y = pd.DataFrame(y)
dis = pd.concat([dis, y], axis=1)

dis.columns = ['pdb', 'res', 'chain', 'mod']

dis2 = pd.read_csv('ds15to222.csv')
dis2 = dis2.iloc[:,1:4]
dis2 = dis2.drop_duplicates()
dis2 = dis2.reset_index()
dis2 = dis2.iloc[:,1:]

y = []
for i in range(len(dis2)):
	y.append('disulphide')

y = pd.DataFrame(y)
dis2 = pd.concat([dis2, y], axis=1)

dis2.columns = ['pdb', 'res', 'chain', 'mod']
# print(dis2)
# print(dis)

non_thio_need_features = pd.concat([metal, sulph, dis, dis2], axis=0)
# print(non_thio_need_features)

non_thio_need_features.to_csv('get_three_features.csv')
