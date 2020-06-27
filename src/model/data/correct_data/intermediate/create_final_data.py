import pandas as pd 
import pickle

m_s_d = pd.read_csv('get_three_features.csv')
m_s_d = m_s_d.iloc[:,1:]
# print(m_s_d)
m_s_d_X = m_s_d.iloc[:,0:3]
m_s_d_y = m_s_d.iloc[:,3:4]

infile = open('pka.pickle','rb')
pka = pickle.load(infile)
infile.close()
pka = pd.DataFrame(pka)

infile = open('bf.pickle','rb')
bf = pickle.load(infile)
infile.close()
bf = pd.DataFrame(bf)

infile = open('rhpy.pickle','rb')
rhpy = pickle.load(infile)
infile.close()
rhpy = pd.DataFrame(rhpy)

# print(pka, bf, rhpy)

m_s_d_X = pd.concat([m_s_d_X, pka, bf, rhpy], axis=1)
# print(m_s_d_X)
m_s_d = pd.concat([m_s_d_X, m_s_d_y], axis=1)
m_s_d.columns = ['pdb', 'res', 'chain', 'pka', 'bf', 'rhpy', 'mod']
print(m_s_d)

thio = pd.read_excel('Thioether_1.5-2.0.xlsx')
thio = thio.iloc[:,1:]
thio.columns = ['pdb', 'res', 'chain', 'pka', 'bf', 'rhpy', 'mod']
print(thio)

ds = pd.concat([m_s_d, thio], axis=0)
ds = ds.reset_index()
ds = ds.iloc[:,1:]
print(ds)

ds.to_csv('dataset.csv')