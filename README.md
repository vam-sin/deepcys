[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Deepcys Structure
A complete Deep Learning solution to predicting the behavior of a given cysteine. The predictions are made using the features from the high resolution protein crystal structures.

# Requirements

```python3
pip3 install -r requirements.txt
sudo apt-get install dssp
sudo apt-get install -y ncbi-entrez-direct
sudo apt-get install -y libgfortran5
sudo apt-get install csh
```

# Download Data

Install the icedtea plugin.

```
sudo apt-get install icedtea-netx
```

In the downloads folder, run the following commands. Store the pdb files in the /pdb folder.

```java
javaws pdb.jnlp
```

# Feature Generation

In /src/model/features.

## Secondary Structure Folds

Folder SSF. 
Create a folder named feature.
Set the window size to 7 and run the command.

```python3
python3 ssf_gen.py
```

## Amino Acid Signatures in Interaction Shells

Folder AASIS. 
Create a folder named feature.

```python3
python3 aasis_gen.py
```

## Enzyme Class

Folder EC. 
Create a folder named feature.

```python3
python3 ec_gen.py
```

## Motifs

Folder Motifs. 
Create a folder named feature.
Run the following command for the window sizes [3, 5, 7, 9, 11, 13] 

```python3
python3 motif_gen.py
```

# Model Training

## Structure Based Prediction

Folder /src/model/nn.

The folder cd-hit had the results from the CD-HIT server which contain the representative sequences for 30% and 100% redundancy.

ann.py is the training file for the entire training dataset.
ann_redundancy.py is the training file for 30% non-redundant training dataset.
ann_identity.py is the training file for 100% non-redundant training dataset.

```python3
python3 ann.py
python3 ann_redundancy.py
python3 ann_identity.py
```

The results mentioned in the paper can be replicated by using the pre-trained models available [here](https://drive.google.com/drive/folders/1Hu_P80OheLdSDpRqPcECbdua1kTjv8Vx?usp=sharing).

# Batch Prediction

## Structure Based Prediction

Folder /src/batch_pred_struct/menv_server

Change the file location in the mentioned lines for all these files accordingly.

new-psfgen-script.tcl : 7, 12.  <br/>
newrenumber.sh : 10, 33, 68, 99.  <br/>
psf-pdb-gen-standalone.tcl : 3.  <br/>
script3 : 14 (twice), 16, 17, 20, 21, 24, 29.  <br/>
temp-psfgen-scripts.tcl : 8.  <br/>

Folder /src/batch_pred_struct.
In the file named "data.txt", mention the tuple of [PDB ID, Residue Number, Chain]. One in a line. Run the following command.

```python3
python3 batch_pred.py
```

In a file named "results.txt", the results from the batch prediction will be stored.

An example residue is already mentioned in the data.txt for reference. 
