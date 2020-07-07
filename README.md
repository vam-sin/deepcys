[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# deepcys
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

## NF1 

Folder NF1. 
Create a folder named feature.
Set the window size to 9 and run the command.

```python3
python3 NF1_gen.py
```

## NF2

Folder NF2. 
Create a folder named feature.

```python3
python3 NF2_gen.py
```

## NF3

Folder NF3. 
Create a folder named feature.

```python3
python3 NF3_gen.py
```

## NF4

Folder NF4. 
Create a folder named feature.
Run the following command for the window sizes [3, 5, 7, 9, 11, 13] 

```python3
python3 NF4_gen.py
```

# Model Training

## Structure Based Prediction

Folder /src/model/nn.

```python3
python3 ann.py
```

## Sequence Based Prediction

Folder /src/model/nn.

```python3
python3 gru.py
```

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

## Sequence Based Prediction

Folder /src/batch_pred_seq.
In the file named "data.txt", mention the tuple of [PDB ID, Residue Number, Chain]. One in a line. Run the following command.

```python3
python3 batch_pred_seq.py
```

In a file named "results.txt", the results from the batch prediction will be stored.
