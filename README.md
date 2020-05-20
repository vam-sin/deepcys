[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# deepcys
A complete Deep Learning solution to predicting the behavior of a given cysteine. The predictions are made using the features from the high resolution protein crystal structures.

# Requirements

```python3
pip3 install -r requirements.txt
sudo apt-get install dssp
```

# Download Data

Install the icedtea plugin.

```
sudo apt-get install icedtea-netx
```

In the downloads folder, run the following commands. Store the fasta files in /fasta folder and the pdb files in the /pdb folder.

```java
javaws fasta_1.jnlp
javaws fasta_2.jnlp
javaws pdb_1.jnlp
javaws pdb_2.jnlp
```

# Feature Generation

In /src/model/features.

## NF1 

Folder NF1. Set the window size to 13 and run the command.

```python3
python3 NF1_gen.py
```

## NF2

Folder NF2. 

```python3
python3 NF2_gen.py
```

## NF3

Folder NF3. 

```python3
python3 NF3_gen.py
```

## NF4

Folder NF4. Run the following command for the window sizes [3, 5, 7, 9, 11, 13] 

```python3
python3 NF4_gen.py
```

## NF5

Folder NF5. Set the window size to 13 and run the command.

```python3
python3 NF5_gen.py
```

# Model Training

Folder /src/model/nn.

```python3
python3 ann.py
```

# Batch Prediction

## Structure Based Prediction

Folder /src/batch_pred_struct/menv_server

Change the file location in the mentioned lines for all these files accordingly.

new-psfgen-script.tcl : 7, 12.
newrenumber.sh : 10, 33, 68, 99.
psf-pdb-gen-standalone.tcl : 3.
script3 : 14 (twice), 16, 17, 20, 21, 24, 29.
temp-psfgen-scripts.tcl : 8.

Folder /src/batch_pred_struct.
In the file named "data.txt", mention the tuple of [PDB ID, Residue Number, Chain]. One in a line. Run the following command.

```python3
python3 batch_pred_struct.py
```

In a file named "results.txt", the results from the batch prediction will be stored.

## Sequence Based Prediction

Folder /src/batch_pred_seq.
In the file named "data.txt", mention the tuple of [PDB ID, Residue Number, Chain]. One in a line. Run the following command.

```python3
python3 batch_pred_seq.py
```

In a file named "results.txt", the results from the batch prediction will be stored.
