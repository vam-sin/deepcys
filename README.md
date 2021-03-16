[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Deepcys Structure <img width="30" height="30" src="downloads/logo.png">
A complete Deep Learning solution to predicting the behavior of a given cysteine. The predictions are made using the features from the high resolution protein crystal structures.

This model can also be accessed via a web server which can be found [here](http://deepcys.herokuapp.com).

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

# Paper

```
@article{https://doi.org/10.1002/prot.26056,
author = {Nallapareddy, Vamsi and Bogam, Shubham and Devarakonda, Himaja and Paliwal, Shubham and Bandyopadhyay, Debashree},
title = {DeepCys: Structure-based multiple cysteine function prediction method trained on deep neural network: Case study on domains of unknown functions belonging to COX2 domains},
journal = {Proteins: Structure, Function, and Bioinformatics},
volume = {n/a},
number = {n/a},
pages = {},
keywords = {deep neural network, multiple cysteine function prediction, protein structure and sequence feature},
doi = {https://doi.org/10.1002/prot.26056},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/prot.26056},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/prot.26056},
abstract = {Abstract Cysteine (Cys) is the most reactive amino acid participating in a wide range of biological functions. In-silico predictions complement the experiments to meet the need of functional characterization. Multiple Cys function prediction algorithm is scarce, in contrast to specific function prediction algorithms. Here we present a deep neural network-based multiple Cys function prediction, available on web-server (DeepCys) (https://deepcys.herokuapp.com/). DeepCys model was trained and tested on two independent datasets curated from protein crystal structures. This prediction method requires three inputs, namely, PDB identifier (ID), chain ID and residue ID for a given Cys and outputs the probabilities of four cysteine functions, namely, disulphide, metal-binding, thioether and sulphenylation and predicts the most probable Cys function. The algorithm exploits the local and global protein properties, like, sequence and secondary structure motifs, buried fractions, microenvironments and protein/enzyme class. DeepCys outperformed most of the multiple and specific Cys function algorithms. This method can predict maximum number of cysteine functions. Moreover, for the first time, explicitly predicts thioether function. This tool was used to elucidate the cysteine functions on domains of unknown functions belonging to cytochrome C oxidase subunit-II like transmembrane domains. Apart from the web-server, a standalone program is also available on GitHub (https://github.com/vam-sin/deepcys).}
}
```
