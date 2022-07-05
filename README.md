# DNA-DDA 
## An algorithm for predicting A/B compartments from the human reference sequence based on dynamical ergodicity DDA


## Overview

DNA-DDA has been adapted from the nonlinear time series analysis techniq dynamical-ergodicity delay differential analysis (DE DDA) ( https://doi.org/11.1063/5.0063724 ), in order to predict chromosomal contacts from a reference sequence. This page includes an example for the first 50 non-empty bins of chromosome 1 at 100kbp resolution (chr1:500001:5600001). 

The general steps of DNA-DDA are:  
 1. convert the reference sequence into a random walk DNA representation  
 2. compute ST and CT DDA features of all genomic bins of interest  
 3. combine ST and CT errors computed in step 2. to DE DDA to get DNA-DDA contact map  
 4. obtained DNA-DDA contact maps can be used as HiC-contact maps. To call A/B compartments (eg with HiCExplorer https://github.com/)  

![DNA-DDA procedure](/Figures/DNA_DDA_precedure.svg)


## Requirements
* the DDA executable requires a Linux enviornment
* post and preprocessing scripts are written in MATLAB R2019b

## Usage
 1. clone repository
```
git clone https://github.com/xX3N1A/DNA-DDA
```


