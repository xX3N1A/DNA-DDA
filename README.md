# DNA-DDA 
## An algorithm for predicting A/B compartments from the human reference sequence based on dynamical ergodicity DDA


## Overview

DNA-DDA has been adapted from the nonlinear time series analysis techniq dynamical-ergodicity delay differential analysis (DE DDA) ( https://doi.org/11.1063/5.0063724 ), in order to predict chromosomal contacts from a reference sequence. This page includes an example for the first 50 non-empty bins of chromosome 1 at 100kbp resolution (chr1:500001:5600001). 

The general steps of DNA-DDA are:  
 1. convert the reference sequence into a random walk DNA representation  
 2. compute ST and CT DDA features of all genomic bins of interest  
 3. combine ST and CT errors computed in step 2. to DE DDA to get DNA-DDA contact map  
 4. obtained DNA-DDA contact maps can be used as HiC-contact maps.

![DNA-DDA procedure](/Figures/DNA_DDA_precedure.svg)


## Requirements
* the DDA executable requires a Linux enviornment
* post and preprocessing scripts are written in MATLAB R2019b

## Usage
 1. clone repository
    ```
    git clone https://github.com/xX3N1A/DNA-DDA
    ```
 2. download genome assembly (GRCh38.p13)

 3. bin sequence with perl script
    ```
    bin_sequence.pl <ChrNr> <Chr_Size> <Path_to_RefSeq> <Output_Path>
    ```
 4. load **vars.mat**

    * **BINs** . . . contains bin number and coordinates of first 50 non-zero bins of chromosome 1
    * **BINs\_ALL** . . . contains bin number and coordinates of all bins of chromosome 1 
    * **ChrNr** . . . chromosome number 
    * **ChrSize** . . . length of reference sequence of chromosome 1
    * **First_NonZero_Bin** . . . first nonempty bin (zero padding)
    * **Resolution** . . . 100kbp 
    * **HiCM** . . . HiC-Contact map chromosome 1

 5. generate **1D DNA walk** and save to file <FN_ASCII> 
	 
    * <FN_ASCII> has dimension **Resolution**$\times$**Nr_of_Bins** 
	    
    ![DNA_1DRW](/Figures/DNA_RW.svg)

 6. generate bash script to run DDA

 7. call bash script 
    ```
    ./<date>_runDDAbash.sh &
    ```
	  
 8. generate DNA-DDA contact matrix from DDA outputs and save to **ERGODICITY.mat**

    * delay values of each cell line can be found in publication     

 9. matrix post processing

    * map high to low values
    * fill in diagonal with neighboring values
    * take logarithm    
    * plot

 10. calling A/B compartments (eg with HiCExplorer https://github.com/deeptools/HiCExplorer)

     * convert matrix to **homer** format (http://homer.ucsd.edu/homer/interactions/) 
     * `hicConvertFormat` to convert **homer** to **h5**
     * `hicNormalize --normalize norm_range` to normalize to 0 and 1
     * `hicPCA --method "dist_norm" --ignoreMaskedBins --pearsonMatrix <OUT_FN.h5>` to generate pearson correlation matrix  
     * `hicConvertFormat --outputFormat ginteractions` to generate tsv
     * in MATAB `pca(DNA_DDA_PEARSON)`

     ![DNA-DDA_DNA](/Figures/ContactMaps.svg)

