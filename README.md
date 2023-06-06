# DNA-DDA 
## An algorithm for predicting A/B compartments from the human reference sequence based on dynamical ergodicity DDA



X. Lainscsek, L. Taher. Predicting chromosomal compartments directly from the nucleotide sequence with DNA-DDA. Briefings in Bioinformatics. 2023. 
https://doi.org/10.1093/bib/bbad198

## Overview

DNA-DDA has been adapted from the nonlinear time series analysis technique dynamical-ergodicity delay differential analysis (DE DDA) ( Lainscsek et. al https://doi.org/10.1063/5.0063724 ), in order to predict chromosomal contacts from a reference sequence. This page includes an example for predicting contacts of chromosome 22 at 100kbp resolution in the GM12878 cell line.

The general steps of DNA-DDA are:  
 1. convert the reference sequence into a random walk DNA representation  
 2. compute ST and CT DDA features of all genomic bins of interest  
 3. combine ST and CT errors computed in step 2. to DE DDA to get DNA-DDA contact map  
 4. obtained DNA-DDA contact maps can be used as HiC-contact maps.

![DNA-DDA procedure](/Figures/DNA_DDA_RW_Procedure.svg)


## Requirements
* the DDA executable requires a Linux enviornment
* post and preprocessing scripts are written in MATLAB R2019b
* raw HiC data for GM12878 avilable in Gene Expression Omnibus under the accession **GSE63525** (Rao et al. (2014), Sanborn et al. (2015)).

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

    * **BINs** . . . contains bin number and coordinates of non-empty bins of chromosome 22
    * **ChrNr** . . . chromosome number 
    * **ChrSize** . . . length of reference sequence of chromosome 22
    * **Resolution** . . . 100 kb 
    * **HiCM** . . . HiC-Contact map chromosome 22
    * **HiCP** . . . HiC Pearson correlation matrix chromosome 22
    * **EXCLUDE** . . . omitted regions from matrices Includes centromeres from UCSD genome browser and low coverage bins

 5. generate **1D DNA walk** and save to file <FN_ASCII> 
	 
    * <FN_ASCII> has dimension **Resolution**$\times$**Nr_of_Bins** 
	    
    ![DNA_1DRW](/Figures/DNA_RW.svg)


 6. generate bash script to run DDA

 7. call bash script 
    ```
    ./<date>_runDDAbash.sh &
    ```
	  
 8. generate DNA-DDA contact matrix `DNA_DDA` from DDA outputs and save to **ERGODICITY.mat**

    * delay values of each cell line can be found in publication     


 9. `DNA-DDA` matrix post processing

    * map high to low values
    * fill in diagonal with neighboring values
    * reflect distance dependency of contact probability: each loci pair is multiplied by their genomic distance to the negative power of one 
    * take logarithm    
    * plot

    ![DNA_DDA](/Figures/ContactMaps.svg)


10. calling A/B compartments (eg with HiCExplorer https://github.com/deeptools/HiCExplorer)

     * convert `DNA-DDA` matrix to **homer** format (http://homer.ucsd.edu/homer/interactions/) 
       * `hicConvertFormat` to convert **homer** to **h5**
       * `hicNormalize --normalize norm_range` to normalize to 0 and 1
       * `hicPCA --method "dist_norm" --ignoreMaskedBins --pearsonMatrix <OUT_FN.h5>` to generate pearson correlation matrix  
       * `hicConvertFormat --outputFormat ginteractions` to generate tsv 
       * `hiCtsv_to_MATLABcsv` to load Pearson DNA-DDA matrix `DNA_DDA_P` into MATLAB
       
     * perform PCA with MATAB's `pca()` function to get principal component coefficients of `DNA_DDA_P` and `HiCP` matrices
     * the columns of resulting `PC_DDA` and `PC_HiC` contain coefficients for each PC in descending order of component variance
       * applying moving average filter to `PC_DDA`
       * call `Norm_PC.m` on `PC_DDA` and `PC_HiC` to normalize and determine which of the first three PCs defines A/B compartments
         * MATLAB ’s `filloutliers(PC,’nearest’)` function replaces outliers by nearest nonoutlier value 
         * Determines the PC with highest correlation to H3K4me1, a histone mark which is associated with open chromatin. The H3K4me1 density was derived from ChIPSeq data set under the GEO accession **GSM733772**

    ![DNA_DDA_P](/Figures/Pearson_Matrices.svg)


    ![PCs](/Figures/PCs.svg)















