# scfccoupling-high



# System requirements
This script was tested using MATLAB(R2020a) and RStudio(version 4.2.2). It requires no additional software packages.

# Installation and usage
Download and unzip the directory, open MATLAB and RStudio, and then navigate to wherever you downloaded the files. You can run the script from within that directory.


# Code

*. 'loaddata.m'
Script to read sc,fc,coordinate data and generates the structural eigenmodes and the functional gradient. Runtime is<1s

*. 'decoupling.m' 
Example script to implement the regional R values via multilinear regression models and to generate the main results . Runtime is <30s

*. 'LASSO_sparse.R' 
Example script to implement the LASSO regression. Runtime is <20s for one region in one fold
