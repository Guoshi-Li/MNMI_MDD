Matlab code and data for the Default mode-Salience (DMN-SAL) network. The detailed method with results can be found in the following paper:

Li et al. (2020) Multiscale Neural Modeling of Resting-state fMRI Reveals Executive-Limbic Malfunction as a Core Mechanism in Major Depressive Disorder. medRxiv, doi: https://doi.org/10.1101/2020.04.29.20084855   

The "FC" folder contains the functional connectivity (FC) matrices for 98 normal control (NC) subjects and 96 individuals with major depressive disorder (MDD).

The "Visualization" folder contains the code and data to visualize the estimated (optimized) effective connectivity (EC)

The "SC" mat file stores the average Structure Connectivity matrix from 14 HCP subjects for the DMN-SAL and EXE-LIM networks

DSN_NC001.m: This m-script estimate effective connectivity (EC) for the NC001 subject

DSN_MDD001.m: This m-script estimate EC for the MDD001 subject

Model_NEURAL.m: This m-script function calculates neural activity

Model_HEMO.m: This m-script function converts neural activity into BOLD responses 
