# PO2PLS Bootstrapping Project
This file contains:
- An overview of the project
- Details of the data and privacy considerations
- Step-by-step instructions for replication of results
- Software details
- Contents of this repository

## Overview of the project
The goal of this project is to assess the performance of bootstrapping and permutation testing for the calculation of p-values for the Two-way Orthogonal Partial Least Squares method (PO2PLS). 

### Method
A simulation study was be done to assess the performance of parametric and non-parametric bootstrapping as well as permutation testing. Four scenarios were created using parameters of chromosome 17 from a genetic atherosclerosis dataset and a corresponding histological dataset with 7 outcomes. The genetic atherosclerosis dataset contained 6103 variables after peak removal (see 'Extra analysis'). In two scenarios, the joint components were removed from the dataset (without joint condition). In the other two the scenarios, the highest 200 joint loadings and all histological components were retained (with joint condition). The scenarios were split up further by altering the noise to create a 1:4 signal-to-noise ratio (low noise condition), while in the other two scenarios the original signal-to-noise ratio was retained (high noise condition). Simulation were run with 1000 bootstrapping iterations and 100 PO2PLS EM-steps.

### Main conclusions
Non-parametric bootstrapping was found to be entirely unsuitable due to PO2PLS being identifiable up to the sign and the expected value under the null hypothesis being zero.

The high noise scenario did not yield any significant results for significance levels of 0.01 or lower (unadjusted) and performed around chance level for both the parametric bootstrapping and permutation testing. In the low noise condition, parametric bootstrapping performed well based on the ROC plot, with a sensitivity of 0.93 and an FDR of 0.28 at alpha = 0.001 (unadjusted p-values). Permutation testing was found to be functional but produced very bad results, with a visually smaller AUC and a sensitivity of 1.00 and an FDR of 0.95 at alpha = 0.001.

### Extra analysis
Early analysis also revealed large peaks in the joint components of the original dataset. Near-identical peaks were found in the specific components and in a separately run SVD. As such, they were considered artifacts. Analysis of the peaks revealed the peaks could not be predicted based on the SVD or specific components, as component numbers did not always match (e.g. the peak in the first joint component appearing in the second SVD component. It was further found that removal of these peaks yielded new peaks on a subsequent PO2PLS run. Due to time constraints, these analyses were not continued.

Datasets created with the parameters after the peak was removed were found to be peakless. However, the peak in the parameters may have affected the noise level and height of the joint component.

## Data & Privacy
### Privacy
As the atherosclerosis dataset contains medical information, it is considered a sensitive dataset. Data handled by the owner of this repository was anonimized by a third party and pre-processed by the supervisor. The dataset included no variables that may be linked back to the participant. As the data is neither open source nor owned by the owner of this repository, it could not be included in the repository. Data was only stored locally on well-protected devices and on the server of the High-Performance Computing Facility, which is suitable for sensitive data storage.

Simulated datasets were based on parameters of the atherosclerosis dataset but cannot be linked to the original data or the participants. It could however not be made available directly due to its size. Steps taken to obtain the parameters for this dataset are described in the code.

### Atherosclerosis data 
Raw data was obtained from the AtheroExpress project. A dataset pre-processed by the supervisor was used for parameter generation. The owner of this repository is aware of the pre-processing steps, but did not handle the raw data directly. Pre-processing included the matching of SNPs to their location, removal of duplicates and filtering out low allele frequencies. Finally, SNPs near or on one gene were aggregated into genetic PCs that retained 80% of the variance. The parameter generation itself is detailed in the code of this repository. One plot of the joint components was included (PO2PLS_res.pdf) for illustration of the peak issue.

### Simulated datasets
The simulated data is included in the file Datasets.RData. Each dataset consist of 6103 variables based on the genetic PCs (after peak removal) and 7 histological outcomes for 1000 samples (participants). Additionally, PO2PLS runs and bootstrapping results of the simulated dataset are included which were used directly for the analysis.

## Step-by-step instructions for replication
Of note is that any section starting with #### can be run separately. To do so: 
1. Always start by running the Setup section at the top of the file.
2. Run any section of one's choice.
Due to seeding being done at the start of every section, running entries separately and running the code file as a whole it will yield the same results.

To replicate this study, the files have to be run in order:
1. The Simulate_dataset.R file details the process for creation of the scenarios. As the atherosclerosis data is not included in the repository, the first part cannot be run. To obtain the datasets for the study, run the setup and second part of the code (from the section where the PO2PLS data is loaded.
2. The Nonparaboot.R, Paraboot.R and Permutation.R files contain the commands for the simulation study. The functions themselves are contained in the Functions.R file. These can be run in any order and can be run per scenario. The variables in th setup are the ones used for the simulation, with only n.cores adjusted based on availability of cores and time constraints. Each step was time-consuming, so running this on multiple cores is recommended. To do so, the computer must have multiple cores and n.cores must be adjusted accordingly. The datasets are consequently stored in the Datasets folder.
3. Figure_creation creates a figure of the PO2PLS loadings of chromosome 17 of the atherosclerosis dataset (not possible due to the dataset not being included), a histogram of the p-values of the non-parametric bootstrapping under the high noise, with joint component scenario (unadjusted) and a ROC plot and table of results for the low noise and high noise scenarios for both parametric bootstrapping and permutation testing. These figures are stored in the Figures folder. Functions for the calculation of p-values and creation of ROC-plots can be found in the Functions.R file. Non-parametric bootstrapping results are given as two-sided p-values by default, but can be changed to one-sided p-values based on the absolute values of the measure of interest by using the function nonpara_pval_abs.

## Software details
All code was run on the linux-based HPC server of the UMC Utrecht. The number of cores used varied per run based on availability of the cores. R version 3.5.1 was used.

## Contents
### Code and Figures
This folder contains all R files used for the primary analysis, as described in the step-by-step instructions, along with a separate file which contains all relevant function made by the owner of the repository. It also contains the Figures and Datasets folder
#### Figures
This folder contains ROC plots of every scenario with joint component for both permutation testing and parametric bootstrapping, a histogram of the p-values of the non-parametric bootstrapping under the high noise, with joint component scenario (unadjusted) and a plot of the joint component of the atherosclerosis dataset (PO2PLS_res.pdf).
#### Datasets
This folder contains the PO2PLS results of chromosome 17 as well as the indexes of the genetic variables of which the joint loadings were retained in the scenarios with joint components

### Extra Files
This section contains 4 knitted RMD files of part of the initial analysis that was later deemed to be beyond the scope of this project.
- Data exploration contains an investigation of the effect of dataset size on the speed of PO2PLS. It was found that computation time increases linearly.
- The three Peak_test files contain an investigation of peak effects using different seeds for data simulation. The details of each simulation are not described, but is included in the code. The conclusion of this was that the peaks in the specific components largely resembled the peaks included in the simulation while SVDs could contain peaks from a strong joint component. The peaks in the joint component, however, were prone to jumping between components and that an alteration in the joint loadings could cause such a jump. Additionally, the component in which peaks were located differed greatly between simulations based on the same parameters. This makes it difficult for bootstrapping methods to function properly without combating the peaks.
