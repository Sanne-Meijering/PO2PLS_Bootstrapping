# PO2PLS Bootstrapping Project

## Goal
The goal of this project is to assess the performance of bootstrapping for the calculation of p-values for the Two-way Orthogonal Partial Least Squares method (PO2PLS). To this end, a simulation study will be done to assess the performance of parametric and non-parametric bootstrapping as well as permutation testing. The best performing method will then be used to analyze twp datasets containing histolgical traits of atherosclerosis plaques and a matching dataset containing single nucleotide polymorphisms (SNPs, small genetic mutations).

## Data

Data is not included in this repository due to sensitivity of the data. Simulation data will be included at a later date. 

## Current contents
### Extra files
Folder containing presentations and reports on the current project.

### Code
* Data exploration: an investigation of the effect of dataset size on the speed of PO2PLS. It was found that computation time increases linearly.
* GWAS_run: a quick correlation-based Genome Wide Association Study for comparison.
* PO2PLS_run, PO2PLS_run_ch10: two files containing the actual PO2PLS runs on the original data. Also include an investigation into strange peaks that appeared in the PO2PLS results.
* R project file
