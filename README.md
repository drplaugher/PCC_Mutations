# PCC_Mutations

Attached in this repository are the codes and data used in our publication *Uncovering potential interventions for pancreatic cancer patients via mathematical modeling*. Specifically, you will find:

TCGA pancreatic cancer data: 
- expressions -- gene expression data
- mutation_data_public -- aggregate statistics of expression data, along with mutation attractor information
- name_map -- a mapping of gene names with our their corresponding nodes in our model
- paad_mutation_status_2018july26 -- time-course data
- PCC_expressionsRed -- subset of the original data focusing on four main mutation combinations

Statistical analysis code (Python3)
- KaplanMeier_PCC -- survival curves
- PC_KW_comb -- Kruskal-Wallis test and plotting 

Model functions
- pcc_logical_functions -- functions used in the model
  - use these functions for stable motifs (https://github.com/jgtz/StableMotifs)
    - can be used to derive attractors and controls

Simulation code (MATLAB)
- comb_pcc_DP - cumulative code to run simulations using SDDS
- MatlabToolboxes - folder of all the needed toolboxes to run MATLAB code
