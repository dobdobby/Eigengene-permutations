# Eigengene-permutations
Use permutation testing to ID significant changes in pairwise correlations between Eigengenes (from WGCNA)

This script takes Eigengenes (representative expression values) from consensus WGCNA gene coexpression modules identified across three conditions (EAA: control, short-lived; SY and RAPA: experimental conditions, both long-lived) and across 5 different fruitfly tissues (brain, gut, fat body, ovary, thorax/muscle). These eigengenes give a high-level representation of the structure of the transcriptional network. 

The script identifies module-module correlations in each condition, and the changes observed in correlation coefficients in the comparisons of EAA:SY and EAA:RAPA. The data are permuted x10000 to generate null distributions for the changes. Observed changes are considered sig. when occurring above the 0.975 quantile or below the 0.025 quantile of the null distribution.
