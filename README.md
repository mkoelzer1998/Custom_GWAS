# Custom_GWAS

Current Overview:

For our project, we are building our own custom GWAS using linear regression that will take in a simulated phenotype data file and a genotype vcf file and output a data file containing p-values for each SNP in the genotype vcf file. We have a python executable that takes in two system arguments corresponding to each data file. The output will be a dataframe where each row corresponds to a SNP and the generated p-value associated with it. The dataframe can be used to generate qq plots and manhattan plots to visualize the generated values. This is the extent of what we have created so far for our project.

Remaining Work:

We would like to test GWAS using different models rather than just linear regression. We also still need to compare our implementation with that of plink's. We want to compare our results and runtime to then analyze the performance of our model. 


