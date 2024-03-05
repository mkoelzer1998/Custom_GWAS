# Custom_GWAS

Current Overview:

For our project, we are building our own custom GWAS using linear regression that will take in a simulated phenotype data file and a genotype vcf file and output a data file containing p-values for each SNP in the genotype vcf file. We have a python executable, GWAS.py, that takes in two system arguments corresponding to each data file. The two arguments are vcf_file and phen_file. The vcf_file should store the genotype data for n samples and the phen_file should store the phenotype data for n samples. In our case, we use the simulated phenotype data and the genotype data given to us in problem set 3. This then allows us to compare our custom GWAS with what we did in problem set 3. The executable will read both files in and organize the data so that we can run a linear regression model on the data to obtain p values corresponding to each SNP. We store each value along with its corresponding chromosome, rID, and position into a dataframe. The dataframe is then written to a output.csv file where we can use that data to produce visualizations such as qq plots and manhattan plots as demonstrated in our visuals.ipynb. The following line shows how to run the python executable on the command line given the files gwas.phen and pwas.vcf.gz:

python GWAS.py "gwas.vcf.gz" "gwas.phen"

For help and a description about the python executable, you can run:

python GWAS.py -h

Remaining Work:

We would like to test GWAS using different models rather than just linear regression. We also still need to do a deeper comparison of our implementation with that of plink's. We want to compare our results and runtime to then analyze the performance of our model. Something that we noticed was that the runtime of the excecutable took a long time. There are two main portions that take up most of the time: reading the vcf file and running GWAS on the data. We want to be able to decrease the runtime by finding out ways to optimize these two portions. This could be done by possibly parallelizing any parts that we can.


