# Custom GWAS

## Project Members
Madison Koelzer

Nathan Bui

Group Number: 19

## Current Overview:
For our project, we are building our own custom GWAS using linear regression that will take in a genotype vcf file and a phenotype data file and output a data file containing p-values for each SNP in the genotype vcf file. 

### Methodology
We use the simulated phenotype data and the genotype data given to us in problem set 3. This will then allow us to compare our custom GWAS with what we did in problem set 3. The executable reads both files in and organizes the data so that we can run linear regression models on the data to obtain p values corresponding to each SNP. We store each p-value along with its corresponding chromosome, rID, and position into a dataframe. The dataframe is then written to a ```output.csv``` file where we can use that data to produce visualizations such as qq plots and manhattan plots, as demonstrated in ```visuals.ipynb```.

## Setup
We include the following setup guide to get started with our tool.
### Setting up a Virtual Environmennt
Run the following script to build a virtual environment.
```
bash setup_virtual.bash 
```
To start up this virtual environment, run the following from the same directory as the project:
```
source custom_env/virtualenv/Scripts/activate
```
Now that you have a fresh virtual environment setup, the dependencies can be started with the following script:
```
bash installs.bash
```
### Getting Additional Files
We perform some analysis on our tool, and the necessary files for this analysis can be found at the following google drive:
```
https://drive.google.com/drive/u/1/folders/1GuMk7VcI40INi8544QBJkMPRPchA4eHy
```
Just copy the files from this drive into the same repository as this tool

## Using the tool
Our tool is conveniently packaged into a python script. You can run it using the following structure:
```
python GWAS.py -vcf_file "file_name.vcf.gz" -phen_file "file_name.phen"
```
For help and a description about this python executable, you can run:
```
python GWAS.py -h
```

## Analysis
Analysis was performed on our tool to compare its efficiency to that of Plink's implementation. Note you must have the files from the above [section](#getting-additional-files). To view our tools efficiency, we first run our python script on the test dataset:
```
python GWAS.py "data/gwas.vcf.gz" "data/gwas.phen"
```
This command may take considerable time to run. On our machine, it takes roughly 4 minutes - 1 to read in all the data, and 3 to perform the model fitting. 

Next, we will run our python notebook ```visuals.ipynb```. Note that the final cell will not run without an additional file. For users who are members of the CSE 210 Course, the file needed is ```ps3_gwas.assoc.linear``` from Problem Set 3. If everything has worked correctly, the graphs generated from our tool should be identical to those generated using Plink in Problem Set 3.
