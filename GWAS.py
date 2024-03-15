import numpy as np
import pandas as pd
import vcf
import pandas as pd
from tqdm import tqdm
import statsmodels.api as sm
import sys
import argparse

#Deprecated
class LinearRegressionModel:
    def init(self, alphas, xs, Ys):
        assert(len(alphas) == len(xs) == len(Ys))
        self.alphas = alphas
        self.xs = xs
        self.Ys = Ys
    def fit(self):
        sum_squared_x = np.sum(self.xs**2)
        sum_x_squared = np.sum(self.xs)**2
        sum_xy = np.sum(self.xs*self.Ys)
        sum_x_sum_y = np.sum(self.xs)*np.sum(self.Ys)
        n = len(self.xs)

        Beta = (sum_xy - (sum_x_sum_y / n)) / (sum_squared_x - (sum_x_squared / n))
        return Beta

#Add column of ones before values
def add_column_ones(x):
    return np.append(np.atleast_2d(np.ones(len(x))).T, np.atleast_2d(x).T, axis=1)

# Any model we use needs to be used in a function like this
# This one uses statsmodels
# Must return Beta, PValue
def runLinearModel(x, y):
    X = add_column_ones(x)
    model = sm.OLS(y, X)
    results = model.fit()
    return results.params[1], results.pvalues[1]

class GWAS:
    def init(self, linearFunc, gts, pts):
        self.gts = gts
        self.linearFunc = linearFunc
        self.pts = pts
    def get_p_values(self):
        #row is snp
        #column is sample
        pValues = []
        for snp in tqdm(self.gts):
            beta, pValue = self.linearFunc(self.pts, snp)
            pValues.append(pValue)
        return pValues

def main():
    # Argument Reading
    try:
        parser = argparse.ArgumentParser(description='''A Custom GWAS Implementation.''')
        parser.add_argument("vcf_file", help="A VCF file storing the genotype data for n samples",
                type=str)
        parser.add_argument("phen_file", help="A .phen file storing the phenotype data for n samples",
                type=str)
        args = parser.parse_args()


        phen_file = args.phen_file
        vcf_file = args.vcf_file
    except:
        e = sys.exc_info()[0]
        print(e)
        return -1

    # File Reading & Data Organization
    try:
        print("Reading Input:")
        df = pd.read_table(phen_file, sep='\s+', header=None)
        phen_values = df[2].to_numpy()
        phen_samples = df[0].to_numpy()

        vcf_reader = vcf.Reader(filename=vcf_file)
        vcf_samples = []

        for record in vcf_reader:
            for sample in record.samples:
                vcf_samples.append(sample.sample)
            break
        
        vcf_reader = vcf.Reader(filename=vcf_file)

        #####################################
        print("Building pandas dataframe...")
        vcf_df = pd.read_csv('data/gwas.vcf.gz', sep="\t", comment='#', header=None)
        CHR = vcf_df[0].to_numpy()
        BP = vcf_df[1].to_numpy()
        SNP = vcf_df[2].to_numpy()

        vcf_arr = vcf_df.to_numpy()

        startSamples = 0
        for i in range(len(vcf_df.columns)):
            if type(vcf_arr[0][i]) == str and len(vcf_arr[0][i]) == 3 and '|' in vcf_arr[0][i]:
                startSamples = i
                break

        Xj = []

        print("Processing dataframe")
        for i in tqdm(range(len(vcf_arr))):
            all_count = []
            for allele in vcf_arr[i][startSamples:len(vcf_arr[0])]:
                all_count.append(int(allele[0] != '0') + int(allele[2] != '0'))
            Xj.append(all_count)
        
        ######################################################################
        organized_phen_values = []
        for vSample in vcf_samples:
            for i in range(len(phen_samples)):
                if vSample == phen_samples[i]:
                    organized_phen_values.append(phen_values[i])
        print("Input Processing Finished")
    except:
        print("Unable to read phenotype and/or genotype data from files.")

    #GWAS
    try:
        print("Beginning GWAS")
        gwas = GWAS()
        gwas.init(runLinearModel, Xj, organized_phen_values)
        p_values = gwas.get_p_values()
        df = pd.DataFrame({'CHR':[int(c) for c in CHR],'BP':BP,'P':p_values, 'SNP':SNP})
        df.to_csv('output.csv', index=False)
        print("Gwas finished. Output sent to output.csv")
    except:
        print("Unable to perform GWAS. Is your data complete?")



main()

