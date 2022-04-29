#!/bin/env python

# Meta-analyses script for burdens - calculate lambdas
    # Jan 2022
    # Start script like this:
        # python calclambdas_meta.py \
        #             --test [SkatO, CMCWald] \
        #             --numvar [minimum number of variants per gene you want to keep] \
        #             --group [what groups are being meta-analyzed? String - for information purposes only] \
        #             --input [input .tab file] \
        #             --ncases [number of cases (int)] \
        #             --ncontrols [number of controls (int)] \
        #             --variant_group [variant group as defined by ANNOVAR or SNPEFF - for information purposes only] \
        #             --maf [MAF defined during RVtests - for information purposes only] \
        #             -o [output directory and file name desired]

# Import the necessary packages
import numpy as np
import pandas as pd
import math
import sys
import joblib
import subprocess
import statsmodels.api as sm
from scipy import stats
from functools import reduce
import argparse
import os
from scipy.stats import ncx2

# Initialize parser and add arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Input file name (with suffix) for meta-analysis .tab file")

parser.add_argument("--group", help="Input string with what groups are being analyzed")
parser.add_argument("--variant_group", help="Input string with what variant group was looked at")
parser.add_argument("--maf", help="Input string with what MAF was looked at")
parser.add_argument("--test", choices=["SkatO", "CMCWald"], help="What test was run?")
parser.add_argument("--numvar", metavar='N', type=int, help="Keep genes with variants equal to or greater than this number")

parser.add_argument("--ncases", type=int, help="Input number of cases")
parser.add_argument("--ncontrols", type=int, help="Input number of controls")
parser.add_argument("--output", "-o", help="Desired output name for files")
args = parser.parse_args()

# Read in your file
print("Reading in your dataset...")
input_df = pd.read_csv(args.input, sep="\t")
df = input_df[['META_PVAL']].copy()
df.dropna(inplace=True)

# Make variables
test = args.test 

# Make lambdas function
def calculate_inflation(pval_array, normalize=False, ncases=None, ncontrols=None):
    """Calculate lambda/genomic inflation.
    Normalize option: when set to True normalizes it to 1000 cases and 1000 controls. Also called lambda1000.
    ncases + ncontrols are necessary if calculating lambda1000"""
    num = ncx2.ppf(1-pval_array, 1, nc=0)
    denom = ncx2.ppf(0.5, 1, nc = 0)
    inflation = np.median(num)/denom
    if normalize:
        inflation1000 = 1 + (inflation -1) * (1/ncases+ 1/ncontrols)/(1/1000 + 1/1000)
        inflation = inflation1000
    return(inflation)

# Calculate values 
lambda_val = calculate_inflation(df['META_PVAL'])
lambda1000_val = calculate_inflation(df['META_PVAL'], normalize=True, ncases=args.ncases, ncontrols=args.ncontrols)

#print(f'Lamda: {lambda_val:.3f}')
#print(f'Lambda 1000: {lambda1000_val:.3f}')

results = []
results.append([args.group, args.variant_group, args.maf, args.test, args.numvar, lambda_val, lambda1000_val])

output = pd.DataFrame(results, columns=['META_COHORT', 'VARIANT_GROUP', 'MAF', 'TEST', 'MIN_NUMVAR', 'LAMBDA', 'LAMBDA1000'])

# Now save out the file
print("Saving out file...")
output_file = args.output + ".lambdas." + args.test + ".tab"
output.to_csv(output_file, sep="\t", index=False)

print(f"File has been saved as {output_file}")
