#!/bin/env python

# Meta-analyses script for burdens
    # Sept 2021 
    # Start script like this:
    # Takes in .assoc files generated by RVTests
        # python AMPNIH_UKB_ALL_meta_analysis_numvar.py \
        #             --test [SkatO, CMC] \
        #             --numvar [minimum number of variants per gene you want to keep] \
        #             --ukb_all_cases [filepath to UKB all cases .assoc file - suffix included] \
        #             --amp_nih [filepath to AMP NIH .assoc file - suffix included] \
        #             --variant_group [variant group as defined by ANNOVAR or SNPEFF - for information purposes only] \
        #             --maf [MAF defined during RVtests - for information purposes only] \
        #             --group [what groups are being meta-analyzed? String - for information purposes only] \
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

# Initialize parser and add arguments
parser = argparse.ArgumentParser()
parser.add_argument("--ukb_all_cases", help="Input file name (with suffix) for UKB cases .assoc")
parser.add_argument("--amp_nih", help="Input file name (with suffix) for AMPxNIH .assoc")

parser.add_argument("--variant_group", help="Input string with what variant group was looked at")
parser.add_argument("--maf", help="Input string with what MAF was looked at")
parser.add_argument("--group", help="Input string with what groups are being analyzed")

parser.add_argument("--test", choices=["SkatO", "CMC"], help="What test was run? Options are SkatO or CMC")
parser.add_argument("--numvar", metavar='N', type=int, help="Keep genes with variants equal to or greater than this number")
parser.add_argument("--output", "-o", help="Desired output name for files")
args = parser.parse_args()

# Read in your separated .assoc files
print("Reading in your datasets...")

test = args.test 

if test == "SkatO": 
    ## For SKATO 
    header_text = ['Gene','RANGE','N_INFORMATIVE','NumVar','NumPolyVar','Q','rho','Pvalue']
    ukb_all_cases_raw_df = pd.read_csv(args.ukb_all_cases, delim_whitespace=True, header=0, names=header_text)
    amp_nih_raw_df = pd.read_csv(args.amp_nih, delim_whitespace=True, header=0, names=header_text)

    # Only keep rows equal to or greater than the limit indicated by the user
    print(f"Keeping only genes with variants equal to or greater than {args.numvar}...")
    ukb_all_cases_raw_df = ukb_all_cases_raw_df[ukb_all_cases_raw_df['NumVar'] >= args.numvar]
    amp_nih_raw_df = amp_nih_raw_df[amp_nih_raw_df['NumVar'] >= args.numvar]

    # Save the NumVars to merge with later...
    ukb_numvar = ukb_all_cases_raw_df[['Gene', 'NumVar']].copy()
    ukb_numvar.columns = ['GENE', 'UKB_ALL_CASES_NUMVAR']
    amp_nih_numvar = amp_nih_raw_df[['Gene','NumVar']].copy()
    amp_nih_numvar.columns = ['GENE', 'AMP_NIH_NUMVAR']
    
    # Drop the unncessary columns, leaving just GENE and P-VALUE
    print(f"Dropping unneccesary columns...")
    ukb_all_cases_raw_df.drop(columns=['RANGE','N_INFORMATIVE','NumVar','NumPolyVar','Q','rho'], inplace=True)
    amp_nih_raw_df.drop(columns=['RANGE','N_INFORMATIVE','NumVar','NumPolyVar','Q','rho'], inplace=True)

elif test == "CMC":
    ## For CMC 
    header_text = ['Gene','RANGE','N_INFORMATIVE','NumVar','NumPolyVar','NonRefSite','Pvalue']
    ukb_all_cases_raw_df = pd.read_csv(args.ukb_all_cases, delim_whitespace=True, header=0, names=header_text)
    amp_nih_raw_df = pd.read_csv(args.amp_nih, delim_whitespace=True, header=0, names=header_text)

    # Only keep rows equal to or greater than the limit indicated by the user
    print(f"Keeping only genes with variants equal to or greater than {args.numvar}...")
    ukb_all_cases_raw_df = ukb_all_cases_raw_df[ukb_all_cases_raw_df['NumVar'] >= args.numvar]
    amp_nih_raw_df = amp_nih_raw_df[amp_nih_raw_df['NumVar'] >= args.numvar]

    # Save the NumVars to merge with later...
    ukb_numvar = ukb_all_cases_raw_df[['Gene', 'NumVar']].copy()
    ukb_numvar.columns = ['GENE', 'UKB_ALL_CASES_NUMVAR']
    amp_nih_numvar = amp_nih_raw_df[['Gene','NumVar']].copy()
    amp_nih_numvar.columns = ['GENE', 'AMP_NIH_NUMVAR']

    # Drop the unncessary columns, leaving just GENE and P-VALUE
    print(f"Dropping unneccesary columns...")
    ukb_all_cases_raw_df.drop(columns=['RANGE','N_INFORMATIVE','NumVar','NumPolyVar','NonRefSite'], inplace=True)
    amp_nih_raw_df.drop(columns=['RANGE','N_INFORMATIVE','NumVar','NumPolyVar','NonRefSite'], inplace=True)

else:
    print("This is not a valid test, options are SkatO or CMC. Please re-run")
    exit()

# Save out the P-value columns to refer to later
data_frames = [ukb_all_cases_raw_df, amp_nih_raw_df]
original_pval_df = reduce(lambda left,right: pd.merge(left,right,on=['Gene'], how='outer'), data_frames).fillna('NA')
original_pval_df.columns = ['GENE', 'UKB_ALL_CASES_PVAL', 'AMP_NIH_PVAL']

# Make a gene list from the merged dataframe
final_gene_list = original_pval_df.GENE.unique()

# Now make a total raw df, with all the P-values
total_raw_df = pd.concat(data_frames)

# Loop it using Fisher to combine p-values.
print("Using Fisher to combine p-values...")
results = []
for i in range(len(final_gene_list)):
  this_gene = final_gene_list[i]
  temp_data = total_raw_df[total_raw_df['Gene'] == this_gene]
  combined_gene = stats.combine_pvalues(temp_data['Pvalue'], method='fisher', weights=None)
  test_stat = combined_gene[0]
  p_val = combined_gene[1]
  #print(this_gene, test_stat, p_val)
  results.append([this_gene, test_stat, p_val])
output = pd.DataFrame(results, columns=['GENE', 'TEST_STAT', 'META_PVAL'])

# Sort the output so that the most promising is at the top
print("Sorting output...")
sorted_output = output.sort_values(by=['META_PVAL'], ascending=True)
sorted_output.reset_index(drop=True, inplace=True)
sorted_output.head()

# Now merge with the original P values
print("Merging with original P-values...")
merge_pvals = sorted_output.merge(original_pval_df, how='left', on=['GENE'])

# Merge with the NumVars 
numvar_data_frames = ukb_numvar.merge(amp_nih_numvar, how='outer', on=['GENE']).fillna('NA')
final_dataframe = merge_pvals.merge(numvar_data_frames, how='outer', on=['GENE']).fillna('NA')

# Add columns to help with merging full results later that indicate variant group, MAF, and meta-analysis it's from 
print("Adding a few extra columns of information...")
final_dataframe['TEST'] = test
final_dataframe['VARIANT_GROUP'] = args.variant_group
final_dataframe['MAF'] = args.maf 
final_dataframe['MIN_NUMVAR'] = args.numvar 
final_dataframe['META_ANALYSIS_GROUP'] = args.group

# Now save out the file
print("Saving out file...")
output_file = args.output + ".combined_Ps." + args.test + ".tab"
final_dataframe.to_csv(output_file, sep="\t", index=False)

print(f"File has been saved as {output_file}")