# Import the necessary packages
import numpy as np
import pandas as pd
import math
import sys
import joblib
import subprocess
import statsmodels.api as sm
from scipy import stats
from functools import partial,reduce
import argparse
import os
import warnings

warnings.filterwarnings('ignore')
print(sys.version) ## Need python 3.10 to run!

# Initialize parser and add arguments
parser = argparse.ArgumentParser()
parser.add_argument("--amp_nih", default="amp_nih", help="Input file name (with suffix) for AMPxNIH .assoc")
parser.add_argument("--ukb_all", default="ukb_all", help="Input file name (with suffix) for UKB all PD .assoc")
parser.add_argument("--ukb_cases", default="ukb_cases", help="Input file name (with suffix) for UKB cases .assoc")
parser.add_argument("--ukb_sib", default="ukb_sib", help="Input file name (with suffix) for UKB proxy siblings .assoc")
parser.add_argument("--ukb_parent", default="ukb_parent", help="Input file name (with suffix) for UKB proxy parents .assoc")
parser.add_argument("--gne", default="gne", help="Input file name (with suffix) for GNE .assoc")

parser.add_argument("--ntotal", type=int, default=1000, help="Input number of total cases+controls")
parser.add_argument("--variant_group", help="Input string with what variant group was looked at")
parser.add_argument("--maf", help="Input string with what MAF was looked at")
parser.add_argument("--group", help="Input string with what groups are being analyzed")

parser.add_argument("--test", choices=["SkatO","CMCWald"], help="What test was run? Options are SkatO or CMCWald")
parser.add_argument("--numvar", metavar='N', type=int, help="Keep genes with variants equal to or greater than this number")
parser.add_argument("--output", "-o", help="Desired prefix output name for files")
args = parser.parse_args()

# Initialize the variables reading in user input
amp_nih = args.amp_nih
gne = args.gne
ukb_all = args.ukb_all
ukb_case = args.ukb_cases
ukb_sib = args.ukb_sib
ukb_parent = args.ukb_parent
    ## Variables needed for analysis
numvar = args.numvar
test = args.test
ntotal = args.ntotal
    ## Variables needed for final output
variant_group = args.variant_group
maf = args.maf
group = args.group #meta-group
output = args.output

# Read in refFlat file
refFlat = "/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/RECESSIVE_SEARCH/refFlat-gene-range-noDupGenes.txt"

# Define function to subset by minimum number of variants
def keep_numvar(df, min_var=numvar):
    if not df.empty:
        df = df[df['NumVar'] >= numvar]
    return df

## Define a function for opening files from RVtests (AMP, UKB) depending on test specified 
if test == "SkatO": 
    def read_rvtests(filepath):
        try:
            with open(filepath) as f:
                print
                header_text = ['GENE','RANGE','N_INFORMATIVE','NumVar','NumPolyVar','Q','rho','Pvalue']
                raw_df = pd.read_csv(filepath, delim_whitespace=True, header=0, names=header_text)
                raw_df['GENE'] = raw_df['GENE'].astype(str)
                raw_df = keep_numvar(raw_df, min_var=numvar)

        except FileNotFoundError:
            print(f'File either not found or not passed through for this analysis: {filepath}')
            raw_df = pd.DataFrame() # Make empty dataframe 

        return raw_df

    ## Make a function for opening files (GNE)
    def read_gne(filepath):
        try:
            with open(filepath) as f:
                print
                header_text = ['GENE', 'NumVar', 'Pvalue']
                raw_df = pd.read_csv(filepath, sep="\t", header=0, names=header_text)
                raw_df['GENE'] = raw_df['GENE'].astype(str)
                raw_df['Pvalue'] = raw_df['Pvalue'].astype(float)
                raw_df = keep_numvar(raw_df, min_var=numvar)

        except FileNotFoundError:
            print(f'File either not found or not passed through for this analysis: {filepath}')
            raw_df = pd.DataFrame() # Make empty dataframe 

        return raw_df
    
    def keep_cols_rvtests(df):
        df.drop(columns=['RANGE','N_INFORMATIVE','NumVar','NumPolyVar','Q','rho'], inplace=True)
        return df
    
    def keep_cols_gne(df):
        df.drop(columns=['NumVar'], inplace=True)
        return df
    
elif test == "CMCWald": 
    def read_rvtests(filepath):
        try:
            with open(filepath) as f:
                print
                header_text = ['GENE','RANGE','N_INFORMATIVE','NumVar','NumPolyVar','NonRefSite','BETA','SE','Pvalue']
                raw_df = pd.read_csv(filepath, delim_whitespace=True, header=0, names=header_text)
                raw_df['GENE'] = raw_df['GENE'].astype(str)
                raw_df = keep_numvar(raw_df, min_var=numvar)

        except FileNotFoundError:
            print(f'File either not found or not passed through for this analysis: {filepath}')
            raw_df = pd.DataFrame() # Make empty dataframe 

        return raw_df

    ## Make a function for opening files (GNE)
    def read_gne(filepath):
        try:
            with open(filepath) as f:
                print
                header_text = ['GENE', 'NumVar', 'BETA', 'SE', 'Pvalue']
                raw_df = pd.read_csv(filepath, sep="\t", header=0, names=header_text)
                raw_df['GENE'] = raw_df['GENE'].astype(str)
                raw_df['Pvalue'] = raw_df['Pvalue'].astype(float)
                raw_df = keep_numvar(raw_df, min_var=numvar)

        except FileNotFoundError:
            print(f'File either not found or not passed through for this analysis: {filepath}')
            raw_df = pd.DataFrame() # Make empty dataframe 

        return raw_df
    
    def keep_cols_rvtests(df):
        df.drop(columns=['RANGE','N_INFORMATIVE','NumVar','NumPolyVar','NonRefSite'], inplace=True)
        return df
    
    def keep_cols_gne(df):
        df.drop(columns=['NumVar'], inplace=True)
        return df
else:
    print("This is not a valid test, options are SkatO or CMCWald. Please re-run")
    exit()

## Define both meta-analysis approaches and sort output by smallest p-val
def meta_analysis(total_raw_df, final_gene_list, approach="combinedP", ntotal=0):
    match approach:
        case "combinedP":
            print("Using Fisher to combine p-values...")
            results = []
            for i in range(len(final_gene_list)):
                this_gene = final_gene_list[i]
                temp_data = total_raw_df[total_raw_df['GENE'] == this_gene]
                combined_gene = stats.combine_pvalues(temp_data['Pvalue'], method='fisher', weights=None)
                test_stat = combined_gene[0]
                p_val = combined_gene[1]
                results.append([this_gene, test_stat, p_val])
                output = pd.DataFrame(results, columns=['GENE', 'TEST_STAT', 'META_PVAL'])
            
            sorted_output = output.sort_values(by=['META_PVAL'], ascending=True)
            sorted_output.reset_index(drop=True, inplace=True)
            sorted_output['-log10'] = np.log10(sorted_output['META_PVAL']) * -1
            
            return sorted_output
        
        case "weightedZ":
            print("Using betas, se, and total n to make weighted Z score and p-values...")
            results = []
            for i in range(len(final_gene_list)):
                this_gene = final_gene_list[i]
                temp_data = total_raw_df[total_raw_df['GENE'] == this_gene].copy()
                temp_data['Z'] = temp_data['BETA']/temp_data['SE']
                temp_data['weighted_Z'] = temp_data['Z']*ntotal
                temp_data['weight_sq'] = ntotal*ntotal
                summed_Z = (temp_data['weighted_Z'].sum())/(np.sqrt(temp_data['weight_sq'].sum()))
                test_stat = summed_Z
                p_val = stats.norm.sf(abs(summed_Z))*2
                results.append([this_gene, test_stat, p_val])
                output = pd.DataFrame(results, columns=['GENE', 'TEST_STAT_WEIGHTEDZ', 'META_PVAL'])
            
            sorted_output = output.sort_values(by=['META_PVAL'], ascending=True)
            sorted_output.reset_index(drop=True, inplace=True)
            sorted_output['-log10'] = np.log10(sorted_output['META_PVAL']) * -1
            
            return sorted_output
        
        case _:
            return "Not a valid approach"

def make_final_df(df):
    merging_dfs = [refflat_genes, df, original_pvals, original_numvar]
    merged_df = reduce(lambda left,right: pd.merge(left, right, on=['GENE'], how='outer'), merging_dfs)
    sorted_merged_df = merged_df.sort_values(by=['META_PVAL'], ascending=True)
    sorted_merged_df.reset_index(drop=True, inplace=True)

    ## Add additional information in columns 
    sorted_merged_df['TEST'] = test
    sorted_merged_df['VARIANT_GROUP'] = variant_group
    sorted_merged_df['MAF'] = maf
    sorted_merged_df['MIN_NUMVAR'] = numvar
    sorted_merged_df['META_ANALYSIS_GROUP'] = group
    
    return sorted_merged_df

def save_sig(df):
    if not df.empty:
        sig_df = df[np.logical_or(df['META_PVAL'] < 0.000001, np.isnan(df['META_PVAL']))].copy()
        sig_df.dropna(subset=['META_PVAL'], how='all', inplace=True)
    return sig_df
    
## Read in files
amp_nih_raw_df = read_rvtests(amp_nih)

ukb_all_raw_df = read_rvtests(ukb_all)
ukb_case_raw_df = read_rvtests(ukb_case)
ukb_sib_raw_df = read_rvtests(ukb_sib)
ukb_parent_raw_df = read_rvtests(ukb_parent)

gne_raw_df = read_gne(gne)

## Read in refFlat
header_text = ['CHR', 'BP', 'END_BP', 'GENE']
refflat_df = pd.read_csv(refFlat, sep=" ", header=0, names=header_text)
refflat = refflat_df[refflat_df[['CHR']].apply(lambda x: x[0].isdigit(), axis=1)].copy()
refflat['CHR'] = refflat['CHR'].astype('Int64')
refflat['GENE'] = refflat['GENE'].astype('str')

refflat_genes = refflat[['GENE']].copy()

## Keeping numvar and pval columns 
if not amp_nih_raw_df.empty:
    amp_nih_numvar = amp_nih_raw_df[['GENE','NumVar']].copy()
    amp_nih_numvar.columns = ['GENE', 'AMP_NIH_NUMVAR']
    amp_nih_pvals = amp_nih_raw_df[['GENE','Pvalue']].copy()
    amp_nih_pvals.columns = ['GENE', 'AMP_NIH_PVAL']
else:
    amp_nih_numvar = pd.DataFrame()
    amp_nih_pvals = pd.DataFrame()
 
if not ukb_all_raw_df.empty:
    ukb_all_numvar = ukb_all_raw_df[['GENE', 'NumVar']].copy()
    ukb_all_numvar.columns = ['GENE', 'UKB_ALL_NUMVAR']
    ukb_all_pvals = ukb_all_raw_df[['GENE','Pvalue']].copy()
    ukb_all_pvals.columns = ['GENE', 'UKB_ALL_PVAL']
else:
    ukb_all_numvar = pd.DataFrame()
    ukb_all_pvals = pd.DataFrame()
    
if not ukb_case_raw_df.empty:
    ukb_numvar = ukb_case_raw_df[['GENE', 'NumVar']].copy()
    ukb_numvar.columns = ['GENE', 'UKB_CASE_NUMVAR']
    ukb_pvals = ukb_case_raw_df[['GENE','Pvalue']].copy()
    ukb_pvals.columns = ['GENE', 'UKB_CASE_PVAL']
else:
    ukb_numvar = pd.DataFrame()
    ukb_pvals = pd.DataFrame()

if not ukb_sib_raw_df.empty:
    ukb_sib_numvar = ukb_sib_raw_df[['GENE', 'NumVar']].copy()
    ukb_sib_numvar.columns = ['GENE', 'UKB_SIBLING_PROXY_NUMVAR']
    ukb_sib_pvals = ukb_sib_raw_df[['GENE','Pvalue']].copy()
    ukb_sib_pvals.columns = ['GENE', 'UKB_SIBLING_PROXY_PVAL']
else:
    ukb_sib_numvar = pd.DataFrame()
    ukb_sib_pvals = pd.DataFrame()
    
if not ukb_parent_raw_df.empty:
    ukb_parent_numvar = ukb_parent_raw_df[['GENE', 'NumVar']].copy()
    ukb_parent_numvar.columns = ['GENE', 'UKB_PARENT_PROXY_NUMVAR']
    ukb_parent_pvals = ukb_parent_raw_df[['GENE','Pvalue']].copy()
    ukb_parent_pvals.columns = ['GENE', 'UKB_PARENT_PROXY_PVAL']
else:
    ukb_parent_numvar = pd.DataFrame()
    ukb_parent_pvals = pd.DataFrame()

if not gne_raw_df.empty:
    gne_numvar = gne_raw_df[['GENE', 'NumVar']].copy()
    gne_numvar.columns = ['GENE', 'GNE_NUMVAR']
    gne_pvals = gne_raw_df[['GENE','Pvalue']].copy()
    gne_pvals.columns = ['GENE', 'GNE_PVAL']
else:
    gne_numvar = pd.DataFrame()
    gne_pvals = pd.DataFrame()
 
## Make dataframe lists to merge later 
pval_dfs = [refflat_genes, amp_nih_pvals, gne_pvals, ukb_all_pvals, ukb_pvals, ukb_sib_pvals, ukb_parent_pvals]
pval_dfs = [df for df in pval_dfs if not df.empty] ## Keep only non-empty dataframes 

numvar_dfs = [refflat_genes, amp_nih_numvar, gne_numvar, ukb_all_numvar, ukb_numvar, ukb_sib_numvar, ukb_parent_numvar]
numvar_dfs = [df for df in numvar_dfs if not df.empty] ## Keep only non-empty dataframes 

## Make dataframes 
original_pvals = reduce(lambda left,right: pd.merge(left, right, on=['GENE'], how='outer'), pval_dfs)
original_numvar = reduce(lambda left,right: pd.merge(left, right, on=['GENE'], how='outer'), numvar_dfs)

## Keep only the necessary columns 
if not amp_nih_raw_df.empty:
    keep_cols_rvtests(amp_nih_raw_df)

if not ukb_all_raw_df.empty:
    keep_cols_rvtests(ukb_all_raw_df)
    
if not ukb_case_raw_df.empty:
    keep_cols_rvtests(ukb_case_raw_df)
    
if not ukb_sib_raw_df.empty:
    keep_cols_rvtests(ukb_sib_raw_df)
    
if not ukb_parent_raw_df.empty:
    keep_cols_rvtests(ukb_parent_raw_df)

if not gne_raw_df.empty:
    keep_cols_gne(gne_raw_df)
    
## Make a total raw df with all the P-values for meta-analysis (using Mike's approach)
data_frames = [amp_nih_raw_df, ukb_all_raw_df, ukb_case_raw_df, ukb_sib_raw_df, ukb_parent_raw_df, gne_raw_df]
data_frames = [df for df in data_frames if not df.empty] ## Keep only non-empty dataframes 
total_raw_df = pd.concat(data_frames)

## Make final gene list 
final_gene_list = original_pvals.GENE.unique()

## Run meta-analysis; save out significant hits, and tell the user where the files are
if test == "SkatO":
    print("Analyzing SkatO input")
    skato_combinedP = meta_analysis(total_raw_df, final_gene_list, approach="combinedP")
    
    results_skato_combinedP = make_final_df(skato_combinedP)
    
    outfile_results_skato_combinedP = output + "." + str(variant_group) + ".freqUpper" + str(maf) +".numvar" + str(numvar) + ".combined_Ps.SkatO.tab"
    results_skato_combinedP.to_csv(outfile_results_skato_combinedP, sep="\t", index=False, na_rep='NA')
    
    print(f"Your final combined P meta-analysis file was created here: {outfile_results_skato_combinedP}")

    sig_results_skato_combinedP = save_sig(results_skato_combinedP)
    if not sig_results_skato_combinedP.empty:
        outfile_sig_results_skato_combinedP = output + "." + str(variant_group) + ".freqUpper" + str(maf) +".numvar" + str(numvar) + ".combined_Ps.SkatO.1e-6.tab"
        sig_results_skato_combinedP.to_csv(outfile_sig_results_skato_combinedP, sep="\t", index=False, na_rep='NA')
        print(f"You had {sig_results_skato_combinedP.shape[0]} genes <1E-6, created a summary here: {outfile_sig_results_skato_combinedP}")
    else:
        print("No genes were <1E-6")
        
elif test == "CMCWald":
    print("Analyzing CMCWald input")
    cmcwald_combinedP = meta_analysis(total_raw_df, final_gene_list, approach="combinedP")
    cmcwald_weightedZ = meta_analysis(total_raw_df, final_gene_list, approach="weightedZ", ntotal=ntotal)
    
    results_cmc_combinedP = make_final_df(cmcwald_combinedP)
    results_cmc_weightedZ = make_final_df(cmcwald_weightedZ)
    
    outfile_results_cmc_combinedP = output + "." + str(variant_group) + ".freqUpper" + str(maf) +".numvar" + str(numvar) + ".combined_Ps.CMCWald.tab"
    results_cmc_combinedP.to_csv(outfile_results_cmc_combinedP, sep="\t", index=False, na_rep='NA')

    outfile_results_cmc_weightedZ = output + "." + str(variant_group) + ".freqUpper" + str(maf) +".numvar" + str(numvar) + ".weightedZ.CMCWald.tab"
    results_cmc_weightedZ.to_csv(outfile_results_cmc_weightedZ, sep="\t", index=False, na_rep='NA')

    print(f"Your final combined P meta-analysis file was created here: {outfile_results_cmc_combinedP}")
    print(f"Your final weighted Z meta-analysis file was created here: {outfile_results_cmc_weightedZ}")
    
    sig_results_cmc_combinedP = save_sig(results_cmc_combinedP)
    if not sig_results_cmc_combinedP.empty:
        outfile_sig_results_cmc_combinedP = output + "." + str(variant_group) + ".freqUpper" + str(maf) +".numvar" + str(numvar) + ".combined_Ps.CMCWald.1e-6.tab"
        sig_results_cmc_combinedP.to_csv(outfile_sig_results_cmc_combinedP, sep="\t", index=False, na_rep='NA')
        print(f"You had {sig_results_cmc_combinedP.shape[0]} gene(s) <1E-6, created a summary here: {outfile_sig_results_cmc_combinedP}")
    else:
        print("No genes were <1E-6")
            
    sig_results_cmc_weightedZ = save_sig(results_cmc_weightedZ)
    if not sig_results_cmc_weightedZ.empty:
        outfile_sig_results_cmc_weightedZ = output + "." + str(variant_group) + ".freqUpper" + str(maf) +".numvar" + str(numvar) + ".weightedZ.CMCWald.1e-6.tab"
        sig_results_cmc_weightedZ.to_csv(outfile_sig_results_cmc_weightedZ, sep="\t", index=False, na_rep='NA')
        print(f"You had {sig_results_cmc_weightedZ.shape[0]} gene(s) <1E-6, created a summary here: {outfile_sig_results_cmc_weightedZ}")
    else:
        print("No genes were <1E-6")
