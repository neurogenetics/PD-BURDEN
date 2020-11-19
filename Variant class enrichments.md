## Variant class enrichments


 => sandbox now variant class enrichments per disease group
 
 
 #### Code Mike below
 Two parts:
 - part1 :
 - part2 :
 
```
# part1
# Get on an interactive biowulf node.

sinteractive --mem 120g

# Set your working directory for the fun to begin.

cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/mike_scratch

# Make one list of SNPs of interest. This permits you to make one smaller genotype file and keep it for all PRS analyses. Also, make one smaple file.

filepath=/data/CARD/UKBIOBANK/EXOME_DATA_200K/annotation_of_plink_files
cp $filepath/all_frameshift.txt .
cp $filepath/all_missense.txt .
cp $filepath/all_nonframeshift.txt .
cp $filepath/all_splice_normal.txt .
cp $filepath/all_stopgain.txt .
cp $filepath/all_stoploss.txt .
cp $filepath/ALL_CADD_10.txt .
cp $filepath/ALL_CADD_20.txt .
cp $filepath/ALL_LOF.txt .
cp $filepath/ALL_MISSENSE_and_LOF.txt .
cat *.txt | grep ':' > variant_list

filepath=/data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN
cp $filepath/UKB_EXOM_AD_CASE_CONTROL_with_PC.txt .
cp $filepath/UKB_EXOM_AD_PARENT_CONTROL_with_PC.txt .
cp $filepath/UKB_EXOM_PD_CASE_CONTROL_with_PC.txt .
cp $filepath/UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt .
cat UKB_EXOM_AD_CASE_CONTROL_with_PC.txt UKB_EXOM_AD_PARENT_CONTROL_with_PC.txt | cut -f -2 > AD_sample_list
cat UKB_EXOM_PD_CASE_CONTROL_with_PC.txt UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt | cut -f -2 > PD_sample_list

# Now we make the genotypes for the AD and PD series. Remember, we always do AD first ;-).
# These need to eventually get to one file.

module load plink/2.0_alpha_1_final

genopath=/data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files

for CHRNUM in {1..23}
do
  echo "Working on chromosome $CHRNUM in the ADRD samples" 
  plink2 --bed $genopath/ukb23155_c${CHRNUM}_b0_v1.bed --bim $genopath/UKBexomeOQFE_chr${CHRNUM}.bim --fam $genopath/ukb23155_c1_b0_v1_s200632.fam --extract variant_list --keep AD_sample_list --make-bed --out AD_chr$CHRNUM
  echo "Working on chromosome $CHRNUM in the PD samples" 
  plink2 --bed $genopath/ukb23155_c${CHRNUM}_b0_v1.bed --bim $genopath/UKBexomeOQFE_chr${CHRNUM}.bim --fam $genopath/ukb23155_c1_b0_v1_s200632.fam --extract variant_list --keep PD_sample_list --make-bed --out PD_chr$CHRNUM
done


# Merge and then delete. Note the block below has been commented out as the files are too big. Looks like we'll have to do this the hard way.

module load plink/1.9

# plink --bfile AD_chr1 --merge-list AD_merge_list --make-bed --out AD_merged

# plink --bfile PD_chr1 --merge-list PD_merge_list --make-bed --out PD_merged

# First we need the minor alleles per variant set per chromosome. AD dataset is large enough to get good estimates of MAF.
for CHRNUM in {1..23}
do
  echo "Working on chromosome $CHRNUM in the ADRD samples" 
  plink --bfile AD_chr$CHRNUM --extract all_frameshift.txt --freq --out frameshift_chr$CHRNUM
  awk '{print $2"\t"$3"\t""1"}' frameshift_chr$CHRNUM.frq > frameshift.chr$CHRNUM.to_score
  plink --bfile AD_chr$CHRNUM --extract all_missense.txt --freq --out missense_chr$CHRNUM
  awk '{print $2"\t"$3"\t""1"}' missense_chr$CHRNUM.frq > missense.chr$CHRNUM.to_score
  plink --bfile AD_chr$CHRNUM --extract all_nonframeshift.txt --freq --out nonframeshift_chr$CHRNUM
  awk '{print $2"\t"$3"\t""1"}' nonframeshift_chr$CHRNUM.frq > nonframeshift.chr$CHRNUM.to_score
  plink --bfile AD_chr$CHRNUM --extract all_splice_normal.txt --freq --out splice_normal_chr$CHRNUM
  awk '{print $2"\t"$3"\t""1"}' splice_normal_chr$CHRNUM.frq > splice_normal.chr$CHRNUM.to_score
  plink --bfile AD_chr$CHRNUM --extract all_stopgain.txt --freq --out stopgain_chr$CHRNUM
  awk '{print $2"\t"$3"\t""1"}' stopgain_chr$CHRNUM.frq > stopgain.chr$CHRNUM.to_score
  plink --bfile AD_chr$CHRNUM --extract all_stoploss.txt --freq --out stoploss_chr$CHRNUM
  awk '{print $2"\t"$3"\t""1"}' stoploss_chr$CHRNUM.frq > stoploss.chr$CHRNUM.to_score
  plink --bfile AD_chr$CHRNUM --extract ALL_CADD_10.txt --freq --out CADD_10_chr$CHRNUM
  awk '{print $2"\t"$3"\t""1"}' CADD_10_chr$CHRNUM.frq > CADD_10.chr$CHRNUM.to_score
  plink --bfile AD_chr$CHRNUM --extract ALL_CADD_20.txt --freq --out CADD_20_chr$CHRNUM
  awk '{print $2"\t"$3"\t""1"}' CADD_20_chr$CHRNUM.frq > CADD_20.chr$CHRNUM.to_score
  plink --bfile AD_chr$CHRNUM --extract ALL_LOF.txt --freq --out LOF_chr$CHRNUM
  awk '{print $2"\t"$3"\t""1"}' LOF_chr$CHRNUM.frq > LOF.chr$CHRNUM.to_score
  plink --bfile AD_chr$CHRNUM --extract ALL_MISSENSE_and_LOF.txt --freq --out missense_and_LOF_chr$CHRNUM
  awk '{print $2"\t"$3"\t""1"}' missense_and_LOF_chr$CHRNUM.frq > missense_and_LOF.chr$CHRNUM.to_score
done

# Now its the home stretch, calculate the cummulative dosages per class using PLINK.
for CHRNUM in {1..23}
# for CHRNUM in {1..1}
do
  for var_type in frameshift missense nonframeshift splice_normal stopgain stoploss CADD_10 CADD_20 LOF missense_and_LOF
  do
    echo "Working on aggregating chromosome $CHRNUM for the $var_type variants in the ADRD samples"
    plink --bfile AD_chr$CHRNUM --score $var_type.chr$CHRNUM.to_score 1 2 header sum --out $var_type.chr$CHRNUM.AD_scored 
    echo "Working on aggregating chromosome $CHRNUM for the $var_type variants in the PD samples"
    plink --bfile PD_chr$CHRNUM --score $var_type.chr$CHRNUM.to_score 1 2 header sum --out $var_type.chr$CHRNUM.PD_scored 
  done
done

# Get everything together now.
for var_type in frameshift missense nonframeshift splice_normal stopgain stoploss CADD_10 CADD_20 LOF missense_and_LOF
do
  paste $var_type.chr1.AD_scored.profile $var_type.chr2.AD_scored.profile $var_type.chr3.AD_scored.profile $var_type.chr4.AD_scored.profile $var_type.chr5.AD_scored.profile $var_type.chr6.AD_scored.profile $var_type.chr7.AD_scored.profile $var_type.chr8.AD_scored.profile $var_type.chr9.AD_scored.profile $var_type.chr10.AD_scored.profile $var_type.chr11.AD_scored.profile $var_type.chr12.AD_scored.profile $var_type.chr13.AD_scored.profile $var_type.chr14.AD_scored.profile $var_type.chr15.AD_scored.profile $var_type.chr16.AD_scored.profile $var_type.chr17.AD_scored.profile $var_type.chr18.AD_scored.profile $var_type.chr19.AD_scored.profile $var_type.chr20.AD_scored.profile $var_type.chr21.AD_scored.profile $var_type.chr22.AD_scored.profile $var_type.chr23.AD_scored.profile > $var_type.AD_scored.aggregate_profile
  paste $var_type.chr1.PD_scored.profile $var_type.chr2.PD_scored.profile $var_type.chr3.PD_scored.profile $var_type.chr4.PD_scored.profile $var_type.chr5.PD_scored.profile $var_type.chr6.PD_scored.profile $var_type.chr7.PD_scored.profile $var_type.chr8.PD_scored.profile $var_type.chr9.PD_scored.profile $var_type.chr10.PD_scored.profile $var_type.chr11.PD_scored.profile $var_type.chr12.PD_scored.profile $var_type.chr13.PD_scored.profile $var_type.chr14.PD_scored.profile $var_type.chr15.PD_scored.profile $var_type.chr16.PD_scored.profile $var_type.chr17.PD_scored.profile $var_type.chr18.PD_scored.profile $var_type.chr19.PD_scored.profile $var_type.chr20.PD_scored.profile $var_type.chr21.PD_scored.profile $var_type.chr22.PD_scored.profile $var_type.chr23.PD_scored.profile > $var_type.PD_scored.aggregate_profile
done

# rm AD_chr*.bed
# rm AD_chr*.bim
# rm AD_chr*.fam
# rm PD_chr*.bed
# rm PD_chr*.bim
# rm PD_chr*.fam

# From here on, its just some basic python work (see separate file for that code)!

# In that code, we just sum the CNT2 columns and run some regressions (both including and excluding CHR23).
```

```
# part2
# Imports and prep.
import h5py
import numpy as np
import pandas as pd
import math
import sys
import joblib
import subprocess
import statsmodels.api as sm
from scipy import stats

# Set the big header.
a_list = ['FID','IID','PHENO','CNT','CNT2','SCORESUM']
b_list = ['CHR1','CHR2','CHR3','CHR4','CHR5','CHR6','CHR7','CHR8','CHR9','CHR10','CHR11','CHR12','CHR13','CHR14','CHR15','CHR16','CHR17','CHR18','CHR19','CHR20','CHR21','CHR22','CHR23']
header_text = ['{}_{}'.format(a, b) for b in b_list for a in a_list]

# Now read in the files and test one.

results = []

var_type_list = ['frameshift','missense','nonframeshift','splice_normal','stopgain','stoploss','CADD_10','CADD_20','LOF','missense_and_LOF']
disease_type_list = ['AD','PD']
gwax_switch_list = ['CASE','PARENT']

for x in range(len(var_type_list)):
  for y in range(len(disease_type_list)):
    for z in range(len(gwax_switch_list)):
      var_type = var_type_list[x] # Set the variables.
      disease_type = disease_type_list[y]
      gwax_switch = gwax_switch_list[z]
      score_file = var_type + '.' + disease_type + "_scored.aggregate_profile" # Set the file paths.
      metadata_file = "UKB_EXOM_" + disease_type + "_" +  gwax_switch + "_CONTROL_with_PC.txt"
      scores_df = pd.read_csv(score_file, delim_whitespace=True, header=0, names=header_text) # Read in the data.
      metadata_df = pd.read_csv(metadata_file, delim_whitespace=True)
      data_df = metadata_df.merge(scores_df, left_on='IID', right_on='IID_CHR1', how='inner') # Start data management.
      data_df['all_CHRs_CNT2'] = data_df['CNT2_CHR23'] + data_df['CNT2_CHR1'] + data_df['CNT2_CHR2'] + data_df['CNT2_CHR3'] + data_df['CNT2_CHR4'] + data_df['CNT2_CHR5'] + data_df['CNT2_CHR6'] + data_df['CNT2_CHR7'] + data_df['CNT2_CHR8'] + data_df['CNT2_CHR9'] + data_df['CNT2_CHR10'] + data_df['CNT2_CHR11'] + data_df['CNT2_CHR12'] + data_df['CNT2_CHR13'] + data_df['CNT2_CHR14'] + data_df['CNT2_CHR15'] + data_df['CNT2_CHR16'] + data_df['CNT2_CHR17'] + data_df['CNT2_CHR18'] + data_df['CNT2_CHR19'] + data_df['CNT2_CHR20'] + data_df['CNT2_CHR21'] + data_df['CNT2_CHR22']
      data_df['autosomal_CHRs_CNT2'] = data_df['CNT2_CHR1'] + data_df['CNT2_CHR2'] + data_df['CNT2_CHR3'] + data_df['CNT2_CHR4'] + data_df['CNT2_CHR5'] + data_df['CNT2_CHR6'] + data_df['CNT2_CHR7'] + data_df['CNT2_CHR8'] + data_df['CNT2_CHR9'] + data_df['CNT2_CHR10'] + data_df['CNT2_CHR11'] + data_df['CNT2_CHR12'] + data_df['CNT2_CHR13'] + data_df['CNT2_CHR14'] + data_df['CNT2_CHR15'] + data_df['CNT2_CHR16'] + data_df['CNT2_CHR17'] + data_df['CNT2_CHR18'] + data_df['CNT2_CHR19'] + data_df['CNT2_CHR20'] + data_df['CNT2_CHR21'] + data_df['CNT2_CHR22']
      data_df['sex_CHR_CNT2'] = data_df['CNT2_CHR23']
      data_df['OUTCOME'] = data_df['PHENO'] - 1 # From here down is jsut regressions and exports.
      all_CHRs_fitted = sm.formula.glm(formula="OUTCOME ~ all_CHRs_CNT2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + TOWNSEND + GENETIC_SEX + AGE_OF_RECRUIT", family=sm.families.Binomial(), data=data_df).fit()
      all_CHRs_beta_coef  = all_CHRs_fitted.params.loc['all_CHRs_CNT2']
      all_CHRs_beta_se  = all_CHRs_fitted.bse.loc['all_CHRs_CNT2']
      all_CHRs_p_val = all_CHRs_fitted.pvalues.loc['all_CHRs_CNT2']
      autosomal_CHRs_fitted = sm.formula.glm(formula="OUTCOME ~ autosomal_CHRs_CNT2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + TOWNSEND + GENETIC_SEX + AGE_OF_RECRUIT", family=sm.families.Binomial(), data=data_df).fit()
      autosomal_CHRs_beta_coef  = autosomal_CHRs_fitted.params.loc['autosomal_CHRs_CNT2']
      autosomal_CHRs_beta_se  = autosomal_CHRs_fitted.bse.loc['autosomal_CHRs_CNT2']
      autosomal_CHRs_p_val = autosomal_CHRs_fitted.pvalues.loc['autosomal_CHRs_CNT2']
      sex_CHR_fitted = sm.formula.glm(formula="OUTCOME ~ sex_CHR_CNT2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + TOWNSEND + GENETIC_SEX + AGE_OF_RECRUIT", family=sm.families.Binomial(), data=data_df).fit()
      sex_CHR_beta_coef  = sex_CHR_fitted.params.loc['sex_CHR_CNT2']
      sex_CHR_beta_se  = sex_CHR_fitted.bse.loc['sex_CHR_CNT2']
      sex_CHR_p_val = sex_CHR_fitted.pvalues.loc['sex_CHR_CNT2']
      print(var_type, disease_type, gwax_switch, all_CHRs_beta_coef, all_CHRs_beta_se, all_CHRs_p_val, autosomal_CHRs_beta_coef, autosomal_CHRs_beta_se, autosomal_CHRs_p_val, sex_CHR_beta_coef, sex_CHR_beta_se, sex_CHR_p_val)
      results.append((var_type, disease_type, gwax_switch, all_CHRs_beta_coef, all_CHRs_beta_se, all_CHRs_p_val, autosomal_CHRs_beta_coef, autosomal_CHRs_beta_se, autosomal_CHRs_p_val, sex_CHR_beta_coef, sex_CHR_beta_se, sex_CHR_p_val))

output = pd.DataFrame(results, columns=('var_type', 'disease_type', 'gwax_or_gwas', 'all_CHRs_BETA_COEF', 'all_CHRs_BETA_SE','all_CHRs_P_VAL', 'autosomal_CHRs_BETA_COEF', 'autosomal_CHRs_BETA_SE','autosomal_CHRs_P_VAL', 'sex_CHR_BETA_COEF', 'sex_CHR_BETA_SE','sex_CHR_P_VAL'))
output.to_csv("results_11032020.csv", index=False)

```

#### Results first pass...

```

```


