# PD Burdens Meta-analysis
	# Mary B Makarious
	# Meta-analyze by combining pval (AMPxNIH + UKB cohorts + Genentech)
	# Began JAN 2022 (Last updated 03-FEB-2022)

# Workflow 
	## GETTING STARTED
		# 0. Summary and Notes
		# 1. Getting Started 

	## META-ANALYZE BY COMBINING P-VALUES (SKAT-O AND CMC WALD)tree
		# 2. (Combined Ps) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD_ALL meta-analysis
			# ALL defined as case-control + proxies 
		# 3. (Combined Ps) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD case-control meta-analysis
		# 4. (Combined Ps) Meta-analyze: GNE + AMPxNIH + UK PD case-control + UKB siblings + UKB parents meta-analysis
		# 5. (Combined Ps) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD case-control + UKB parents meta-analysis
		# 6. Summary of Results 

	## META-ANALYZE BY SUMMED Z APPROACH (CMC WALD) 	
		# 7. (Summed Z) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD_ALL meta-analysis
			# ALL defined as case-control + both parent and sibling proxies 
		# 8. (Summed Z) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD case-control meta-analysis
		# 9. (Summed Z) Meta-analyze: GNE + AMPxNIH + UK PD case-control + UKB siblings + UKB parents meta-analysis
		# 10. (Summed Z) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD case-control + UKB parents meta-analysis
		# 11. Summary of Results 

	## META-ANALYZE BY WEIGHTED Z APPROACH (CMC WALD) 	
		# 12. (Weighted Z) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD_ALL meta-analysis
			# ALL defined as case-control + proxies 
		# 13. (Weighted Z) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD case-control meta-analysis
		# 14. (Weighted Z) Meta-analyze: GNE + AMPxNIH + UK PD case-control + UKB siblings + UKB parents meta-analysis
		# 15. (Weighted Z) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD case-control + UKB parents meta-analysis
		# 16. Summary of Results 

	## CALCULATE LAMBDAS 
		# 17. (Combined Ps) Calculate lambdas per meta-analyzed cohort
		# 18. (Summed Z) Calculate lambdas per meta-analyzed cohort
		# 19. (Weighted Z) Calculate lambdas per meta-analyzed cohort
		# 20. Calculate lambdas per individual cohort


# README by Mary Makarious; Last Updated 03-FEB-2022
	# PD Burdens meta-analysis with AMPxNIH, Genentech, and UKB cohorts 
		# Working directory on Biowulf: /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/
		# ├── AMPxNIH/
		# ├── geneburden_meta_genentech_data.sh
		# ├── GENENTECH/
		# ├── LAMBDA_SUMMARIES/
		# ├── README.txt
		# ├── RESULTS_COMBINED_P/
		# ├── RESULTS_SUMMEDZ_CMCWALD/
		# ├── RESULTS_WEIGHTEDZ_CMCWALD/
		# ├── scripts/
		# ├── SWARMS/
		# ├── UKB_EXOM_ALL_PD/
		# ├── UKB_EXOM_PD_CASE/
		# ├── UKB_EXOM_PD_PARENT/
		# └── UKB_EXOM_PD_SIBLING/

	# Results broken by minimum number of variants 
		# Combined p-values method (Fisher) for SkatO and CMC Wald: /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/RESULTS_COMBINED_P/SUMMARIES/
			# LAMBDAS/ directory included where lambda and lambda 1000 values were calculated 
		# Summed Z for CMC Wald: /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/RESULTS_SUMMEDZ_CMCWALD/SUMMARIES
			# LAMBDAS/ directory included where lambda and lambda 1000 values were calculated 
		# Weighted Z for CMC Wald: /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/RESULTS_WEIGHTEDZ_CMCWALD/SUMMARIES
			# LAMBDAS/ directory included where lambda and lambda 1000 values were calculated 

	# Cohorts
		# GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
			# Genentech vs. AMPxNIH WGS vs. UKB PD case-control+proxies 
		# GNE_AMP_NIH_UKB_CASE_CONTROL_META
			# Genentech vs. AMPxNIH WGS vs. UKB PD case-control
		# GNE_AMP_NIH_UKB_CASE_PROXIES_META
			# Genentech vs. AMPxNIH WGS vs. UKB PD case-control vs. UKB parent proxies vs. UKB sibling proxies 
		# GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META
			# Genentech vs. AMPxNIH WGS vs. UKB PD case-control vs. UKB parent proxies (no siblings)
	
	# Variant Groups
		# ALL_MISSENSE_SNPEFF
			# missense variants as defined by SnpEff
		# ALL_LOF_HC_LOFTEE
			# LOF variants assigned "HIGH CONFIDENCE" as defined by LOFTEE
		# ALL_LOF_HC_LOFTEE_and_CADD_20_VEP
			# variants that have a CADD PHRED >20 or are LOF variants that are "HIGH CONFIDENCE" as defined by LOFTEE
		# ALL_MODERATE_HIGH_IMPACT_SNPEFF
			# variants with an “IMPACT” assigned “MODERATE” or “HIGH” by SnpEff (this includes LoF, missense, indels and others)
	
	# MAF Cut-offs 
		# 0.001 (Variants with a MAF <0.1%)
		# 0.01 (Variants with a MAF <1%)

	# Differences between meta-analysis approaches
		# Combined P approach: 
			# Uses the Fisher test, can be used on any P-values
			# Takes in SkatO and CMC Wald results
		
		# Summed Z approach: 
			# Where Z-score is computed, summed across per gene, and the P-value is derived from the Z-score
			# TEST_STAT is this summed Z-score per gene (un-weighted)
			# Takes in CMC Wald results
			# Takes into account directionality
		
		# Weighted Z approach: 
			# Where Z-score is computed, multiplied by the total size of cohort to make a weighted Z-score, summed per gene by weighted Z-score divided by weight squared, and P-value is derived from weighted Z-score 
			# TEST_STAT is this summed Z-score per gene (weighted)
			# Takes in CMC Wald results
			# Takes into account directionality and total cohort size 

	## Cohort Numbers 

	# |                 Cohort                 	|  Cases 	| Controls 	|  TOTAL 	|
	# |:--------------------------------------:	|:------:	|:--------:	|:------:	|
	# |                                AMPxNIH 	|  3376 	|    4610 	|  7986 	|
	# |                                UKB ALL 	|   7806 	|    38051 	|  45857 	|
	# |                       UKB case-control 	|   1105 	|     5643 	|   6748 	|
	# |                           UKB siblings 	|    668 	|     3463 	|   4131 	|
	# |                            UKB parents 	|   6033 	|    28945 	|  34978 	|
	# |                                    GNE 	|   2700 	|     9000 	|  11700 	|
	# |                                        	|        	|          	|        	|
	# |      GNE_AMP_NIH_UKB_ALL_PD_PHENO_META 	| 13882 	|   51661 	| 65543 	|
	# |      GNE_AMP_NIH_UKB_CASE_CONTROL_META 	|  7181 	|   19253 	| 26434 	|
	# |      GNE_AMP_NIH_UKB_CASE_PROXIES_META 	| 13882 	|   51661 	| 65543 	|
	# | GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META 	| 13214 	|   48198 	| 61412 	|


##########################################################################################################################################
##########################################################################################################################################
##### 0. SUMMARY AND NOTES ###############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## Goal
	# Meta-analysis at MAF 0.1 and 0.001 for SkatO and CMC wald tests for 4 variant groups: 
		# ALL_MISSENSE_SNPEFF
		# ALL_LOF_HC_LOFTEE
		# ALL_MODERATE_HIGH_IMPACT_SNPEFF
		# ALL_LOF_HC_LOFTEE_and_CADD_20_VEP

## Data 
### from Genentech @ /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/PREP_DATA/fromGenentech/11182021
	# ├── gne.maf0.1pct.cadd20_lofhc.CMCWALD.tsv
	# ├── gne.maf0.1pct.cadd20_lofhc.SKATO.tsv
	# ├── gne.maf0.1pct.himod.CMCWALD.tsv
	# ├── gne.maf0.1pct.himod.SKATO.tsv
	# ├── gne.maf0.1pct.lofhc.CMCWALD.tsv
	# ├── gne.maf0.1pct.lofhc.SKATO.tsv
	# ├── gne.maf0.1pct.missense.CMCWALD.tsv
	# ├── gne.maf0.1pct.missense.SKATO.tsv
	# ├── gne.maf1pct.cadd20_lofhc.CMCWALD.tsv
	# ├── gne.maf1pct.cadd20_lofhc.SKATO.tsv
	# ├── gne.maf1pct.himod.CMCWALD.tsv
	# ├── gne.maf1pct.himod.SKATO.tsv
	# ├── gne.maf1pct.lofhc.CMCWALD.tsv
	# ├── gne.maf1pct.lofhc.SKATO.tsv
	# ├── gne.maf1pct.missense.CMCWALD.tsv
	# ├── gne.maf1pct.missense.SKATO.tsv
	# └── README.txt

# maf1pct: <1%
# maf0.1pct: <0.1%

# missense: missense (ALL_MISSENSE_SNPEFF)
# lofhc: LOF (ALL_LOF_HC_LOFTEE)
# cadd20_lofhc: CADD >20 or LOF (ALL_LOF_HC_LOFTEE_and_CADD_20_VEP)
# himod: All moderate or high impact events (ALL_MODERATE_HIGH_IMPACT_SNPEFF)

# Format: COHORT_${variants}.freqUpper${cutoff}.${test}.assoc

##########################################################################################################################################
##########################################################################################################################################
##### 1. GETTING STARTED #################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

module load python/3.7 #3.7

WORK_DIR="/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech"
SCRIPTS="/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/scripts"

# Set up arrays to loop through

variant_group=(
	ALL_MISSENSE_SNPEFF
	ALL_LOF_HC_LOFTEE
	ALL_MODERATE_HIGH_IMPACT_SNPEFF
	ALL_LOF_HC_LOFTEE_and_CADD_20_VEP)

test_type=(
	SkatO
	CMCWald)

MAF=(
	0.001
	0.01)

number_variants=(
	1
	2
	3
	4)

metas=(
	GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
	GNE_AMP_NIH_UKB_CASE_CONTROL_META
	GNE_AMP_NIH_UKB_CASE_PROXIES_META
	GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META)

## Rename gne files
cd /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/PREP_DATA/fromGenentech/11182021
scp *.tsv $WORK_DIR/GENENTECH

cd $WORK_DIR/GENENTECH

## Rename genentech files 
mv gne.maf0.1pct.cadd20_lofhc.CMCWALD.tsv gne.ALL_LOF_HC_LOFTEE_and_CADD_20_VEP.freqUpper0.001.CMCWald.noCovEffect.assoc
mv gne.maf0.1pct.cadd20_lofhc.SKATO.tsv gne.ALL_LOF_HC_LOFTEE_and_CADD_20_VEP.freqUpper0.001.SkatO.assoc 
mv gne.maf0.1pct.himod.CMCWALD.tsv gne.ALL_MODERATE_HIGH_IMPACT_SNPEFF.freqUpper0.001.CMCWald.noCovEffect.assoc
mv gne.maf0.1pct.himod.SKATO.tsv gne.ALL_MODERATE_HIGH_IMPACT_SNPEFF.freqUpper0.001.SkatO.assoc
mv gne.maf0.1pct.lofhc.CMCWALD.tsv gne.ALL_LOF_HC_LOFTEE.freqUpper0.001.CMCWald.noCovEffect.assoc
mv gne.maf0.1pct.lofhc.SKATO.tsv gne.ALL_LOF_HC_LOFTEE.freqUpper0.001.SkatO.assoc
mv gne.maf0.1pct.missense.CMCWALD.tsv gne.ALL_MISSENSE_SNPEFF.freqUpper0.001.CMCWald.noCovEffect.assoc
mv gne.maf0.1pct.missense.SKATO.tsv gne.ALL_MISSENSE_SNPEFF.freqUpper0.001.SkatO.assoc
mv gne.maf1pct.cadd20_lofhc.CMCWALD.tsv gne.ALL_LOF_HC_LOFTEE_and_CADD_20_VEP.freqUpper0.01.CMCWald.noCovEffect.assoc
mv gne.maf1pct.cadd20_lofhc.SKATO.tsv gne.ALL_LOF_HC_LOFTEE_and_CADD_20_VEP.freqUpper0.01.SkatO.assoc
mv gne.maf1pct.himod.CMCWALD.tsv gne.ALL_MODERATE_HIGH_IMPACT_SNPEFF.freqUpper0.01.CMCWald.noCovEffect.assoc
mv gne.maf1pct.himod.SKATO.tsv gne.ALL_MODERATE_HIGH_IMPACT_SNPEFF.freqUpper0.01.SkatO.assoc
mv gne.maf1pct.lofhc.CMCWALD.tsv gne.ALL_LOF_HC_LOFTEE.freqUpper0.01.CMCWald.noCovEffect.assoc
mv gne.maf1pct.lofhc.SKATO.tsv gne.ALL_LOF_HC_LOFTEE.freqUpper0.01.SkatO.assoc
mv gne.maf1pct.missense.CMCWALD.tsv gne.ALL_MISSENSE_SNPEFF.freqUpper0.01.CMCWald.noCovEffect.assoc
mv gne.maf1pct.missense.SKATO.tsv gne.ALL_MISSENSE_SNPEFF.freqUpper0.01.SkatO.assoc

## Might need to rename columns? 
	# SKATO: gene    numvars P
	# CMC Wald: gene    numvars beta    se      P

## Copy over cohort files to one place 
cd /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/PREP_DATA/gene-burden

scp -r AMPxNIH $WORK_DIR/
scp -r UKB_EXOM_ALL_PD $WORK_DIR/
scp -r UKB_EXOM_PD_CASE $WORK_DIR/
scp -r UKB_EXOM_PD_PARENT $WORK_DIR/
scp -r UKB_EXOM_PD_SIBLING $WORK_DIR/

## Rename CMCWald.noCovEffect to just CMCWald for downstream meta-analysis 
cd $WORK_DIR/GENENTECH
for file in *; do mv "${file}" "${file//\.noCovEffect/}"; done

cd $WORK_DIR/AMPxNIH
for file in *; do mv "${file}" "${file//\.noCovEffect/}"; done

cd $WORK_DIR/UKB_EXOM_ALL_PD
for file in *; do mv "${file}" "${file//\.noCovEffect/}"; done

cd $WORK_DIR/UKB_EXOM_PD_CASE
for file in *; do mv "${file}" "${file//\.noCovEffect/}"; done

cd $WORK_DIR/UKB_EXOM_PD_PARENT
for file in *; do mv "${file}" "${file//\.noCovEffect/}"; done

cd $WORK_DIR/UKB_EXOM_PD_SIBLING
for file in *; do mv "${file}" "${file//\.noCovEffect/}"; done

## Great, all data should be here

## Make scripts 
## CMC WALD HAS DIFFERENT HEADERS !! Scripts have been re-written to reflect this change 
#mkdir $WORK_DIR/scripts 

##########################################################################################################################################
##########################################################################################################################################
##### 2. (COMBINED PS) META-ANALYZE: GNE + AMPXNIH (WGS) + UKB PD_ALL META-ANALYSIS ######################################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UKB PD ALL

# Make directories 
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for test in "${test_type[@]}"
do
    for variants in "${variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/gne_AMPNIH_UKB_ALL_meta_analysis_numvar.py \
    			--test ${test} \
	    		--numvar ${numvar} \
				--ukb_all_cases $WORK_DIR/UKB_EXOM_ALL_PD/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
				--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021 \
				-o $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> GNE_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm
    		done
    	done
    done
done

swarm -f GNE_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 64 commands run in 8 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 8 cores and 16 cpus
# 30605188 - done

## Checks
# 2 tests * 4 variant groups * 2 MAFs = 16 files per numvar directory
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_1 | wc -l # 16 
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2 | wc -l # 16 
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3 | wc -l # 16 
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4 | wc -l # 16 


##########################################################################################################################################
##########################################################################################################################################
##### 3. (COMBINED PS) META-ANALYZE: GNE + AMPXNIH (WGS) + UKB PD CASE-CONTROL META-ANALYSIS #############################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UKB PD CASE-CONTROL

# Make directories 
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_4


cd $WORK_DIR/SWARMS

for test in "${test_type[@]}"
do
    for variants in "${variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/gne_AMPNIH_UKB_meta_analysis_numvar.py \
    			--test ${test} \
	    		--numvar ${numvar} \
				--ukb_cases $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
				--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021 \
				-o $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}" >> GNE_AMP_NIH_UKB_CASE_CONTROL_META.swarm
    		done
    	done
    done
done

swarm -f GNE_AMP_NIH_UKB_CASE_CONTROL_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 64 commands run in 8 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 8 cores and 16 cpus
# 30649846 - done

## Checks
# 2 tests * 4 variant groups * 2 MAFs = 16 files per numvar directory
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_1 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_2 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_3 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_4 | wc -l # 16

##########################################################################################################################################
##########################################################################################################################################
##### 4. (COMBINED PS) META-ANALYZE: GNE + AMPXNIH + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS  ####################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS 

# Make directories 
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for test in "${test_type[@]}"
do
    for variants in "${variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/gne_AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py \
    			--test ${test} \
	    		--numvar ${numvar} \
				--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
				--ukb_cases $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--ukb_sib $WORK_DIR/UKB_EXOM_PD_SIBLING/UKB_EXOM_PD_SIBLING_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--ukb_parent $WORK_DIR/UKB_EXOM_PD_PARENT/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021 \
				-o $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021_${variants}_freqUpper${cutoff}" >> GNE_AMP_NIH_UKB_CASE_PROXIES_META.swarm
    		done
    	done
    done
done

swarm -f GNE_AMP_NIH_UKB_CASE_PROXIES_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 64 commands run in 8 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 8 cores and 16 cpus
# 30655270 - done

## Checks
# 2 tests * 4 variant groups * 2 MAFs = 16 files per numvar directory
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_1 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_2 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_3 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_4 | wc -l # 16

##########################################################################################################################################
##########################################################################################################################################
##### 5. (COMBINED PS) META-ANALYZE: GNE + AMPXNIH + UK PD CASE-CONTROL + UKB PARENTS  ###################################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UK PD CASE-CONTROL + UKB PARENTS 

# Make directories 
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for test in "${test_type[@]}"
do
    for variants in "${variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/gne_AMPNIH_UKB_CASES_PROXIES_noSIB_meta_analysis_numvar.py \
    			--test ${test} \
	    		--numvar ${numvar} \
				--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
				--ukb_cases $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--ukb_parent $WORK_DIR/UKB_EXOM_PD_PARENT/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PARENT_PROXY_2021 \
				-o $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021_${variants}_freqUpper${cutoff}" >> GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META.swarm
    		done
    	done
    done
done

swarm -f GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 64 commands run in 8 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 8 cores and 16 cpus
# 31614538 - done 

## Checks
# 2 tests * 4 variant groups * 2 MAFs = 16 files per numvar directory
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_1 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_2 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_3 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_4 | wc -l # 16

##########################################################################################################################################
##########################################################################################################################################
##### 6. SUMMARY OF RESULTS ###############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

mkdir $WORK_DIR/RESULTS_COMBINED_P/SUMMARIES 

for meta in "${metas[@]}"
do
	for variants in "${variant_group[@]}"
	do
		for cutoff in "${MAF[@]}"
		do
			for test in "${test_type[@]}"
			do
				for numvar in "${number_variants[@]}"
		    	do
					echo "Now looking at cohort: ${meta}; variant group ${variants} at MAF ${cutoff} for ${test} at minimum numvar ${numvar}"
					echo "COHORT: ${meta}; VARIANT GROUP: ${variants}; MAF: ${cutoff}; NUMVAR >= ${numvar} for ${test}" >> $WORK_DIR/RESULTS_COMBINED_P/SUMMARIES/combinedPs_top_25_lines_${meta}_minvar${numvar}.txt
					cat $WORK_DIR/RESULTS_COMBINED_P/${meta}/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}*${test}* | head -25 >> $WORK_DIR/RESULTS_COMBINED_P/SUMMARIES/combinedPs_top_25_lines_${meta}_minvar${numvar}.txt
					echo " " >> $WORK_DIR/RESULTS_COMBINED_P/SUMMARIES/combinedPs_top_25_lines_${meta}_minvar${numvar}.txt
				done
			done
		done
	done
done

# 4 cohorts * 4 minimum numvar = 16 summary files (each summary has each type of variant group, test, and MAF cut-off)
# done

##########################################################################################################################################
##########################################################################################################################################
##### 7. (SUMMED Z) META-ANALYZE: GNE + AMPXNIH (WGS) + UKB PD_ALL META-ANALYSIS #######################################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UKB PD ALL

# Make directories 
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
		do
			echo "python $SCRIPTS/cmcwald_summedZ_gne_AMPNIH_UKB_ALL_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
			--ukb_all_cases $WORK_DIR/UKB_EXOM_ALL_PD/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021 \
			-o $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> cmc_summedZ_GNE_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm
		done
	done
done


swarm -f cmc_summedZ_GNE_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 31608018 - done

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_1 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4 | wc -l # 8

##########################################################################################################################################
##########################################################################################################################################
##### 8. (SUMMED Z) META-ANALYZE: GNE + AMPXNIH (WGS) + UKB PD CASE-CONTROL META-ANALYSIS ##############################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UKB PD CASE-CONTROL

# Make directories 
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
		do
			echo "python $SCRIPTS/cmcwald_summedZ_gne_AMPNIH_UKB_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
			--ukb_cases $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021 \
			-o $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}" >> cmc_summedZ_GNE_AMP_NIH_UKB_CASE_CONTROL_META.swarm
		done
	done
done

swarm -f cmc_summedZ_GNE_AMP_NIH_UKB_CASE_CONTROL_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 31608260 - done

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_1 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_4 | wc -l # 8


##########################################################################################################################################
##########################################################################################################################################
##### 9. (SUMMED Z) META-ANALYZE: GNE + AMPXNIH + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS  #####################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS 

# Make directories 
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
		do
			echo "python $SCRIPTS/cmcwald_summedZ_gne_AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--ukb_cases $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--ukb_sib $WORK_DIR/UKB_EXOM_PD_SIBLING/UKB_EXOM_PD_SIBLING_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--ukb_parent $WORK_DIR/UKB_EXOM_PD_PARENT/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021 \
			-o $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021_${variants}_freqUpper${cutoff}" >> cmc_summedZ_GNE_AMP_NIH_UKB_CASE_PROXIES_META.swarm
		done
	done
done

swarm -f cmc_summedZ_GNE_AMP_NIH_UKB_CASE_PROXIES_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 31609075 - done

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_1 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_4 | wc -l # 8

##########################################################################################################################################
##########################################################################################################################################
##### 10. (SUMMED Z) META-ANALYZE: GNE + AMPXNIH + UK PD CASE-CONTROL + UKB PARENTS  #####################################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UK PD CASE-CONTROL + UKB PARENTS 

# Make directories 
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
		do
			echo "python $SCRIPTS/cmcwald_summedZ_gne_AMPNIH_UKB_CASES_PROXIES_noSIB_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--ukb_cases $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--ukb_parent $WORK_DIR/UKB_EXOM_PD_PARENT/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PARENT_PROXY_2021 \
			-o $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021_${variants}_freqUpper${cutoff}" >> cmc_summedZ_GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META.swarm
		done
	done
done

swarm -f cmc_summedZ_GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 31615533 - done

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_1 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_4 | wc -l # 8


##########################################################################################################################################
##########################################################################################################################################
##### 11. SUMMARY OF RESULTS ##############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/SUMMARIES/

for meta in "${metas[@]}"
do
	for variants in "${variant_group[@]}"
	do
		for cutoff in "${MAF[@]}"
		do
			for numvar in "${number_variants[@]}"
	    	do
				echo "Now looking at cohort: ${meta}; variant group ${variants} at MAF ${cutoff} at minimum numvar ${numvar}"
				echo "COHORT: ${meta}; VARIANT GROUP: ${variants}; MAF: ${cutoff}; NUMVAR >= ${numvar}" >> $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/SUMMARIES/cmcwald_summedZ_top_25_lines_${meta}_minvar${numvar}.txt
				cat $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/${meta}/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}* | head -25 >> $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/SUMMARIES/cmcwald_summedZ_top_25_lines_${meta}_minvar${numvar}.txt
				echo " " >> $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/SUMMARIES/cmcwald_summedZ_top_25_lines_${meta}_minvar${numvar}.txt
			done
		done
	done
done

# 4 cohorts * 4 minimum numvar = 16 summary files (each summary has each type of variant group and MAF cut-off)
# done

##########################################################################################################################################
##########################################################################################################################################
##### 12. (WEIGHTED Z) META-ANALYZE: GNE + AMPXNIH (WGS) + UKB PD_ALL META-ANALYSIS ######################################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UKB PD ALL
	## Cohort Numbers 

	# |                 Cohort                 	|  Cases 	| Controls 	|  TOTAL 	|
	# |:--------------------------------------:	|:------:	|:--------:	|:------:	|
	# |                                AMPxNIH 	|  3376 	|    4610 	|  7986 	|
	# |                                UKB ALL 	|   7806 	|    38051 	|  45857 	|
	# |                       UKB case-control 	|   1105 	|     5643 	|   6748 	|
	# |                           UKB siblings 	|    668 	|     3463 	|   4131 	|
	# |                            UKB parents 	|   6033 	|    28945 	|  34978 	|
	# |                                    GNE 	|   2700 	|     9000 	|  11700 	|
	# |                                        	|        	|          	|        	|
	# |      GNE_AMP_NIH_UKB_ALL_PD_PHENO_META 	| 13882 	|   51661 	| 65543 	|
	# |      GNE_AMP_NIH_UKB_CASE_CONTROL_META 	|  7181 	|   19253 	| 26434 	|
	# |      GNE_AMP_NIH_UKB_CASE_PROXIES_META 	| 13882 	|   51661 	| 65543 	|
	# | GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META 	| 13214 	|   48198 	| 61412 	|

# Make directories 
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
		do
			echo "python $SCRIPTS/cmcwald_weightedZ_gne_AMPNIH_UKB_ALL_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
    		--ntotal 65543 \
			--ukb_all_cases $WORK_DIR/UKB_EXOM_ALL_PD/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021 \
			-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> cmc_weightedZ_GNE_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm
		done
	done
done


swarm -f cmc_weightedZ_GNE_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 31669797 - done

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_1 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4 | wc -l # 8

##########################################################################################################################################
##########################################################################################################################################
##### 13. (WEIGHTED Z) META-ANALYZE: GNE + AMPXNIH (WGS) + UKB PD CASE-CONTROL META-ANALYSIS #############################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UKB PD CASE-CONTROL

# Make directories 
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
		do
			echo "python $SCRIPTS/cmcwald_weightedZ_gne_AMPNIH_UKB_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
    		--ntotal 26434 \
			--ukb_cases $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021 \
			-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}" >> cmc_weightedZ_GNE_AMP_NIH_UKB_CASE_CONTROL_META.swarm
		done
	done
done

swarm -f cmc_weightedZ_GNE_AMP_NIH_UKB_CASE_CONTROL_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 31672730 - done 

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_1 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_4 | wc -l # 8


##########################################################################################################################################
##########################################################################################################################################
##### 14. (WEIGHTED Z) META-ANALYZE: GNE + AMPXNIH + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS  ####################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS 

# Make directories 
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
		do
			echo "python $SCRIPTS/cmcwald_weightedZ_gne_AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
    		--ntotal 65543 \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--ukb_cases $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--ukb_sib $WORK_DIR/UKB_EXOM_PD_SIBLING/UKB_EXOM_PD_SIBLING_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--ukb_parent $WORK_DIR/UKB_EXOM_PD_PARENT/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021 \
			-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021_${variants}_freqUpper${cutoff}" >> cmc_weightedZ_GNE_AMP_NIH_UKB_CASE_PROXIES_META.swarm
		done
	done
done

swarm -f cmc_weightedZ_GNE_AMP_NIH_UKB_CASE_PROXIES_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 31672004 - done

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_1 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_4 | wc -l # 8

##########################################################################################################################################
##########################################################################################################################################
##### 15. (WEIGHTED Z) META-ANALYZE: GNE + AMPXNIH + UK PD CASE-CONTROL + UKB PARENTS  ###################################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UK PD CASE-CONTROL + UKB PARENTS 

# Make directories 
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
		do
			echo "python $SCRIPTS/cmcwald_weightedZ_gne_AMPNIH_UKB_CASES_PROXIES_noSIB_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
    		--ntotal 61412 \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--ukb_cases $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--ukb_parent $WORK_DIR/UKB_EXOM_PD_PARENT/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PARENT_PROXY_2021 \
			-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021_${variants}_freqUpper${cutoff}" >> cmc_weightedZ_GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META.swarm
		done
	done
done

swarm -f cmc_weightedZ_GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 31672720 - done 

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_1 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_4 | wc -l # 8


##########################################################################################################################################
##########################################################################################################################################
##### 16. SUMMARY OF RESULTS #############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/SUMMARIES/

for meta in "${metas[@]}"
do
	for variants in "${variant_group[@]}"
	do
		for cutoff in "${MAF[@]}"
		do
			for numvar in "${number_variants[@]}"
	    	do
				echo "Now looking at cohort: ${meta}; variant group ${variants} at MAF ${cutoff} at minimum numvar ${numvar}"
				echo "COHORT: ${meta}; VARIANT GROUP: ${variants}; MAF: ${cutoff}; NUMVAR >= ${numvar}" >> $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/SUMMARIES/cmcwald_weightedz_top_25_lines_${meta}_minvar${numvar}.txt
				cat $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/${meta}/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}* | head -25 >> $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/SUMMARIES/cmcwald_weightedz_top_25_lines_${meta}_minvar${numvar}.txt
				echo " " >> $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/SUMMARIES/cmcwald_weightedz_top_25_lines_${meta}_minvar${numvar}.txt
			done
		done
	done
done

# 4 cohorts * 4 minimum numvar = 12 summary files (each summary has each type of variant group and MAF cut-off)
# done

##########################################################################################################################################
##########################################################################################################################################
##### 17. (COMBINED PS) CALCULATE LAMBDAS PER META COHORT ################################################################################
##########################################################################################################################################
##########################################################################################################################################

	# |                 Cohort                 	|  Cases 	| Controls 	|  TOTAL 	|
	# |:--------------------------------------:	|:------:	|:--------:	|:------:	|
	# |                                AMPxNIH 	|  3376 	|    4610 	|  7986 	|
	# |                                UKB ALL 	|   7806 	|    38051 	|  45857 	|
	# |                       UKB case-control 	|   1105 	|     5643 	|   6748 	|
	# |                           UKB siblings 	|    668 	|     3463 	|   4131 	|
	# |                            UKB parents 	|   6033 	|    28945 	|  34978 	|
	# |                                    GNE 	|   2700 	|     9000 	|  11700 	|
	# |                                        	|        	|          	|        	|
	# |      GNE_AMP_NIH_UKB_ALL_PD_PHENO_META 	| 13882 	|   51661 	| 65543 	|
	# |      GNE_AMP_NIH_UKB_CASE_CONTROL_META 	|  7181 	|   19253 	| 26434 	|
	# |      GNE_AMP_NIH_UKB_CASE_PROXIES_META 	| 13882 	|   51661 	| 65543 	|
	# | GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META 	| 13214 	|   48198 	| 61412 	|

cd $WORK_DIR/
mkdir $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS

cd $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for test in "${test_type[@]}"
		do
			for numvar in "${number_variants[@]}"
	    	do
			python $SCRIPTS/calclambdas_meta.py \
			--test ${test} \
    		--numvar ${numvar} \
			--input $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.combined_Ps.${test}.tab \
			--variant_group ${variants} \
			--ncases 13882 \
			--ncontrols 51661 \
			--maf ${cutoff} \
			--group GNE_AMP_NIH_UKB_ALL_PD_PHENO_META \
			-o $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META_${variants}_freqUpper${cutoff}_minvar${numvar}.combined_Ps
			done
		done
	done
done

cd $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for test in "${test_type[@]}"
		do
			for numvar in "${number_variants[@]}"
	    	do
			python $SCRIPTS/calclambdas_meta.py \
			--test ${test} \
    		--numvar ${numvar} \
			--input $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.combined_Ps.${test}.tab \
			--variant_group ${variants} \
			--ncases 7181 \
			--ncontrols 19253 \
			--maf ${cutoff} \
			--group GNE_AMP_NIH_UKB_CASE_CONTROL_META \
			-o $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS/GNE_AMP_NIH_UKB_CASE_CONTROL_META_${variants}_freqUpper${cutoff}_minvar${numvar}.combined_Ps
			done
		done
	done
done

cd $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for test in "${test_type[@]}"
		do
			for numvar in "${number_variants[@]}"
	    	do
			python $SCRIPTS/calclambdas_meta.py \
			--test ${test} \
    		--numvar ${numvar} \
			--input $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.combined_Ps.${test}.tab \
			--variant_group ${variants} \
			--ncases 13882 \
			--ncontrols 51661 \
			--maf ${cutoff} \
			--group GNE_AMP_NIH_UKB_CASE_PROXIES_META \
			-o $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS/GNE_AMP_NIH_UKB_CASE_PROXIES_META_${variants}_freqUpper${cutoff}_minvar${numvar}.combined_Ps
			done
		done
	done
done

cd $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for test in "${test_type[@]}"
		do
			for numvar in "${number_variants[@]}"
	    	do
			python $SCRIPTS/calclambdas_meta.py \
			--test ${test} \
    		--numvar ${numvar} \
			--input $WORK_DIR/RESULTS_COMBINED_P/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.combined_Ps.${test}.tab \
			--variant_group ${variants} \
			--ncases 13214 \
			--ncontrols 48198 \
			--maf ${cutoff} \
			--group GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META \
			-o $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META_${variants}_freqUpper${cutoff}_minvar${numvar}.combined_Ps
			done
		done
	done
done

## Summarize 
cd $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS/
head -1 GNE_AMP_NIH_UKB_CASE_CONTROL_META_ALL_LOF_HC_LOFTEE_freqUpper0.001_minvar1.combined_Ps.lambdas.CMCWald.tab >> lambdas_combinedPs_summary.txt
for f in *.tab; do sed -n '2p' $f >> lambdas_combinedPs_summary.txt ; done

## Copy to working directory
scp lambdas_combinedPs_summary.txt $WORK_DIR/LAMBDA_SUMMARIES
# done

##########################################################################################################################################
##########################################################################################################################################
##### 18. (SUMMED ZS) CALCULATE LAMBDAS PER META COHORT ################################################################################
##########################################################################################################################################
##########################################################################################################################################

	# |                 Cohort                 	|  Cases 	| Controls 	|  TOTAL 	|
	# |:--------------------------------------:	|:------:	|:--------:	|:------:	|
	# |                                AMPxNIH 	|  3376 	|    4610 	|  7986 	|
	# |                                UKB ALL 	|   7806 	|    38051 	|  45857 	|
	# |                       UKB case-control 	|   1105 	|     5643 	|   6748 	|
	# |                           UKB siblings 	|    668 	|     3463 	|   4131 	|
	# |                            UKB parents 	|   6033 	|    28945 	|  34978 	|
	# |                                    GNE 	|   2700 	|     9000 	|  11700 	|
	# |                                        	|        	|          	|        	|
	# |      GNE_AMP_NIH_UKB_ALL_PD_PHENO_META 	| 13882 	|   51661 	| 65543 	|
	# |      GNE_AMP_NIH_UKB_CASE_CONTROL_META 	|  7181 	|   19253 	| 26434 	|
	# |      GNE_AMP_NIH_UKB_CASE_PROXIES_META 	| 13882 	|   51661 	| 65543 	|
	# | GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META 	| 13214 	|   48198 	| 61412 	|

cd $WORK_DIR/
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/LAMBDAS

cd $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
    	do
		python $SCRIPTS/calclambdas_meta.py \
		--test CMCWald \
		--numvar ${numvar} \
		--input $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.summedZ.CMCWald.tab \
		--variant_group ${variants} \
		--ncases 13882 \
		--ncontrols 51661 \
		--maf ${cutoff} \
		--group GNE_AMP_NIH_UKB_ALL_PD_PHENO_META \
		-o $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/LAMBDAS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META_${variants}_freqUpper${cutoff}_minvar${numvar}.summedZ
		done
	done
done

cd $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
    	do
		python $SCRIPTS/calclambdas_meta.py \
		--test CMCWald \
		--numvar ${numvar} \
		--input $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.summedZ.CMCWald.tab \
		--variant_group ${variants} \
		--ncases 7181 \
		--ncontrols 19253 \
		--maf ${cutoff} \
		--group GNE_AMP_NIH_UKB_CASE_CONTROL_META \
		-o $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/LAMBDAS/GNE_AMP_NIH_UKB_CASE_CONTROL_META_${variants}_freqUpper${cutoff}_minvar${numvar}.summedZ
		done
	done
done

cd $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
    	do
		python $SCRIPTS/calclambdas_meta.py \
		--test CMCWald \
		--numvar ${numvar} \
		--input $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.summedZ.CMCWald.tab \
		--variant_group ${variants} \
		--ncases 13882 \
		--ncontrols 51661 \
		--maf ${cutoff} \
		--group GNE_AMP_NIH_UKB_CASE_PROXIES_META \
		-o $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/LAMBDAS/GNE_AMP_NIH_UKB_CASE_PROXIES_META_${variants}_freqUpper${cutoff}_minvar${numvar}.summedZ
		done
	done
done

cd $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
    	do
		python $SCRIPTS/calclambdas_meta.py \
		--test CMCWald \
		--numvar ${numvar} \
		--input $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.summedZ.CMCWald.tab \
		--variant_group ${variants} \
		--ncases 13214 \
		--ncontrols 48198 \
		--maf ${cutoff} \
		--group GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META \
		-o $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/LAMBDAS/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META_${variants}_freqUpper${cutoff}_minvar${numvar}.summedZ
		done
	done
done

## Summarize 
cd $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/LAMBDAS/
head -1 GNE_AMP_NIH_UKB_ALL_PD_PHENO_META_ALL_LOF_HC_LOFTEE_freqUpper0.01_minvar4.summedZ.lambdas.CMCWald.tab >> lambdas_cmcwald_summedZ_summary.txt
for f in *.tab; do sed -n '2p' $f >> lambdas_cmcwald_summedZ_summary.txt ; done

## Copy to working directory
scp lambdas_cmcwald_summedZ_summary.txt $WORK_DIR/LAMBDA_SUMMARIES
# done 

##########################################################################################################################################
##########################################################################################################################################
##### 19. (WEIGHTED ZS) CALCULATE LAMBDAS PER META COHORT ################################################################################
##########################################################################################################################################
##########################################################################################################################################

	# |                 Cohort                 	|  Cases 	| Controls 	|  TOTAL 	|
	# |:--------------------------------------:	|:------:	|:--------:	|:------:	|
	# |                                AMPxNIH 	|  3376 	|    4610 	|  7986 	|
	# |                                UKB ALL 	|   7806 	|    38051 	|  45857 	|
	# |                       UKB case-control 	|   1105 	|     5643 	|   6748 	|
	# |                           UKB siblings 	|    668 	|     3463 	|   4131 	|
	# |                            UKB parents 	|   6033 	|    28945 	|  34978 	|
	# |                                    GNE 	|   2700 	|     9000 	|  11700 	|
	# |                                        	|        	|          	|        	|
	# |      GNE_AMP_NIH_UKB_ALL_PD_PHENO_META 	| 13882 	|   51661 	| 65543 	|
	# |      GNE_AMP_NIH_UKB_CASE_CONTROL_META 	|  7181 	|   19253 	| 26434 	|
	# |      GNE_AMP_NIH_UKB_CASE_PROXIES_META 	| 13882 	|   51661 	| 65543 	|
	# | GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META 	| 13214 	|   48198 	| 61412 	|

cd $WORK_DIR/
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/LAMBDAS

cd $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
    	do
		python $SCRIPTS/calclambdas_meta.py \
		--test CMCWald \
		--numvar ${numvar} \
		--input $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.weightedZ.CMCWald.tab \
		--variant_group ${variants} \
		--ncases 13882 \
		--ncontrols 51661 \
		--maf ${cutoff} \
		--group GNE_AMP_NIH_UKB_ALL_PD_PHENO_META \
		-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/LAMBDAS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META_${variants}_freqUpper${cutoff}_minvar${numvar}.weightedZ
		done
	done
done

cd $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
    	do
		python $SCRIPTS/calclambdas_meta.py \
		--test CMCWald \
		--numvar ${numvar} \
		--input $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.weightedZ.CMCWald.tab \
		--variant_group ${variants} \
		--ncases 7181 \
		--ncontrols 19253 \
		--maf ${cutoff} \
		--group GNE_AMP_NIH_UKB_CASE_CONTROL_META \
		-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/LAMBDAS/GNE_AMP_NIH_UKB_CASE_CONTROL_META_${variants}_freqUpper${cutoff}_minvar${numvar}.weightedZ
		done
	done
done

cd $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
    	do
		python $SCRIPTS/calclambdas_meta.py \
		--test CMCWald \
		--numvar ${numvar} \
		--input $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.weightedZ.CMCWald.tab \
		--variant_group ${variants} \
		--ncases 13882 \
		--ncontrols 51661 \
		--maf ${cutoff} \
		--group GNE_AMP_NIH_UKB_CASE_PROXIES_META \
		-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/LAMBDAS/GNE_AMP_NIH_UKB_CASE_PROXIES_META_${variants}_freqUpper${cutoff}_minvar${numvar}.weightedZ
		done
	done
done

cd $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
    	do
		python $SCRIPTS/calclambdas_meta.py \
		--test CMCWald \
		--numvar ${numvar} \
		--input $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.weightedZ.CMCWald.tab \
		--variant_group ${variants} \
		--ncases 13214 \
		--ncontrols 48198 \
		--maf ${cutoff} \
		--group GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META \
		-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/LAMBDAS/GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META_${variants}_freqUpper${cutoff}_minvar${numvar}.weightedZ
		done
	done
done

## Summarize 
cd $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/LAMBDAS/
head -1 GNE_AMP_NIH_UKB_ALL_PD_PHENO_META_ALL_LOF_HC_LOFTEE_freqUpper0.01_minvar4.weightedZ.lambdas.CMCWald.tab >> lambdas_cmcwald_weightedZ_summary.txt
for f in *.tab; do sed -n '2p' $f >> lambdas_cmcwald_weightedZ_summary.txt ; done

## Copy to working directory
scp lambdas_cmcwald_weightedZ_summary.txt $WORK_DIR/LAMBDA_SUMMARIES
# done

##########################################################################################################################################
##########################################################################################################################################
##### 20. CALCULATE LAMBDAS PER INDIVIDUAL COHORT ########################################################################################
##########################################################################################################################################
##########################################################################################################################################

	# |                 Cohort                 	|  Cases 	| Controls 	|  TOTAL 	|
	# |:--------------------------------------:	|:------:	|:--------:	|:------:	|
	# |                                AMPxNIH 	|  3376 	|    4610 	|  7986 	|
	# |                                UKB ALL 	|   7806 	|    38051 	|  45857 	|
	# |                       UKB case-control 	|   1105 	|     5643 	|   6748 	|
	# |                           UKB siblings 	|    668 	|     3463 	|   4131 	|
	# |                            UKB parents 	|   6033 	|    28945 	|  34978 	|
	# |                                    GNE 	|   2700 	|     9000 	|  11700 	|
	# |                                        	|        	|          	|        	|
	# |      GNE_AMP_NIH_UKB_ALL_PD_PHENO_META 	| 13882 	|   51661 	| 65543 	|
	# |      GNE_AMP_NIH_UKB_CASE_CONTROL_META 	|  7181 	|   19253 	| 26434 	|
	# |      GNE_AMP_NIH_UKB_CASE_PROXIES_META 	| 13882 	|   51661 	| 65543 	|
	# | GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META 	| 13214 	|   48198 	| 61412 	|


cd $WORK_DIR/
mkdir $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS_PER_COHORT

## AMPxNIH

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for test in "${test_type[@]}"
		do
			for numvar in "${number_variants[@]}"
			do
			python $SCRIPTS/calclambdas_percohort.py \
			--test ${test} \
			--input $WORK_DIR/AMPxNIH/*${variants}_freqUpper${cutoff}.${test}.assoc \
			--pval Pvalue \
			--variant_group ${variants} \
			--ncases 3376 \
			--ncontrols 4610 \
			--maf ${cutoff} \
			--numvar_col NumVar \
			--numvar ${numvar} \
			--group AMP_NIH \
			-o $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS_PER_COHORT/AMP_NIH_${variants}_freqUpper${cutoff}_minvar${numvar}
			done
		done
	done
done

# UKB ALL 

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for test in "${test_type[@]}"
		do
			for numvar in "${number_variants[@]}"
			do
			python $SCRIPTS/calclambdas_percohort.py \
			--test ${test} \
			--input $WORK_DIR/UKB_EXOM_ALL_PD/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--pval Pvalue \
			--variant_group ${variants} \
			--ncases 7806 \
			--ncontrols 38051 \
			--maf ${cutoff} \
			--numvar_col NumVar \
			--numvar ${numvar} \
			--group UKB_EXOM_ALL_PD \
			-o $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS_PER_COHORT/UKB_EXOM_ALL_PD_${variants}_freqUpper${cutoff}_minvar${numvar}
		done
		done
	done
done

## done

# UKB case-control

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for test in "${test_type[@]}"
		do
			for numvar in "${number_variants[@]}"
			do
			python $SCRIPTS/calclambdas_percohort.py \
			--test ${test} \
			--input $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--pval Pvalue \
			--variant_group ${variants} \
			--ncases 1105 \
			--ncontrols 5643 \
			--maf ${cutoff} \
			--numvar_col NumVar \
			--numvar ${numvar} \
			--group UKB_EXOM_PD_CASE \
			-o $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS_PER_COHORT/UKB_EXOM_PD_CASE_${variants}_freqUpper${cutoff}_minvar${numvar}
		done
		done
	done
done


# UKB siblings 

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for test in "${test_type[@]}"
		do
			for numvar in "${number_variants[@]}"
			do
			python $SCRIPTS/calclambdas_percohort.py \
			--test ${test} \
			--input $WORK_DIR/UKB_EXOM_PD_SIBLING/UKB_EXOM_PD_SIBLING_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--pval Pvalue \
			--variant_group ${variants} \
			--ncases 668 \
			--ncontrols 3463 \
			--maf ${cutoff} \
			--numvar_col NumVar \
			--numvar ${numvar} \
			--group UKB_EXOM_PD_SIBLING \
			-o $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS_PER_COHORT/UKB_EXOM_PD_SIBLING_${variants}_freqUpper${cutoff}_minvar${numvar}
		done
		done
	done
done

# UKB parents

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for test in "${test_type[@]}"
		do
			for numvar in "${number_variants[@]}"
			do
			python $SCRIPTS/calclambdas_percohort.py \
			--test ${test} \
			--input $WORK_DIR/UKB_EXOM_PD_PARENT/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--pval Pvalue \
			--variant_group ${variants} \
			--ncases 6033 \
			--ncontrols 28945 \
			--maf ${cutoff} \
			--numvar_col NumVar \
			--numvar ${numvar} \
			--group UKB_EXOM_PD_PARENT \
			-o $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS_PER_COHORT/UKB_EXOM_PD_PARENT_${variants}_freqUpper${cutoff}_minvar${numvar}
		done
		done
	done
done

# GNE
for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for test in "${test_type[@]}"
		do
			for numvar in "${number_variants[@]}"
			do
			python $SCRIPTS/calclambdas_percohort.py \
			--test ${test} \
			--input $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.${test}.assoc \
			--pval P \
			--variant_group ${variants} \
			--ncases 2700 \
			--ncontrols 9000 \
			--maf ${cutoff} \
			--numvar_col numvars \
			--numvar ${numvar} \
			--group GENENTECH \
			-o $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS_PER_COHORT/GENENTECH_${variants}_freqUpper${cutoff}_minvar${numvar}
		done
		done
	done
done

## Summarize 
cd $WORK_DIR/RESULTS_COMBINED_P/LAMBDAS_PER_COHORT/
head -1 AMP_NIH_ALL_LOF_HC_LOFTEE_freqUpper0.001_minvar1.lambdas.CMCWald.tab >> lambdas_per_cohort_summary.txt
for f in *.tab; do sed -n '2p' $f >> lambdas_per_cohort_summary.txt ; done

## Copy to working directory
scp lambdas_per_cohort_summary.txt $WORK_DIR/LAMBDA_SUMMARIES
# done 



