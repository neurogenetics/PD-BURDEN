# PD Burdens Meta-analysis
	# Mary B Makarious
	# Meta-analyze by combining pval (AMPxNIH + UKB cohorts + Genentech)
	# JAN 2022

# Workflow 
	# 0. Summary and Notes
	# 1. Getting Started 
## META-ANALYZE BY COMBINING P-VALUES (SKAT-O AND CMC WALD)
	# 2. (Combined Ps) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD_ALL + Genentech meta-analysis
		# ALL defined as case-control + proxies 
	# 3. (Combined Ps) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD case-control + Genentech meta-analysis
	# 4. (Combined Ps) Meta-analyze: GNE + AMPxNIH + UK PD case-control + UKB siblings + UKB parents + Genentech meta-analysis
	# 5. Summary of Results 
## META-ANALYZE BY WEIGHTED Z APPROACH (CMC WALD) 	
	# 6. (Weighted Z) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD_ALL + Genentech meta-analysis
		# ALL defined as case-control + proxies 
	# 7. (Weighted Z) Meta-analyze: GNE + AMPxNIH (WGS) + UKB PD case-control + Genentech meta-analysis
	# 8. (Weighted Z) Meta-analyze: GNE + AMPxNIH + UK PD case-control + UKB siblings + UKB parents + Genentech meta-analysis
	# 9. Summary of Results 

# README by Mary Makarious; Last Updated 24-JAN-2022
	# PD Burdens meta-analysis with AMPxNIH, Genentech, and UKB cohorts 

	# Results broken by minimum number of variants 
		# Combined p-values method (Fisher) for SkatO and CMC Wald: /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/RESULTS/SUMMARIES/
		# Weighted Z for CMC Wald: /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/RESULTS_WEIGHTEDZ_CMCWALD/SUMMARIES

	# Cohorts
		# GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
			# Genentech vs. AMPxNIH WGS vs. UKB PD case-control+proxies 
		# GNE_AMP_NIH_UKB_CASE_CONTROL_META
			# Genentech vs. AMPxNIH WGS vs. UKB PD case-control
		# GNE_AMP_NIH_UKB_CASE_PROXIES_META
			# Genentech vs. AMPxNIH WGS vs. UKB PD case-control vs. UKB parent proxies vs. UKB sibling proxies 
	
	# Variant Groups
		# ALL_MISSENSE_SNPEFF
			# missense variants as defined by SnpEff
		# ALL_LOF_HC_LOFTEE
			# LOF variants assigned "HIGH CONFIDENCE" as defined by LOFTEE
		# ALL_LOF_HC_LOFTEE_and_CADD_20_VEP
			# variants that have a CADD PHRED >20 or are LOF variants that are "HIGH CONFIDENCE" as defined by LOFTEE
		# ALL_MODERATE_HIGH_IMPACT_SNPEFF
			# variants with an “IMPACT” assigned “MODERATE” or “HIGH” by SnpEff (this includes LoF, missense, indels and others)
	
	# MAF
		# 0.001
		# 0.01

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
	GNE_AMP_NIH_UKB_CASE_PROXIES_META)

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
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4

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
				-o $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> GNE_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm
    		done
    	done
    done
done

swarm -f GNE_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 64 commands run in 8 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 8 cores and 16 cpus
# 30605188 - done

## Checks
# 2 tests * 4 variant groups * 2 MAFs = 16 files per numvar directory
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_1 | wc -l # 16 
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2 | wc -l # 16 
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3 | wc -l # 16 
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4 | wc -l # 16 


##########################################################################################################################################
##########################################################################################################################################
##### 3. (COMBINED PS) META-ANALYZE: GNE + AMPXNIH (WGS) + UKB PD CASE-CONTROL META-ANALYSIS #############################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UKB PD CASE-CONTROL

# Make directories 
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_CONTROL_META
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_4


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
				-o $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}" >> GNE_AMP_NIH_UKB_CASE_CONTROL_META.swarm
    		done
    	done
    done
done

swarm -f GNE_AMP_NIH_UKB_CASE_CONTROL_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 64 commands run in 8 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 8 cores and 16 cpus
# 30649846 - done

## Checks
# 2 tests * 4 variant groups * 2 MAFs = 16 files per numvar directory
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_1 | wc -l # 16
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_2 | wc -l # 16
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_3 | wc -l # 16
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_4 | wc -l # 16

##########################################################################################################################################
##########################################################################################################################################
##### 4. (COMBINED PS) META-ANALYZE: GNE + AMPXNIH + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS  ####################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS 

# Make directories 
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_PROXIES_META
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_1
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_4

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
				-o $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021_${variants}_freqUpper${cutoff}" >> GNE_AMP_NIH_UKB_CASE_PROXIES_META.swarm
    		done
    	done
    done
done

swarm -f GNE_AMP_NIH_UKB_CASE_PROXIES_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 64 commands run in 8 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 8 cores and 16 cpus
# 30655270 - done

## Checks
# 2 tests * 4 variant groups * 2 MAFs = 16 files per numvar directory
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_1 | wc -l # 16
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_2 | wc -l # 16
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_3 | wc -l # 16
ls $WORK_DIR/RESULTS/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_4 | wc -l # 16


##########################################################################################################################################
##########################################################################################################################################
##### 5. SUMMARY OF RESULTS ###############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

mkdir $WORK_DIR/RESULTS/SUMMARIES 

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
					echo "COHORT: ${meta}; VARIANT GROUP: ${variants}; MAF: ${cutoff}; NUMVAR >= ${numvar} for ${test}" >> $WORK_DIR/RESULTS/SUMMARIES/combinedPs_top_25_lines_${meta}_minvar${numvar}.txt
					cat $WORK_DIR/RESULTS/${meta}/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}*${test}* | head -25 >> $WORK_DIR/RESULTS/SUMMARIES/combinedPs_top_25_lines_${meta}_minvar${numvar}.txt
					echo " " >> $WORK_DIR/RESULTS/SUMMARIES/combinedPs_top_25_lines_${meta}_minvar${numvar}.txt
				done
			done
		done
	done
done

# 3 cohorts * 4 minimum numvar = 12 summary files (each summary has each type of variant group, test, and MAF cut-off)

##########################################################################################################################################
##########################################################################################################################################
##### 6. (WEIGHTED Z) META-ANALYZE: GNE + AMPXNIH (WGS) + UKB PD_ALL META-ANALYSIS #######################################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UKB PD ALL

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
			echo "python $SCRIPTS/cmcwald_gne_AMPNIH_UKB_ALL_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
			--ukb_all_cases $WORK_DIR/UKB_EXOM_ALL_PD/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021 \
			-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> cmc_GNE_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm
		done
	done
done


swarm -f cmc_GNE_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 30741705 - done

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_1 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4 | wc -l # 8

##########################################################################################################################################
##########################################################################################################################################
##### 7. (WEIGHTED Z) META-ANALYZE: GNE + AMPXNIH (WGS) + UKB PD CASE-CONTROL META-ANALYSIS #############################################
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
			echo "python $SCRIPTS/cmcwald_gne_AMPNIH_UKB_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
			--ukb_cases $WORK_DIR/UKB_EXOM_PD_CASE/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--gne $WORK_DIR/GENENTECH/gne.${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021 \
			-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}" >> cmc_GNE_AMP_NIH_UKB_CASE_CONTROL_META.swarm
		done
	done
done

swarm -f cmc_GNE_AMP_NIH_UKB_CASE_CONTROL_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 30744152 - running

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_1 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_CONTROL_META/minimum_numvar_4 | wc -l # 8


##########################################################################################################################################
##########################################################################################################################################
##### 8. (WEIGHTED Z) META-ANALYZE: GNE + AMPXNIH + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS  ####################################
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
			echo "python $SCRIPTS/cmcwald_gne_AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py \
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
			-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_${numvar}/meta_GNE_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_PROXIES_2021_${variants}_freqUpper${cutoff}" >> cmc_GNE_AMP_NIH_UKB_CASE_PROXIES_META.swarm
		done
	done
done

swarm -f cmc_GNE_AMP_NIH_UKB_CASE_PROXIES_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 32 commands run in 4 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 4 cores and 8 cpus
# 30770156

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_1 | wc -l # 
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_2 | wc -l # 
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_3 | wc -l # 
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/GNE_AMP_NIH_UKB_CASE_PROXIES_META/minimum_numvar_4 | wc -l # 


##########################################################################################################################################
##########################################################################################################################################
##### 9. SUMMARY OF RESULTS ###############################################################################################################
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

# 3 cohorts * 4 minimum numvar = 12 summary files (each summary has each type of variant group and MAF cut-off)



