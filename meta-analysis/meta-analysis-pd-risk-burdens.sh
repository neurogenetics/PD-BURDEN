# !/bin/env bash

# PD Burdens - PD Risk Meta-analysis
    # Sept 2021
    # Mary B. Makarious, Julie Lake, and Cornelis Blauwendraat (LNG/NIA/NIH)

# Project Summary
    # Running a meta-analysis on the outputs of RVtests for UKB and AMPxNIH datasets 

# Workflow 
	# 0. Summary and Notes
	# 1. Getting Started 
	# 2. Meta-analyze: AMPxNIH (WGS) + UKB PD_ALL meta-analysis
	# 3. Meta-analyze: AMPxNIH (WGS) + UKB PD case-control meta-analysis
	# 4. Meta-analyze: AMPxNIH + UK PD case-control + UKB siblings + UKB parents
	# 5. Summary of Results 


##########################################################################################################################################
##########################################################################################################################################
##### 0. SUMMARY AND NOTES ###############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

# final UKB:
	# /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/meta_prep_annovar
	# /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/meta_prep_snpeff_loftee

# final AMPxNIH:
	# /data/CARD/PD/AMP_NIH/no_relateds/burden_annovar
	# /data/CARD/PD/AMP_NIH/no_relateds/burden_snpEff_loftee

# Meta-analysis for PD risk 
	# Per ANNOVAR and SnpEff/LOFTEE group
	# Per MAF and per variant group 
	# For SkatO and CMC separately 
		# WGS + UKB PD (@ AMPNIH_UKB_meta_analysis_SkatO.py and AMPNIH_UKB_meta_analysis_CMC.py) --ukb_cases, --amp_nih, -o 
		# WGS + UKB PD + PROXY1 + PROXY 2 (@ AMPNIH_UKB_CASES_PROXIES_meta_analysis_SkatO.py and AMPNIH_UKB_CASES_PROXIES_meta_analysis_CMC.py )
		# WGS + UKB PD_ALL (@ AMPNIH_UKB_ALL_meta_analysis_SkatO.py and AMPNIH_UKB_ALL_meta_analysis_CMC.py )

# Structure 
WORK_DIR="/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis"

$WORK_DIR
├── AMP_NIH_UKB_ALL_PD_PHENO_META/
│			├── ANNOVAR/
│			└── SNPEFF_LOFTEE/
├── AMP_NIH_UKB_CASE_CONTROL_META/
│			├── ANNOVAR/
│			└── SNPEFF_LOFTEE/
├── AMP_NIH_UKB_CASE_PROXIES_META/
│			├── ANNOVAR/
│			└── SNPEFF_LOFTEE/
└── etc... (Scripts and Swarm Files)

# AMP_NIH_UKB_ALL_PD_PHENO_META/
	# AMPxNIH (WGS) + UKB PD_ALL meta-analysis
	# Used AMPNIH_UKB_ALL_meta_analysis_SkatO.py and AMPNIH_UKB_ALL_meta_analysis_CMC.py

# AMP_NIH_UKB_CASE_CONTROL_META/
	# AMPxNIH (WGS) + UKB PD case-control meta-analysis
	# Used AMPNIH_UKB_meta_analysis_SkatO.py and AMPNIH_UKB_meta_analysis_CMC.py

# AMP_NIH_UKB_CASE_PROXIES_META/
	# AMPxNIH + UK PD case-control + UKB siblings + UKB parents
	# Used AMPNIH_UKB_CASES_PROXIES_meta_analysis_SkatO.py and AMPNIH_UKB_CASES_PROXIES_meta_analysis_CMC.py

# 48 files per subdirectory expected
	# 2 tests * 6 variant groups * 4 MAFs = 48 files per annotation type (ANNOVAR or SnpEff/LOFTEE)

##########################################################################################################################################
##########################################################################################################################################
##### 1. GETTING STARTED #################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

# Setting up global variables 
WORK_DIR="/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis"
UKB_ANNOVAR="/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/meta_prep_annovar"
UKB_SNPEFF="/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/meta_prep_snpeff_loftee"
AMPNIH_ANNOVAR="/data/CARD/PD/AMP_NIH/no_relateds/burden_annovar"
AMPNIH_SNPEFF="/data/CARD/PD/AMP_NIH/no_relateds/burden_snpEff_loftee"

# Loading the necessary modules 
module load python 

# Creating the working directory 
cd /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis
mkdir meta_risk_analysis
cd meta_risk_analysis

# Summarizing the file names etc to properly loop through later 
# AMP file 
	# AMP_NIH_noRelateds

# UKB disease groups
	# UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021
	# UKB_EXOM_PD_CASE_CONTROL_2021
	# UKB_EXOM_PD_PARENT_CONTROL_2021
	# UKB_EXOM_PD_SIBLING_CONTROL_2021

# Annovar variant groups 
	# ALL_MISSENSE
	# ALL_LOF
	# ALL_CADD_20
	# ALL_CADD_10
	# ALL_MISSENSE_and_LOF
	# ALL_CADD_20_and_LOF

# SnpEff/LOFTEE variant groups 
	# ALL_MISSENSE_and_LOF_SNPEFF
	# ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE
	# ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE
	# ALL_HIGH_IMPACT_and_LOF_SNPEFF
	# ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE
	# ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE

# Frequency groups
	# 0.001
	# 0.005
	# 0.01
	# 0.05

# Format: 
	# {disease_group}_{variant_group}.freqUpper{MAF}.{test_type}.assoc

# Setting up bash arrays to loop through 
test_type=(
	SkatO
	CMC)

MAF=(
	0.001
	0.005
	0.01
	0.05)

annovar_variant_group=(
	ALL_MISSENSE 
	ALL_LOF 
	ALL_CADD_20 
	ALL_CADD_10 
	ALL_MISSENSE_and_LOF 
	ALL_CADD_20_and_LOF)

snpeff_variant_group=(
	ALL_MISSENSE_and_LOF_SNPEFF
	ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE
	ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE
	ALL_HIGH_IMPACT_and_LOF_SNPEFF
	ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE
	ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE)


##########################################################################################################################################
##########################################################################################################################################
##### 2. META-ANALYZE: AMPXNIH (WGS) + UKB PD_ALL META-ANALYSIS ##########################################################################
##########################################################################################################################################
##########################################################################################################################################


## WGS + UKB PD CASE-CONTROL
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE

### ANNOVAR
# 2 tests * 6 variant groups * 4 MAFs = 48 files 

for test in "${test_type[@]}"
do
    for variants in "${annovar_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		echo "python $WORK_DIR/AMPNIH_UKB_meta_analysis_${test}.py \
			--ukb_cases $UKB_ANNOVAR/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--amp_nih $AMPNIH_ANNOVAR/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
			-o $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}" >> $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META.annovar.swarm
    	done
    done
done

swarm -f $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META.annovar.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 4
# 48 commands run in 12 subjobs, each command requiring 1.5 gb and 1 thread, running 4 processes serially per subjob, allocating 12 cores and 24 cpus
# 22982137


### SNPEFF/LOFTEE
# 2 tests * 6 variant groups * 4 MAFs = 48 files 
for test in "${test_type[@]}"
do
    for variants in "${snpeff_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		echo "python $WORK_DIR/AMPNIH_UKB_meta_analysis_${test}.py \
			--ukb_cases $UKB_SNPEFF/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--amp_nih $AMPNIH_SNPEFF/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
			-o $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}" >> $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META.loftee.swarm
    	done
    done
done

swarm -f $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META.loftee.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 4
# 48 commands run in 12 subjobs, each command requiring 1.5 gb and 1 thread, running 4 processes serially per subjob, allocating 12 cores and 24 cpus
# 22982146

ls $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/ | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/ | wc -l #48 



##########################################################################################################################################
##########################################################################################################################################
##### 3. META-ANALYZE: AMPXNIH (WGS) + UKB PD CASE-CONTROL META-ANALYSIS #################################################################
##########################################################################################################################################
##########################################################################################################################################


## WGS + UKB PD + PROXY1 + PROXY2 
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE

# ANNOVAR
for test in "${test_type[@]}"
do
    for variants in "${annovar_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		echo "python $WORK_DIR/AMPNIH_UKB_CASES_PROXIES_meta_analysis_${test}.py \
			--ukb_cases $UKB_ANNOVAR/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--amp_nih $AMPNIH_ANNOVAR/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
			--ukb_sib $UKB_ANNOVAR/UKB_EXOM_PD_SIBLING_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--ukb_parent $UKB_ANNOVAR/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			-o $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META.annovar.swarm
    	done
    done
done

swarm -f $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META.annovar.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 4
# 48 commands run in 12 subjobs, each command requiring 1.5 gb and 1 thread, running 4 processes serially per subjob, allocating 12 cores and 24 cpus
# 22982431

# SNPEFF/LOFTEE
for test in "${test_type[@]}"
do
    for variants in "${snpeff_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		echo "python $WORK_DIR/AMPNIH_UKB_CASES_PROXIES_meta_analysis_${test}.py \
			--ukb_cases $UKB_SNPEFF/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--amp_nih $AMPNIH_SNPEFF/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
			--ukb_sib $UKB_SNPEFF/UKB_EXOM_PD_SIBLING_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--ukb_parent $UKB_SNPEFF/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			-o $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META.loftee.swarm
    	done
    done
done

swarm -f $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META.loftee.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 4
# 48 commands run in 12 subjobs, each command requiring 1.5 gb and 1 thread, running 4 processes serially per subjob, allocating 12 cores and 24 cpus
# 22982515

ls $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/ | wc -l 
ls $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/ | wc -l 

##########################################################################################################################################
##########################################################################################################################################
##### 4. META-ANALYZE: AMPXNIH + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS #########################################################
##########################################################################################################################################
##########################################################################################################################################


## WGS + UKB PD ALL 
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE

# ANNOVAR
for test in "${test_type[@]}"
do
    for variants in "${annovar_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		echo "python AMPNIH_UKB_ALL_meta_analysis_${test}.py \
			--ukb_all_cases $UKB_ANNOVAR/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--amp_nih $AMPNIH_ANNOVAR/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
			-o $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> AMP_NIH_UKB_ALL_PD_PHENO_META_annovar.swarm
    	done
    done
done

swarm -f AMP_NIH_UKB_ALL_PD_PHENO_META_annovar.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 4
# 48 commands run in 12 subjobs, each command requiring 1.5 gb and 1 thread, running 4 processes serially per subjob, allocating 12 cores and 24 cpus
# 22981359

# SNPEFF/LOFTEE
for test in "${test_type[@]}"
do
    for variants in "${snpeff_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		echo "python AMPNIH_UKB_ALL_meta_analysis_${test}.py \
			--ukb_all_cases $UKB_SNPEFF/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--amp_nih $AMPNIH_SNPEFF/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
			-o $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> AMP_NIH_UKB_ALL_PD_PHENO_META_loftee.swarm
    	done
    done
done

swarm -f AMP_NIH_UKB_ALL_PD_PHENO_META_loftee.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 4
# 48 commands run in 12 subjobs, each command requiring 1.5 gb and 1 thread, running 4 processes serially per subjob, allocating 12 cores and 24 cpus
# 22981369

ls $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/ | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/ | wc -l #48

##########################################################################################################################################
##########################################################################################################################################
##### 5. SUMMARY OF RESULTS ##############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

