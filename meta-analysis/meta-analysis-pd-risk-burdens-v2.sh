# !/bin/env bash

# PD Burdens - PD Risk Meta-analysis
    # Sept 2021
    # Mary B. Makarious, Julie Lake, and Cornelis Blauwendraat (LNG/NIA/NIH)

# Project Summary
    # Running a meta-analysis on the outputs of RVtests for UKB and AMPxNIH datasets 

# Workflow 
	# 0. Summary and Notes
	# 1. Getting Started 
	# 2. Meta-analyze: AMPxNIH (WGS) + UKB PD case-control meta-analysis
	# 3. Meta-analyze: AMPxNIH + UK PD case-control + UKB siblings + UKB parents
	# 4. Meta-analyze: AMPxNIH (WGS) + UKB PD_ALL meta-analysis 
	# 5. Summary of Results 
	# 6. LRRK2 Conditional Burden: Getting Started 
	# 7. LRRK2 Conditional Burden: AMPxNIH (WGS) + UKB PD case-control meta-analysis
	# 8. LRRK2 Conditional Burden: Meta-analyze: AMPxNIH + UK PD case-control + UKB siblings + UKB parents
	# 9. LRRK2 Conditional Burden: Meta-analyze: AMPxNIH (WGS) + UKB PD_ALL meta-analysis 
	# 10. LRRK2 Conditional Burden: Summary of Results 


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
		# WGS + UKB PD 
		# WGS + UKB PD + PROXY1 + PROXY 2 
		# WGS + UKB PD_ALL 

# Structure 
WORK_DIR="/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis"

# $WORK_DIR
# ├── AMP_NIH_UKB_ALL_PD_PHENO_META/
# │			├── ANNOVAR/
# │			│	├── minimum_numvar_1/
# │			│	├── minimum_numvar_2/
# │			│	├── minimum_numvar_3/
# │			│	└── minimum_numvar_4/
# │			└── SNPEFF_LOFTEE/
# │				├── minimum_numvar_1/
# │				├── minimum_numvar_2/
# │				├── minimum_numvar_3/
# │				└── minimum_numvar_4/
# ├── AMP_NIH_UKB_CASE_CONTROL_META/
# │			├── ANNOVAR/
# │			│	├── minimum_numvar_1/
# │			│	├── minimum_numvar_2/
# │			│	├── minimum_numvar_3/
# │			│	└── minimum_numvar_4/
# │			└── SNPEFF_LOFTEE/
# │				├── minimum_numvar_1/
# │				├── minimum_numvar_2/
# │				├── minimum_numvar_3/
# │				└── minimum_numvar_4/
# ├── AMP_NIH_UKB_CASE_PROXIES_META/
# │			├── ANNOVAR/
# │			│	├── minimum_numvar_1/
# │			│	├── minimum_numvar_2/
# │			│	├── minimum_numvar_3/
# │			│	└── minimum_numvar_4/
# │			└── SNPEFF_LOFTEE/
# │				├── minimum_numvar_1/
# │				├── minimum_numvar_2/
# │				├── minimum_numvar_3/
# │				└── minimum_numvar_4/
# ├── *.swarm files 
# ├── swarm_logs/
# └── scripts/

# AMP_NIH_UKB_ALL_PD_PHENO_META/
	# AMPxNIH (WGS) + UKB PD_ALL meta-analysis
	# Used AMPNIH_UKB_ALL_meta_analysis_numvar.py

# AMP_NIH_UKB_CASE_CONTROL_META/
	# AMPxNIH (WGS) + UKB PD case-control meta-analysis
	# Used AMPNIH_UKB_meta_analysis_numvar.py

# AMP_NIH_UKB_CASE_PROXIES_META/
	# AMPxNIH + UK PD case-control + UKB siblings + UKB parents
	# Used AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py 

# 48 files per ANNOVAR subdirectory expected
	# 2 tests * 6 variant groups * 4 MAFs = 48 files per numvariant cut-off 
# 96 files per SnpEff/LOFTEE subdirectory expected
	# 2 tests * 12 variant groups * 4 MAFs = 96 files per numvar cut-off 

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
SCRIPTS="/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/scripts"

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

# Number of variants cut-offs
	# 1
	# 2
	# 3 
	# 4

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
	ALL_MODERATE_HIGH_IMPACT_SNPEFF
	ALL_HIGH_IMPACT_SNPEFF
	ALL_MISSENSE_and_LOF_SNPEFF
	ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE
	ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE
	ALL_MISSENSE_SNPEFF
	ALL_LOF_SNPEFF
	ALL_LOF_HC_LOFTEE
	ALL_LOF_and_HC_LOFTEE
	ALL_HIGH_IMPACT_and_LOF_SNPEFF
	ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE
	ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE)

number_variants=(
	1
	2
	3
	4)


##########################################################################################################################################
##########################################################################################################################################
##### 2. META-ANALYZE: AMPXNIH (WGS) + UKB PD CASE-CONTROL META-ANALYSIS #################################################################
##########################################################################################################################################
##########################################################################################################################################

## WGS + UKB PD CASE-CONTROL

# Make directories 
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_1
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_2
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_3
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_4

mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_1
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_2
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_3
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_4

### ANNOVAR
# 2 tests * 6 variant groups * 4 MAFs = 48 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${annovar_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/AMPNIH_UKB_meta_analysis_numvar.py \
	    		--test ${test} \
				--ukb_cases $UKB_ANNOVAR/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--amp_nih $AMPNIH_ANNOVAR/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021 \
				--numvar ${numvar} \
				-o $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}" >> $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META.annovar.swarm
    		done
    	done
    done
done


swarm -f $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META.annovar.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 192 commands run in 24 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 24 cores and 48 cpus
# 23478404

### SNPEFF/LOFTEE
# 2 tests * 12 variant groups * 4 MAFs = 96 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${snpeff_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/AMPNIH_UKB_meta_analysis_numvar.py \
	    		--test ${test} \
				--ukb_cases $UKB_SNPEFF/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--amp_nih $AMPNIH_SNPEFF/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021 \
				--numvar ${numvar} \
				-o $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}" >> $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META.loftee.swarm
    		done
    	done
    done
done

swarm -f $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META.loftee.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 384 commands run in 48 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 48 cores and 96 cpus
# 23566179

## Checks 
ls $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_1 | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_2 | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_3 | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_4 | wc -l #48

ls $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_1 | wc -l #96
ls $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_2 | wc -l #96
ls $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_3 | wc -l #96
ls $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_4 | wc -l #96 

##########################################################################################################################################
##########################################################################################################################################
##### 3. META-ANALYZE: AMPXNIH + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS #########################################################
##########################################################################################################################################
##########################################################################################################################################

## WGS + UKB PD + PROXY1 + PROXY2 

# Make directories 
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_1
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_2
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_3
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_4

mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_1
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_2
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_3
mkdir $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_4

### ANNOVAR
# 2 tests * 6 variant groups * 4 MAFs = 48 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${annovar_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
	    		echo "python $SCRIPTS/AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py \
	    		--test ${test} \
	    		--numvar ${numvar} \
				--ukb_cases $UKB_ANNOVAR/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--amp_nih $AMPNIH_ANNOVAR/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
				--ukb_sib $UKB_ANNOVAR/UKB_EXOM_PD_SIBLING_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--ukb_parent $UKB_ANNOVAR/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021 \
				-o $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META.annovar.swarm
    		done
    	done
    done
done

swarm -f $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META.annovar.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 192 commands run in 24 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 24 cores and 48 cpus
# 23479243

### SNPEFF/LOFTEE
# 2 tests * 12 variant groups * 4 MAFs = 96 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${snpeff_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
	    		echo "python $SCRIPTS/AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py \
	    		--test ${test} \
	    		--numvar ${numvar} \
				--ukb_cases $UKB_SNPEFF/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--amp_nih $AMPNIH_SNPEFF/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
				--ukb_sib $UKB_SNPEFF/UKB_EXOM_PD_SIBLING_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--ukb_parent $UKB_SNPEFF/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021 \
				-o $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META.loftee.swarm
    		done
    	done
    done
done

swarm -f $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META.loftee.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 384 commands run in 48 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 48 cores and 96 cpus
# 23566278

## Checks 
ls $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_1 | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_2 | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_3 | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_4 | wc -l #48

ls $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_1 | wc -l #96
ls $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_2 | wc -l #96
ls $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_3 | wc -l #96
ls $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_4 | wc -l #96


##########################################################################################################################################
##########################################################################################################################################
##### 4. META-ANALYZE: AMPXNIH (WGS) + UKB PD_ALL META-ANALYSIS ##########################################################################
##########################################################################################################################################
##########################################################################################################################################

## WGS + UKB PD ALL

# Make directories 
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_1
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_2
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_3
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_4

mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_1
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_2
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_3
mkdir $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_4

### ANNOVAR
# 2 tests * 6 variant groups * 4 MAFs = 48 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${annovar_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/AMPNIH_UKB_ALL_meta_analysis_numvar.py \
    			--test ${test} \
	    		--numvar ${numvar} \
				--ukb_all_cases $UKB_ANNOVAR/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--amp_nih $AMPNIH_ANNOVAR/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021 \
				-o $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> AMP_NIH_UKB_ALL_PD_PHENO_META_annovar.swarm
    		done
    	done
    done
done

swarm -f $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META_annovar.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 192 commands run in 24 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 24 cores and 48 cpus
# 23477432

### SNPEFF/LOFTEE
# 2 tests * 12 variant groups * 4 MAFs = 96 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${snpeff_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
	    	echo "python $SCRIPTS/AMPNIH_UKB_ALL_meta_analysis_numvar.py \
	    	--test ${test} \
	    	--numvar ${numvar} \
			--ukb_all_cases $UKB_SNPEFF/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
			--amp_nih $AMPNIH_SNPEFF/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021 \
			-o $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> AMP_NIH_UKB_ALL_PD_PHENO_META_loftee.swarm
    		done
    	done
    done
done

swarm -f $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META_loftee.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 384 commands run in 48 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 48 cores and 96 cpus
# 23566348

## Checks 
ls $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_1 | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_2 | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_3 | wc -l #48
ls $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_4 | wc -l #48

ls $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_1 | wc -l #96
ls $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_2 | wc -l #96
ls $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_3 | wc -l #96
ls $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_4 | wc -l #96 

## Quick clean-up
mv swarm_* swarm_logs/

##########################################################################################################################################
##########################################################################################################################################
##### 5. SUMMARY OF RESULTS ##############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## Getting unique gene IDs and their smallest P-value from top 10 lines 

# AMP_NIH_UKB_CASE_CONTROL_META/
	# AMPxNIH (WGS) + UKB PD case-control meta-analysis
	# Used AMPNIH_UKB_meta_analysis_numvar.py

cd $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_1
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1.txt
done

cd $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_1
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1.txt
done


## Min var 2
cd $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_2
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2.txt
done

cd $WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_2
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2.txt
done


# AMP_NIH_UKB_ALL_PD_PHENO_META/
	# AMPxNIH (WGS) + UKB PD_ALL meta-analysis
	# Used AMPNIH_UKB_ALL_meta_analysis_numvar.py

cd $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_1
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1.txt
done

cd $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_1
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1.txt
done

## min var 2
cd $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_2
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2.txt
done

cd $WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_2
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2.txt
done


# AMP_NIH_UKB_CASE_PROXIES_META/
	# AMPxNIH + UK PD case-control + UKB siblings + UKB parents
	# Used AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py 

cd $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_1
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1.txt
done

cd $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_1
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1.txt
done

## min var 2
cd $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_2
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2.txt
done

cd $WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_2
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2.txt
	cat *${variants}_freqUpper*  | head -1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2.txt
	echo " " >> $WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2.txt
done


##########################################################################################################################################
##########################################################################################################################################
##### 6. LRRK2 CONDITIONAL BURDEN: GETTING STARTED #######################################################################################
##########################################################################################################################################
##########################################################################################################################################

# final UKB - LRRK2 Conditional RVTests Results:
	# /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/LRRK2_conditional_burden/burden_annovar
	# /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/LRRK2_conditional_burden/burden_snpeff_loftee

# final AMPxNIH - LRRK2 Conditional RVTests Results:
	# /data/CARD/PD/AMP_NIH/no_relateds/LRRK2_conditional_burden/burden_annovar
	# /data/CARD/PD/AMP_NIH/no_relateds/LRRK2_conditional_burden/burden_snpeff_loftee

# Meta-analysis for PD risk 
	# Per ANNOVAR and SnpEff/LOFTEE group
	# Per MAF and per variant group 
	# For SkatO and CMC separately 
		# WGS + UKB PD 
		# WGS + UKB PD + PROXY1 + PROXY 2 
		# WGS + UKB PD_ALL 

# Setting up global variables 
LRRK2_WORK_DIR="/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/LRRK2_CONDITIONAL_META"
LRRK2_UKB_ANNOVAR="/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/LRRK2_conditional_burden/burden_annovar"
LRRK2_UKB_SNPEFF="/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/LRRK2_conditional_burden/burden_snpeff_loftee"
LRRK2_AMPNIH_ANNOVAR="/data/CARD/PD/AMP_NIH/no_relateds/LRRK2_conditional_burden/burden_annovar"
LRRK2_AMPNIH_SNPEFF="/data/CARD/PD/AMP_NIH/no_relateds/LRRK2_conditional_burden/burden_snpeff_loftee"
SCRIPTS="/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/scripts"

# Loading the necessary modules 
module load python 

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
	ALL_MODERATE_HIGH_IMPACT_SNPEFF
	ALL_HIGH_IMPACT_SNPEFF
	ALL_MISSENSE_and_LOF_SNPEFF
	ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE
	ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE
	ALL_MISSENSE_SNPEFF
	ALL_LOF_SNPEFF
	ALL_LOF_HC_LOFTEE
	ALL_LOF_and_HC_LOFTEE
	ALL_HIGH_IMPACT_and_LOF_SNPEFF
	ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE
	ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE)

number_variants=(
	1
	2
	3
	4)

# Setting up directories 
mkdir $LRRK2_WORK_DIR

##########################################################################################################################################
##########################################################################################################################################
##### 7. LRRK2 CONDITIONAL BURDEN: AMPXNIH (WGS) + UKB PD CASE-CONTROL META-ANALYSIS #####################################################
##########################################################################################################################################
##########################################################################################################################################

## WGS + UKB PD CASE-CONTROL

# Make directories 
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_1
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_2
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_3
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_4

mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_1
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_2
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_3
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_4


### ANNOVAR
# 2 tests * 6 variant groups * 4 MAFs = 48 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${annovar_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/AMPNIH_UKB_meta_analysis_numvar.py \
	    		--test ${test} \
				--ukb_cases $LRRK2_UKB_ANNOVAR/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.chr12_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--amp_nih $LRRK2_AMPNIH_ANNOVAR/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_LRRK2 \
				--numvar ${numvar} \
				-o $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}_LRRK2" >> $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META_LRRK2.annovar.swarm
    		done
    	done
    done
done


swarm -f $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META_LRRK2.annovar.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 192 commands run in 24 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 24 cores and 48 cpus
# 23575926

### SNPEFF/LOFTEE
# 2 tests * 12 variant groups * 4 MAFs = 96 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${snpeff_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/AMPNIH_UKB_meta_analysis_numvar.py \
	    		--test ${test} \
				--ukb_cases $LRRK2_UKB_SNPEFF/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.chr12_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--amp_nih $LRRK2_AMPNIH_SNPEFF/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_LRRK2 \
				--numvar ${numvar} \
				-o $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021_${variants}_freqUpper${cutoff}_LRRK2" >> $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META_LRRK2.loftee.swarm
    		done
    	done
    done
done

swarm -f $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META_LRRK2.loftee.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 384 commands run in 48 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 48 cores and 96 cpus
# 23575933

## Checks 
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_1 | wc -l #48
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_2 | wc -l #48
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_3 | wc -l #48
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_4 | wc -l #48

ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_1 | wc -l #96
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_2 | wc -l #96
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_3 | wc -l #96
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_4 | wc -l #96


##########################################################################################################################################
##########################################################################################################################################
###### 8. LRRK2 CONDITIONAL BURDEN: META-ANALYZE: AMPXNIH + UK PD CASE-CONTROL + UKB SIBLINGS + UKB PARENTS ##############################
##########################################################################################################################################
##########################################################################################################################################

## WGS + UKB PD + PROXY1 + PROXY2 

# Make directories 
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_1
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_2
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_3
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_4

mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_1
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_2
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_3
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_4

### ANNOVAR
# 2 tests * 6 variant groups * 4 MAFs = 48 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${annovar_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
	    		echo "python $SCRIPTS/AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py \
	    		--test ${test} \
	    		--numvar ${numvar} \
				--ukb_cases $LRRK2_UKB_ANNOVAR/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.chr12_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--amp_nih $LRRK2_AMPNIH_ANNOVAR/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--ukb_sib $LRRK2_UKB_ANNOVAR/UKB_EXOM_PD_SIBLING_CONTROL_2021_${variants}.chr12_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--ukb_parent $LRRK2_UKB_ANNOVAR/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.chr12_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021_LRRK2 \
				-o $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021_${variants}_freqUpper${cutoff}_LRRK2" >> $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META_LRRK2.annovar.swarm
    		done
    	done
    done
done

swarm -f $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META_LRRK2.annovar.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 192 commands run in 24 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 24 cores and 48 cpus
# 23576384

### SNPEFF/LOFTEE
# 2 tests * 12 variant groups * 4 MAFs = 96 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${snpeff_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
	    		echo "python $SCRIPTS/AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py \
	    		--test ${test} \
	    		--numvar ${numvar} \
				--ukb_cases $LRRK2_UKB_SNPEFF/UKB_EXOM_PD_CASE_CONTROL_2021_${variants}.chr12_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--amp_nih $LRRK2_AMPNIH_SNPEFF/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--ukb_sib $LRRK2_UKB_SNPEFF/UKB_EXOM_PD_SIBLING_CONTROL_2021_${variants}.chr12_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--ukb_parent $LRRK2_UKB_SNPEFF/UKB_EXOM_PD_PARENT_CONTROL_2021_${variants}.chr12_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021_LRRK2 \
				-o $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021_${variants}_freqUpper${cutoff}_LRRK2" >> $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META_LRRK2.loftee.swarm
    		done
    	done
    done
done

swarm -f $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META_LRRK2.loftee.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 384 commands run in 48 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 48 cores and 96 cpus
# 23576798

## Checks 
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_1 | wc -l #48
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_2 | wc -l #48
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_3 | wc -l #48
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_4 | wc -l #48

ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_1 | wc -l #96
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_2 | wc -l #96
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_3 | wc -l #96
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_4 | wc -l #96

 
##########################################################################################################################################
##########################################################################################################################################
##### 9. LRRK2 CONDITIONAL BURDEN: META-ANALYZE: AMPXNIH (WGS) + UKB PD_ALL META-ANALYSIS  ###############################################
##########################################################################################################################################
##########################################################################################################################################

## WGS + UKB PD ALL

# Make directories 
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_1
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_2
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_3
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_4

mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_1
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_2
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_3
mkdir $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_4

### ANNOVAR
# 2 tests * 6 variant groups * 4 MAFs = 48 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${annovar_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/AMPNIH_UKB_ALL_meta_analysis_numvar.py \
    			--test ${test} \
	    		--numvar ${numvar} \
				--ukb_all_cases $LRRK2_UKB_ANNOVAR/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.chr12_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--amp_nih $LRRK2_AMPNIH_ANNOVAR/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}_LRRK2.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_LRRK2 \
				-o $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}_LRRK2" >> AMP_NIH_UKB_ALL_PD_PHENO_META_LRRK2.annovar.swarm
    		done
    	done
    done
done

swarm -f $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META_LRRK2.annovar.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 192 commands run in 24 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 24 cores and 48 cpus
# 23581273

### SNPEFF/LOFTEE
# 2 tests * 12 variant groups * 4 MAFs = 96 files per numvar directory

for test in "${test_type[@]}"
do
    for variants in "${snpeff_variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
	    	echo "python $SCRIPTS/AMPNIH_UKB_ALL_meta_analysis_numvar.py \
	    	--test ${test} \
	    	--numvar ${numvar} \
			--ukb_all_cases $LRRK2_UKB_SNPEFF/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.chr12_freqUpper${cutoff}_LRRK2.${test}.assoc \
			--amp_nih $LRRK2_AMPNIH_SNPEFF/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}_LRRK2.${test}.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_LRRK2 \
			-o $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}_LRRK2" >> AMP_NIH_UKB_ALL_PD_PHENO_META_LRRK2.loftee.swarm
    		done
    	done
    done
done

swarm -f $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META_LRRK2.loftee.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 384 commands run in 48 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 48 cores and 96 cpus
# 23581597

## Checks 
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_1 | wc -l #48
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_2 | wc -l #48
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_3 | wc -l #48
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_4 | wc -l #48

ls $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_1 | wc -l #96
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_2 | wc -l #96
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_3 | wc -l #96
ls $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_4 | wc -l #96

## Quick clean-up
mv swarm_* swarm_logs/

##########################################################################################################################################
##########################################################################################################################################
##### 10. LRRK2 CONDITIONAL BURDEN: SUMMARY OF RESULTS ###################################################################################
##########################################################################################################################################
##########################################################################################################################################

## Getting unique gene IDs and their smallest P-value from top 10 lines 

# AMP_NIH_UKB_CASE_CONTROL_META/
	# AMPxNIH (WGS) + UKB PD case-control meta-analysis
	# Used AMPNIH_UKB_meta_analysis_numvar.py

cd $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_1
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1_LRRK2.txt
done

cd $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_1
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar1_LRRK2.txt
done


## Min var 2
cd $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/ANNOVAR/minimum_numvar_2
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2_LRRK2.txt
done

cd $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_CONTROL_META/SNPEFF_LOFTEE/minimum_numvar_2
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_CONTROL_META_minvar2_LRRK2.txt
done


# AMP_NIH_UKB_ALL_PD_PHENO_META/
	# AMPxNIH (WGS) + UKB PD_ALL meta-analysis
	# Used AMPNIH_UKB_ALL_meta_analysis_numvar.py

cd $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_1
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1_LRRK2.txt
done

cd $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_1
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar1_LRRK2.txt
done

## min var 2
cd $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/ANNOVAR/minimum_numvar_2
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2_LRRK2.txt
done

cd $LRRK2_WORK_DIR/AMP_NIH_UKB_ALL_PD_PHENO_META/SNPEFF_LOFTEE/minimum_numvar_2
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_ALL_PD_PHENO_META_minvar2_LRRK2.txt
done


# AMP_NIH_UKB_CASE_PROXIES_META/
	# AMPxNIH + UK PD case-control + UKB siblings + UKB parents
	# Used AMPNIH_UKB_CASES_PROXIES_meta_analysis_numvar.py 

cd $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_1
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1_LRRK2.txt
done

cd $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_1
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar1_LRRK2.txt
done

## min var 2
cd $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/ANNOVAR/minimum_numvar_2
for variants in "${annovar_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2_LRRK2.txt
done

cd $LRRK2_WORK_DIR/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_2
for variants in "${snpeff_variant_group[@]}"
do
	echo "Now looking at variant group ${variants}"
	echo "${variants}" >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2_LRRK2.txt
	cat *${variants}_freqUpper*  | sort -k2 -nr | head -10 | sort -u -k1,1 >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2_LRRK2.txt
	echo " " >> $LRRK2_WORK_DIR/top_unique_genes_top_10_lines_AMP_NIH_UKB_CASE_PROXIES_META_minvar2_LRRK2.txt
done

