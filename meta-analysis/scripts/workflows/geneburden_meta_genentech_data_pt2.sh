# PD Burdens Meta-analysis - Part 2 
	# Mary B Makarious
	# Meta-analyze by combining pval (AMPxNIH + UKB cohorts + Genentech)
	# Began FEB 2022 (Last updated 15-FEB-2022)

# README by Mary Makarious; Last Updated 15-FEB-2022
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

# Workflow 
	## GETTING STARTED
		# 0. Summary and Notes
		# 1. Getting Started 

	## META-ANALYZE BY COMBINING P-VALUES (SKAT-O AND CMC WALD)tree
		# 2. (Combined Ps) Meta-analyze: AMPxNIH (WGS) + UKB PD_ALL meta-analysis
			# ALL defined as case-control + proxies 
		# 3. Summary of Results 

	## META-ANALYZE BY SUMMED Z APPROACH (CMC WALD) 	
		# 4. (Summed Z) Meta-analyze: AMPxNIH (WGS) + UKB PD_ALL meta-analysis
			# ALL defined as case-control + both parent and sibling proxies 
		# 5. Summary of Results 

	## META-ANALYZE BY WEIGHTED Z APPROACH (CMC WALD) 	
		# 6. (Weighted Z) Meta-analyze: AMPxNIH (WGS) + UKB PD_ALL meta-analysis
			# ALL defined as case-control + proxies 
		# 7. Summary of Results 

	## CALCULATE LAMBDAS 
		# 8. (Combined Ps) Calculate lambdas per meta-analyzed cohort
		# 9. (Summed Z) Calculate lambdas per meta-analyzed cohort
		# 10. (Weighted Z) Calculate lambdas per meta-analyzed cohort

	## MANHATTAN PLOTS
		# 11. (Combined Ps) Manhattan Plots: AMPxNIH (WGS) + UKB PD_ALL meta-analysis
		# 12. (Summed Z) Manhattan Plots: AMPxNIH (WGS) + UKB PD_ALL meta-analysis
		# 13. (Weighted Z) Manhattan Plots: AMPxNIH (WGS) + UKB PD_ALL meta-analysis

	## SUMMARIZE + INVESTIGATE 

##########################################################################################################################################
##########################################################################################################################################
##### 0. SUMMARY AND NOTES ###############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

## Goal
	# - Use minvar 2 and above 
	# - Make sure p-value available in at least 2 datasets 
	# - Only focus on UKB all (cases-controls; sibling+parent proxies)
	# UKB x AMP meta using 1x10-6 as significant, check for gene P-values in Genentech
	# UKB x AMP x GNE as meta and check genes 1x10-6, check for gene P-values


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
	2
	3
	4)

gne_metas=(
	GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
	GNE_AMP_NIH_UKB_CASE_CONTROL_META
	GNE_AMP_NIH_UKB_CASE_PROXIES_META
	GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META)

metas=(
	AMP_NIH_UKB_ALL_PD_PHENO_META
	)

## Made 3 new scripts to not include Genentech
	# AMPNIH_UKB_ALL_meta_analysis_numvar.py
	# cmcwald_summedZ_AMPNIH_UKB_ALL_meta_analysis_numvar.py
	# cmcwald_weightedZ_AMPNIH_UKB_ALL_meta_analysis_numvar.py

##########################################################################################################################################
##########################################################################################################################################
##### 2. (COMBINED PS) META-ANALYZE: AMPXNIH (WGS) + UKB PD_ALL META-ANALYSIS ############################################################
##########################################################################################################################################
##########################################################################################################################################

# Make directories 
mkdir $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META
mkdir $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for test in "${test_type[@]}"
do
    for variants in "${variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "python $SCRIPTS/AMPNIH_UKB_ALL_meta_analysis_numvar.py \
    			--test ${test} \
	    		--numvar ${numvar} \
				--ukb_all_cases $WORK_DIR/UKB_EXOM_ALL_PD/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.${test}.assoc \
				--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.${test}.assoc \
				--variant_group ${variants} \
				--maf ${cutoff} \
				--group meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021 \
				-o $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> AMP_NIH_UKB_ALL_PD_PHENO_META.swarm
    		done
    	done
    done
done

swarm -f AMP_NIH_UKB_ALL_PD_PHENO_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 48 commands run in 6 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 6 cores and 12 cpus
# 32502760 - done 

## Checks
# 2 tests * 4 variant groups * 2 MAFs = 16 files per numvar directory
ls $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3 | wc -l # 16
ls $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4 | wc -l # 16 
# done

##########################################################################################################################################
##########################################################################################################################################
##### 3. SUMMARY OF RESULTS ##############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

mkdir $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_SUMMARIES 

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
					echo "COHORT: ${meta}; VARIANT GROUP: ${variants}; MAF: ${cutoff}; NUMVAR >= ${numvar} for ${test}" >> $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_SUMMARIES/combinedPs_top_25_lines_${meta}_minvar${numvar}.txt
					cat $WORK_DIR/RESULTS_COMBINED_P/${meta}/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}*${test}* | head -25 >> $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_SUMMARIES/combinedPs_top_25_lines_${meta}_minvar${numvar}.txt
					echo " " >> $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_SUMMARIES/combinedPs_top_25_lines_${meta}_minvar${numvar}.txt
				done
			done
		done
	done
done

# done 

##########################################################################################################################################
##########################################################################################################################################
##### 4. (SUMMED Z) META-ANALYZE: AMPXNIH (WGS) + UKB PD_ALL META-ANALYSIS ###############################################################
##########################################################################################################################################
##########################################################################################################################################

## GNE + AMPxNIH WGS + UKB PD ALL

# Make directories 
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
		do
			echo "python $SCRIPTS/cmcwald_summedZ_AMPNIH_UKB_ALL_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
			--ukb_all_cases $WORK_DIR/UKB_EXOM_ALL_PD/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021 \
			-o $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> cmc_summedZ_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm
		done
	done
done


swarm -f cmc_summedZ_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 24 commands run in 3 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 3 cores and 6 cpus
# 32503182 - done 

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3 | wc -l # 8
ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4 | wc -l # 8
# done 

##########################################################################################################################################
##########################################################################################################################################
##### 5. SUMMARY OF RESULTS ##############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_SUMMARIES/

for meta in "${metas[@]}"
do
	for variants in "${variant_group[@]}"
	do
		for cutoff in "${MAF[@]}"
		do
			for numvar in "${number_variants[@]}"
	    	do
				echo "Now looking at cohort: ${meta}; variant group ${variants} at MAF ${cutoff} at minimum numvar ${numvar}"
				echo "COHORT: ${meta}; VARIANT GROUP: ${variants}; MAF: ${cutoff}; NUMVAR >= ${numvar}" >> $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_SUMMARIES/cmcwald_summedZ_top_25_lines_${meta}_minvar${numvar}.txt
				cat $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/${meta}/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}* | head -25 >> $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_SUMMARIES/cmcwald_summedZ_top_25_lines_${meta}_minvar${numvar}.txt
				echo " " >> $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_SUMMARIES/cmcwald_summedZ_top_25_lines_${meta}_minvar${numvar}.txt
			done
		done
	done
done

# 1 cohorts * 3 minimum numvar = 3 summary files (each summary has each type of variant group and MAF cut-off)
# done 

##########################################################################################################################################
##########################################################################################################################################
##### 6. (WEIGHTED Z) META-ANALYZE: AMPXNIH (WGS) + UKB PD_ALL META-ANALYSIS #############################################################
##########################################################################################################################################
##########################################################################################################################################

	## Cohort Numbers 

	# |                 Cohort                 	| Cases 	| Controls 	| TOTAL 	|
	# |:--------------------------------------:	|:-----:	|:--------:	|:-----:	|
	# |                                AMPxNIH 	|  3376 	|     4610 	|  7986 	|
	# |                                UKB ALL 	|  7806 	|    38051 	| 45857 	|
	# |                       UKB case-control 	|  1105 	|     5643 	|  6748 	|
	# |                           UKB siblings 	|   668 	|     3463 	|  4131 	|
	# |                            UKB parents 	|  6033 	|    28945 	| 34978 	|
	# |                                    GNE 	|  2700 	|     9000 	| 11700 	|
	# |                                        	|       	|          	|       	|
	# |      GNE_AMP_NIH_UKB_ALL_PD_PHENO_META 	| 13882 	|    51661 	| 65543 	|
	# |      GNE_AMP_NIH_UKB_CASE_CONTROL_META 	|  7181 	|    19253 	| 26434 	|
	# |      GNE_AMP_NIH_UKB_CASE_PROXIES_META 	| 13882 	|    51661 	| 65543 	|
	# | GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META 	| 13214 	|    48198 	| 61412 	|
	# |          AMP_NIH_UKB_ALL_PD_PHENO_META 	| 11182 	|    42661 	| 53843 	|


# Make directories 
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4

cd $WORK_DIR/SWARMS

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
		do
			echo "python $SCRIPTS/cmcwald_weightedZ_AMPNIH_UKB_ALL_meta_analysis_numvar.py \
			--test CMCWald \
    		--numvar ${numvar} \
    		--ntotal 53843 \
			--ukb_all_cases $WORK_DIR/UKB_EXOM_ALL_PD/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}.freqUpper${cutoff}.CMCWald.assoc \
			--amp_nih $WORK_DIR/AMPxNIH/AMP_NIH_noRelateds_${variants}_freqUpper${cutoff}.CMCWald.assoc \
			--variant_group ${variants} \
			--maf ${cutoff} \
			--group meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021 \
			-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}" >> cmc_weightedZ_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm
		done
	done
done

swarm -f cmc_weightedZ_AMP_NIH_UKB_ALL_PD_PHENO_META.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module python --bundle 8
# 24 commands run in 3 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 3 cores and 6 cpus
# 32504301 - done

## Checks
# 1 test * 4 variant groups * 2 MAFs = 8 files per numvar directory
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_2 | wc -l # 
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_3 | wc -l # 
ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_4 | wc -l # 

# done 

##########################################################################################################################################
##########################################################################################################################################
##### 7. SUMMARY OF RESULTS #############################################################################################################
##########################################################################################################################################
##########################################################################################################################################

mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_SUMMARIES/

for meta in "${metas[@]}"
do
	for variants in "${variant_group[@]}"
	do
		for cutoff in "${MAF[@]}"
		do
			for numvar in "${number_variants[@]}"
	    	do
				echo "Now looking at cohort: ${meta}; variant group ${variants} at MAF ${cutoff} at minimum numvar ${numvar}"
				echo "COHORT: ${meta}; VARIANT GROUP: ${variants}; MAF: ${cutoff}; NUMVAR >= ${numvar}" >> $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_SUMMARIES/cmcwald_weightedz_top_25_lines_${meta}_minvar${numvar}.txt
				cat $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/${meta}/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}* | head -25 >> $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_SUMMARIES/cmcwald_weightedz_top_25_lines_${meta}_minvar${numvar}.txt
				echo " " >> $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_SUMMARIES/cmcwald_weightedz_top_25_lines_${meta}_minvar${numvar}.txt
			done
		done
	done
done

# 1 cohorts * 3 minimum numvar = 3 summary files (each summary has each type of variant group and MAF cut-off)
# done


##########################################################################################################################################
##########################################################################################################################################
##### 8. (COMBINED PS) CALCULATE LAMBDAS PER META COHORT ################################################################################
##########################################################################################################################################
##########################################################################################################################################

	## Cohort Numbers 

	# |                 Cohort                 	| Cases 	| Controls 	| TOTAL 	|
	# |:--------------------------------------:	|:-----:	|:--------:	|:-----:	|
	# |                                AMPxNIH 	|  3376 	|     4610 	|  7986 	|
	# |                                UKB ALL 	|  7806 	|    38051 	| 45857 	|
	# |                       UKB case-control 	|  1105 	|     5643 	|  6748 	|
	# |                           UKB siblings 	|   668 	|     3463 	|  4131 	|
	# |                            UKB parents 	|  6033 	|    28945 	| 34978 	|
	# |                                    GNE 	|  2700 	|     9000 	| 11700 	|
	# |                                        	|       	|          	|       	|
	# |      GNE_AMP_NIH_UKB_ALL_PD_PHENO_META 	| 13882 	|    51661 	| 65543 	|
	# |      GNE_AMP_NIH_UKB_CASE_CONTROL_META 	|  7181 	|    19253 	| 26434 	|
	# |      GNE_AMP_NIH_UKB_CASE_PROXIES_META 	| 13882 	|    51661 	| 65543 	|
	# | GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META 	| 13214 	|    48198 	| 61412 	|
	# |          AMP_NIH_UKB_ALL_PD_PHENO_META 	| 11182 	|    42661 	| 53843 	|
cd $WORK_DIR/
mkdir $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_LAMBDAS

cd $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META

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
			--input $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.combined_Ps.${test}.tab \
			--variant_group ${variants} \
			--ncases 11182 \
			--ncontrols 42661 \
			--maf ${cutoff} \
			--group AMP_NIH_UKB_ALL_PD_PHENO_META \
			-o $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_LAMBDAS/AMP_NIH_UKB_ALL_PD_PHENO_META_${variants}_freqUpper${cutoff}_minvar${numvar}.combined_Ps
			done
		done
	done
done

# done 

## Summarize 
cd $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_LAMBDAS/
head -1 AMP_NIH_UKB_ALL_PD_PHENO_META_ALL_LOF_HC_LOFTEE_and_CADD_20_VEP_freqUpper0.001_minvar2.combined_Ps.lambdas.CMCWald.tab >> lambdas_combinedPs_summary.txt
for f in *.tab; do sed -n '2p' $f >> lambdas_combinedPs_summary.txt ; done

## Copy to working directory
mkdir $WORK_DIR/AMP_NIH_UKB_LAMBDA_SUMMARIES
scp lambdas_combinedPs_summary.txt $WORK_DIR/AMP_NIH_UKB_LAMBDA_SUMMARIES
# done


##########################################################################################################################################
##########################################################################################################################################
##### 9. (SUMMED ZS) CALCULATE LAMBDAS PER META COHORT ################################################################################
##########################################################################################################################################
##########################################################################################################################################

cd $WORK_DIR/
mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_LAMBDAS

cd $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
    	do
		python $SCRIPTS/calclambdas_meta.py \
		--test CMCWald \
		--numvar ${numvar} \
		--input $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.summedZ.CMCWald.tab \
		--variant_group ${variants} \
		--ncases 11182 \
		--ncontrols 42661 \
		--maf ${cutoff} \
		--group AMP_NIH_UKB_ALL_PD_PHENO_META \
		-o $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_LAMBDAS/AMP_NIH_UKB_ALL_PD_PHENO_META_${variants}_freqUpper${cutoff}_minvar${numvar}.summedZ
		done
	done
done

# done

## Summarize 
cd $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_LAMBDAS/
head -1 AMP_NIH_UKB_ALL_PD_PHENO_META_ALL_MODERATE_HIGH_IMPACT_SNPEFF_freqUpper0.01_minvar3.summedZ.lambdas.CMCWald.tab >> lambdas_summedZ_summary.txt
for f in *.tab; do sed -n '2p' $f >> lambdas_summedZ_summary.txt ; done

## Copy to working directory
scp lambdas_summedZ_summary.txt $WORK_DIR/AMP_NIH_UKB_LAMBDA_SUMMARIES
# done

##########################################################################################################################################
##########################################################################################################################################
##### 10. (WEIGHTED ZS) CALCULATE LAMBDAS PER META COHORT ################################################################################
##########################################################################################################################################
##########################################################################################################################################

cd $WORK_DIR/
mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_LAMBDAS

cd $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META

for variants in "${variant_group[@]}"
do
	for cutoff in "${MAF[@]}"
	do
		for numvar in "${number_variants[@]}"
    	do
		python $SCRIPTS/calclambdas_meta.py \
		--test CMCWald \
		--numvar ${numvar} \
		--input $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/*${variants}_freqUpper${cutoff}.weightedZ.CMCWald.tab \
		--variant_group ${variants} \
		--ncases 11182 \
		--ncontrols 42661 \
		--maf ${cutoff} \
		--group AMP_NIH_UKB_ALL_PD_PHENO_META \
		-o $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_LAMBDAS/AMP_NIH_UKB_ALL_PD_PHENO_META_${variants}_freqUpper${cutoff}_minvar${numvar}.weightedZ
		done
	done
done

# done

## Summarize 
cd $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_LAMBDAS/
head -1 AMP_NIH_UKB_ALL_PD_PHENO_META_ALL_LOF_HC_LOFTEE_freqUpper0.001_minvar2.summedZ.lambdas.CMCWald.tab >> lambdas_weightedZ_summary.txt
for f in *.tab; do sed -n '2p' $f >> lambdas_weightedZ_summary.txt ; done

## Copy to working directory
scp lambdas_weightedZ_summary.txt $WORK_DIR/AMP_NIH_UKB_LAMBDA_SUMMARIES
# done

##########################################################################################################################################
##########################################################################################################################################
##### 11. (COMBINED PS) META-ANALYZE: AMPXNIH (WGS) + UKB PD_ALL MANHATTANS ##############################################################
##########################################################################################################################################
##########################################################################################################################################

mkdir $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS

cd $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS

for test in "${test_type[@]}"
do
    for variants in "${variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "Rscript --vanilla ./../../manhattan_burden_1e-6.R \
    			$WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}.combined_Ps.${test}.tab" >> $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS/manhattans.swarm
    		done
    	done
    done
done

swarm -f $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS/manhattans.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module R --bundle 8
# 48 commands run in 6 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 6 cores and 12 cpus
# 32506376 - done 

ls $WORK_DIR/RESULTS_COMBINED_P/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS/*.pdf | wc -l # 48
# done 

##########################################################################################################################################
##########################################################################################################################################
##### 12. (SUMMED Z) META-ANALYZE: AMPXNIH (WGS) + UKB PD_ALL MANHATTANS #################################################################
##########################################################################################################################################
##########################################################################################################################################

mkdir $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS

cd $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS

for test in "${test_type[@]}"
do
    for variants in "${variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "Rscript --vanilla ./../../manhattan_burden_1e-6.R \
    			$WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}.summedZ.${test}.tab" >> $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS/manhattans.swarm
    		done
    	done
    done
done

swarm -f $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS/manhattans.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module R --bundle 8
# 48 commands run in 6 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 6 cores and 12 cpus
# 32540196 - done 

ls $WORK_DIR/RESULTS_SUMMEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS/*.pdf | wc -l # 24 
# done

##########################################################################################################################################
##########################################################################################################################################
##### 12. (WEIGHTED Z) META-ANALYZE: AMPXNIH (WGS) + UKB PD_ALL MANHATTANS #################################################################
##########################################################################################################################################
##########################################################################################################################################

mkdir $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS

cd $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS

for test in "${test_type[@]}"
do
    for variants in "${variant_group[@]}"
    do
    	for cutoff in "${MAF[@]}"
    	do
    		for numvar in "${number_variants[@]}"
    		do
    			echo "Rscript --vanilla ./../../manhattan_burden_1e-6.R \
    			$WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_META/minimum_numvar_${numvar}/meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_${variants}_freqUpper${cutoff}.weightedZ.${test}.tab" >> $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS/manhattans.swarm
    		done
    	done
    done
done

swarm -f $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS/manhattans.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module R --bundle 8
# 48 commands run in 6 subjobs, each command requiring 1.5 gb and 1 thread, running 8 processes serially per subjob, allocating 6 cores and 12 cpus

ls $WORK_DIR/RESULTS_WEIGHTEDZ_CMCWALD/AMP_NIH_UKB_ALL_PD_PHENO_MANHATTANS/*.pdf | wc -l # 24 
# done
