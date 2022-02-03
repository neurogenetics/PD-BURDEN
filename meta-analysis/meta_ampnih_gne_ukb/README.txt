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
		# Weighted Z for CMC Wald: /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/RESULTS_SUMMEDZ_CMCWALD/SUMMARIES
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