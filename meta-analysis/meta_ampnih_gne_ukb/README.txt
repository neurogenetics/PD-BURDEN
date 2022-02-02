# README by Mary Makarious; Last Updated 02-FEB-2022
	# PD Burdens meta-analysis with AMPxNIH, Genentech, and UKB cohorts 

	# Results broken by minimum number of variants 
		# Combined p-values method (Fisher) for SkatO and CMC Wald: /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/RESULTS/SUMMARIES/
			# LAMBDAS/ directory included where lambda and lambda 1000 values were calculated 
		# Summed Z for CMC Wald: /data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/geneburden_risk_wGenentech/RESULTS_WEIGHTEDZ_CMCWALD/SUMMARIES
			# LAMBDAS/ directory included where lambda and lambda 1000 values were calculated 
			# Directory naming currently misleading; will be corrected
		# Weighted Z for CMC Wald: Coming soon...
			# Will include a LAMBDAS/ directory as well 

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

	# Cohort                                   	| Cases  | Controls
	# ---------------------------------------- 	| ------ | --------
	# AMPxNIH                                  	| 3,376  | 4,610   
	# UKB ALL                                  	| 7806   | 38051   
	# UKB case-control                         	| 1105   | 5643    
	# UKB siblings                             	| 668    | 3463    
	# UKB parents                              	| 6033   | 28945   
	# GNE                                      	| 2700   | 9000    
	#                                          	|        |         
	# GNE_AMP_NIH_UKB_ALL_PD_PHENO_META  		| 13,882 | 51,661  
	# GNE_AMP_NIH_UKB_CASE_CONTROL_META  		| 7,181  | 19,253  
	# GNE_AMP_NIH_UKB_CASE_PROXIES_META  		| 13,882 | 51,661  
	# GNE_AMP_NIH_UKB_CASE_PARENT_PROXY_META	| 13,214 | 48,198