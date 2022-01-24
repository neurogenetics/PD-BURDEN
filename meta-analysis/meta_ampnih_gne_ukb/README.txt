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