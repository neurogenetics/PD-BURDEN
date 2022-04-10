NIA Data Overview - PD Burdens 
	- Prepared by: Mary Makarious, Julie Lake, Cornelis Blauwendraat 
	- Last Updated: 25 MAR 2022

Directory Structure Overview 
.
├── README.txt
├── Robak_setFile_vs_geneFile
│		├── 4 *.CMC.assoc files 
│		├── 4 *.log files 
│		├── 4 *.SkatO.assoc files
│		├── refFlat_HG38_ROBAK_LSD_no_GBA.txt
│		└── refFlat_HG38_ROBAK_LSD_with_GBA.txt
├── gene-burden-meta
│	├── MANHATTANS
│	│	├── COMBINED_P
│	│	├── GROUP1-GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
│	│	├── GROUP2-GNE_AMP_NIH_UKB_CASE_CONTROL_META
│	│	├── GROUP3-GNE_AMP_NIH_UKB_CASE_PROXIES_META
│	│	├── INITIAL-DISCOVERY
│	│	└── WEIGHTED_Z
│	├── PER_COHORT
│	│	├── AMPxNIH
│	│	├── GENENTECH
│	│	├── UKB_EXOM_ALL_PD
│	│	├── UKB_EXOM_PD_CASE
│	│	├── UKB_EXOM_PD_PARENT
│	│	└── UKB_EXOM_PD_SIBLING
│	├── RESULTS
│	│	├── RESULTS_COMBINED_P
│	│	│	├── GROUP1-GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
│	│	│	├── GROUP2-GNE_AMP_NIH_UKB_CASE_CONTROL_META
│	│	│	├── GROUP3-GNE_AMP_NIH_UKB_CASE_PROXIES_META
│	│	│	└── INITIAL-DISCOVERY-AMP_NIH_UKB_ALL_PD_PHENO_META
│	│	└── RESULTS_WEIGHTEDZ_CMCWALD
│	│	    ├── GROUP1-GNE_AMP_NIH_UKB_ALL_PD_PHENO_META
│	│	    ├── GROUP2-GNE_AMP_NIH_UKB_CASE_CONTROL_META
│	│	    ├── GROUP3-GNE_AMP_NIH_UKB_CASE_PROXIES_META
│	│	    └── INITIAL-DISCOVERY-AMP_NIH_UKB_ALL_PD_PHENO_META
│	└── SUMMARIES
│	    ├── COMBINEDP_SUMMARIES
│	    ├── LAMBDA_SUMMARIES
│	    └── WEIGHTEDZ_SUMMARIES
└── pathway-burden
	└── complete_c2_P_0.05
		├── forReplication_allConditions.txt
		├── c2.cp.v7.4.symbols.setFile.Oct52021.final.all_pathways.noX.FINAL.txt (C2 Canonical v7.4 MSigDB)
		├── ROBAK_SETLIST_hg38_noX.txt (Robak et al., 2017; Lysosomal)
		├── 8 files AMP_NIH_*_ROBAK_LSD.SkatO.assoc (Robak et al., 2017; Lysosomal)
		└── 8 files *PATHWAYS.SkatO.assoc (C2 Canonical v7.4 MSigDB)


Cohorts 
	- UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL
	- UKB_EXOM_PD_CASE_CONTROL
	- UKB_EXOM_PD_PARENT_CONTROL
	- UKB_EXOM_PD_SIBLING_CONTROL
	- AMP_NIH_noRelateds

Variant Groups 
	- ALL_MISSENSE_SNPEFF
	- ALL_LOF_HC_LOFTEE
	- ALL_LOF_HC_LOFTEE_and_CADD_20_VEP
	- ALL_MODERATE_HIGH_IMPACT_SNPEFF

MAF Cut-Offs
	- 0.01
	- 0.001

Rare Variant Burden Analyses 
	- Gene burden 
	- Gene-set pathway burden

Tests via RVTests
	- SkatO (gene; pathways)
	- CMCWald (gene)

Overview of Cohorts 
	- UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL; PD cases, parent proxies and sibling proxies and controls (Cases = 7,806; Controls = 38,051)
	- UKB_EXOM_PD_CASE_CONTROL; PD cases and controls (Cases = 1,105; Controls = 5,643)
	- UKB_EXOM_PD_PARENT_CONTROL; PD parent proxy-cases and controls (Proxy-cases = 6,033; Controls = 28,945)
	- UKB_EXOM_PD_SIBLING_CONTROL; PD sibling proxy-cases and controls (Proxy-cases = 668; Controls = 3,463)
	- AMP_NIH_noRelateds; combined NIH clinical center PD WGS with AMP-PD v2.5 data (Cases = 3,376; Controls = 4,610)

Overview of Variant Groups
	- ALL_MISSENSE_SNPEFF; where all missense variants as defined by SnpEff
	- ALL_LOF_HC_LOFTEE; where all LOF variants assigned "HIGH CONFIDENCE" as defined by LOFTEE
	- ALL_LOF_HC_LOFTEE_and_CADD_20_VEP; where all variants that have a CADD PHRED >20 or are LOF variants that are "HIGH CONFIDENCE" as defined by LOFTEE
	- ALL_MODERATE_HIGH_IMPACT_SNPEFF; where all variants with an “IMPACT” assigned “MODERATE” or “HIGH” by SnpEff (this includes LoF, missense, indels and others)

Overview of Meta-analysis Groups
	- INITIAL-DISCOVERY: AMP-PDxNIH and UKB using Genentech as a replication cohort. Stratified by numvar minimum, burden test, variant class, and MAF
	- GROUP1-GNE_AMP_NIH_UKB_ALL_PD_PHENO_META: Genentech vs. AMPxNIH vs. UKB all PD phenotypes (PD cases+1st degree proxies). Stratified by numvar minimum, burden test, variant class, and MAF 
	- GROUP2-GNE_AMP_NIH_UKB_CASE_CONTROL_META: Genentech vs. AMPxNIH vs. UKB PD cases. Stratified by numvar minimum, burden test, variant class, and MAF 
	- GROUP3-GNE_AMP_NIH_UKB_CASE_PROXIES_META: Genentech vs. AMPxNIH vs. UKB PD cases vs. UKB PD parents vs. UKB PD sibling proxies. Stratified by numvar minimum, burden test, variant class, and MAF 

Overview of Meta-analysis Approaches
	- COMBINEDP: Combined P-value approach, from CMC Wald or SKAT-O results. Usesd Fisher's test.
	- WEIGHTEDZ: Weighted Z-score, from CMC Wald results. Takes into account directionality of effect and size of cohort.

Notes
	- General 
		- Data and analysis all in GRCh38 (hg38)
		- Naming convention: ${COHORT}_${VARIANT_GROUP}.freqUpper${MAF}.${TEST}.assoc
		- Software: VEP v104; RVTests v2.1.0
		- European ancestry only, not including related individuals (0.125 PIHAT cut-off)
		- UKB cohorts: Age-matched controls chosen, different controls per group
		- Phenotypes in AMP v2.5 data defined as: 1=Healthy control with no neurological disorders and 2=Positive PD case not in a genetic enrichment study, prodromal, or a SWEDD
		- Covariates for AMPxNIH Cohort: SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5; Covariates for UKB Cohorts: SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5
		- Total Gene Burden .assoc Files Count: 5 cohorts * 4 variant groups * 2 MAF cut-offs * 2 tests = 80 files (16 files per cohort)
		- Total Pathway Burden .assoc Files Count: 2 pathway sets * 1 cohort * 4 variant groups * 2 MAF cut-offs * 1 test = 16 files
	
	- Gene Burden Analysis 
		- First result reported in CMC Wald .assoc files were kept, p-value accounts for all covariates (suffix: *.CMCWald.noCovEffect.assoc)
		- Meta-analysis to be performed by combining P-values using Fisher's method (custom Python scripts)
	
Pathways Notes 
	Overview of Gene-Set Pathway Files with p<0.05 in AMP_NIH_noRelateds Cohort (C2 canonical pathways from v7.4 MSigDB)
		- Pathways with GBA and LRRK2 were also run without GBA and LRRK2, denoted with suffixes (ex. *_no_GBA; *_with_GBA)
		- 908 forReplication_allConditions.txt
			- 908 combinations (pathway/MAF/variant group) had p<0.05 in AMPxNIH cohort with no related individuals 
		- Gene-Set Pathway Burden Analysis
			- Pathways investigated currently include 
				- MSigDB v7.4 C2 canonical pathways (do not include X chromosome)
				- Robak et al., 2017 manuscript (DOI: 10.1093/brain/awx285) 
			- Total of 908 pathway combinations to be run; overview @ forReplication_allConditions.txt
			- setFiles used with RVTests (hg38) @ 
				- c2.cp.v7.4.symbols.setFile.Oct52021.final.all_pathways.noX.txt (MSigDB; C2) 
				- ROBAK_SETLIST_hg38_noX.txt (Robak et al., 2017; Lysosomal)

CHANGELOG 
	- 30-MAR-2022: Updated RESULTS/ directory with the full results files from meta-analysis and updated Manht
	- 25-MAR-2022: Added gene burden meta-analysis results, lambda calculations, Manhattan plots, and summaries. Restructured directory + updated README
	- 02-DEC-2021: Added both setFile and geneFile results from Robak 2017 to Robak_setFile_vs_geneFile/ directory for troubleshooting 
	- 01-DEC-2021: Added Robak results + added setFile + updated README
	- 30-NOV-2021: Prepared p<0.05 SKAT-O .assoc files following gene-set pathway burden analysis in AMP_NIH_noRelateds Cohort + added setFile (hg38)
	- 29-NOV-2021: Added covariate information + Replaced CMC Wald files with updated results (accounting for all covariates)	
	- 23-NOV-2021: Replacing CMC files with CMC Wald files 
	- 19-NOV-2021: Preparation of PD burdens data + overviews 