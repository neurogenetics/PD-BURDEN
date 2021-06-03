## TODO 
- [ ] QC the AMP v2.5+NIH data 
	- In Terra, recreate 2.5 covariate file and download
	- Remove duplicates in cohorts separately, merge, and remove duplicates 
	- Re-create covariate file (Mary has notebook for this)
	- Still waiting on the release, but have a full GWAS-flavored script written for data on Biowulf 
- [ ] Split up the AMP+NIH data into cases, parents, siblings, and all (like UKB) 
- [ ] Annotate one sample 
- [ ] Subset the groups of variants (missense, LOF, splicing, CADD>10, CADD>20) 
- [ ] Run frequencies in PLINK2
- [ ] Run burdens in RVTESTs ==(run PD_allCASEandPROXY_CONTROL ?)==
	- AD_CASE_CONTROL
	- AD_PARENT_CONTROL
	- PD_CASE_CONTROL
	- PD_PARENT_CONTROL
	- ==5th group here?==

	For the following variant levels 
	-  ALL_MISSENSE_and_LOF
	- ALL_CADD_20
	- ALL_CADD_10
	- ALL_LOF
	- ALL_MISSENSE
	- ALL_CADD_20_and_LOF

	For the following frequencies 
	- 0.05
	- 0.01
	- 0.005
	- 0.001
- [ ] Re-run all LOF and all CADD 20 for UKB (re-submitted today with more memory wondering why those 2 chromosomes didn't work?) 
- [ ] Prep for meta-analysis (outlined in UKB burden README) 
- [ ] Calculate cumulative frequencies in R (outlined in UKB burden README) 
- [ ] Merge results and frequencies 
	??
- [ ] 

## Cleaned READMEs
- [ ] GitHub README
- [ ] UKBiobank_Burden (adapted from v2)

---
# Burden Checklist - UKB 
- [ ] AD_CASE_CONTROL
	- [ ] ALL_MISSENSE_and_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] ALL_CADD_20
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] ALL_CADD_10
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_MISSENSE
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_CADD_20_and_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001

- [ ] AD_PARENT_CONTROL
	- [ ] ALL_MISSENSE_and_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] ALL_CADD_20
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] ALL_CADD_10
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_MISSENSE
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_CADD_20_and_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	
- [ ] PD_CASE_CONTROL
	- [ ] ALL_MISSENSE_and_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] ALL_CADD_20
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] ALL_CADD_10
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_MISSENSE
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_CADD_20_and_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
		
- [ ] PD_PARENT_CONTROL
	- [ ] ALL_MISSENSE_and_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] ALL_CADD_20
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] ALL_CADD_10
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_MISSENSE
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
	- [ ] 	ALL_CADD_20_and_LOF
		- [ ] 0.05
		- [ ] 0.01
		- [ ] 0.005
		- [ ] 0.001
---
### Question 1: Why did you set a random seed here?
Subsetting the controls into different groups depending on the case/proxy ?
Why are we splitting up siblings and parents (are we assuming the person is healthy and their parent/sibling isnt?)

```R
# From UKBiobank_Burden 

### Generate 4 lists + EUROPEAN == 1
	# 1) All PD cases (2998) vs control (N=14372, ratio -> 0.15)
	# 2) All PD parent (15010) vs control (N=72817, ratio -> 0.76)
	# 3) All PD sibling (1669) vs control (N=8623, ratio -> 0.09)
	# 4) Something with PD (19677) vs control
		# total controls => 95,812

PD <- subset(PD_mergedv2, PHENO=="PD" & EUROPEAN==1)
PARENT <- subset(PD_mergedv2, PHENO=="parent" & EUROPEAN==1)
SIBLING <- subset(PD_mergedv2, PHENO=="sibling" & EUROPEAN==1)

set.seed(123458)
RANDOM_order <- sample(1:3, size=95812, replace=TRUE, prob=c(.15,.76,.09))
CONTROLSv3 <- cbind(CONTROLSv2, RANDOM_order)
PD_control <- subset(CONTROLSv3, RANDOM_order==1)
PARENT_control <- subset(CONTROLSv3, RANDOM_order==2)
SIBLING_control <- subset(CONTROLSv3, RANDOM_order==3)
```

---

