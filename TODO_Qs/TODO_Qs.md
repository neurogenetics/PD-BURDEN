

## AMPxNIH Dataset 
### Done 
1. Generate NIH PLINK2 merged all chromosomes file
2. Generate AMP PLINK2 merged all chromosomes file 
3. Reformat the NIH data to chr:bp:ref:alt and convert to PLINK1.9
4. Split+normalize+reformat the AMP data to chr:bp:ref:alt and convert to PLINK1.9
5. QC the AMP v2.5 data (ancestry+relatedness)
6. KING analysis on AMP v2.5 and remove duplicates (none found)
7. KING analysis on NIH data and remove duplicates 
8. Merge AMP and NIH data in PLINK v1.9 
9. QC the merged autosomes data (ancestry+relatedness; remove any non-Europeans)
10. Make list of all duplicates and remove from the NIH data 
11. Re-run the QC (ancestry+relatedness; no duplicates found)
12. Create covariate file for merged data 
13. Merge the Euro PCs to the merged data
14. Make list of samples to keep removing non-Euros and duplicates
15. Make separate covariate file with just PD_EXTRA_PHENO 
16. Make PLINK1.9 files using cases and controls from PD_EXTRA_PHENO
17. Prune the PD_EXTRA_PHENO 
18. Re-create PCs with just the pruned PD_EXTRA_PHENO people
19. Remake merged PCs+cov file


### Upcoming 

* --keep the proper Euro, no dups, PD_EXTRA_PHENO samples from AMP PLINK2 files
* --keep the proper Euro, no dups, PD_EXTRA_PHENO samples from NIH PLINK2 files
* Pick out random samples from AMP and convert to VCF
* Pick out random samples from NIH and convert to VCF
* Annotate in ANNOVAR the subset AMP VCF
* Annotate in ANNOVAR the subset NIH VCF
* ...?
* Eventually clean directory...

# Files
QC done here: `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5`

- Final PD_EXTRA AMPxNIH files to use (Euro, no dups, with relatives):
    - Covariate file with PD_EXTRA PCs: `/data/CARD/PD/AMP_NIH/EuroONLY_PCs_noDups_withRelatives_PDEXTRA_amp_nih_merged_covariateFile_July2021.txt`
    - PLINKv1.9: `/data/CARD/PD/AMP_NIH/data/CARD/PD/AMP_NIH/PD_extra_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_withRelatives_wPheno_wSex`
    - PLINKv2 (AMP v2.5 Only):
    - PLINKv2 (NIH Only): 

- Final Full AMPxNIH PLINKv1.9 files to use: 
    - PLINKv1.9 European only, no duplicates (w/ relatives): `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/NIH_AMPv2.5_sampleQC_EURO_noDups.* `
    - Covariate File: `/data/CARD/PD/AMP_NIH/EuroONLY_PCs_noDups_withRelatives_full_amp_nih_merged_covariateFile_July2021.txt`
    - Both sets were formatted to be chr:bp:ref:alt hg38 
    - Within NIH only: Split+normalize multi-allelics, KING analysis with duplicates removed
    - AMP was checked for in-cohort duplicates, none found
    - QC was done on merged AMPxNIH (ancestry+relatedness), but 3846 duplicates were found in Europeans
    - Duplicates were removed, and QC (ancestry+relatedness) was re-run on Europeans 
        - 0 duplicates found
        - 667 1st degree relatives
        - 206 2nd degree relatives 
        - Total of 12,616 European samples to use
    - Rest of QC files (eigens, ancestry files, Scree plot): `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5`
    - Set of IDs to keep (prioritizing AMP; removing NIH duplicates; Euro only; keeping relatives): `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_withRelatives.txt`
    - Counts 
        - PD_EXTRA_PHENO (1=Healthy control with no neurological disorders, 2=Positive PD case not in a genetic enrichment study, but includes prodromal or SWEDDs; -9 for all other phenotypes):
            - control=1; 4731 samples
            - other dx=-9; 4411 samples 
            - case=2; 3474 samples
        - PD_PHENO (1=Healthy control with no neurological disorders and 2=Positive PD case not in a genetic enrichment study, prodromal, or a SWEDD): 
            - control=1; 4731 samples
            - other dx=-9; 4488 samples 
            - case=2; 3397 samples


- [FULL COHORT - EURO] Relatives files:
    - With within-cohort duplicates removed, KING full analysis output: `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/king_all_chr_noDups.kin0`
    - Text file with within-cohort duplicates removed, all duplicates between AMP and NIH (Euro only): `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/duplicate_samples_between_NIH_AMP_v25.txt`
    - Text file with all duplicates removed, all 1st degree relatives between AMP and NIH (Euro only): `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/first_degree_samples_between_NIH_AMP_v25.txt`
    - Text file with all duplicates removed, all 2nd degree relatives between AMP and NIH (Euro only): `data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/second_degree_samples_between_NIH_AMP_v25.txt`
    - Text file with all duplicates removed, all 1st+2nd degree relatives between AMP and NIH (Euro only): `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/relative_samples_between_NIH_AMP_v25.txt`
    - NIH duplicates to remove from merged set, in PLINK format: `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/NIH_duplicate_samples_toRemove.txt`

- [FULL COHORT - EURO] PCA files: 
    - PLINKv1.9 European only, no duplicates (w/ relatives) Eigenvectors: `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/euro_pca.eigenvec`
    - PLINKv1.9 European only, no duplicates (w/ relatives) Eigenvalues: `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/euro_pca.eigenval`
    - Scree Plot of European only, no duplicates (w/ relatives): `/data/CARD/PD/AMP_NIH/euro_king_pca_AMP_NIH_v2.5/Scree-NIH-AMPv2.5-noDuplicates-EURO.png
`



## BURDEN TODO
- [ ] Subset the groups of variants (missense, LOF, splicing, CADD>10, CADD>20; see below) 
- [ ] Run frequencies in PLINK2
- [ ] Run burdens in RVTESTs (case-control and AAO)
	- PD_CASE_CONTROL
	- PD_PARENT_CONTROL (check if sibs are in there)

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
- [ ] Prep for meta-analysis (outlined in UKB burden README) 
- [ ] Calculate cumulative frequencies in R (outlined in UKB burden README) 
- [ ] Merge results and frequencies 
	??
- [ ] 

## Cleaned READMEs
- [x] GitHub README
- [x] UKBiobank_Burden (adapted from v2)

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
# Questions
1. What are the proxies exactly (just parents or both parents and siblings)
2. What are some downstream analyses 
---
