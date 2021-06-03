# Case-Control Parkinson's Disease Burden Testing 

`LNG ❤️ Open Science 😍`

 - **Project:** PD Burdens in WGS (NIH, AMP v2.5, and UKB exomes) 
 - **Authors:** Mary B Makarious, Cornelis Blauwendraat, Ben Jacobs, Julie Lake, Hampton Leonard, Mike Nalls, Alastair Noyce, Andrew Singleton 
 - **Date Last Updated:** June 2021 
    - **Update Description:** Edits to README


---
### Quick Description: 
Performing genome-wide burden testing in all exome/genome Parkinson's disease data of European ancestry available in AMP, UKB, and the NIH clinical center

### Motivation/Goals:
To perform case-control burden tests in large-scale PD cohorts to look at both known and candidate genes associated with each of the phenotypes in our cohort and 1) confirm the known genes 2) assess the burden associated with candidate genes.

Analysis will be split as the following:
1. Positive PD cases vs controls
2. Positive PD proxies (parents and siblings with PD) vs controls 
3. All PD cases and proxies vs controls 

### Background:
There are a number of known genes implicated in Parkinson's disease (PD). Missense and loss-of-function mutations (defined as stop-gain, stop-loss, frameshift, and splicing variants) have been shown to cause autosomal dominant Parkinson's disease. **Gene burden tests are conducted to find the association of rare variants on a particular phenotype.** Generally, information is aggregated across several variant sites within a gene to calculate association signals and to reduce the penalty of multiple testing. Burden tests generate a burden score by taking a weighted linear combination of the mutations found within a gene or indicating whether there is any mutation within a gene. 

Here we screen the UK Biobank cohort, the AMP v2.5 cohort, and the NIH clinical center cohort for pathogenic mutations and to look at both known and candidate genes associated with each of the phenotypes in our cohort and 1) confirm the known genes 2) assess the burden associated with candidate genes.




We have identified xx


### Link to Manuscript:

### Data 
- PD WGS (Combination of AMP-PD data and internal LNG WGS data)
	- xx controls, xx cases 
	- Controls defined as healthy control with no neurological disorders
	- Cases defined as a positive PD case not in a genetic enrichment study, but includes prodromal or SWEDDs
- UKB exomes (Combination of PD cases, PD proxies [parent or sibling with PD] and controls)
	-  ~50K controls, ~600 cases, and ~6K proxies 


## Structure of README:
### [0. Repository Information](#0)
This section describes the directory and files available in this GitHub repository 

### [1. Data Pre-Processing](#1)
This section goes through:
- Downloading and splitting the refFlat hg38 file into chromosomes 
- QC the UKB exomes
- QC the AMP v2.5 data
-  Explanation of data in the covariate file 

### [2. Individual Level Data Filtering](#2)
 This section goes through: 
- xx

### [3. Annotation of Variants](#3)
 This section goes through: 
- xx

### [4. Subset Variant Classes](#4)
 This section goes through: 
- xx

### [5. Subset Genetic Data](#5)
 This section goes through: 
- xx
- 
### [6. Run Gene Burdens in RVTests](#6)
 This section goes through: 
- xx

### [7. Filter Burdens for Number of Variants (min. of 3)](#7)
 This section goes through: 
- xx

### [8. Meta-analyze Burden Tests](#8)
 This section goes through: 
- xx

### [9. Create Cumulative Case-Control Frequency Files](#9)
 This section goes through: 
- xx

### [10. Create Final Results File for Interpretation](#10)
 This section goes through: 
- xx

---
<a id="0"></a>
# 0. Repository Information
#### Files Available on GitHub

- `UKBiobank_Burden.md` : describes workflow of UK Biobank data PD WGS
- `Burden.md` => describes workflow of PD WGS data Meta analyses of
- `Burden.md` => describes the meta-analyses of burden results PD WGS Age
- `at onset Burden.md` => describes workflow of PD WGS age at onset
- `burden testing Pathways testing.md` => sandbox now for ideas on
- `pathways Variant class enrichments.md` => sandbox now variant class
   enrichments per disease group

<a id="1"></a>
# 1. Data Pre-Processing

#### Downloading and splitting the refFlat hg38 file into chromosomes 
```bash
# Go to the proper directory
cd ${WORK_DIR}/UKBIOBANK/EXOME_DATA_200K/REFFLAT/

# Download the RefFlat file from UCSC 
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz

# Separate per chromosome and sort
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
	zless refFlat.txt.gz | awk '$3 == "chr'$chnum'"' | sort -nk5 > refFlat_HG38_chr"$chnum".txt
done

# Separate one for chromosome X
zless refFlat.txt.gz | awk '$3 == "chrX"' | sort -nk5 > refFlat_HG38_chr23.txt
```


<a id="2"></a>
# 2. Individual Level Filtering 
Samples that remained in the analysis were 
- European ancestry (as determined by overlapping SNPs with HapMap)
- Non-related (PIHAT <0.125)

Additional control filtering for UK Biobank 
- Controls have age of recruit >60 and no AD or PD parent and no PD-ism, no dementia diagnosis


<a id="3"></a>
# 3. Annotation of Variants

#### Breakdown of Files:

 - `ALL_MISSENSE.txt` : All missense variants at 5, 1, 0.5, and 0.1%
   frequency
 - `ALL_LOF.txt` : All loss-of-function variants at (splicing,
   frameshift, stop-gain, stop-loss) at 5, 1, 0.5 and 0.1% frequency
 - `ALL_CADD_20.txt` : All variants with a CADD score above 20 at 5, 1,
   0.5 and 0.1% frequency 
 - `ALL_CADD_10.txt` : All variants with a CADD score above 10 at 5, 1, 0.5 and 0.1% frequency
 - `ALL_MISSENSE_and_LOF.txt` : All missense, splicing, frameshift,
   stop-gain, stop-loss variants at 5, 1, 0.5 and 0.1% frequency


<a id="4"></a>
# 4. Subset Variant Classes
xx

<a id="5"></a>
# 5. Subset Genetic Data
xx

<a id="6"></a>
# 6. Run Gene Burdens in RVTests

**Analysis:** The methodology described here calculates gene burdens using several different statistical tests. 
Adapted from [RVTests GitHub](https://github.com/zhanxw/rvtests):

Burden Tests | Traits* | Covariates | Un/Related | Description
|:--------------:|:------:|:----------:|:-------------------:|:-----------
**CMC**             |  B, Q  |     Y      |         U           | Collapsing and combine rare variants by Bingshan Li
**Zeggini**         |  B, Q  |     Y      |         U           | Aggregate counts of rare variants by Morris Zeggini
**Madsen-Browning** | B     |     Y      |         U           | Up-weight rare variant using inverse frequency from controls by Madsen
**Fp**              |  B     |     Y      |         U           | Up-weight rare variant using inverse frequency from controls by Danyu Lin
**CMC Wald**        |  B, Q  |     Y      |         U           | Collapsing and combine rare variants, then perform Wald test
**Wald  Test**    |  B, Q  |     Y      |         U           | *(Single-Variant Model)* Only fit alternative model, and effect size will be estimated
**SKAT**     |  B, Q  |     Y      |         U           | *(Kernel)* Sequencing kernel association test by Shawn Lee.
**SKAT-O**     | B, Q  |     Y      |         U           | *(Kernel)* Optimal sequencing kernel association test (SKAT-O) by Shawn Lee. Based on optimal combination of burden and variance statistics 

\* Traits are (B)inary or (Q)uantitative 
