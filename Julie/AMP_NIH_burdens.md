# AMP x NIH burdens for PD

---

## Structure of README:

### [1. Remove related samples](#1-remove-related-samples-1)
### [2. Make covariate file with PCs](#2-make-covariate-file-with-pcs-1)
### [3. Annotate with snpEff and loftee](#3-annotate-with-snpeff-and-loftee-1)
### [4. Subset variant classes](#4-subset-variant-classes-1)
### [5. Subset genetic data](#5-subset-genetic-data-1)
### [6. Run burdens with annovar annotations](#6-run-burdens-with-annovar-annotations-1)
### [7. Run burdens with snpEff/loftee annotations](#7-run-burdens-with-snpeffloftee-annotations-1)

---
### Here are some important files from this analysis:

- Working directory: /data/CARD/PD/AMP_NIH/no_relateds
- Covariate file for burdens: /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt
- Plink file merged AMP x NIH with no relateds: /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.*
- Annovar annotations: /data/CARD/PD/AMP_NIH/PD_PHENO/VCF_annotation/; /data/CARD/PD/AMP_NIH/PD_PHENO/pathoscores/
- snpEff/loftee annotations: /data/CARD/PD/AMP_NIH/no_relateds/snpEff_loftee/*annotated.txt
- Annovar burden results: /data/CARD/PD/AMP_NIH/no_relateds/burden_annovar
- snpEff/loftee burden results: /data/CARD/PD/AMP_NIH/no_relateds/burden_snpEff_loftee

---

## 1. Remove related samples
```
cd /data/CARD/PD/AMP_NIH/
mkdir no_relateds
cd no_relateds

## First remove relateds from the merged AMPxNIH data
sbatch --cpus-per-task=20 --mem=150g --mail-type=ALL --time=24:00:00 AMP_NIH_relatedness.sh

#!/bin/bash
# sh AMP_NIH_relatedness.sh
module load plink
module load GCTA
plink --bfile /data/CARD/PD/AMP_NIH/PD_PHENO/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_withRelatives_wPheno_wSex \
 --geno 0.01 --maf 0.05 --out input_for_prune --make-bed
plink --bfile input_for_prune --indep-pairwise 1000 10 0.02 --out TEMP_pruning
plink --bfile input_for_prune --extract TEMP_pruning.prune.in --make-bed --out TEMP_pruned_data
gcta64 --bfile TEMP_pruned_data --make-grm --out GRM_matrix --autosome --maf 0.05 --threads 20
gcta64 --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0.125 --make-grm --threads 20
plink --bfile /data/CARD/PD/AMP_NIH/PD_PHENO/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_withRelatives_wPheno_wSex \
 --keep GRM_matrix_0.125.grm.id --make-bed --out PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins

# Count number of people before and after relatedness of AMPxNIH
wc -l /data/CARD/PD/AMP_NIH/PD_PHENO/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_withRelatives_wPheno_wSex.fam
# 8128
wc -l PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.fam
# 7986
```

## 2. Make covariate file with PCs
```
cd /data/CARD/PD/AMP_NIH/no_relateds
mkdir PCA
cd PCA

sbatch --cpus-per-task=10 --mem=25g --mail-type=ALL --time=24:00:00 make_PCs_AMP_NIH.sh

#!/bin/sh
# sh make_PCs_AMP_NIH.sh
module load flashpca
module load plink
# PLINK1.9 - Prune data for flashPCA
plink --bfile /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins \
--maf 0.01 \
--geno 0.01 \
--hwe 1e-10 \
--exclude range /data/CARD/GENERAL/longLD_hg38.txt \
--make-bed --out prune_step1
plink --bfile prune_step1 \
--indep-pairwise 50 5 0.01
plink --bfile prune_step1 \
--extract plink.prune.in \
--make-bed \
--out PRUNED_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins
# PLINK1.9 - flashPCA - make 10 PCs 
flashpca --bfile PRUNED_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins \
--suffix _PD_PHENO_EURO_noDups_no_cousins.txt --numthreads 19
```
```
R
require(data.table)
require(dplyr)

# Load the full covariates file
cov <- fread("/data/CARD/PD/AMP_NIH/PD_PHENO/EuroONLY_PCs_noDups_withRelatives_PD_amp_nih_merged_covariateFile_July2021.txt",header=T)

# Remove the PCs since these were calculated with relateds
cov2 <- cov %>% select(1:25)

pc <- fread(file="pcs_PD_PHENO_EURO_noDups_no_cousins.txt",header=T)
pc$IID <- NULL
Mrg <- merge(cov2,pc,by="FID")
Mrg2 <- Mrg %>% select(FID, IID, SEX, AGE_ANALYSIS, PC1, PC2, PC3, PC4, PC5, PD_PHENO)

write.table(Mrg2, file="../COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt", quote=FALSE,row.names=F,sep="\t")
q()
n
```

## 3. Annotate with snpEff and loftee
```
cd /data/CARD/PD/AMP_NIH/no_relateds
mkdir snpEff_loftee
cd snpEff_loftee
```

### snpEff/snpSift

- Biowulf documentation:  https://hpc.nih.gov/apps/snpEff.html
- See here for all the fields you can extract: https://pcingola.github.io/SnpEff/ss_extractfields/

```
echo "CHROM
POS
ID
REF
ALT
FILTER
AF
AC
DP
MQ
ANN[*].ALLELE
ANN[*].EFFECT
ANN[*].IMPACT
ANN[*].GENE
ANN[*].GENEID
ANN[*].FEATURE
ANN[*].FEATUREID
ANN[*].BIOTYPE
ANN[*].RANK
ANN[*].HGVS_C
ANN[*].HGVS_P
ANN[*].CDNA_POS
ANN[*].CDNA_LEN
ANN[*].CDS_POS
ANN[*].CDS_LEN
ANN[*].AA_POS
ANN[*].AA_LEN
ANN[*].DISTANCE
ANN[*].ERRORS
LOF[*].GENE
LOF[*].GENEID
LOF[*].NUMTR
LOF[*].PERC
NMD[*].GENE
NMD[*].GENEID
NMD[*].NUMTR
NMD[*].PERC" > snpSift_extract_fields.txt

module load snpEff

## Annotate NIH VCFs
for chnum in {1..22} X Y;
do
    echo "java -Xmx20g -jar $SNPEFF_JAR -v GRCh38.86 /data/CARD/PD/AMP_NIH/PD_PHENO/subset_NIH_ONLY_forANNOVAR_chr${chnum}.vcf > subset_NIH_ONLY_forSNPEFF_chr${chnum}.eff.vcf
    java -jar $SNPSIFT_JAR dbnsfp -v -db /fdb/dbNSFP2/dbNSFP3.2a.txt.gz subset_NIH_ONLY_forSNPEFF_chr${chnum}.eff.vcf > subset_NIH_ONLY_forSNPEFF_chr${chnum}.annotated.vcf
    java -jar $SNPSIFT_JAR extractFields -s "," -e "." subset_NIH_ONLY_forSNPEFF_chr${chnum}.annotated.vcf $(tr '\n' ' ' < snpSift_extract_fields.txt) > subset_NIH_ONLY_forSNPEFF_chr${chnum}.annotated.txt" >> nih_snpeff.swarm
done

## Annotate AMP VCFs
for chnum in {1..22} X Y;
do
    echo "java -Xmx20g -jar $SNPEFF_JAR -v GRCh38.86 /data/CARD/PD/AMP_NIH/PD_PHENO/subset_AMP_ONLY_forANNOVAR_chr${chnum}.vcf > subset_AMP_ONLY_forSNPEFF_chr${chnum}.eff.vcf
    java -jar $SNPSIFT_JAR dbnsfp -v -db /fdb/dbNSFP2/dbNSFP3.2a.txt.gz subset_AMP_ONLY_forSNPEFF_chr${chnum}.eff.vcf > subset_AMP_ONLY_forSNPEFF_chr${chnum}.annotated.vcf
    java -jar $SNPSIFT_JAR extractFields -s "," -e "." subset_AMP_ONLY_forSNPEFF_chr${chnum}.annotated.vcf $(tr '\n' ' ' < snpSift_extract_fields.txt) > subset_AMP_ONLY_forSNPEFF_chr${chnum}.annotated.txt" >> amp_snpeff.swarm
done

## Run all swarms for snpEff annotation and snpSift field extraction
swarm -f nih_snpeff.swarm  --bundle 3 --verbose 1 --module snpEff -g 20 -t 3
# 72 commands run in 24 subjobs, each command requiring 20 gb and 3 threads, running 3 processes serially per subjob, allocating 48 cores and 96 cpus

swarm -f amp_snpeff.swarm  --bundle 3 --verbose 1 --module snpEff -g 20 -t 3
# 72 commands run in 24 subjobs, each command requiring 20 gb and 3 threads, running 3 processes serially per subjob, allocating 48 cores and 96 cpus
```

### LOFTEE

- GitHub:  https://github.com/konradjk/loftee
- Biowulf documentation:  https://hpc.nih.gov/apps/VEP.html

```
cd /data/CARD/PD/AMP_NIH/no_relateds/snpEff_loftee

# To test in interactive, load:
ml VEP/101
export PERL5LIB=$PERL5LIB:${VEPCACHEDIR}/Plugins/loftee_GRCh38

## Annotate NIH VCFs
for chnum in {1..22} X Y;
do
	echo "export PERL5LIB=$PERL5LIB:${VEPCACHEDIR}/Plugins/loftee_GRCh38
	vep \
	--offline --cache --dir_cache ${VEPCACHEDIR} \
	--input_file /data/CARD/PD/AMP_NIH/PD_PHENO/subset_NIH_ONLY_forANNOVAR_chr${chnum}.vcf \
	--species human --assembly GRCh38 --fasta ${VEPCACHEDIR}/GRCh38.fa \
	--output_file subset_NIH_ONLY_forLOFTEE_chr${chnum}.loftee --stats_file subset_NIH_ONLY_forLOFTEE_chr${chnum}_summary.loftee.txt \
	--plugin LoF,loftee_path:${VEPCACHEDIR}/Plugins/loftee_GRCh38,\
	human_ancestor_fa:${VEPCACHEDIR}/Plugins/loftee_GRCh38/human_ancestor.fa.gz,\
	conservation_file:${VEPCACHEDIR}/Plugins/loftee_GRCh38/loftee.sql,\
	gerp_bigwig:${VEPCACHEDIR}/Plugins/loftee_GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
	grep 'LoF=HC' subset_NIH_ONLY_forLOFTEE_chr${chnum}.loftee | uniq > subset_NIH_ONLY_forLOFTEE.chr${chnum}.refined_LoF.all
	grep 'LoF=HC' subset_NIH_ONLY_forLOFTEE_chr${chnum}.loftee | cut -f1 | uniq > subset_NIH_ONLY_forLOFTEE.chr${chnum}.refined_LoF.snps" >> nih_loftee.swarm
done

## Annotate AMP VCFs
for chnum in {1..22} X Y;
do
	echo "export PERL5LIB=$PERL5LIB:${VEPCACHEDIR}/Plugins/loftee_GRCh38
	vep \
	--offline --cache --dir_cache ${VEPCACHEDIR} \
	--input_file /data/CARD/PD/AMP_NIH/PD_PHENO/subset_AMP_ONLY_forANNOVAR_chr${chnum}.vcf \
	--species human --assembly GRCh38 --fasta ${VEPCACHEDIR}/GRCh38.fa \
	--output_file subset_AMP_ONLY_forLOFTEE_chr${chnum}.loftee --stats_file subset_AMP_ONLY_forLOFTEE_chr${chnum}_summary.loftee.txt \
	--plugin LoF,loftee_path:${VEPCACHEDIR}/Plugins/loftee_GRCh38,\
	human_ancestor_fa:${VEPCACHEDIR}/Plugins/loftee_GRCh38/human_ancestor.fa.gz,\
	conservation_file:${VEPCACHEDIR}/Plugins/loftee_GRCh38/loftee.sql,\
	gerp_bigwig:${VEPCACHEDIR}/Plugins/loftee_GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
	grep 'LoF=HC' subset_AMP_ONLY_forLOFTEE_chr${chnum}.loftee | uniq > subset_AMP_ONLY_forLOFTEE.chr${chnum}.refined_LoF.all
	grep 'LoF=HC' subset_AMP_ONLY_forLOFTEE_chr${chnum}.loftee | cut -f1 | uniq > subset_AMP_ONLY_forLOFTEE.chr${chnum}.refined_LoF.snps" >> amp_loftee.swarm
done

# Test with chromosome 22
# Need to grep the export statement, all chromosome 22 lines from the swarm file
grep -e export -e 22 nih_loftee.swarm | uniq | head -4  > test.swarm
swarm -f test.swarm  --bundle 4 --verbose 1 --module VEP/101 -g 5 -t 3
# 4 commands run in 1 subjob, each command requiring 5 gb and 3 threads, running 4 processes serially per subjob, allocating 2 cores and 4 cpus

# Remove before proceeding
rm *LoF* 
rm *loftee

## Run all swarms for loftee annotation
swarm -f nih_loftee.swarm  --bundle 4 --verbose 1 --module VEP/101 -g 5 -t 3
# 96 commands run in 24 subjobs, each command requiring 5 gb and 3 threads, running 4 processes serially per subjob, allocating 48 cores and 96 cpus

swarm -f amp_loftee.swarm  --bundle 4 --verbose 1 --module VEP/101 -g 5 -t 3
# 96 commands run in 24 subjobs, each command requiring 5 gb and 3 threads, running 4 processes serially per subjob, allocating 48 cores and 96 cpus
```

## 4. Subset variant classes

### 4a. Subset variant classes using annovar annotations
```
cd /data/CARD/PD/AMP_NIH/no_relateds
mkdir subset_variant_categories
cd subset_variant_categories
```

#### Reformat the annotations

i.e. combine chromosomes, merge AMP and NIH annotations, filter for unique variants

```
# First subset the headers only
head -1 /data/CARD/PD/AMP_NIH/PD_PHENO/VCF_annotation/annovar_subsetAMPonly_PD_chr1.hg38_multianno.txt > header_refGene.txt
head -1 /data/CARD/PD/AMP_NIH/PD_PHENO/pathoscores/annovar_subsetAMPonly_PD_chr1.hg38_multianno.txt > header_patho.txt

# Concat all of the annotations across chromosomes for the AMP dataset
# grep -v Otherinfo6 to remove the header line from each of the files since adding back only once
cat /data/CARD/PD/AMP_NIH/PD_PHENO/VCF_annotation/annovar_subsetAMPonly_PD_chr*.hg38_multianno.txt | grep -v Otherinfo6 > temp1
cat header_refGene.txt temp1 > annovar_subsetAMPonly_PD_all_chr.hg38_multianno.refGene.txt

cat /data/CARD/PD/AMP_NIH/PD_PHENO/pathoscores/annovar_subsetAMPonly_PD_chr*.hg38_multianno.txt | grep -v Otherinfo6 > temp2
cat header_patho.txt temp2 > annovar_subsetAMPonly_PD_all_chr.hg38_multianno.patho.txt

# Concat all of the annotations across chromosomes for the NIH dataset
cat /data/CARD/PD/AMP_NIH/PD_PHENO/VCF_annotation/annovar_subsetNIHonly_PD_chr*.hg38_multianno.txt | grep -v Otherinfo6 > temp3
cat header_refGene.txt temp3 > annovar_subsetNIHonly_PD_all_chr.hg38_multianno.refGene.txt

cat /data/CARD/PD/AMP_NIH/PD_PHENO/pathoscores/annovar_subsetNIHonly_PD_chr*.hg38_multianno.txt | grep -v Otherinfo6 > temp4
cat header_patho.txt temp4 > annovar_subsetNIHonly_PD_all_chr.hg38_multianno.patho.txt

wc -l annovar_subsetAMPonly_PD_all_chr.hg38_multianno.refGene.txt
# 3236381
wc -l annovar_subsetAMPonly_PD_all_chr.hg38_multianno.patho.txt
# 3236381
wc -l annovar_subsetNIHonly_PD_all_chr.hg38_multianno.refGene.txt
# 2800334
wc -l annovar_subsetNIHonly_PD_all_chr.hg38_multianno.patho.txt
# 2800334

## Concat the annotations from the NIH and AMP datasets and filter for unique variants (based on the Otherinfo6 column)
cat annovar_subsetAMPonly_PD_all_chr.hg38_multianno.refGene.txt annovar_subsetNIHonly_PD_all_chr.hg38_multianno.refGene.txt | awk 'NR == 1; NR > 1 {print $0 | "sort -u -k22,22"}' > annovar_AMPandNIH_PD_all_chr.hg38_multianno.refGene.txt
cat annovar_subsetAMPonly_PD_all_chr.hg38_multianno.patho.txt annovar_subsetNIHonly_PD_all_chr.hg38_multianno.patho.txt | awk 'NR == 1; NR > 1 {print $0 | "sort -u -k62,62"}' > annovar_AMPandNIH_PD_all_chr.hg38_multianno.patho.txt

wc -l annovar_AMPandNIH_PD_all_chr.hg38_multianno.refGene.txt
# 3743294
wc -l annovar_AMPandNIH_PD_all_chr.hg38_multianno.patho.txt
# 3814645
```

#### Subset "groups" of variants
```
# Missense
grep exonic annovar_AMPandNIH_PD_all_chr.hg38_multianno.refGene.txt | grep nonsynonymous | cut -f22 > ALL_MISSENSE.txt
# n=1196826

# LOF (stop, frame)
grep stopgain annovar_AMPandNIH_PD_all_chr.hg38_multianno.refGene.txt | cut -f22 > all_stopgain.txt
# n=35187
grep stoploss annovar_AMPandNIH_PD_all_chr.hg38_multianno.refGene.txt | cut -f22 > all_stoploss.txt
# n=1685
grep nonframeshift annovar_AMPandNIH_PD_all_chr.hg38_multianno.refGene.txt | cut -f22 > all_nonframeshift.txt
# n=27926
grep frame annovar_AMPandNIH_PD_all_chr.hg38_multianno.refGene.txt | grep -v nonframeshift | grep -v nonsynonymous | cut -f22 > all_frameshift.txt
# n=58053

# Splicing 
grep splicing annovar_AMPandNIH_PD_all_chr.hg38_multianno.refGene.txt | grep -v ncRNA | cut -f22 > all_splice.txt
# n=16759
# CADD > 10
awk '{ if($31 >= 10) { print }}' annovar_AMPandNIH_PD_all_chr.hg38_multianno.patho.txt | cut -f62 > ALL_CADD_10.txt
# n=1054994

# CADD > 20
awk '{ if($31 >= 20) { print }}' annovar_AMPandNIH_PD_all_chr.hg38_multianno.patho.txt | cut -f62 > ALL_CADD_20.txt
# n=755639

## Prepping final files:
cat all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_LOF.txt
cat ALL_MISSENSE.txt all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_MISSENSE_and_LOF.txt
cat ALL_CADD_20.txt all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_CADD_20_and_LOF.txt

wc -l ALL_LOF.txt
# 111684
wc -l ALL_MISSENSE_and_LOF.txt
# 1308510
wc -l ALL_CADD_20_and_LOF.txt
# 867323
```

#### Use these variant categories for burdens:
```
1. ALL_CADD_10.txt
2. ALL_CADD_20.txt
3. ALL_MISSENSE_and_LOF.txt
4. ALL_MISSENSE.txt
5. ALL_LOF.txt
6. ALL_CADD_20_and_LOF.txt
```

### 4b. Subset variant classes using snpeff/loftee

#### Reformat the annotations
```
cd /data/CARD/PD/AMP_NIH/no_relateds/snpEff_loftee

# First subset the header only
head -1 subset_NIH_ONLY_forSNPEFF_chr1.annotated.txt > header_snpEff.txt

# Concat all of the annotations across chromosomes for the AMP dataset
cat subset_NIH_ONLY_forSNPEFF_chr*.annotated.txt | grep -v CHROM > temp1
cat header_snpEff.txt temp1 > subset_NIH_ONLY_forSNPEFF_all_chr.annotated.txt
cat subset_NIH_ONLY_forLOFTEE.chr*.refined_LoF.snps > subset_NIH_ONLY_forLOFTEE.all_chr.refined_LoF.snps

# Concat all of the annotations across chromosomes for the NIH dataset
cat subset_AMP_ONLY_forSNPEFF_chr*.annotated.txt | grep -v CHROM > temp2
cat header_snpEff.txt temp2 > subset_AMP_ONLY_forSNPEFF_all_chr.annotated.txt
cat subset_AMP_ONLY_forLOFTEE.chr*.refined_LoF.snps > subset_AMP_ONLY_forLOFTEE.all_chr.refined_LoF.snps

wc -l subset_NIH_ONLY_forSNPEFF_all_chr.annotated.txt
# 2800334
wc -l subset_AMP_ONLY_forSNPEFF_all_chr.annotated.txt
# 3236381
wc -l subset_NIH_ONLY_forLOFTEE.all_chr.refined_LoF.snps
# 60021
wc -l subset_AMP_ONLY_forLOFTEE.all_chr.refined_LoF.snps
# 80697

## Concat the annotations from the NIH and AMP datasets and filter for unique variants (based on the ID column)
cat subset_NIH_ONLY_forSNPEFF_all_chr.annotated.txt subset_AMP_ONLY_forSNPEFF_all_chr.annotated.txt | awk 'NR == 1; NR > 1 {print $0 | "sort -u -k3,3"}' > subset_AMPandNIH_forSNPEFF_all_chr.annotated.txt
cat subset_NIH_ONLY_forLOFTEE.all_chr.refined_LoF.snps subset_AMP_ONLY_forLOFTEE.all_chr.refined_LoF.snps | uniq > ALL_LOF_HC_LOFTEE.txt

wc -l subset_AMPandNIH_forSNPEFF_all_chr.annotated.txt
# 3814645
wc -l ALL_LOF_HC_LOFTEE.txt
# 140718
```

#### Subset "groups" of variants
```
# Missense
grep missense_variant subset_AMPandNIH_forSNPEFF_all_chr.annotated.txt | cut -f3 > ALL_MISSENSE_SNPEFF.txt
# n=1303382

# LOF (stop, frame)
grep stop_gained subset_AMPandNIH_forSNPEFF_all_chr.annotated.txt | cut -f3 > all_stopgain.txt
# n=37273

grep stop_lost subset_AMPandNIH_forSNPEFF_all_chr.annotated.txt | cut -f3 > all_stoploss.txt
# n=3186

grep frameshift_variant subset_AMPandNIH_forSNPEFF_all_chr.annotated.txt | cut -f3 > all_frameshift.txt
# n=57489

# Splicing 
grep splice_region_variant subset_AMPandNIH_forSNPEFF_all_chr.annotated.txt | cut -f3 > all_splice.txt
# n=216347

# Moderate or high impact
grep -e MODERATE -e HIGH subset_AMPandNIH_forSNPEFF_all_chr.annotated.txt | cut -f3 > ALL_MODERATE_HIGH_IMPACT_SNPEFF.txt
# n=1488572

# High impact
grep HIGH subset_AMPandNIH_forSNPEFF_all_chr.annotated.txt | cut -f3 > ALL_HIGH_IMPACT_SNPEFF.txt
# n=162558

## Prepping final files:
cat all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_LOF_SNPEFF.txt

# How many LOF determined by snpEff
wc -l ALL_LOF_SNPEFF.txt
# 314295

# How many LOF determined by loftee with high confidence
wc -l ALL_LOF_HC_LOFTEE.txt
# 140718

# How many SNPs are in both the loftee high confidence and SnpEff LOF
comm -12 <(sort ALL_LOF_SNPEFF.txt) <(sort ALL_LOF_HC_LOFTEE.txt) > ALL_LOF_and_HC_LOFTEE.txt
# n=79473

# How many SNPs were identified by SnpEff but weren't in high confidence loftee
comm -23 <(sort ALL_LOF_SNPEFF.txt) <(sort ALL_LOF_HC_LOFTEE.txt) > ALL_LOF_not_HC_LOFTEE.txt
# n=234822

# How many SNPs identified by loftee high confidence aren't in the SnpEff LOF 
comm -13 <(sort ALL_LOF_SNPEFF.txt) <(sort ALL_LOF_HC_LOFTEE.txt) > HC_LOFTEE_not_ALL_LOF.txt
# n=61245

# Concat LOF with missense
cat ALL_MISSENSE_SNPEFF.txt ALL_LOF_SNPEFF.txt > ALL_MISSENSE_and_LOF_SNPEFF.txt
cat ALL_MISSENSE_SNPEFF.txt ALL_LOF_HC_LOFTEE.txt> ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE.txt
cat ALL_MISSENSE_SNPEFF.txt ALL_LOF_and_HC_LOFTEE.txt > ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE.txt

# Concat LOF with high impact variants
cat ALL_HIGH_IMPACT_SNPEFF.txt ALL_LOF_SNPEFF.txt > ALL_HIGH_IMPACT_and_LOF_SNPEFF.txt
cat ALL_HIGH_IMPACT_SNPEFF.txt ALL_LOF_HC_LOFTEE.txt> ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE.txt
cat ALL_HIGH_IMPACT_SNPEFF.txt ALL_LOF_and_HC_LOFTEE.txt > ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE.txt

wc -l ALL_MISSENSE_and_LOF_SNPEFF.txt
# 1617677
wc -l ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE.txt
# 1444100
wc -l ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE.txt
# 1382855
wc -l ALL_HIGH_IMPACT_and_LOF_SNPEFF.txt
# 476853
wc -l ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE.txt
# 303276
wc -l ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE.txt
# 242031
```

#### Use these variant categories for burdens:
```
1. ALL_MODERATE_HIGH_IMPACT_SNPEFF.txt
2. ALL_HIGH_IMPACT_SNPEFF.txt
3. ALL_MISSENSE_and_LOF_SNPEFF.txt
4. ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE.txt
5. ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE.txt
6. ALL_MISSENSE_SNPEFF.txt
7. ALL_LOF_SNPEFF.txt
8. ALL_LOF_HC_LOFTEE.txt
9. ALL_LOF_and_HC_LOFTEE.txt
10. ALL_HIGH_IMPACT_and_LOF_SNPEFF.txt
11. ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE.txt
12. ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE.txt
```

## 5. Subset genetic data

### 5a. Subset genetic data using annovar annotations
```
cd /data/CARD/PD/AMP_NIH/no_relateds
mkdir subset_genetic_data_annovar
cd subset_genetic_data_annovar

# Make vcf files for input for rvtests:
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden.sh ALL_MISSENSE.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden.sh ALL_LOF.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden.sh ALL_CADD_20.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden.sh ALL_CADD_10.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden.sh ALL_MISSENSE_and_LOF.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden.sh ALL_CADD_20_and_LOF.txt

#!/bin/sh
# sh make_vcf_for_burden.sh ALL_MISSENSE.txt
module load plink/2.0-dev-20191128
module load samtools
VARIANT_FILE=$1
OUTNAME=${VARIANT_FILE/".txt"/""}
plink2 --bfile /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins \
--extract /data/CARD/PD/AMP_NIH/no_relateds/subset_variant_categories/$VARIANT_FILE \
--export vcf bgz id-paste=iid --out AMP_NIH_noRelateds_${OUTNAME} --mac 1
tabix -p vcf  AMP_NIH_noRelateds_${OUTNAME}.vcf.gz
```

### 5b. Subset genetic data using snpEff/loftee
```
cd /data/CARD/PD/AMP_NIH/no_relateds
mkdir subset_genetic_data_snpEff_loftee
cd subset_genetic_data_snpEff_loftee

# Make vcf files for input for rvtests:
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_MODERATE_HIGH_IMPACT_SNPEFF.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_HIGH_IMPACT_SNPEFF.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_MISSENSE_and_LOF_SNPEFF.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_MISSENSE_SNPEFF.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_LOF_SNPEFF.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_LOF_HC_LOFTEE.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_LOF_and_HC_LOFTEE.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_HIGH_IMPACT_and_LOF_SNPEFF.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE.txt
sbatch --cpus-per-task=10 --mem=20g --mail-type=ALL --time=24:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE.txt

#!/bin/sh
# sh make_vcf_for_burden_snpEff_loftee.sh ALL_MISSENSE.txt
module load plink/2.0-dev-20191128
module load samtools
VARIANT_FILE=$1
OUTNAME=${VARIANT_FILE/".txt"/""}
plink2 --bfile /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins \
--extract /data/CARD/PD/AMP_NIH/no_relateds/snpEff_loftee/$VARIANT_FILE \
--export vcf bgz id-paste=iid --out AMP_NIH_noRelateds_${OUTNAME} --mac 1
tabix -p vcf  AMP_NIH_noRelateds_${OUTNAME}.vcf.gz
```

## 6. Run burdens with annovar annotations

### 6a. Sanity checks with GBA and LRRK2
```
cd /data/CARD/PD/AMP_NIH/no_relateds
mkdir burden_annovar
cd burden_annovar

module load rvtests

### PD risk

# GBA
rvtest --noweb --hide-covar --out GBA --burden cmc --kernel skato \
--inVcf /data/CARD/PD/AMP_NIH/no_relateds/subset_genetic_data_annovar/AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF.vcf.gz \
--pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--pheno-name PD_PHENO \
--covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--freqUpper 0.05 \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# CMC: 2.29221e-13
# SKATO: 2.1519e-10

# LRRK2
rvtest --noweb --hide-covar --out LRRK2 --burden cmc --kernel skato \
--inVcf /data/CARD/PD/AMP_NIH/no_relateds/subset_genetic_data_annovar/AMP_NIH_noRelateds_ALL_MISSENSE.vcf.gz \
--pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--pheno-name PD_PHENO \
--covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--freqUpper 0.01 \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
# CMC: 0.0033884
# SKATO: 1.95963e-07
```
*[Will add the sanity checks for AAO later]*

### 6b. Run burden (risk)
```
cd /data/CARD/PD/AMP_NIH/no_relateds/burden_annovar

## First make a new gene file with all chromosomes
R
require(dplyr)
require(data.table)
## Combine the genefiles for all chromosomes into one
genefile.names <- dir("/data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT",pattern="refFlat_HG38_chr.*.txt",full.names=TRUE)
names <- c("geneName","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")
# Read in the gene files per chromosome
genefile.dfs <- lapply(genefile.names, fread, col.names=names)
# Combine the gene files of all chromosomes
combined_genefile <- bind_rows(genefile.dfs)
write.table(combined_genefile, file="refFlat_HG38_all_chr.txt", quote=FALSE,row.names=F,sep="\t",col.names=FALSE)
q()
n

ls /data/CARD/PD/AMP_NIH/no_relateds/subset_genetic_data_annovar/AMP*vcf.gz | sed 's@.*/@@' > vcf_files_for_annovar_burden.txt

cat vcf_files_for_annovar_burden.txt
# AMP_NIH_noRelateds_ALL_CADD_10.vcf.gz
# AMP_NIH_noRelateds_ALL_CADD_20_and_LOF.vcf.gz
# AMP_NIH_noRelateds_ALL_CADD_20.vcf.gz
# AMP_NIH_noRelateds_ALL_LOF.vcf.gz
# AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF.vcf.gz
# AMP_NIH_noRelateds_ALL_MISSENSE.vcf.gz

cat vcf_files_for_annovar_burden.txt | while read line
do
    sbatch --cpus-per-task=4 --mem=10g --mail-type=FAIL --time=24:00:00 run_burden_annovar.sh $line 0.05
    sbatch --cpus-per-task=4 --mem=10g --mail-type=FAIL --time=24:00:00 run_burden_annovar.sh $line 0.01
    sbatch --cpus-per-task=4 --mem=10g --mail-type=FAIL --time=24:00:00 run_burden_annovar.sh $line 0.005
    sbatch --cpus-per-task=4 --mem=10g --mail-type=FAIL --time=24:00:00 run_burden_annovar.sh $line 0.001
done

#!/bin/sh
# sh run_burden_annovar.sh AMP_NIH_noRelateds_ALL_LOF.vcf.gz 0.05
module load rvtests
FILENAME=$1
MAF=$2
OUTNAME=${FILENAME/".vcf.gz"/""}
rvtest --noweb --hide-covar --out ${OUTNAME}_freqUpper${MAF} --burden cmc --kernel skato \
--inVcf /data/CARD/PD/AMP_NIH/no_relateds/subset_genetic_data_annovar/${FILENAME} \
--pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--pheno-name PD_PHENO \
--covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--freqUpper $MAF \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
--geneFile refFlat_HG38_all_chr.txt
```
```
## Analyze results

ls AMP*SkatO.assoc | wc -l
# 24
 --> 24 / 4 frequency cutoffs = 6 (makes sense)

ls AMP*SkatO.assoc | while read line
do
    sed "s/$/ $line/" $line >> all_SkatO.txt
done

awk '$8 < 5e-8'  all_SkatO.txt
# LRRK2	12:40224996-40369284	7974	52	52	168787	0	2.1557e-10 AMP_NIH_noRelateds_ALL_CADD_10_freqUpper0.005.SkatO.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	17	17	657326	0.3	2.15194e-10 AMP_NIH_noRelateds_ALL_CADD_10_freqUpper0.05.SkatO.assoc
# LRRK2	12:40224996-40369284	7974	33	33	164825	0	2.15428e-10 AMP_NIH_noRelateds_ALL_CADD_20_and_LOF_freqUpper0.005.SkatO.assoc
# LRRK2	12:40224996-40369284	7974	33	33	164825	0	2.15428e-10 AMP_NIH_noRelateds_ALL_CADD_20_freqUpper0.005.SkatO.assoc
# LRRK2	12:40224996-40369284	7974	52	52	166746	0	2.15788e-10 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_freqUpper0.005.SkatO.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	562974	0.2	2.1519e-10 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_freqUpper0.05.SkatO.assoc
# LRRK2	12:40224996-40369284	7974	51	51	166618	0	2.15752e-10 AMP_NIH_noRelateds_ALL_MISSENSE_freqUpper0.005.SkatO.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	562974	0.2	2.1519e-10 AMP_NIH_noRelateds_ALL_MISSENSE_freqUpper0.05.SkatO.assoc

ls AMP*CMC.assoc | while read line
do
    sed "s/$/ $line/" $line >> all_CMC.txt
done

awk '$7 < 5e-8' all_CMC.txt
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	17	17	524	4.86596e-14 AMP_NIH_noRelateds_ALL_CADD_10_freqUpper0.05.CMC.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	530	2.29221e-13 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_freqUpper0.05.CMC.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	530	2.29221e-13 AMP_NIH_noRelateds_ALL_MISSENSE_freqUpper0.05.CMC.assoc
```

### 6c. Run burden (AAO) [to be completed later]

## 7. Run burdens with snpEff/loftee annotations

### 7a. Sanity checks with GBA and LRRK2
```
cd /data/CARD/PD/AMP_NIH/no_relateds
mkdir burden_snpEff_loftee
cd burden_snpEff_loftee

module load rvtests

### PD risk

# GBA
rvtest --noweb --hide-covar --out GBA --burden cmc --kernel skato \
--inVcf /data/CARD/PD/AMP_NIH/no_relateds/subset_genetic_data_snpEff_loftee/AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE.vcf.gz \
--pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--pheno-name PD_PHENO \
--covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--freqUpper 0.05 \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# CMC: 2.29221e-13
# SKATO: 2.1519e-10

# LRRK2
rvtest --noweb --hide-covar --out LRRK2 --burden cmc --kernel skato \
--inVcf /data/CARD/PD/AMP_NIH/no_relateds/subset_genetic_data_snpEff_loftee/AMP_NIH_noRelateds_ALL_MISSENSE_SNPEFF.vcf.gz \
--pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--pheno-name PD_PHENO \
--covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--freqUpper 0.01 \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
# CMC: 0.00243298
# SKATO: 1.96111e-07
```
*[Will add the sanity checks for AAO later]*

### 7b. Run burden (risk)
```
cd /data/CARD/PD/AMP_NIH/no_relateds/burden_snpEff_loftee

ls /data/CARD/PD/AMP_NIH/no_relateds/subset_genetic_data_snpEff_loftee/AMP*vcf.gz | sed 's@.*/@@' > vcf_files_for_snpEff_loftee_burden.txt

cat vcf_files_for_snpEff_loftee_burden.txt
# AMP_NIH_noRelateds_ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE.vcf.gz
# AMP_NIH_noRelateds_ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE.vcf.gz
# AMP_NIH_noRelateds_ALL_HIGH_IMPACT_and_LOF_SNPEFF.vcf.gz
# AMP_NIH_noRelateds_ALL_HIGH_IMPACT_SNPEFF.vcf.gz
# AMP_NIH_noRelateds_ALL_LOF_and_HC_LOFTEE.vcf.gz
# AMP_NIH_noRelateds_ALL_LOF_HC_LOFTEE.vcf.gz
# AMP_NIH_noRelateds_ALL_LOF_SNPEFF.vcf.gz
# AMP_NIH_noRelateds_ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE.vcf.gz
# AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE.vcf.gz
# AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF.vcf.gz
# AMP_NIH_noRelateds_ALL_MISSENSE_SNPEFF.vcf.gz
# AMP_NIH_noRelateds_ALL_MODERATE_HIGH_IMPACT_SNPEFF.vcf.gz

cat vcf_files_for_snpEff_loftee_burden.txt | while read line
do
    sbatch --cpus-per-task=4 --mem=10g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.05
    sbatch --cpus-per-task=4 --mem=10g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.01
    sbatch --cpus-per-task=4 --mem=10g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.005
    sbatch --cpus-per-task=4 --mem=10g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.001
done

#!/bin/sh
# sh run_burden_snpeff_loftee.sh AMP_NIH_noRelateds_ALL_LOF_SNPEFF.vcf.gz 0.05
module load rvtests
FILENAME=$1
MAF=$2
OUTNAME=${FILENAME/".vcf.gz"/""}
rvtest --noweb --hide-covar --out ${OUTNAME}_freqUpper${MAF} --burden cmc --kernel skato \
--inVcf /data/CARD/PD/AMP_NIH/no_relateds/subset_genetic_data_snpEff_loftee/${FILENAME} \
--pheno /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--pheno-name PD_PHENO \
--covar /data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt \
--freqUpper $MAF \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/PD/AMP_NIH/no_relateds/burden_annovar/refFlat_HG38_all_chr.txt
```
```
## Analyze results

ls AMP*SkatO.assoc | wc -l
# 48
 --> 48 / 4 frequency cutoffs = 12 (makes sense)

ls AMP*SkatO.assoc | while read line
do
    sed "s/$/ $line/" $line >> all_SkatO.txt
done

awk '$8 < 5e-8'  all_SkatO.txt
# SLC50A1	1:155135374-155138853,1:155135374-155138853,1:155135374-155138853,1:155135811-155138858,1:155135811-155138858,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853	7974	2	2	287131	0	1.69865e-09 AMP_NIH_noRelateds_ALL_HIGH_IMPACT_and_LOF_SNPEFF_freqUpper0.05.SkatO.assoc
# SLC50A1	1:155135374-155138853,1:155135374-155138853,1:155135374-155138853,1:155135811-155138858,1:155135811-155138858,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853	7974	1	1	287116	0	1.89958e-09 AMP_NIH_noRelateds_ALL_LOF_SNPEFF_freqUpper0.05.SkatO.assoc
# LRRK2	12:40224996-40369284	7974	56	56	169164	0	2.15529e-10 AMP_NIH_noRelateds_ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE_freqUpper0.005.SkatO.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	562974	0.2	2.1519e-10 AMP_NIH_noRelateds_ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE_freqUpper0.05.SkatO.assoc
# LRRK2	12:40224996-40369284	7974	56	56	169164	0	2.15529e-10 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE_freqUpper0.005.SkatO.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	562974	0.2	2.1519e-10 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE_freqUpper0.05.SkatO.assoc
# LRRK2	12:40224996-40369284	7974	65	65	172557	0	2.15673e-10 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF_freqUpper0.005.SkatO.assoc
# SLC50A1	1:155135374-155138853,1:155135374-155138853,1:155135374-155138853,1:155135811-155138858,1:155135811-155138858,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853	7974	4	4	295417	0.2	1.69191e-09 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF_freqUpper0.05.SkatO.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	21	21	581071	0.2	2.15192e-10 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF_freqUpper0.05.SkatO.assoc
# LRRK2	12:40224996-40369284	7974	55	55	169037	0	2.15543e-10 AMP_NIH_noRelateds_ALL_MISSENSE_SNPEFF_freqUpper0.005.SkatO.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	562974	0.2	2.1519e-10 AMP_NIH_noRelateds_ALL_MISSENSE_SNPEFF_freqUpper0.05.SkatO.assoc
# LRRK2	12:40224996-40369284	7974	57	57	170841	0	2.15934e-10 AMP_NIH_noRelateds_ALL_MODERATE_HIGH_IMPACT_SNPEFF_freqUpper0.005.SkatO.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	562974	0.2	2.1519e-10 AMP_NIH_noRelateds_ALL_MODERATE_HIGH_IMPACT_SNPEFF_freqUpper0.05.SkatO.assoc

ls AMP*CMC.assoc | while read line
do
    sed "s/$/ $line/" $line >> all_CMC.txt
done

awk '$7 < 5e-8' all_CMC.txt
# SLC50A1	1:155135374-155138853,1:155135374-155138853,1:155135374-155138853,1:155135811-155138858,1:155135811-155138858,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853	7974	2	2	241	2.90449e-09 AMP_NIH_noRelateds_ALL_HIGH_IMPACT_and_LOF_SNPEFF_freqUpper0.05.CMC.assoc
# SLC50A1	1:155135374-155138853,1:155135374-155138853,1:155135374-155138853,1:155135811-155138858,1:155135811-155138858,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853	7974	1	1	240	2.27067e-09 AMP_NIH_noRelateds_ALL_LOF_SNPEFF_freqUpper0.05.CMC.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	530	2.29221e-13 AMP_NIH_noRelateds_ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE_freqUpper0.05.CMC.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	530	2.29221e-13 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE_freqUpper0.05.CMC.assoc
# SLC50A1	1:155135374-155138853,1:155135374-155138853,1:155135374-155138853,1:155135811-155138858,1:155135811-155138858,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853,1:155135851-155138853	7974	4	4	268	7.77553e-09 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF_freqUpper0.05.CMC.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	21	21	539	8.06281e-14 AMP_NIH_noRelateds_ALL_MISSENSE_and_LOF_SNPEFF_freqUpper0.05.CMC.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	530	2.29221e-13 AMP_NIH_noRelateds_ALL_MISSENSE_SNPEFF_freqUpper0.05.CMC.assoc
# GBA	1:155234451-155241249,1:155234451-155244627,1:155234451-155244627,1:155234451-155244627,1:155234451-155241249	7974	18	18	530	2.29221e-13 AMP_NIH_noRelateds_ALL_MODERATE_HIGH_IMPACT_SNPEFF_freqUpper0.05.CMC.assoc
```

### 7c. Run burden (AAO) [to be completed later]
