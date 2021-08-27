# UKB burdens for PD

---

## Structure of README:

### [1. Make covariate files with PCs for UKB Exomes](#1-make-covariate-files-with-pcs-for-ukb-exomes-1)
### [2. Run dbnsfp41a with annovar to get CADD scores](#2-run-dbnsfp41a-with-annovar-to-get-cadd-scores-1)
### [3. Subset variant classes](#3-subset-variant-classes-1)
### [4. Subset genetic data](#4-subset-genetic-data-1)
### [5. Run burdens with annovar annotations](#5-run-burdens-with-annovar-annotations-1)
### [6. Run burdens with snpEff/loftee annotations](#6-run-burdens-with-snpeffloftee-annotations-1)

---
### Here are some important files from this analysis:

- Working directory: /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis
- Covariate files for burdens: 
  - /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt
  - /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt
  - /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt
  - /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt
- Plink file merged with no relateds: /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass.*
- Annovar annotations:
  - /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar.txt
  - /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar_Y.txt
- snpEff/loftee annotations: /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/snpEff_loftee/*annotated.txt
- Annovar burden results: /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar
- snpEff/loftee burden results: /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_snpEff_loftee

---

## 1. Make covariate files with PCs for UKB Exomes

### Calculate PCs
``` 
cd /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups
ls UKB_EXOM*PD*2021.txt > need_PCs_EXOM.txt
 
cat need_PCs_EXOM.txt 
# UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021.txt
# UKB_EXOM_PD_CASE_CONTROL_2021.txt
# UKB_EXOM_PD_PARENT_CONTROL_2021.txt
# UKB_EXOM_PD_SIBLING_CONTROL_2021.txt
 
cat need_PCs_EXOM.txt | while read line
do
    sbatch --cpus-per-task=2 --mem=20g --mail-type=ALL --time=24:00:00 make_PC_flash_UKB_EXOM_2021.sh $line
done
 
#!/bin/sh
# sh make_PC_flash_UKB_EXOM_2021.sh UKB_EXOM_PD_CASE_CONTROL_2021.txt
module load flashpca
module load plink/2.0-dev-20191128
FILENAME=$1
OUTNAME=${FILENAME%%.*}_PCA
tmpfile=$(mktemp)
cut -f1 $FILENAME | tail -n+2  > ${tmpfile}
plink2 --pfile /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass --memory 195000 \
--maf 0.05 --geno 0.01 --hwe 5e-6 --autosome --keep ${tmpfile} \
--pheno-name PHENO --make-bed --out ${OUTNAME}_2 
module load plink
plink --bfile ${OUTNAME}_2 --indep-pairwise 1000 10 0.02 --autosome --out ${OUTNAME}_pruned_data --memory 195000
plink --bfile ${OUTNAME}_2 --extract ${OUTNAME}_pruned_data.prune.in --make-bed --out ${OUTNAME}_3 --memory 195000
flashpca --bfile ${OUTNAME}_3 --suffix _"$OUTNAME" --numthreads 19
```

### Merge PCs with phenotype information to make final covariate files
```
cd /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/
module load R
R
# Import PCs
EXOM_ALL_PD_pcs <- read.table("pcs_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_PCA",header=T)
EXOM_PD_CASE_pcs <- read.table("pcs_UKB_EXOM_PD_CASE_CONTROL_2021_PCA",header=T)
EXOM_PD_PARENT_pcs <- read.table("pcs_UKB_EXOM_PD_PARENT_CONTROL_2021_PCA",header=T)
EXOM_PD_SIBLING_pcs <- read.table("pcs_UKB_EXOM_PD_SIBLING_CONTROL_2021_PCA",header=T)
 
# Remove FID so not duplicated when merged
EXOM_ALL_PD_pcs$FID <- NULL
EXOM_PD_CASE_pcs$FID <- NULL
EXOM_PD_PARENT_pcs$FID <- NULL
EXOM_PD_SIBLING_pcs$FID <- NULL
 
# Import other covariates
EXOM_ALL_PD_cov <- read.table("UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021.txt",header=T)
EXOM_PD_CASE_cov <- read.table("UKB_EXOM_PD_CASE_CONTROL_2021.txt",header=T)
EXOM_PD_PARENT_cov <- read.table("UKB_EXOM_PD_PARENT_CONTROL_2021.txt",header=T)
EXOM_PD_SIBLING_cov <- read.table("UKB_EXOM_PD_SIBLING_CONTROL_2021.txt",header=T)
 
# Merge PCs with other covariates
MM1 <- merge(EXOM_ALL_PD_cov, EXOM_ALL_PD_pcs, by.x="FID", by.y="IID")
MM2 <- merge(EXOM_PD_CASE_cov, EXOM_PD_CASE_pcs,  by.x="FID", by.y="IID")
MM3 <- merge(EXOM_PD_PARENT_cov, EXOM_PD_PARENT_pcs,  by.x="FID", by.y="IID")
MM4 <- merge(EXOM_PD_SIBLING_cov, EXOM_PD_SIBLING_pcs,  by.x="FID", by.y="IID")
 
# Make sure no duplicates (there are none)
MM11 <- MM1[!duplicated(MM1), ]
MM22 <- MM2[!duplicated(MM2), ]
MM33 <- MM3[!duplicated(MM3), ]
MM44 <- MM4[!duplicated(MM4), ]
 
# Check PHENO column
unique(MM11$PHENO)
# 1 2
 
# Save
write.table(MM11, file="UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt", quote=FALSE,row.names=F,sep="\t")
write.table(MM22, file="UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt", quote=FALSE,row.names=F,sep="\t")
write.table(MM33, file="UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt", quote=FALSE,row.names=F,sep="\t")
write.table(MM44, file="UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt", quote=FALSE,row.names=F,sep="\t")
q()
n
```

## 2. Run dbnsfp41a with annovar to get CADD scores
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/
mkdir VCF_annotation_pathoscores
cd VCF_annotation_pathoscores

cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/PVCF_TO_LOOP_OVER.txt | while read line
do
   sbatch --cpus-per-task=10 --mem=10g --time=1:00:00 annotate_UKB_200K_pathoscores_July2021.sh $line
done

ls patho_ukb23156_c*v1.hg38_multianno.txt | wc -l
# 977

cat patho_ukb23156_c*v1.hg38_multianno.txt > MERGED_annotation.txt

wc -l MERGED_annotation.txt
# 17982874

R
require(data.table)
anno <- fread("MERGED_annotation.txt",header=T)
anno_Y <- fread("patho_ukb23156_cY_b0_v1.hg38_multianno.txt",header=T)
plink <- fread("/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_PLINK_with_header.txt",header=T)
plink_Y <- fread("/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_PLINK_Y_with_header.txt",header=T)
MM <- merge(plink,anno,by.x="ID",by.y="Otherinfo6")
MM_Y <- merge(plink_Y,anno_Y,by.x="ID",by.y="Otherinfo6")
write.table(MM, file="MERGED_annotation_added_pvar.txt",quote=F,row.names=F,sep="\t")
write.table(MM_Y, file="MERGED_annotation_added_pvar_Y.txt",quote=F,row.names=F,sep="\t")
q()
n

#!/bin/bash
# sbatch --cpus-per-task=10 --mem=10g --time=1:00:00 annotate_UKB_200K_pathoscores_July2021.sh ukb23156_c10_b15_v1
# load packages
module load samtools
module load annovar
module load bcftools/1.9
# variables
FILE=$1
echo "this is"
echo $FILE 
echo "starting annovar"
table_annovar.pl /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/filtered_"$FILE".for.anno.names.vcf $ANNOVAR_DATA/hg38 --thread 16 -buildver hg38 \
-out patho_"$FILE" -remove --otherinfo -polish -protocol dbnsfp41a \
-operation f -nastring . -vcfinput
echo "done"
```

## 3. Subset variant classes

### 3a. Subset variant classes using annovar annotations
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES
mkdir burden_analysis
cd burden_analysis

## Subset variant categories for burden analysis
mkdir subset_variant_categories_annovar
cd subset_variant_categories_annovar

# Merged annotations here
/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar.txt
/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar_Y.txt

wc -l /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar.txt
# 16645959
wc -l /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar_Y.txt
# 6634
```

### Subset "groups" of variants
```
# Missense
grep exonic /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar.txt | grep nonsynonymous | cut -f1 > ALL_MISSENSE.txt
# n=4605401

# LOF (stop, frame)
grep stopgain /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar.txt | cut -f1 > all_stopgain.txt
# n=158558
grep stoploss /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar.txt | cut -f1 > all_stoploss.txt
# n=6895
grep nonframeshift /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar.txt | cut -f1 > all_nonframeshift.txt
# n=86648
grep frame /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar.txt | grep -v nonframeshift | grep -v nonsynonymous | cut -f1 > all_frameshift.txt
# n=223981

# Splicing 
grep splicing /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/MERGED_annotation_added_pvar.txt | grep -v ncRNA | cut -f1 > all_splice.txt
# n=91192

# CADD > 10
awk '{ if($38 >= 10) { print }}' /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation_pathoscores/MERGED_annotation_added_pvar.txt | cut -f1 > ALL_CADD_10.txt
# n=4114318

# CADD > 20
awk '{ if($38 >= 20) { print }}' /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation_pathoscores/MERGED_annotation_added_pvar.txt | cut -f1 > ALL_CADD_20.txt
# n=3060317

#### Prepping final files:
cat all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_LOF.txt
cat ALL_MISSENSE.txt all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_MISSENSE_and_LOF.txt
cat ALL_CADD_20.txt all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_CADD_20_and_LOF.txt

wc -l ALL_LOF.txt
# 480626
wc -l ALL_MISSENSE_and_LOF.txt
# 5086027
wc -l ALL_CADD_20_and_LOF.txt
# 3540943

ls *txt | while read line
do
    sed -i 's/_/:/g' $line
done
```

### Use these variant categories for burdens:
```
1. ALL_CADD_10.txt
2. ALL_CADD_20.txt
3. ALL_MISSENSE_and_LOF.txt
4. ALL_MISSENSE.txt
5. ALL_LOF.txt
6. ALL_CADD_20_and_LOF.txt
```

### 3b. Subset variant classes using snpeff/loftee

### Reformat the annotations
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis
mkdir subset_variant_categories_snpeff_loftee
cd subset_variant_categories_snpeff_loftee

ls /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/snpEff_loftee/subset_UKB_forSNPEFF.filtered_ukb23156_c*v1.annotated.txt | wc -l
# 977

# First subset the header only
head -1 /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/snpEff_loftee/subset_UKB_forSNPEFF.filtered_ukb23156_c6_b1_v1.annotated.txt  > header_snpEff.txt

# Concat all of the SnpEff annotations across chromosomes
cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/snpEff_loftee/subset_UKB_forSNPEFF.filtered_ukb23156_c*v1.annotated.txt | grep -v CHROM > temp1
cat header_snpEff.txt temp1 > subset_UKB_forSNPEFF_all_chr.annotated.txt
wc -l subset_UKB_forSNPEFF_all_chr.annotated.txt
# 17981898

# # Concat all of the Loftee high confidence annotations across chromosomes
cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/snpEff_loftee/subset_UKB_forLOFTEE.filtered_ukb23156_c*v1.refined_LoF.snps > ALL_LOF_HC_LOFTEE.txt
wc -l ALL_LOF_HC_LOFTEE.txt
# 147997
```

### Subset "groups" of variants
```
# Missense
grep missense_variant subset_UKB_forSNPEFF_all_chr.annotated.txt | cut -f3 > ALL_MISSENSE_SNPEFF.txt
# n=4807373

# LOF (stop, frame)
grep stop_gained subset_UKB_forSNPEFF_all_chr.annotated.txt | cut -f3 > all_stopgain.txt
# n=159504

grep stop_lost subset_UKB_forSNPEFF_all_chr.annotated.txt | cut -f3 > all_stoploss.txt
# n=11115

grep frameshift_variant subset_UKB_forSNPEFF_all_chr.annotated.txt | cut -f3 > all_frameshift.txt
# n=227453

# Splicing 
grep splice_region_variant subset_UKB_forSNPEFF_all_chr.annotated.txt | cut -f3 > all_splice.txt
# n=784832

# Moderate or high impact
grep -e MODERATE -e HIGH subset_UKB_forSNPEFF_all_chr.annotated.txt | cut -f3 > ALL_MODERATE_HIGH_IMPACT_SNPEFF.txt
# n=5523869

# High impact
grep HIGH subset_UKB_forSNPEFF_all_chr.annotated.txt | cut -f3 > ALL_HIGH_IMPACT_SNPEFF.txt
# n=660085

#### Prepping final files:
cat all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_LOF_SNPEFF.txt

# How many LOF determined by snpEff
wc -l ALL_LOF_SNPEFF.txt
# 1182904

# How many LOF determined by loftee with high confidence
wc -l ALL_LOF_HC_LOFTEE.txt
# 147997

# How many SNPs are in both the loftee high confidence and SnpEff LOF
comm -12 <(sort ALL_LOF_SNPEFF.txt) <(sort ALL_LOF_HC_LOFTEE.txt) > ALL_LOF_and_HC_LOFTEE.txt
# n=117907

# How many SNPs were identified by SnpEff but weren't in high confidence loftee
comm -23 <(sort ALL_LOF_SNPEFF.txt) <(sort ALL_LOF_HC_LOFTEE.txt) > ALL_LOF_not_HC_LOFTEE.txt
# n=1064997

# How many SNPs identified by loftee high confidence aren't in the SnpEff LOF 
comm -13 <(sort ALL_LOF_SNPEFF.txt) <(sort ALL_LOF_HC_LOFTEE.txt) > HC_LOFTEE_not_ALL_LOF.txt
# n=30090

# Concat LOF with missense
cat ALL_MISSENSE_SNPEFF.txt ALL_LOF_SNPEFF.txt > ALL_MISSENSE_and_LOF_SNPEFF.txt
cat ALL_MISSENSE_SNPEFF.txt ALL_LOF_HC_LOFTEE.txt> ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE.txt
cat ALL_MISSENSE_SNPEFF.txt ALL_LOF_and_HC_LOFTEE.txt > ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE.txt

# Concat LOF with high impact variants
cat ALL_HIGH_IMPACT_SNPEFF.txt ALL_LOF_SNPEFF.txt > ALL_HIGH_IMPACT_and_LOF_SNPEFF.txt
cat ALL_HIGH_IMPACT_SNPEFF.txt ALL_LOF_HC_LOFTEE.txt> ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE.txt
cat ALL_HIGH_IMPACT_SNPEFF.txt ALL_LOF_and_HC_LOFTEE.txt > ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE.txt

wc -l ALL_MISSENSE_and_LOF_SNPEFF.txt
# 5990277
wc -l ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE.txt
# 4955370
wc -l ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE.txt
# 4925280
wc -l ALL_HIGH_IMPACT_and_LOF_SNPEFF.txt
# 1842989
wc -l ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE.txt
# 808082
wc -l ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE.txt
# 777992

ls ALL*txt all*txt | while read line
do
    sed -i 's/_/:/g' $line
done
```

### Use these variant categories for burdens:
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

## 4. Subset genetic data

### 4a. Subset genetic data using annovar annotations
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis
mkdir subset_genetic_data_annovar
cd subset_genetic_data_annovar

# Make vcf files for input for rvtests
ls /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_*_2021_with_PC.txt | sed 's@.*/@@' > EXOM_sample_selection.txt

cat EXOM_sample_selection.txt
# UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt
# UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt
# UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt
# UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt

cat EXOM_sample_selection.txt | while read line
do
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=1:00:00 make_vcf_for_burden.sh ALL_MISSENSE.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=1:00:00 make_vcf_for_burden.sh ALL_LOF.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=1:00:00 make_vcf_for_burden.sh ALL_CADD_20.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=1:00:00 make_vcf_for_burden.sh ALL_CADD_10.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=1:00:00 make_vcf_for_burden.sh ALL_MISSENSE_and_LOF.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=1:00:00 make_vcf_for_burden.sh ALL_CADD_20_and_LOF.txt $line
done

#!/bin/sh
# sh make_vcf_for_burden.sh ALL_MISSENSE.txt UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt
module load plink/2.0-dev-20191128
module load samtools
VARIANT_FILE=$1
SAMPLE_FILE=$2
OUTNAME=${SAMPLE_FILE/"with_PC.txt"/""}${VARIANT_FILE/".txt"/""}
tmpfile=$(mktemp)
cut -f1 /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/$SAMPLE_FILE | tail -n+2  > ${tmpfile}
for chnum in {1..22};
do
    plink2 --pfile /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass \
    --chr $chnum \
    --keep ${tmpfile} \
    --extract /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_variant_categories_annovar/$VARIANT_FILE \
    --export vcf bgz id-paste=iid --out ${OUTNAME}.chr${chnum} --mac 1
    tabix -p vcf  ${OUTNAME}.chr${chnum}.vcf.gz
done
# Also subset the X chromosome
plink2 --pfile /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_chromosome_X \
--keep ${tmpfile} \
--extract /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_variant_categories_annovar/$VARIANT_FILE \
--export vcf bgz id-paste=iid --out ${OUTNAME}.chr23 --mac 1
tabix -p vcf  ${OUTNAME}.chr23.vcf.gz

## We expect 552 vcfs (4 sample selection files X 23 chromosomes X 6 variant selection files)
ls UKB*vcf.gz | wc -l
# 552
```

### 4b. Subset genetic data using snpEff/loftee
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis
mkdir subset_genetic_data_snpeff_loftee
cd subset_genetic_data_snpeff_loftee

# Make vcf files for input for rvtests

cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_annovar/EXOM_sample_selection.txt | while read line
do
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_MODERATE_HIGH_IMPACT_SNPEFF.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_HIGH_IMPACT_SNPEFF.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_MISSENSE_and_LOF_SNPEFF.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_MISSENSE_SNPEFF.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_LOF_SNPEFF.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_LOF_HC_LOFTEE.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_LOF_and_HC_LOFTEE.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_HIGH_IMPACT_and_LOF_SNPEFF.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE.txt $line
    sbatch --cpus-per-task=7 --mem=5g --mail-type=FAIL --time=2:00:00 make_vcf_for_burden_snpEff_loftee.sh ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE.txt $line
done

#!/bin/sh
# sh make_vcf_for_burden_snpEff_loftee.sh ALL_MISSENSE_SNPEFF.txt UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt
module load plink/2.0-dev-20191128
module load samtools
VARIANT_FILE=$1
SAMPLE_FILE=$2
OUTNAME=${SAMPLE_FILE/"with_PC.txt"/""}${VARIANT_FILE/".txt"/""}
tmpfile=$(mktemp)
cut -f1 /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/$SAMPLE_FILE | tail -n+2  > ${tmpfile}
for chnum in {1..22};
do
    plink2 --pfile /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass \
    --chr $chnum \
    --keep ${tmpfile} \
    --extract /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_variant_categories_snpeff_loftee/$VARIANT_FILE \
    --export vcf bgz id-paste=iid --out ${OUTNAME}.chr${chnum} --mac 1
    tabix -p vcf  ${OUTNAME}.chr${chnum}.vcf.gz
done
# Also subset the X chromosome
plink2 --pfile /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_chromosome_X \
--keep ${tmpfile} \
--extract /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_variant_categories_snpeff_loftee/$VARIANT_FILE \
--export vcf bgz id-paste=iid --out ${OUTNAME}.chr23 --mac 1
tabix -p vcf  ${OUTNAME}.chr23.vcf.gz

## We expect 1104 vcfs (4 sample selection files X 23 chromosomes X 12 variant selection files)
ls UKB*gz | wc -l
# 1104
```

## 5. Run burdens with annovar annotations

### 5a. Sanity checks with GBA and LRRK2
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis
mkdir burden_annovar
cd burden_annovar

## Reformat the covariate files so that the pheno column (either PHENO for case-control or AGE_OF_RECRUIT)
ls /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_*_2021_with_PC.txt | sed 's@.*/@@'> reformat_cov_files.txt

cat reformat_cov_files.txt | while read line
do
    awk 'OFS="\t" { print $1,$2,$4,$5,$10,$13,$14,$15,$16,$17,$9 }' /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/$line > reformatted_${line}
    ***[Note: will need to redo this once AGE is changed]*** 
    awk 'NR==1 || $9 == 2'  /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/$line > temp
    awk 'OFS="\t" { print $1,$2,$4,$9,$10,$13,$14,$15,$16,$17,$5 }' temp > reformatted_AAO_${line}
done

module load rvtests

### PD risk

# GBA
rvtest --noweb --hide-covar --out GBA --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_annovar/UKB_EXOM_PD_CASE_CONTROL_2021_ALL_MISSENSE_and_LOF.chr1.vcf.gz \
--pheno reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--pheno-name PHENO \
--covar reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--freqUpper 0.05 \
--covar-name SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# CMC: 2.29815e-08
# SKATO: 7.01163e-08

# LRRK2
rvtest --noweb --hide-covar --out LRRK2 --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_annovar/UKB_EXOM_PD_CASE_CONTROL_2021_ALL_MISSENSE.chr12.vcf.gz \
--pheno reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--pheno-name PHENO \
--covar reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--freqUpper 0.01 \
--covar-name SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
# CMC: 0.237723
# SKATO: 0.414194
```

*[Will add the sanity checks for AAO later]*

### 5b. Run burden (risk)
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar

# Use the new covariate files with the phenotype column at the end, otherwise RVTESTS will fail
ls reformatted_UKB* | while read line
do
>     wc -l $line
done
# 45858 reformatted_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt
# 6749 reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt
# 34979 reformatted_UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt
# 4132 reformatted_UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt

# We expect 552 input VCFs…
6 variant categories x 4 phenotypes x 23 chromosomes = 552
ls /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_annovar/*vcf.gz | sed 's@.*/@@' > vcf_files_for_annovar_burden.txt
wc -l vcf_files_for_annovar_burden.txt
# 552

cat vcf_files_for_annovar_burden.txt | grep SIBLING | while read line
do
    sbatch --cpus-per-task=2 --mem=1g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=1g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=1g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=1g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.001
done

cat vcf_files_for_annovar_burden.txt | grep CASE_CONTROL | while read line
do
    sbatch --cpus-per-task=2 --mem=1g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=1g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=1g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=1g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.001
done

cat vcf_files_for_annovar_burden.txt | grep PARENT | while read line
do
    sbatch --cpus-per-task=2 --mem=8g --mail-type=FAIL --time=6:00:00 run_burden_annovar.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=8g --mail-type=FAIL --time=6:00:00 run_burden_annovar.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=8g --mail-type=FAIL --time=6:00:00 run_burden_annovar.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=8g --mail-type=FAIL --time=6:00:00 run_burden_annovar.sh $line 0.001
done

cat vcf_files_for_annovar_burden.txt | grep ALL_PD_PHENOTYPES | while read line
do
    sbatch --cpus-per-task=2 --mem=9g --mail-type=FAIL --time=10:00:00 run_burden_annovar.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=9g --mail-type=FAIL --time=10:00:00 run_burden_annovar.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=9g --mail-type=FAIL --time=10:00:00 run_burden_annovar.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=9g --mail-type=FAIL --time=10:00:00 run_burden_annovar.sh $line 0.001
done

# for cases and sibs
ls UKB*SkatO.assoc | wc -l
# 2208
# 552 * 4 = 2,208 

#!/bin/sh
# sh run_burden_annovar.sh UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_ALL_CADD_10.chr10.vcf.gz 0.05
module load rvtests
FILENAME=$1
MAF=$2
OUTNAME=${FILENAME/".vcf.gz"/""}
COV_NAME=reformatted_${FILENAME%_ALL*}_with_PC.txt
rvtest --noweb --hide-covar --out ${OUTNAME}_freqUpper${MAF} --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_annovar/${FILENAME} \
--pheno ${COV_NAME} \
--pheno-name PHENO \
--covar ${COV_NAME} \
--freqUpper $MAF \
--covar-name SEX,TOWNSEND,AGE_OF_RECRUIT,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/PD/AMP_NIH/no_relateds/burden_annovar/refFlat_HG38_all_chr.txt
```
```
## Analyze results

# See if any files are empty, meaning failed
ls -p | while read line
do
     [ -s $line ] || echo $line "is empty"
done
# none

***[ Note: add here the significant genes once finished ]***
```

### 5c. Run burden (AAO) [to be completed later]

## 6. Run burdens with snpEff/loftee annotations

### 6a. Sanity checks with GBA and LRRK2
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis
mkdir burden_snpEff_loftee
cd burden_snpEff_loftee

module load rvtests

### PD risk

# GBA
rvtest --noweb --hide-covar --out GBA --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_snpeff_loftee/UKB_EXOM_PD_CASE_CONTROL_2021_ALL_MISSENSE_and_LOF_SNPEFF.chr1.vcf.gz \
--pheno /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--pheno-name PHENO \
--covar /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--freqUpper 0.05 \
--covar-name SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# CMC: 2.12215e-08
# SKATO: 1.28641e-08

# LRRK2
rvtest --noweb --hide-covar --out LRRK2 --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_snpeff_loftee/UKB_EXOM_PD_CASE_CONTROL_2021_ALL_MISSENSE_SNPEFF.chr12.vcf.gz \
--pheno /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--pheno-name PHENO \
--covar /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--freqUpper 0.01 \
--covar-name SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
# CMC: 0.261744
# SKATO: 0.419593
```

*[Will add the sanity checks for AAO later]*

### 6b. Run burden (risk)

[will add here once finished]

### 6c. Run burden (AAO) [to be completed later]

