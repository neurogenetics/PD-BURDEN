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
--covar-name SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# CMC: .17592e-09
# SKATO: 5.58925e-08

# LRRK2
rvtest --noweb --hide-covar --out LRRK2 --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_annovar/UKB_EXOM_PD_CASE_CONTROL_2021_ALL_MISSENSE.chr12.vcf.gz \
--pheno reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--pheno-name PHENO \
--covar reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--freqUpper 0.01 \
--covar-name SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
# CMC: 0.194985
# SKATO: 0.337583

```

*[Will add the sanity checks for AAO later]*

### 5b. Run burden (risk)
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar

# Use the new covariate files with the phenotype column at the end, otherwise RVTESTS will fail
ls reformatted_UKB* | while read line
do
    wc -l $line
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
    sbatch --cpus-per-task=2 --mem=2g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=2g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=2g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=2g --mail-type=FAIL --time=1:00:00 run_burden_annovar.sh $line 0.001
done

cat vcf_files_for_annovar_burden.txt | grep CASE_CONTROL | while read line
do
    sbatch --cpus-per-task=2 --mem=3g --mail-type=FAIL --time=2:00:00 run_burden_annovar.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=3g --mail-type=FAIL --time=2:00:00 run_burden_annovar.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=3g --mail-type=FAIL --time=2:00:00 run_burden_annovar.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=3g --mail-type=FAIL --time=2:00:00 run_burden_annovar.sh $line 0.001
done

cat vcf_files_for_annovar_burden.txt | grep PARENT | while read line
do
    sbatch --cpus-per-task=2 --mem=12g --mail-type=FAIL --time=12:00:00 run_burden_annovar.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=12g --mail-type=FAIL --time=12:00:00 run_burden_annovar.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=12g --mail-type=FAIL --time=12:00:00 run_burden_annovar.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=12g --mail-type=FAIL --time=12:00:00 run_burden_annovar.sh $line 0.001
done

cat vcf_files_for_annovar_burden.txt | grep ALL_PD_PHENOTYPES | while read line
do
    sbatch --cpus-per-task=2 --mem=16g --mail-type=FAIL --time=24:00:00 run_burden_annovar.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=16g --mail-type=FAIL --time=24:00:00 run_burden_annovar.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=16g --mail-type=FAIL --time=24:00:00 run_burden_annovar.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=16g --mail-type=FAIL --time=24:00:00 run_burden_annovar.sh $line 0.001
done

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
--covar-name SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/PD/AMP_NIH/no_relateds/burden_annovar/refFlat_HG38_all_chr.txt
```

### Troubleshoot the gene that is "stuck"

Problem: RVTESTS sometimes fails to move past a gene that should be skipped, causing the burden to run indefinitely

Solution: Loop through each gene for the burdnes that failed to isolate the genes that should have been skipped

The corresponding burden commands that failed are: 
- sh run_burden_annovar.sh UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5.vcf.gz 0.001
- sh run_burden_annovar.sh UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5.vcf.gz 0.05
- sh run_burden_annovar.sh UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5.vcf.gz 0.005
- sh run_burden_annovar.sh UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5.vcf.gz 0.01

```
## Confirm these are the only burdens that failed

# See if any files are empty
ls -p | while read line
do
     [ -s $line ] || echo $line "is empty"
done

# UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.001.CMC.assoc is empty
# UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.001.SkatO.assoc is empty
# UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.005.CMC.assoc is empty
# UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.005.SkatO.assoc is empty
# UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.01.CMC.assoc is empty
# UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.01.SkatO.assoc is empty
# UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.05.CMC.assoc is empty
# UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.05.SkatO.assoc is empty
```

### Gene loop for chromosome 5

Condition we need to run per gene:
- Chromosome 5
- MAF cutoffs 0.05, 0.01, 0.005, 0.001
- UKB case-control 
- Annovar all CADD > 20
  
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis
mkdir gene_loop_chr5
cd gene_loop_chr5

cut -f1 /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr5.txt | sort -u > chr5_gene_list.txt

wc -l chr5_gene_list.txt
# 1314

## Split up the gene list into 2 so that the swarm is only grouped by gene and can therefore identify failed genes (max is 1000 subjobs)
head -n 657 chr5_gene_list.txt > chr5_gene_list_part1.txt
tail -n 657 chr5_gene_list.txt > chr5_gene_list_part2.txt

wc -l chr5_gene_list_part1.txt
# 657
wc -l chr5_gene_list_part2.txt
# 657

cat chr5_gene_list_part1.txt| while read GENE
do
    echo "sh run_burden_by_gene.sh $GENE 0.05" >> burden_by_gene_part1.swarm
    echo "sh run_burden_by_gene.sh $GENE 0.01" >> burden_by_gene_part1.swarm
    echo "sh run_burden_by_gene.sh $GENE 0.005" >> burden_by_gene_part1.swarm
    echo "sh run_burden_by_gene.sh $GENE 0.001" >> burden_by_gene_part1.swarm
done

cat chr5_gene_list_part2.txt| while read GENE
do
    echo "sh run_burden_by_gene.sh $GENE 0.05" >> burden_by_gene_part2.swarm
    echo "sh run_burden_by_gene.sh $GENE 0.01" >> burden_by_gene_part2.swarm
    echo "sh run_burden_by_gene.sh $GENE 0.005" >> burden_by_gene_part2.swarm
    echo "sh run_burden_by_gene.sh $GENE 0.001" >> burden_by_gene_part2.swarm
done

wc -l burden_by_gene_part1.swarm 
# 2628
wc -l burden_by_gene_part2.swarm 
# 2628

swarm -f burden_by_gene_part1.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module rvtests --bundle 4
# 2628 commands run in 657 subjobs, each command requiring 1.5 gb and 1 thread, running 4 processes serially per subjob, allocating 657 cores and 1314 cpus
# 22836805

swarm -f burden_by_gene_part2.swarm --verbose 1 --sbatch "--mail-type=FAIL" --module rvtests --bundle 4
# 2628 commands run in 657 subjobs, each command requiring 1.5 gb and 1 thread, running 4 processes serially per subjob, allocating 657 cores and 1314 cpus
# 22836813

# The only gene that failed --> 22836805_257 --> gene is FAM153A

## Check the genes that got stuck previously
# For example, LINC01574 and MIR4281

cat *LINC01574*assoc
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	NonRefSite	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	Q	rho	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	NonRefSite	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	Q	rho	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	NonRefSite	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	Q	rho	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	NonRefSite	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	Q	rho	Pvalue

cat *MIR4281*assoc
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	NonRefSite	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	Q	rho	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	NonRefSite	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	Q	rho	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	NonRefSite	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	Q	rho	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	NonRefSite	Pvalue
Gene	RANGE	N_INFORMATIVE	NumVar	NumPolyVar	Q	rho	Pvalue

## Double check all genes are completed except FAM153A 
cat chr5_gene_list.txt | while read GENE
do
    [ -s UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.05_${GENE}.log ] || echo $GENE "is empty"
    [ -s UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.01_${GENE}.log ] || echo $GENE "is empty"
    [ -s UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.005_${GENE}.log ] || echo $GENE "is empty"
    [ -s UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.001_${GENE}.log ] || echo $GENE "is empty"
done

# FAM153A is empty
# FAM153A is empty
# FAM153A is empty

ls UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper*.log | wc -l
# 5253 (makes sense --> 2628 * 2 = 5256  - 3 for FAM153A)

#!/bin/sh
# sh run_burden_by_gene.sh LRRK2 0.05
module load rvtests
GENE=$1
MAF=$2
rvtest --noweb --hide-covar --out UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper${MAF}_${GENE} --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_annovar/UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5.vcf.gz \
--pheno /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--pheno-name PHENO \
--covar /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--freqUpper $MAF \
--covar-name SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr5.txt --gene $GENE
```

### Reformat results
```
## Combine the chromosome 5 results per gene

cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/gene_loop_chr5

head -1 UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.05_ZSWIM6.SkatO.assoc > skat_header.txt
head -1 UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.05_ZSWIM6.CMC.assoc  > cmc_header.txt

## Combine genes
cat UKB*freqUpper0.05*SkatO.assoc | grep -v N_INFORMATIVE > temp
cat skat_header.txt temp > UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.05.SkatO.assoc
cat UKB*freqUpper0.01*SkatO.assoc | grep -v N_INFORMATIVE > temp
cat skat_header.txt temp > UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.01.SkatO.assoc
cat UKB*freqUpper0.005*SkatO.assoc | grep -v N_INFORMATIVE > temp
cat skat_header.txt temp > UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.005.SkatO.assoc
cat UKB*freqUpper0.001*SkatO.assoc | grep -v N_INFORMATIVE > temp
cat skat_header.txt temp > UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.001.SkatO.assoc

cat UKB*freqUpper0.05*CMC.assoc | grep -v N_INFORMATIVE > temp
cat cmc_header.txt temp > UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.05.CMC.assoc
cat UKB*freqUpper0.01*CMC.assoc | grep -v N_INFORMATIVE > temp
cat cmc_header.txt temp > UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.01.CMC.assoc
cat UKB*freqUpper0.005*CMC.assoc | grep -v N_INFORMATIVE > temp
cat cmc_header.txt temp > UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.005.CMC.assoc
cat UKB*freqUpper0.001*CMC.assoc | grep -v N_INFORMATIVE > temp
cat cmc_header.txt temp > UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper0.001.CMC.assoc

## Move to the directory with the other burden results
ls UKB_EXOM_PD_CASE_CONTROL_2021_ALL_CADD_20.chr5_freqUpper{0.05,0.01,0.005,0.001}.{SkatO,CMC}.assoc | while read line
do
    cp $line /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/
done

## Now check that nothing is empty in the final burdens directory
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar
ls -p | while read line
do
     [ -s $line ] || echo $line "is empty"
done
# none :)
```
```
## Combine chromosomes for meta-analysis

cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis
mkdir meta_prep_annovar
cd meta_prep_annovar

ls /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/UKB*SkatO.assoc | sed 's/.chr.*//g' | sort -u | sed 's@.*/@@' > file_prefixes_annovar.txt

wc -l file_prefixes_annovar.txt
# 24 (makes sense, 6 categories x 4 phenotypes)

head -1 /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/UKB_EXOM_PD_SIBLING_CONTROL_2021_ALL_MISSENSE.chr5_freqUpper0.05.SkatO.assoc > header_SkatO.txt
head -1 /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/UKB_EXOM_PD_SIBLING_CONTROL_2021_ALL_MISSENSE.chr5_freqUpper0.05.CMC.assoc > header_CMC.txt

cat file_prefixes_annovar.txt | while read line
do
     for type in "SkatO" "CMC";
     do
          cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/${line}.chr{1..23}_freqUpper0.05.${type}.assoc | grep -v N_INFORMATIVE > temp
          cat header_${type}.txt temp > ${line}.freqUpper0.05.${type}.assoc
          cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/${line}.chr{1..23}_freqUpper0.01.${type}.assoc | grep -v N_INFORMATIVE > temp
          cat header_${type}.txt temp > ${line}.freqUpper0.01.${type}.assoc
          cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/${line}.chr{1..23}_freqUpper0.005.${type}.assoc | grep -v N_INFORMATIVE > temp
          cat header_${type}.txt temp > ${line}.freqUpper0.005.${type}.assoc
          cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/${line}.chr{1..23}_freqUpper0.001.${type}.assoc | grep -v N_INFORMATIVE > temp
          cat header_${type}.txt temp > ${line}.freqUpper0.001.${type}.assoc
     done
done
```

### Analyze results
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/meta_prep_annovar

ls UKB*SkatO.assoc | wc -l
# 96 (makes sense, 6 variant categories x 4 maf categories x 4 phenotypes)

ls UKB*SkatO.assoc | while read line
do
    sed "s/$/ $line/" $line >> all_SkatO.txt
done

awk '$8 < 5e-8'  all_SkatO.txt | cut -f1 | sort -u
# AUNIP
# GBA
# HPS4

awk '$8 < 1e-6'  all_SkatO.txt | cut -f1 | sort -u
# ADH5
# AUNIP
# GBA
# HPS4
# KIF9
# OR1G1
# PDZRN4
# SEPT5-GP1BB
# ZNF454
```
```
ls UKB*CMC.assoc | wc -l
# 96

ls UKB*CMC.assoc | while read line
do
    sed "s/$/ $line/" $line >> all_CMC.txt
done

awk '$7 < 5e-8' all_CMC.txt | cut -f1 | sort -u
# GBA

awk '$7 < 1e-6' all_CMC.txt | cut -f1 | sort -u
# ADH5
# GBA
# KIF9
# NOXO1
# SEPT5-GP1BB
# TREML1
# ZNF454
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
--covar-name SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# CMC: 4.27527e-09
# SKATO: 2.68661e-09

# LRRK2
rvtest --noweb --hide-covar --out LRRK2 --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_snpeff_loftee/UKB_EXOM_PD_CASE_CONTROL_2021_ALL_MISSENSE_SNPEFF.chr12.vcf.gz \
--pheno /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--pheno-name PHENO \
--covar /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt \
--freqUpper 0.01 \
--covar-name SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
# CMC: 0.218849
# SKATO: 0.371834
```

*[Will add the sanity checks for AAO later]*

### 6b. Run burden (risk)
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_snpeff_loftee

# Expect 1104 input VCFs…
12 variant categories x 4 phenotypes x 23 chromosomes = 1104

ls /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_snpeff_loftee/*vcf.gz | sed 's@.*/@@' > vcf_files_for_snpeff_loftee_burden.txt
wc -l vcf_files_for_snpeff_loftee_burden.txt
# 1104

23 * 12 = 276 * 4 = 1,104 

cat vcf_files_for_snpeff_loftee_burden.txt | grep SIBLING | while read line
do
    sbatch --cpus-per-task=2 --mem=2g --mail-type=FAIL --time=2:00:00 run_burden_snpeff_loftee.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=2g --mail-type=FAIL --time=2:00:00 run_burden_snpeff_loftee.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=2g --mail-type=FAIL --time=2:00:00 run_burden_snpeff_loftee.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=2g --mail-type=FAIL --time=2:00:00 run_burden_snpeff_loftee.sh $line 0.001
done

cat vcf_files_for_snpeff_loftee_burden.txt | grep CASE_CONTROL | while read line
do
    sbatch --cpus-per-task=2 --mem=4g --mail-type=FAIL --time=5:00:00 run_burden_snpeff_loftee.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=4g --mail-type=FAIL --time=5:00:00 run_burden_snpeff_loftee.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=4g --mail-type=FAIL --time=5:00:00 run_burden_snpeff_loftee.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=4g --mail-type=FAIL --time=5:00:00 run_burden_snpeff_loftee.sh $line 0.001
done

cat vcf_files_for_snpeff_loftee_burden.txt | grep PARENT | while read line
do
    sbatch --cpus-per-task=2 --mem=16g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=16g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=16g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=16g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.001
done

cat vcf_files_for_snpeff_loftee_burden.txt | grep ALL_PD_PHENOTYPES | while read line
do
    sbatch --cpus-per-task=2 --mem=25g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.05
    sbatch --cpus-per-task=2 --mem=25g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.01
    sbatch --cpus-per-task=2 --mem=25g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.005
    sbatch --cpus-per-task=2 --mem=25g --mail-type=FAIL --time=24:00:00 run_burden_snpeff_loftee.sh $line 0.001
done

#!/bin/sh
# sh run_burden_snpeff_loftee.sh UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_ALL_MISSENSE_and_LOF_SNPEFF.chr1.vcf.gz 0.05
module load rvtests
FILENAME=$1
MAF=$2
OUTNAME=${FILENAME/".vcf.gz"/""}
COV_NAME=/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_${FILENAME%_2021*}_2021_with_PC.txt
rvtest --noweb --hide-covar --out ${OUTNAME}_freqUpper${MAF} --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_snpeff_loftee/${FILENAME} \
--pheno ${COV_NAME} \
--pheno-name PHENO \
--covar ${COV_NAME} \
--freqUpper $MAF \
--covar-name SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/PD/AMP_NIH/no_relateds/burden_annovar/refFlat_HG38_all_chr.txt
```

```
## Confirm that none of the burdens failed

ls -p | while read line
do
     [ -s $line ] || echo $line "is empty"
done
# none
```

### Reformat results
```
## Combine chromosomes for meta-analysis

cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis
mkdir meta_prep_snpeff_loftee
cd meta_prep_snpeff_loftee

ls /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_snpeff_loftee/UKB*SkatO.assoc | sed 's/.chr.*//g' | sort -u | sed 's@.*/@@' > file_prefixes_snpeff.txt

wc -l file_prefixes_snpeff.txt
# 48 (makes sense, 12 categories x 4 phenotypes)

head -1 /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_snpeff_loftee/UKB_EXOM_PD_SIBLING_CONTROL_2021_ALL_HIGH_IMPACT_SNPEFF.chr5_freqUpper0.05.SkatO.assoc > header_SkatO.txt
head -1 /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_snpeff_loftee/UKB_EXOM_PD_SIBLING_CONTROL_2021_ALL_HIGH_IMPACT_SNPEFF.chr5_freqUpper0.05.CMC.assoc > header_CMC.txt

cat file_prefixes_snpeff.txt | while read line
do
     for type in "SkatO" "CMC";
     do
          cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_snpeff_loftee/${line}.chr{1..23}_freqUpper0.05.${type}.assoc | grep -v N_INFORMATIVE > temp
          cat header_${type}.txt temp > ${line}.freqUpper0.05.${type}.assoc
          cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_snpeff_loftee/${line}.chr{1..23}_freqUpper0.01.${type}.assoc | grep -v N_INFORMATIVE > temp
          cat header_${type}.txt temp > ${line}.freqUpper0.01.${type}.assoc
          cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_snpeff_loftee/${line}.chr{1..23}_freqUpper0.005.${type}.assoc | grep -v N_INFORMATIVE > temp
          cat header_${type}.txt temp > ${line}.freqUpper0.005.${type}.assoc
          cat /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_snpeff_loftee/${line}.chr{1..23}_freqUpper0.001.${type}.assoc | grep -v N_INFORMATIVE > temp
          cat header_${type}.txt temp > ${line}.freqUpper0.001.${type}.assoc
     done
done
```

### Analyze results
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/meta_prep_snpeff_loftee

ls UKB*SkatO.assoc | wc -l
# 192 (makes sense, 12 variant categories x 4 maf categories x 4 phenotypes)

ls UKB*SkatO.assoc | while read line
do
    sed "s/$/ $line/" $line >> all_SkatO.txt
done

awk '$8 < 5e-8'  all_SkatO.txt | cut -f1 | sort -u
# AUNIP
# DYNC1I1
# GBA
# PSMD11
# SLC26A11

awk '$8 < 1e-6'  all_SkatO.txt | cut -f1 | sort -u
# AUNIP
# DYNC1I1
# GBA
# KIF9
# PDZRN4
# PSMD11
# SEPT5-GP1BB
# SLC26A11
# TUBA1B
```
```
ls UKB*CMC.assoc | wc -l
# 192

ls UKB*CMC.assoc | while read line
do
    sed "s/$/ $line/" $line >> all_CMC.txt
done

awk '$7 < 5e-8' all_CMC.txt | cut -f1 | sort -u
# GBA
# PSMD11
# SLC26A11

awk '$7 < 1e-6' all_CMC.txt | cut -f1 | sort -u
# CLRN1
# DYNC1I1
# GBA
# KIF9
# PSG7
# PSMD11
# SEPT5-GP1BB
# SLC26A11
```

### 6c. Run burden (AAO) [to be completed later]

