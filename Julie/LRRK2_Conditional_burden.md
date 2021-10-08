# LRRK2 Conditional burdens for PD risk

---

## Structure of README:

### [1. AMP x NIH LRRK2 Conditional Burdens](#1-amp-x-nih-lrrk2-conditional-burdens-1)
### [2. UK Biobank LRRK2 Conditional Burdens](#2-uk-biobank-lrrk2-conditional-burdens-1)

---
### Here are some important files from this analysis:
- Working directory: /data/CARD/PD/AMP_NIH/no_relateds/LRRK2_conditional_burden
- Annovar burden results:
  - AMPxNIH: /data/CARD/PD/AMP_NIH/no_relateds/LRRK2_conditional_burden/burden_annovar
  - UKB: /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/LRRK2_conditional_burden/burden_annovar
- snpEff/loftee burden results:
  - AMPxNIH: /data/CARD/PD/AMP_NIH/no_relateds/LRRK2_conditional_burden/burden_snpeff_loftee
  - UKB: /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/LRRK2_conditional_burden/burden_snpeff_loftee

---

## 1. AMP x NIH LRRK2 Conditional Burdens
```
cd /data/CARD/PD/AMP_NIH/no_relateds
mkdir LRRK2_conditional_burden
cd LRRK2_conditional_burden
```

### Recode p.G2019S
```
module load plink
plink --bfile /data/CARD/PD/AMP_NIH/no_relateds/PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins --recode A --snp chr12:40340400:G:A --out G2019S_status

## Make a new covariate file with p.G2019S status
R 
require(dplyr)
require(data.table)
cov <- fread("/data/CARD/PD/AMP_NIH/no_relateds/COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins.txt",header=T)
G2019S <- fread("G2019S_status.raw",header=T) %>% select(FID, "chr12:40340400:G:A_A")
Mrg <- merge(cov, G2019S, by="FID")
Mrg2 <- Mrg %>% rename(G2019S_status = "chr12:40340400:G:A_A")
Mrg3 <- Mrg2 %>% relocate(PD_PHENO, .after="G2019S_status")
write.table(Mrg3, file="COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins_with_G2019S.txt", quote=FALSE,row.names=F,sep="\t")
q()
n
```

### ANNOVAR
```
mkdir burden_annovar

cat /data/CARD/PD/AMP_NIH/no_relateds/burden_annovar/vcf_files_for_annovar_burden.txt | while read line
do
    sh run_burden_annovar_LRRK2_condi.sh $line 0.05
    sh run_burden_annovar_LRRK2_condi.sh $line 0.01
    sh run_burden_annovar_LRRK2_condi.sh $line 0.005
    sh run_burden_annovar_LRRK2_condi.sh $line 0.001
done

#!/bin/sh
# sh run_burden_annovar_LRRK2_condi.sh AMP_NIH_noRelateds_ALL_LOF.vcf.gz 0.05
module load rvtests
FILENAME=$1
MAF=$2
OUTNAME=${FILENAME/".vcf.gz"/""}
rvtest --noweb --hide-covar --out ./burden_annovar/${OUTNAME}_freqUpper${MAF}_LRRK2 --burden cmc --kernel skato \
--inVcf /data/CARD/PD/AMP_NIH/no_relateds/subset_genetic_data_annovar/${FILENAME} \
--pheno COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins_with_G2019S.txt \
--pheno-name PD_PHENO \
--covar COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins_with_G2019S.txt \
--freqUpper $MAF \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5,G2019S_status \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
```

### SnpEff / LOFTEE
```
mkdir burden_snpeff_loftee

cat /data/CARD/PD/AMP_NIH/no_relateds/burden_snpEff_loftee/vcf_files_for_snpEff_loftee_burden.txt | while read line
do
    sh run_burden_snpeff_loftee_LRRK2_condi.sh $line 0.05
    sh run_burden_snpeff_loftee_LRRK2_condi.sh $line 0.01
    sh run_burden_snpeff_loftee_LRRK2_condi.sh $line 0.005
    sh run_burden_snpeff_loftee_LRRK2_condi.sh $line 0.001
done

#!/bin/sh
# sh run_burden_snpeff_loftee_LRRK2_condi.sh AMP_NIH_noRelateds_ALL_LOF_SNPEFF.vcf.gz 0.05
module load rvtests
FILENAME=$1
MAF=$2
OUTNAME=${FILENAME/".vcf.gz"/""}
rvtest --noweb --hide-covar --out ./burden_snpeff_loftee/${OUTNAME}_freqUpper${MAF}_LRRK2 --burden cmc --kernel skato \
--inVcf /data/CARD/PD/AMP_NIH/no_relateds/subset_genetic_data_snpEff_loftee/${FILENAME} \
--pheno COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins_with_G2019S.txt \
--pheno-name PD_PHENO \
--covar COV_PD_NIH_AMPv2.5_samplestoKeep_EuroOnly_noDups_noNIHDups_wPheno_wSex_no_cousins_with_G2019S.txt \
--freqUpper $MAF \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5,G2019S_status \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
```

## 2. UK Biobank LRRK2 Conditional Burdens
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis
mkdir LRRK2_conditional_burden
cd LRRK2_conditional_burden
```

### Recode p.G2019S
```
## Need to manually select "G" as the reference to recode 
echo chr12:40340400:G:A$'\t'A > snplist.txt

module load plink/2.0-dev-20191128
plink2 --pfile /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/MERGED_UKB_first_pass \
--snp chr12:40340400:G:A --make-bed --out G2019S

module load plink
plink --bfile G2019S --recode A --reference-allele snplist.txt --out G2019S_status

## Make a new covariate file with p.G2019S status
R 
require(dplyr)
require(data.table)
all_pheno <- fread("/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC.txt",header=T)
case_control <- fread("/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC.txt",header=T)
parent_control <- fread("/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC.txt",header=T)
sibling_control <- fread("/data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/reformatted_UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC.txt",header=T)

G2019S <- fread("G2019S_status.raw",header=T) %>% select(IID, "chr12:40340400:G:A_A")

all_pheno_Mrg <- merge(all_pheno, G2019S, by.x="FID", by.y="IID")
all_pheno_Mrg2 <- all_pheno_Mrg %>% rename(G2019S_status = "chr12:40340400:G:A_A")
all_pheno_Mrg3 <- all_pheno_Mrg2 %>% relocate(PHENO, .after="G2019S_status")

case_control_Mrg <- merge(case_control, G2019S, by.x="FID", by.y="IID")
case_control_Mrg2 <- case_control_Mrg %>% rename(G2019S_status = "chr12:40340400:G:A_A")
case_control_Mrg3 <- case_control_Mrg2 %>% relocate(PHENO, .after="G2019S_status")

parent_control_Mrg <- merge(parent_control, G2019S, by.x="FID", by.y="IID")
parent_control_Mrg2 <- parent_control_Mrg %>% rename(G2019S_status = "chr12:40340400:G:A_A")
parent_control_Mrg3 <- parent_control_Mrg2 %>% relocate(PHENO, .after="G2019S_status")

sibling_control_Mrg <- merge(sibling_control, G2019S, by.x="FID", by.y="IID")
sibling_control_Mrg2 <- sibling_control_Mrg %>% rename(G2019S_status = "chr12:40340400:G:A_A")
sibling_control_Mrg3 <- sibling_control_Mrg2 %>% relocate(PHENO, .after="G2019S_status")

write.table(all_pheno_Mrg3, file="reformatted_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_with_PC_with_G2019S.txt", quote=FALSE,row.names=F,sep="\t")
write.table(case_control_Mrg3, file="reformatted_UKB_EXOM_PD_CASE_CONTROL_2021_with_PC_with_G2019S.txt", quote=FALSE,row.names=F,sep="\t")
write.table(parent_control_Mrg3, file="reformatted_UKB_EXOM_PD_PARENT_CONTROL_2021_with_PC_with_G2019S.txt", quote=FALSE,row.names=F,sep="\t")
write.table(sibling_control_Mrg3, file="reformatted_UKB_EXOM_PD_SIBLING_CONTROL_2021_with_PC_with_G2019S.txt", quote=FALSE,row.names=F,sep="\t")
q()
n
```

### ANNOVAR
```
mkdir burden_annovar

grep chr12 /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_annovar/vcf_files_for_annovar_burden.txt | while read line
do
    sh run_burden_annovar_LRRK2_condi.sh $line 0.05
    sh run_burden_annovar_LRRK2_condi.sh $line 0.01
    sh run_burden_annovar_LRRK2_condi.sh $line 0.005
    sh run_burden_annovar_LRRK2_condi.sh $line 0.001
done

#!/bin/sh
# sh run_burden_annovar_LRRK2_condi.sh UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_ALL_CADD_10.chr10.vcf.gz 0.05
module load rvtests
FILENAME=$1
MAF=$2
OUTNAME=${FILENAME/".vcf.gz"/""}
COV_NAME=reformatted_${FILENAME%_ALL*}_with_PC_with_G2019S.txt
rvtest --noweb --hide-covar --out ./burden_annovar/${OUTNAME}_freqUpper${MAF}_LRRK2 --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_annovar/${FILENAME} \
--pheno ${COV_NAME} \
--pheno-name PHENO \
--covar ${COV_NAME} \
--freqUpper $MAF \
--covar-name SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5,G2019S_status \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
```

### SnpEff / LOFTEE
```
mkdir burden_snpeff_loftee

grep chr12 /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/burden_snpeff_loftee/vcf_files_for_snpeff_loftee_burden.txt | while read line
do
    sh run_burden_snpeff_loftee_LRRK2_condi.sh $line 0.05
    sh run_burden_snpeff_loftee_LRRK2_condi.sh $line 0.01
    sh run_burden_snpeff_loftee_LRRK2_condi.sh $line 0.005
    sh run_burden_snpeff_loftee_LRRK2_condi.sh $line 0.001
done

#!/bin/sh
# sh run_burden_snpeff_loftee_LRRK2_condi.sh UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021_ALL_MISSENSE_and_LOF_SNPEFF.chr1.vcf.gz 0.05
module load rvtests
FILENAME=$1
MAF=$2
OUTNAME=${FILENAME/".vcf.gz"/""}
COV_NAME=reformatted_${FILENAME%_2021*}_2021_with_PC_with_G2019S.txt
rvtest --noweb --hide-covar --out ./burden_snpeff_loftee/${OUTNAME}_freqUpper${MAF}_LRRK2 --burden cmc --kernel skato \
--inVcf /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/burden_analysis/subset_genetic_data_snpeff_loftee/${FILENAME} \
--pheno ${COV_NAME} \
--pheno-name PHENO \
--covar ${COV_NAME} \
--freqUpper $MAF \
--covar-name SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5,G2019S_status \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
```

