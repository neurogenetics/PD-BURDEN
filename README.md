# PD-BURDEN


### UK biobank work

##### Variant selection

```
all_missense.txt => 5, 1, 0.5 and 0.1%
ALL_LOF.txt => (splicing, frameshift, stopgain, stoploss) => 5, 1, 0.5 and 0.1%
ALL_CADD_20.txt => 5, 1, 0.5 and 0.1%
ALL_CADD_10.txt => 5, 1, 0.5 and 0.1%
ALL_MISSENSE_and_LOF.txt => All missense + splicing, frameshift, stopgain, stoploss => 5, 1, 0.5 and 0.1%
```

PDF file for more information
https://www.ukbiobank.ac.uk/wp-content/uploads/2020/11/UK-Biobank-Exome-Release-FAQ_v8-October-2020.pdf




previously done:

- annotated
- filtered for European
- filtered for relatedness


to do now:
- Make final phenotype sheet
all cases
all proxy => no cases

use controls => no PD parent, no AD parent, no PD-ism, no dementia and over >60
```
cd /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/

module load R
R
cov <- read.table("../covariates_phenome_to_use.txt",header=T)
# parkinson
PD <- read.table("PD.txt",header=T)
PD$PD <- 1
PD <- PD[,c(1,5)]
# PD parent
PD_parent <- read.table("PD_parent_no_PD.txt",header=F)
PD_parent$PD_parent <- 1
names(PD_parent)[1] <- "eid"
PD_parent$V2 <- NULL
# AD parent
AD_parent <- read.table("ALL_indi_with_AD_parent.txt",header=F)
AD_parent$AD_parent <- 1
names(AD_parent)[1] <- "eid"
AD_parent$V2 <- NULL
# dementia
Dementia <- read.table("Dementia.txt",header=T)
Dementia$AD <- 1
Dementia <- Dementia[,c(1,5)]
# parkinsonism
PD_ism <- read.table("Parkinsonism.txt",header=T)
PD_ism$PDism <- 1
PD_ism <- PD_ism[,c(1,5)]
# merge all
MM <- merge(cov,PD,by.x="FID",by.y="eid",all=TRUE)
MM1 <- merge(MM,PD_parent,by.x="FID",by.y="eid",all=TRUE)
MM2 <- merge(MM1,AD_parent,by.x="FID",by.y="eid",all=TRUE)
MM3 <- merge(MM2,Dementia,by.x="FID",by.y="eid",all=TRUE)
MM4 <- merge(MM3,PD_ism,by.x="FID",by.y="eid",all=TRUE)
# remove samples no consent...
remove <- read.csv("/data/CARD/UKBIOBANK/SAMPLES_WHO_WITHDRAWAL.csv",header=F)
names(remove)[1] <- "FID"
library(dplyr)
MM5 <- anti_join(MM4, remove, by="FID",)
# replace NA with 0
MM5[is.na(MM5)] <- 0
### start making final lists.... (PRE EXOME FILTER)
CONTROLS <- subset(MM5, EUROPEAN==1 & AGE_OF_RECRUIT >=60 & PD==0 & PD_parent==0 & AD_parent==0 & AD==0 & PDism==0)
# 148948
PD <- subset(MM5, EUROPEAN==1 & PD==1)
# 1735
PDPARENT <- subset(MM5, EUROPEAN==1 & PD_parent==1 & PD==0 )
# 15331
AD <- subset(MM5, EUROPEAN==1 & AD==1)
# 2227
ADPARENT <- subset(MM5, EUROPEAN==1 & AD==0 & AD_parent==1)
# 53460
#### Merge with (filtered) exome data for relatedness...
unrelateds <- read.table("/data/CARD/UKBIOBANK/EXOME_DATA_200K/genotype_data_of_exome_people/genotype_data_of_exome_people_N200469_no_cousins.fam", header=F)
unrelateds <- unrelateds[,c(1,5)]
names(unrelateds)[1] <- "eid"
names(unrelateds)[2] <- "SEX"
MM6 <- merge(MM5,unrelateds,by.x="FID",by.y="eid")
# MM6 => 159253 
# unrelateds => 158317 

### start making final lists.... (POST EXOME FILTER)
CONTROLS <- subset(MM6, EUROPEAN==1 & AGE_OF_RECRUIT >=60 & PD==0 & PD_parent==0 & AD_parent==0 & AD==0 & PDism==0)
# 56311
PD <- subset(MM6, EUROPEAN==1 & PD==1)
# 633
PDPARENT <- subset(MM6, EUROPEAN==1 & PD_parent==1 & PD==0 )
# 6174
AD <- subset(MM6, EUROPEAN==1 & AD==1)
# 794
ADPARENT <- subset(MM6, EUROPEAN==1 & AD==0 & AD_parent==1)
# 22421

#### now randomly assign controls...
# PD first
# PD, PDPARENT, CONTROLS
# split controls 
n <- 56311
ind <- sample(c(TRUE, FALSE), n, replace=TRUE, prob=c(0.1, 0.9))
control_for_PD <- CONTROLS[ind, ]
# 5562
control_for_PDparents <- CONTROLS[!ind, ]
# 50749
PD_CASE_CONTROL <- rbind(PD, control_for_PD)
# 6195 <= 5562 + 633
PD_PARENT_CONTROL <- rbind(PDPARENT, control_for_PDparents)
# 56923 <= 50749 + 6174
write.table(PD_CASE_CONTROL, file="UKB_EXOM_PD_CASE_CONTROL.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PD_PARENT_CONTROL, file="UKB_EXOM_PD_PARENT_CONTROL.txt", quote=FALSE,row.names=F,sep="\t")

# then AD
# AD, ADPARENT, CONTROLS
# split controls 
n <- 56311
ind <- sample(c(TRUE, FALSE), n, replace=TRUE, prob=c(0.1, 0.9))
control_for_AD <- CONTROLS[ind, ]
# 5582
control_for_ADparents <- CONTROLS[!ind, ]
# 50729
AD_CASE_CONTROL <- rbind(AD, control_for_AD)
# 6376 <= 5562 + 794
AD_PARENT_CONTROL <- rbind(ADPARENT, control_for_ADparents)
# 73150 <= 50729 + 22421
write.table(AD_CASE_CONTROL, file="UKB_EXOM_AD_CASE_CONTROL.txt", quote=FALSE,row.names=F,sep="\t")
write.table(AD_PARENT_CONTROL, file="UKB_EXOM_AD_PARENT_CONTROL.txt", quote=FALSE,row.names=F,sep="\t")

# DONE

```

Now create PC's from genotypes....

```
/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups

UKB_EXOM_PD_CASE_CONTROL.txt
UKB_EXOM_AD_PARENT_CONTROL.txt
UKB_EXOM_AD_CASE_CONTROL.txt
UKB_EXOM_PD_PARENT_CONTROL.txt

module load flashpca
module load plink

cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/genotype_data_of_exome_people/

plink --bfile genotype_data_of_exome_people_N200469_no_cousins \
--maf 0.05 --geno 0.01 --hwe 5e-6 --autosome \
--exclude /data/CARD/GENERAL/exclusion_regions_hg19.txt --make-bed --out FILENAME_2 \
--keep /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/UKB_EXOM_PD_CASE_CONTROL.txt
plink --bfile FILENAME_2 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data
plink --bfile FILENAME_2 --extract pruned_data.prune.in --make-bed --out FILENAME_3 
flashpca --bfile FILENAME_3 --suffix _UKB_EXOM_PD_CASE_CONTROL.txt --numthreads 19
# for all 4 subsets...

# merge with phenotype file

cd /data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/

module load R
R
AD_case <- read.table("/data/CARD/UKBIOBANK/EXOME_DATA_200K/genotype_data_of_exome_people/pcs_UKB_EXOM_AD_CASE_CONTROL.txt",header=T)
AD_parent <- read.table("/data/CARD/UKBIOBANK/EXOME_DATA_200K/genotype_data_of_exome_people/pcs_UKB_EXOM_AD_PARENT_CONTROL.txt",header=T)
PD_case <- read.table("/data/CARD/UKBIOBANK/EXOME_DATA_200K/genotype_data_of_exome_people/pcs_UKB_EXOM_PD_CASE_CONTROL.txt",header=T)
PD_parent <- read.table("/data/CARD/UKBIOBANK/EXOME_DATA_200K/genotype_data_of_exome_people/pcs_UKB_EXOM_PD_PARENT_CONTROL.txt",header=T)
AD_case$IID <- NULL
AD_parent$IID <- NULL
PD_case$IID <- NULL
PD_parent$IID <- NULL
covAD <- read.table("UKB_EXOM_AD_CASE_CONTROL.txt",header=T)
covADpat <- read.table("UKB_EXOM_AD_PARENT_CONTROL.txt",header=T)
covPD <- read.table("UKB_EXOM_PD_CASE_CONTROL.txt",header=T)
covPDpat <- read.table("UKB_EXOM_PD_PARENT_CONTROL.txt",header=T)
MM1 <- merge(covAD,AD_case,by="FID")
MM2 <- merge(covADpat,AD_parent,by="FID")
MM3 <- merge(covPD,PD_case,by="FID")
MM4 <- merge(covPDpat,PD_parent,by="FID")
# make sure no duplicates are present...
MM11 <- MM1[!duplicated(MM1), ]
MM22 <- MM2[!duplicated(MM2), ]
MM33 <- MM3[!duplicated(MM3), ]
MM44 <- MM4[!duplicated(MM4), ]
# add in additional pheno column
MM11$PHENO <- MM11$AD+1
MM22$PHENO <- MM22$AD_parent+1
MM33$PHENO <- MM33$PD+1
MM44$PHENO <- MM44$PD_parent+1
# save
write.table(MM11, file="UKB_EXOM_AD_CASE_CONTROL_with_PC.txt", quote=FALSE,row.names=F,sep="\t")
write.table(MM22, file="UKB_EXOM_AD_PARENT_CONTROL_with_PC.txt", quote=FALSE,row.names=F,sep="\t")
write.table(MM33, file="UKB_EXOM_PD_CASE_CONTROL_with_PC.txt", quote=FALSE,row.names=F,sep="\t")
write.table(MM44, file="UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt", quote=FALSE,row.names=F,sep="\t")
q()
n

```

### prep for burden

```
module load plink/2.0-dev-20191128
module load samtools #1.11
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/

# update X to 23...

scp ukb23155_cX_b0_v1.bed ukb23155_c23_b0_v1.bed
scp UKBexomeOQFE_chrX.bim UKBexomeOQFE_chr23.bim

# make vcf files for input for rvtests....
# need to make the following combinations:
all_missense.txt
ALL_LOF.txt
ALL_CADD_20.txt
ALL_CADD_10.txt
ALL_MISSENSE_and_LOF.txt

# make subfolders
mkdir PD_CASE_CONTROL
mkdir PD_PARENT_CONTROL
mkdir AD_CASE_CONTROL
mkdir AD_PARENT_CONTROL

### big loop

module load plink/2.0-dev-20191128
module load samtools #1.11
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/

for CHNUM in {1..23};
do
plink2 --bed ukb23155_c"$CHNUM"_b0_v1.bed --bim UKBexomeOQFE_chr$CHNUM.bim --fam ukb23155_c1_b0_v1_s200632.fam \
--extract ../annotation_of_plink_files/all_missense.txt --keep ../BURDEN/UKB_EXOM_PD_CASE_CONTROL_with_PC.txt \
--export vcf bgz id-paste=iid --out PD_CASE_CONTROL/UKB_EXOM_PD_CASE_CONTROL_chr$CHNUM --mac 1
tabix -p vcf PD_CASE_CONTROL/UKB_EXOM_PD_CASE_CONTROL_chr$CHNUM.vcf.gz
done

and do for all 4 groups:
PD_CASE_CONTROL
PD_PARENT_CONTROL
AD_CASE_CONTROL
AD_PARENT_CONTROL

### examples below for just 1 chromosome....
# PD - control
plink2 --bed ukb23155_c1_b0_v1.bed --bim UKBexomeOQFE_chr1.bim --fam ukb23155_c1_b0_v1_s200632.fam \
--extract ../annotation_of_plink_files/all_missense.txt --keep ../BURDEN/UKB_EXOM_PD_CASE_CONTROL_with_PC.txt \
--export vcf bgz id-paste=iid --out UKB_EXOM_PD_CASE_CONTROL_chr1 --mac 1 
tabix -p vcf UKB_EXOM_PD_CASE_CONTROL_chr1.vcf.gz

# PD parent - control
plink2 --bed ukb23155_c1_b0_v1.bed --bim UKBexomeOQFE_chr1.bim --fam ukb23155_c1_b0_v1_s200632.fam \
--extract ../annotation_of_plink_files/all_missense.txt --keep ../BURDEN/UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt \
--export vcf bgz id-paste=iid --out UKB_EXOM_PD_PARENT_CONTROL_chr1 --mac 1
tabix -p vcf UKB_EXOM_PD_PARENT_CONTROL_chr1.vcf.gz

# AD - control
plink2 --bed ukb23155_c1_b0_v1.bed --bim UKBexomeOQFE_chr1.bim --fam ukb23155_c1_b0_v1_s200632.fam \
--extract ../annotation_of_plink_files/all_missense.txt --keep ../BURDEN/UKB_EXOM_AD_CASE_CONTROL_with_PC.txt \
--export vcf bgz id-paste=iid --out UKB_EXOM_AD_CASE_CONTROL_chr1 --mac 1 
tabix -p vcf UKB_EXOM_AD_CASE_CONTROL_chr1.vcf.gz

# AD parent - control
plink2 --bed ukb23155_c1_b0_v1.bed --bim UKBexomeOQFE_chr1.bim --fam ukb23155_c1_b0_v1_s200632.fam \
--extract ../annotation_of_plink_files/all_missense.txt --keep ../BURDEN/UKB_EXOM_AD_PARENT_CONTROL_with_PC.txt \
--export vcf bgz id-paste=iid --out UKB_EXOM_AD_PARENT_CONTROL_chr1 --mac 1
tabix -p vcf UKB_EXOM_AD_PARENT_CONTROL_chr1.vcf.gz

```

### run burden

Using RVTEST (http://zhanxw.github.io/rvtests/) publication => https://pubmed.ncbi.nlm.nih.gov/27153000/
=> with algorithms:
SKATO, CMC, MB
More info here => http://zhanxw.github.io/rvtests/#burden-tests
SKATO => --kernel skato
CMC => --burden cmc
MB => --burden mb

```
Output folder:
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN

Download rvtests latest version (version: 20190205):
wget https://github.com/zhanxw/rvtests/releases/download/v2.1.0/rvtests_linux64.tar.gz
tar -xf rvtests_linux64.tar.gz

INPUT FORMATS:
UKB_EXOM_AD_CASE_CONTROL_ALL_MISSENSE_and_LOF_chr
UKB_EXOM_AD_CASE_CONTROL_ALL_CADD_20_chr
UKB_EXOM_AD_CASE_CONTROL_ALL_CADD_10_chr
UKB_EXOM_AD_CASE_CONTROL_ALL_LOF_chr
UKB_EXOM_AD_CASE_CONTROL_ALL_MISSENSE_chr

module load rvtests

# sanity checks...
# GBA
rvtest --noweb --hide-covar --out UKB_EXOM_PD_CASE_CONTROL_chr1 --burden cmc  \
--inVcf ../PLINK_files/PD_CASE_CONTROL/UKB_EXOM_PD_CASE_CONTROL_ALL_MISSENSE_and_LOF_chr1.vcf.gz \
--pheno UKB_EXOM_PD_CASE_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_PD_CASE_CONTROL_with_PC.txt --freqUpper 0.05 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --geneFile ../refFlat_HG38.txt --gene GBA
# 1.57857e-05
rvtest --noweb --hide-covar --out GBA --burden cmc  \
--inVcf ../PLINK_files/PD_PARENT_CONTROL/UKB_EXOM_PD_PARENT_CONTROL_ALL_MISSENSE_and_LOF_chr1.vcf.gz \
--pheno UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt --freqUpper 0.05 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --geneFile ../REFFLAT/refFlat_HG38.txt --gene PCSK9
# 1.8764e-08
# GBA all good 
# LRRK2
rvtest --noweb --hide-covar --out UKB_EXOM_PD_CASE_CONTROL_chr12 --burden cmc  \
--inVcf ../PLINK_files/PD_CASE_CONTROL/UKB_EXOM_PD_CASE_CONTROL_ALL_MISSENSE_chr12.vcf.gz \
--pheno UKB_EXOM_PD_CASE_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_PD_CASE_CONTROL_with_PC.txt --freqUpper 0.01 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --geneFile ../refFlat_HG38.txt --gene LRRK2
# 0.0757585
rvtest --noweb --hide-covar --out UKB_EXOM_PD_PARENT_CONTROL_chr12 --burden cmc  \
--inVcf ../PLINK_files/PD_PARENT_CONTROL/UKB_EXOM_PD_PARENT_CONTROL_ALL_MISSENSE_chr12.vcf.gz \
--pheno UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt --freqUpper 0.01 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --geneFile ../refFlat_HG38.txt --gene LRRK2
# 0.586234 => probably too rare...
# TREM2
rvtest --noweb --hide-covar --out UKB_EXOM_AD_CASE_CONTROL_chr6 --burden cmc  \
--inVcf ../PLINK_files/AD_CASE_CONTROL/UKB_EXOM_AD_CASE_CONTROL_ALL_MISSENSE_and_LOF_chr6.vcf.gz \
--pheno UKB_EXOM_AD_CASE_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_AD_CASE_CONTROL_with_PC.txt --freqUpper 0.05 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --geneFile ../refFlat_HG38.txt --gene TREM2
# 0.33655, probably too rare here
rvtest --noweb --hide-covar --out UKB_EXOM_AD_PARENT_CONTROL_chr6 --burden cmc  \
--inVcf ../PLINK_files/AD_PARENT_CONTROL/UKB_EXOM_AD_PARENT_CONTROL_ALL_MISSENSE_and_LOF_chr6.vcf.gz \
--pheno UKB_EXOM_AD_PARENT_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_AD_PARENT_CONTROL_with_PC.txt --freqUpper 0.05 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --geneFile ../refFlat_HG38.txt --gene TREM2
# 7.26046e-06

```

Script format:
```
#!/bin/bash
# sbatch --cpus-per-task=10 --mem=5g --mail-type=END --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020.sh disease variant frequency
# sbatch --cpus-per-task=10 --mem=5g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020.sh PD_CASE_CONTROL ALL_MISSENSE_and_LOF 0.05
DISEASE=$1
VARIANT=$2
FREQUENCY=$3
###
echo "this is"
echo $DISEASE 
echo "whit variant type"
echo $VARIANT
echo "using frequency cut off of"
echo $FREQUENCY
###
# input examples disease level:
# AD_CASE_CONTROL
# AD_PARENT_CONTROL
# PD_CASE_CONTROL
# PD_PARENT_CONTROL
###
# input examples variant level:
# ALL_MISSENSE_and_LOF
# ALL_CADD_20
# ALL_CADD_10
# ALL_LOF
# ALL_MISSENSE
###
# input examples frequency level:
# 0.05
# 0.01
# 0.005
# 0.001
module load rvtests
for CHNUM in {1..22};
do
rvtest --noweb --hide-covar --out RESULTS/$DISEASE/"$VARIANT"_chr$CHNUM --burden cmc --kernel skato \
--inVcf ../PLINK_files/$DISEASE/UKB_EXOM_"$DISEASE"_"$VARIANT"_chr$CHNUM.vcf.gz \
--pheno UKB_EXOM_"$DISEASE"_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_"$DISEASE"_with_PC.txt --freqUpper $FREQUENCY \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --geneFile ../REFFLAT/refFlat_HG38_chr$CHNUM.txt
done

```

```
## checking results files:
# PD_CASE_CONTROL
ls | grep ALL_CADD_10 | wc -l #264
ls | grep ALL_CADD_20 | wc -l #264
ls | grep ALL_MISSENSE_and_LOF | wc -l  #264
ls | grep ALL_LOF | wc -l  #12 
ls | grep ALL_MISSENSE_0 | wc -l #264
---> redo ALL_LOF!!!!
# 264 => 3 files, 4 frequency, 22 chromosomes => 264

# PD_PARENT_CONTROL
ls | grep ALL_CADD_10 | wc -l #264
ls | grep ALL_CADD_20 | wc -l #264
ls | grep ALL_MISSENSE_and_LOF | wc -l #264
ls | grep ALL_LOF | wc -l #264 
ls | grep ALL_MISSENSE_0 | wc -l #264
# 264 => 3 files, 4 frequency, 22 chromosomes => 264




```


