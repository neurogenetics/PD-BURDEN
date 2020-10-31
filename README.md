# PD-BURDEN


### UK biobank work

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



```







- run burden
SKATO, CMC, MB

