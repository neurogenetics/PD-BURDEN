## PD UK Biobank exome data burden testing

```
November 2020
Cornelis

For more information on UK Biobank exome data please visit: 
https://www.ukbiobank.ac.uk/wp-content/uploads/2020/11/UK-Biobank-Exome-Release-FAQ_v8-October-2020.pdf


Steps:
- sort out phenotype file including QC (relatedness, ancestry, PC's)
- annotation 
- subset data
- burden testing
- post burden testing file prepping for meta-analyses
- create cumulative frequency list for each group

```

### Phenotype file

- filtering for Europeans and relatedness

```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/genotype_data_of_exome_people
module load plink/2.0-dev-20191128

# subset samples
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
 	plink2 --bed /data/CARD/UKBIOBANK/GENOTYPE_DATA/ukb_cal_chr"$chnum"_v2.bed \
	--bim /data/CARD/UKBIOBANK/GENOTYPE_DATA/ukb_snp_chr"$chnum"_v2.bim \
	--fam /data/CARD/UKBIOBANK/GENOTYPE_DATA/ukb33601_cal_chr1_v2_s488363.fam \
	--keep /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/ukb23155_c1_b0_v1_s200632.fam \
	--make-bed --out rawgeno"$chnum"
done

# merge all
ls | grep fam | sed 's/.fam//g' > mergelist.txt
module load plink
plink --merge-list mergelist.txt --make-bed --out genotype_data_of_exome_people_N200469
# Performing single-pass merge (200469 people, 784256 variants).

#!/bin/bash
# sbatch --cpus-per-task=20 --mem=150g --mail-type=END --time=4:00:00 GCTA_genotypes.sh

module load plink/2.0-dev-20191128
module load GCTA

plink --bfile genotype_data_of_exome_people_N200469 --geno 0.01 --maf 0.05 \
--keep /data/CARD/UKBIOBANK/PHENOTYPE_DATA/EUROPEANS.txt --out input_for_prune --make-bed

plink --bfile input_for_prune --indep-pairwise 1000 10 0.02 --out TEMP_pruning
plink --bfile input_for_prune --extract TEMP_pruning.prune.in --make-bed --out TEMP_pruned_data

gcta64 --bfile TEMP_pruned_data --make-grm --out GRM_matrix --autosome --maf 0.05 --threads 20
gcta64 --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0.125 --make-grm --threads 20
plink --bfile genotype_data_of_exome_people_N200469 --keep GRM_matrix_0.125.grm.id --make-bed --out genotype_data_of_exome_people_N200469_no_cousins
# 158,317 people remaining....
```

- Make final phenotype sheet
```
all cases
all proxies => no cases
using controls => no PD parent, no AD parent, no PD-ism, no dementia and over >60

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

- Make PC's from genotypes....

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

### Annotation

```
# create frequency files
module load plink/2.0-dev-20191128
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,};
  do
  	plink --bim UKBexomeOQFE_chr"$chnum".bim --fam ukb23155_c1_b0_v1_s200632.fam \
	--bed ukb23155_c"$chnum"_b0_v1.bed --freq --out FREQchr"$chnum"
done

# selecting one test sample to make annotation faster...
echo 4250267 4250267 > textsample.txt

# subset 1 one sample
module load plink/2.0-dev-20191128
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  	plink --bim UKBexomeOQFE_chr"$chnum".bim --fam ukb23155_c1_b0_v1_s200632.fam \
	--bed ukb23155_c"$chnum"_b0_v1.bed --export vcf id-paste=iid --out anno"$chnum" --keep textsample.txt
done

# annotate all
# note optional.... => -arg '-splicing 15',,,

-------------------
#!/bin/bash
# sbatch --cpus-per-task=10 --mem=200g --time=24:00:00 annotate_UKBv2.sh

cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/

module load annovar
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  	table_annovar.pl anno"$chnum".vcf $ANNOVAR_DATA/hg38 --thread 16 -buildver hg38 \
	-out UKB_exomes_200K_v2_chr"$chnum" -remove -protocol refGene,avsnp150,clinvar_20200316,dbnsfp41a \
	-operation g,f,f,f -nastring . -vcfinput
done

table_annovar.pl anno23.vcf $ANNOVAR_DATA/hg38 --thread 16 -buildver hg38 \
-out UKB_exomes_200K_v2_chr23 -remove -protocol refGene,avsnp150,clinvar_20200316,dbnsfp41a \
-operation g,f,f,f -nastring . -vcfinput

table_annovar.pl anno24.vcf $ANNOVAR_DATA/hg38 --thread 16 -buildver hg38 \
-out UKB_exomes_200K_v2_chr24 -remove -protocol refGene,avsnp150,clinvar_20200316,dbnsfp41a \
-operation g,f,f,f -nastring . -vcfinput
-------------------


# move data to annotation folder

mv UKB_exomes_200K_v2_chr* ../annotation_of_plink_files/
mv FREQchr*afreq ../annotation_of_plink_files/

# clean up folder
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/annotation_of_plink_files/
mkdir inputfiles
mv *.avinput inputfiles/
rm *.vcf 

# merge data with frequency file...
mkdir freq_files
mv *.afreq freq_files/

for chnum in {1..24};
  do
	paste UKB_exomes_200K_v2_chr"$chnum".hg38_multianno.txt freq_files/FREQchr"$chnum".afreq > UKB_exomes_200K_v2_chr"$chnum".hg38_multianno.withafreq.txt
done

# clean up folder

mkdir anno_no_freq
mv *.hg38_multianno.txt anno_no_freq/

### subset "groups" of variants...

# missense
grep exonic UKB_exomes_200K_v2_chr*.hg38_multianno.withafreq.txt | grep nonsynonymous | cut -f 92 > ALL_MISSENSE.txt
# n=4672717
grep exonic UKB_exomes_200K_v2_chr*.hg38_multianno.withafreq.txt | grep nonsynonymous | cut -f 73 > ALL_MISSENSE.txt

# LOF (stop, frame)
grep stopgain UKB_exomes_200K_v2_chr*.hg38_multianno.withafreq.txt | cut -f 73 > all_stopgain.txt
# n=159722
grep stoploss UKB_exomes_200K_v2_chr*.hg38_multianno.withafreq.txt | cut -f 73 > all_stoploss.txt
# n=7051
grep nonframeshift UKB_exomes_200K_v2_chr*.hg38_multianno.withafreq.txt | cut -f 73 > all_nonframeshift.txt
# n=107177
grep frame UKB_exomes_200K_v2_chr*.hg38_multianno.withafreq.txt | grep -v nonframeshift | grep -v nonsynonymous | cut -f 73 > all_frameshift.txt
# n=203946

# splicing 
grep splicing UKB_exomes_200K_v2_chr*.hg38_multianno.withafreq.txt | \
grep -v ncRNA | cut -f 6,73 | grep splicing | cut -f 2 > all_splice.txt
# n=93257

# CADD <10
awk '{ if($42 >= 10) { print }}' UKB_exomes_200K_v2_chr*.hg38_multianno.withafreq.txt | cut -f 73 > ALL_CADD_10.txt
# n=276502

# CADD <20
awk '{ if($42 >= 20) { print }}' UKB_exomes_200K_v2_chr*.hg38_multianno.withafreq.txt | cut -f 73 > ALL_CADD_20.txt
# n=244442


#### Prepping final files:

cat all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_LOF.txt

cat all_missense.txt all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_MISSENSE_and_LOF.txt

cat ALL_CADD_20.txt all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice.txt > ALL_CADD_20_and_LOF.txt

wc -l ALL_CADD_10.txt # 276502
wc -l ALL_CADD_20.txt # 244442
wc -l ALL_MISSENSE_and_LOF.txt # 5136354
wc -l ALL_LOF.txt # 463976
wc -l ALL_MISSENSE.txt # 4672717
wc -l ALL_CADD_20_and_LOF.txt # 708418

```


### Subset data

```
module load plink/2.0-dev-20191128
module load samtools #1.11
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/

# update X to 23...

scp ukb23155_cX_b0_v1.bed ukb23155_c23_b0_v1.bed
scp UKBexomeOQFE_chrX.bim UKBexomeOQFE_chr23.bim

# make vcf files for input for rvtests....
# need to make the following combinations:
ALL_MISSENSE.txt
ALL_LOF.txt
ALL_CADD_20.txt
ALL_CADD_10.txt
ALL_MISSENSE_and_LOF.txt
ALL_CADD_20_and_LOF.txt

# make subfolders
mkdir PD_CASE_CONTROL
mkdir PD_PARENT_CONTROL
mkdir AD_CASE_CONTROL
mkdir AD_PARENT_CONTROL

### big loop # check => subset_data_prior_burden.sh

module load plink/2.0-dev-20191128
module load samtools #1.11
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/

for CHNUM in {1..23};
do
plink2 --bed ukb23155_c"$CHNUM"_b0_v1.bed --bim UKBexomeOQFE_chr$CHNUM.bim --fam ukb23155_c1_b0_v1_s200632.fam \
--extract ../annotation_of_plink_files/ALL_MISSENSE.txt --keep ../BURDEN/UKB_EXOM_PD_CASE_CONTROL_with_PC.txt \
--export vcf bgz id-paste=iid --out PD_CASE_CONTROL/UKB_EXOM_PD_CASE_CONTROL_ALL_MISSENSE_chr$CHNUM --mac 1
tabix -p vcf PD_CASE_CONTROL/UKB_EXOM_PD_CASE_CONTROL_ALL_MISSENSE_chr$CHNUM.vcf.gz
done

for CHNUM in {1..23};
do
plink2 --bed ukb23155_c"$CHNUM"_b0_v1.bed --bim UKBexomeOQFE_chr$CHNUM.bim --fam ukb23155_c1_b0_v1_s200632.fam \
--extract ../annotation_of_plink_files/ALL_MISSENSE.txt --keep ../BURDEN/UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt \
--export vcf bgz id-paste=iid --out PD_PARENT_CONTROL/UKB_EXOM_PD_PARENT_CONTROL_ALL_MISSENSE_chr$CHNUM --mac 1
tabix -p vcf PD_PARENT_CONTROL/UKB_EXOM_PD_PARENT_CONTROL_ALL_MISSENSE_chr$CHNUM.vcf.gz
done

for CHNUM in {1..23};
do
plink2 --bed ukb23155_c"$CHNUM"_b0_v1.bed --bim UKBexomeOQFE_chr$CHNUM.bim --fam ukb23155_c1_b0_v1_s200632.fam \
--extract ../annotation_of_plink_files/ALL_MISSENSE.txt --keep ../BURDEN/UKB_EXOM_AD_CASE_CONTROL_with_PC.txt \
--export vcf bgz id-paste=iid --out AD_CASE_CONTROL/UKB_EXOM_AD_CASE_CONTROL_ALL_MISSENSE_chr$CHNUM --mac 1
tabix -p vcf AD_CASE_CONTROL/UKB_EXOM_AD_CASE_CONTROL_ALL_MISSENSE_chr$CHNUM.vcf.gz
done

for CHNUM in {1..23};
do
plink2 --bed ukb23155_c"$CHNUM"_b0_v1.bed --bim UKBexomeOQFE_chr$CHNUM.bim --fam ukb23155_c1_b0_v1_s200632.fam \
--extract ../annotation_of_plink_files/ALL_MISSENSE.txt --keep ../BURDEN/UKB_EXOM_AD_PARENT_CONTROL_with_PC.txt \
--export vcf bgz id-paste=iid --out AD_PARENT_CONTROL/UKB_EXOM_AD_PARENT_CONTROL_ALL_MISSENSE_chr$CHNUM --mac 1
tabix -p vcf AD_PARENT_CONTROL/UKB_EXOM_AD_PARENT_CONTROL_ALL_MISSENSE_chr$CHNUM.vcf.gz
done

and do for all 6 groups:
ALL_MISSENSE.txt
ALL_LOF.txt
ALL_CADD_20.txt
ALL_CADD_10.txt
ALL_MISSENSE_and_LOF.txt
ALL_CADD_20_and_LOF.txt


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

### Run burden

Using RVTEST (http://zhanxw.github.io/rvtests/) publication => https://pubmed.ncbi.nlm.nih.gov/27153000/
=> with algorithms:
SKATO and CMC
More info here => http://zhanxw.github.io/rvtests/#burden-tests
SKATO => --kernel skato
CMC => --burden cmc

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
UKB_EXOM_AD_CASE_CONTROL_ALL_MISSENSE_chr
UKB_EXOM_AD_CASE_CONTROL_ALL_CADD_20_and_LOF_chr

module load rvtests # 2.1.0

# sanity checks...
# GBA
rvtest --noweb --hide-covar --out GBA --burden cmc  \
--inVcf ../PLINK_files/PD_CASE_CONTROL/UKB_EXOM_PD_CASE_CONTROL_ALL_MISSENSE_and_LOF_chr1.vcf.gz \
--pheno UKB_EXOM_PD_CASE_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_PD_CASE_CONTROL_with_PC.txt --freqUpper 0.05 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# 1.57857e-05
rvtest --noweb --hide-covar --out GBA_parent --burden cmc  \
--inVcf ../PLINK_files/PD_PARENT_CONTROL/UKB_EXOM_PD_PARENT_CONTROL_ALL_MISSENSE_and_LOF_chr1.vcf.gz \
--pheno UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt --freqUpper 0.05 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# 1.8764e-08
rvtest --noweb --hide-covar --out GBA_parent2 --burden cmc  \
--inVcf ../PLINK_files/PD_PARENT_CONTROL/UKB_EXOM_PD_PARENT_CONTROL_ALL_CADD_20_and_LOF_chr1.vcf.gz \
--pheno UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt --freqUpper 0.05 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# 0.0606762
# GBA all good 
# LRRK2
rvtest --noweb --hide-covar --out LRRK2 --burden cmc  \
--inVcf ../PLINK_files/PD_CASE_CONTROL/UKB_EXOM_PD_CASE_CONTROL_ALL_MISSENSE_chr12.vcf.gz \
--pheno UKB_EXOM_PD_CASE_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_PD_CASE_CONTROL_with_PC.txt --freqUpper 0.01 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
# 0.0757585
rvtest --noweb --hide-covar --out LRRK2_parent --burden cmc  \
--inVcf ../PLINK_files/PD_PARENT_CONTROL/UKB_EXOM_PD_PARENT_CONTROL_ALL_MISSENSE_chr12.vcf.gz \
--pheno UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_PD_PARENT_CONTROL_with_PC.txt --freqUpper 0.01 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
# 0.586234 => probably too rare...
# TREM2
rvtest --noweb --hide-covar --out TREM2 --burden cmc  \
--inVcf ../PLINK_files/AD_CASE_CONTROL/UKB_EXOM_AD_CASE_CONTROL_ALL_MISSENSE_and_LOF_chr6.vcf.gz \
--pheno UKB_EXOM_AD_CASE_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_AD_CASE_CONTROL_with_PC.txt --freqUpper 0.05 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr6.txt --gene TREM2
# 0.33655, probably too rare here
rvtest --noweb --hide-covar --out TREM2_parent --burden cmc  \
--inVcf ../PLINK_files/AD_PARENT_CONTROL/UKB_EXOM_AD_PARENT_CONTROL_ALL_MISSENSE_and_LOF_chr6.vcf.gz \
--pheno UKB_EXOM_AD_PARENT_CONTROL_with_PC.txt --pheno-name PHENO \
--covar UKB_EXOM_AD_PARENT_CONTROL_with_PC.txt --freqUpper 0.05 \
--covar-name GENETIC_SEX,AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 \
--geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr6.txt --gene TREM2
# 7.26046e-06

```

Script format:
```
This is mainly for small data => AD case-control and PD case-control
BURDEN_TESTING_CMC_SKAT_2020.sh # general script to run per group etc
LAUNCH_BURDEN_TESTING_CMC_SKAT_2020.sh # launches above
sh LAUNCH_BURDEN_TESTING_CMC_SKAT_2020.sh

This is mainly for big data => AD proxy-control and PD proxy-control
BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh # general script to run per chromosome etc
LAUNCH_BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh # launches above
sh LAUNCH_BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh
```

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
# ALL_CADD_20_and_LOF
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

# DOUBLE CHECK HERE
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

# confirm filesizes....
# PD_CASE_CONTROL
wc -l *.assoc | sort -nk1 | head
-> all look good except for chr1 files of LOF files ??? times 4 (frequency levels)
# PD_PARENT_CONTROL
wc -l *.assoc | sort -nk1 | head
-> all look good except for chr2 files of ALL_MISSENSE and ALL_MISSENSE_and_LOF ????? both times 4 (frequency levels)
wc -l *.assoc | sort -nk1 | head

# trying to understand what goes wrong....
module load plink/2.0-dev-20191128

PD_CASE_CONTROL chr1 files of LOF files

UKB_EXOM_PD_CASE_CONTROL_ALL_LOF_chr1.log
UKB_EXOM_PD_CASE_CONTROL_ALL_LOF_chr1.vcf.gz
UKB_EXOM_PD_CASE_CONTROL_ALL_LOF_chr1.vcf.gz.tbi
"4738 variants remaining after main filters."
"6193 samples (3182 females, 3011 males; 6193 founders)"
according to plink:
6193 samples and 4738 variants
# rerun.....
sbatch --cpus-per-task=10 --mem=25g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_CASE_CONTROL ALL_LOF 0.05 1
sbatch --cpus-per-task=10 --mem=25g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_CASE_CONTROL ALL_LOF 0.01 1
sbatch --cpus-per-task=10 --mem=25g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_CASE_CONTROL ALL_LOF 0.005 1
sbatch --cpus-per-task=10 --mem=25g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_CASE_CONTROL ALL_LOF 0.001 1

# if still doesnt work?
Recreate vcf file?


PD_PARENT_CONTROL chr2 files of ALL_MISSENSE and ALL_MISSENSE_and_LOF
# rerun.....

sbatch --cpus-per-task=15 --mem=50g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_PARENT_CONTROL ALL_MISSENSE 0.05 2
sbatch --cpus-per-task=15 --mem=50g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_PARENT_CONTROL ALL_MISSENSE 0.01 2
sbatch --cpus-per-task=15 --mem=50g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_PARENT_CONTROL ALL_MISSENSE 0.005 2
sbatch --cpus-per-task=15 --mem=50g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_PARENT_CONTROL ALL_MISSENSE 0.001 2

sbatch --cpus-per-task=15 --mem=50g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_PARENT_CONTROL ALL_MISSENSE_and_LOF 0.05 2
sbatch --cpus-per-task=15 --mem=50g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_PARENT_CONTROL ALL_MISSENSE_and_LOF 0.01 2
sbatch --cpus-per-task=15 --mem=50g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_PARENT_CONTROL ALL_MISSENSE_and_LOF 0.005 2
sbatch --cpus-per-task=15 --mem=50g --time=24:00:00 BURDEN_TESTING_CMC_SKAT_2020_per_chr.sh PD_PARENT_CONTROL ALL_MISSENSE_and_LOF 0.001 2


```

###  Post burden testing file prepping for meta-analyses...

```
Steps:
- merge per chromosome
- sort based on P-value

## PD_CASE_CONTROL
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN/RESULTS/PD_CASE_CONTROL

cat ../variant_types.txt | while read line
do 
	head -1 ALL_MISSENSE_0.01_chr20.CMC.assoc > CMC_header.txt
	head -1 ALL_MISSENSE_0.01_chr20.SkatO.assoc > SkatO_header.txt
	# SKATO
	cat "$line"_0.05_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../PD_CASE_CONTROL_"$line"_0.05.SkatO.assoc
	cat "$line"_0.01_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../PD_CASE_CONTROL_"$line"_0.01.SkatO.assoc
	cat "$line"_0.005_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../PD_CASE_CONTROL_"$line"_0.005.SkatO.assoc
	cat "$line"_0.001_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../PD_CASE_CONTROL_"$line"_0.001.SkatO.assoc
	# CMC
	cat "$line"_0.05_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../PD_CASE_CONTROL_"$line"_0.05.CMC.assoc
	cat "$line"_0.01_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../PD_CASE_CONTROL_"$line"_0.01.CMC.assoc
	cat "$line"_0.005_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../PD_CASE_CONTROL_"$line"_0.005.CMC.assoc
	cat "$line"_0.001_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../PD_CASE_CONTROL_"$line"_0.001.CMC.assoc
done

## PD_PARENT_CONTROL

cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN/RESULTS/PD_PARENT_CONTROL

cat ../variant_types.txt | while read line
do 
	head -1 ALL_MISSENSE_0.01_chr20.CMC.assoc > CMC_header.txt
	head -1 ALL_MISSENSE_0.01_chr20.SkatO.assoc > SkatO_header.txt
	# SKATO
	cat "$line"_0.05_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../PD_PARENT_CONTROL_"$line"_0.05.SkatO.assoc
	cat "$line"_0.01_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../PD_PARENT_CONTROL_"$line"_0.01.SkatO.assoc
	cat "$line"_0.005_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../PD_PARENT_CONTROL_"$line"_0.005.SkatO.assoc
	cat "$line"_0.001_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../PD_PARENT_CONTROL_"$line"_0.001.SkatO.assoc
	# CMC
	cat "$line"_0.05_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../PD_PARENT_CONTROL_"$line"_0.05.CMC.assoc
	cat "$line"_0.01_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../PD_PARENT_CONTROL_"$line"_0.01.CMC.assoc
	cat "$line"_0.005_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../PD_PARENT_CONTROL_"$line"_0.005.CMC.assoc
	cat "$line"_0.001_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../PD_PARENT_CONTROL_"$line"_0.001.CMC.assoc
done

```
# DOUBLE CHECK HERE

```
# move files to results file
mkdir RESULTS
mv *.assoc RESULTS/

40 files => 2 tests, 4 variant thresholds, 5 variant groups => 40
```


### Create cumulative frequency list for each group...

```
----------------------------------
#!/bin/bash
# sbatch --cpus-per-task=10 --mem=20g --time=24:00:00 Cumulative_MAF_prep.sh PD_CASE_CONTROL ALL_LOF 0.05
module load plink
# mkdir freq_case_control_per_variant_class
# run quick assoc to create F_A and F_U formats.
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
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/
# lets start
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
plink --bed ukb23155_c"$chnum"_b0_v1.bed \
--fam ukb23155_c1_b0_v1_s200632.fam \
--bim UKBexomeOQFE_chr"$chnum".bim \
--assoc --out freq_case_control_per_variant_class/UKB_EXOM_"$DISEASE"_"$VARIANT"_"$FREQUENCY"_"$chnum" \
--pheno ../BURDEN/UKB_EXOM_"$DISEASE"_with_PC.txt \
--keep ../BURDEN/UKB_EXOM_"$DISEASE"_with_PC.txt \
--extract ../annotation_of_plink_files/"$VARIANT".txt \
--pheno-name PHENO --mac 1 --max-maf "$FREQUENCY"
done

# prep .assoc files
cd freq_case_control_per_variant_class
cat UKB_EXOM_"$DISEASE"_"$VARIANT"_"$FREQUENCY"_*.assoc > UKB_EXOM_"$DISEASE"_"$VARIANT"_"$FREQUENCY"_ALL.assoc

# R session...
module load R
Rscript --vanilla Cumulative_MAF_step2.R "$DISEASE" "$VARIANT" "$FREQUENCY"
# Rscript --vanilla Cumulative_MAF_step2.R PD_CASE_CONTROL ALL_LOF 0.05

# done
echo "bye bye"
----------------------------------
```

```
# prep annotation...
# cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/annotation_of_plink_files/
# cut -f 1-11,92 UKB_exomes_200K_chr*.hg38_multianno.withafreq.txt > SHORT_annotation.txt

# then merge in R and filter....
```

```
----------------------------------
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# start like this
# Rscript --vanilla Cumulative_MAF_step2.R DISEASE VARIANT FREQUENCY
# Rscript --vanilla Cumulative_MAF_step2.R PD_CASE_CONTROL ALL_LOF 0.05
DISEASE = args[1]
print(args[1])
print(DISEASE)
VARIANT = args[2]
print(args[2])
print(VARIANT)
FREQUENCY = args[3]
print(args[3])
print(FREQUENCY)
# R session...
# module load R
# Rscript --vanilla Cumulative_MAF_step2.R DISEASE VARIANT FREQUENCY

# module load R
# R
require("data.table")
assoc <- fread(paste("UKB_EXOM_",DISEASE,"_",VARIANT,"_",FREQUENCY,"_ALL.assoc",sep=""),header=T)
# assoc <- fread("UKB_EXOM_PD_CASE_CONTROL_ALL_LOF_0.05_ALL.assoc",header=T)
# assoc <- assoc[c("CHR","SNP","F_A","F_U")]
annotation <- fread("../../annotation_of_plink_files/SHORT_annotation.txt",header=T)
dim(assoc)
dim(annotation)
MM <- merge(assoc,annotation,by.x="SNP",by.y="Otherinfo6")
dim(MM)
# storing frequencies as numbers? kinda nonsense but lets continue
MM$F_A <- as.numeric(MM$F_A)
MM$F_U <- as.numeric(MM$F_U)
# now calculate cumulative case-control frequency
case_feq <- aggregate(x = MM$F_A, by = list(MM$Gene.refGene),FUN = sum)   
control_feq <- aggregate(x = MM$F_U, by = list(MM$Gene.refGene),FUN = sum)
# update header
names(case_feq)[1] <- "GENE"
names(case_feq)[2] <- "case_freq"
names(control_feq)[1] <- "GENE"
names(control_feq)[2] <- "control_freq"
MM2 <- merge(case_feq,control_feq,by="GENE")
write.table(MM2, file=paste("UKB_EXOM_",DISEASE,"_",VARIANT,"_",FREQUENCY,"_cumulative_MAF.txt",sep=""), quote=FALSE,row.names=F,sep="\t")## output format:
#    GENE case_freq control_freq
#   ABCA1 0.0005379    0.0006090
#   ABCA2 0.0005379    0.0001218
q()
n
----------------------------------

## note major downside here is that if annovar annotates the variants with this format: AMY1A;AMY1B;AMY1C then these variants do net get counted towards any of these three genes... Probably genome-wide not a real problem, but something worth double checking for genes of interest.
```

```
## merge files with result files...

work in progress....


```





