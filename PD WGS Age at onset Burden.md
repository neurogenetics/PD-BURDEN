### PD WGS Age at onset Burden

```
November 2020
Cornelis

Steps:
(Done at PD WGS Burden.md)
- sort out phenotype file
- rename variants
- annotation
New steps:
- QC (relatedness, ancestry, PC's)
- subset data
- burden testing
- post burden testing file prepping for meta-analyses

```

### QC (relatedness, ancestry)

```
cd /data/CARD/PD/WGS/june2019
module load plink/2.0-dev-20191128

# subset samples
Using cases only => PD_cases_only_with_AGE.txt
subset from PHENO_FOR_GWAS_v1_november11_with_PC.txt using cases only and with AGE available...

mkdir PCA_AAO
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
 	plink2 --pgen pd.june2019.chr$chnum.freeze9.sqc.pgen --psam pd.june2019.chr$chnum.freeze9.sqc.pgen \
  --pvar PVAR_files/NEW$chnum.pvar \
	--keep PD_cases_only_with_AGE.txt --maf 0.01 \
	--make-bed --out PCA_AAO/geno"$chnum"
done

```
# Continue here
```

# failes due to multi-allelic's
cd PCA
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
 	plink2 --bfile geno"$chnum" \
	--exclude PD_WGS-merge.missnp --maf 0.01 \
	--make-bed --out merge"$chnum"
done

# now merge....

ls | grep fam | grep "merge" | grep -v PD_WGS | sed 's/.fam//g' > mergelist.txt
module load plink
plink --merge-list mergelist.txt --make-bed --out PD_WGS_PCA_input
# Performing single-pass merge (8931 people, 8352699 variants).

# now pruning and relatedness...

#!/bin/bash
# sbatch --cpus-per-task=20 --mem=150g --mail-type=END --time=4:00:00 GCTA_genotypes.sh

module load plink/2.0-dev-20191128
module load GCTA
module load flashpca

plink --bfile PD_WGS_PCA_input --geno 0.01 --maf 0.05 --out input_for_prune --make-bed
plink --bfile input_for_prune --indep-pairwise 500 10 0.02 --out TEMP_pruning
plink --bfile input_for_prune --extract TEMP_pruning.prune.in --make-bed --out TEMP_pruned_data
gcta64 --bfile TEMP_pruned_data --make-grm --out GRM_matrix --autosome --maf 0.05 --threads 20
gcta64 --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0.125 --make-grm --threads 20
plink --bfile PD_WGS_PCA_input --keep GRM_matrix_0.125.grm.id --make-bed --out PD_WGS_no_relatedness
# 8441 samples (3959 females, 4482 males; 8441 founders) remaining after main filters.

# make PC's of no relatedness data
plink --bfile PD_WGS_no_relatedness --geno 0.01 --maf 0.05 --out input_for_prune --make-bed
plink --bfile input_for_prune --indep-pairwise 500 10 0.02 --out TEMP_pruning
plink --bfile input_for_prune --extract TEMP_pruning.prune.in --make-bed --out TEMP_pruned_data
flashpca --bfile TEMP_pruned_data --suffix _PD_WGS_PCs_no_relatedness.txt --numthreads 19

### plotting looks good no visible outliers....

# make PC's and remove relatedness samples from only keeping GWAS viable samples...

plink --bfile PD_WGS_PCA_input --geno 0.01 --maf 0.05 --out input_for_prune --make-bed \
--keep ../PHENO_FOR_GWAS_v1_november11.txt
plink --bfile input_for_prune --indep-pairwise 500 10 0.02 --out TEMP_pruning
plink --bfile input_for_prune --extract TEMP_pruning.prune.in --make-bed --out TEMP_pruned_data
gcta64 --bfile TEMP_pruned_data --make-grm --out GRM_matrix --autosome --maf 0.05 --threads 20
gcta64 --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0.125 --make-grm --threads 20
plink --bfile PD_WGS_PCA_input --keep GRM_matrix_0.125.grm.id --make-bed --out PD_WGS_no_relatedness_GWAS_viable
# 6893 samples (3095 females, 3798 males; 6893 founders) remaining after main filters.

# make PC's of no relatedness data + only to keep for GWAS....
plink --bfile PD_WGS_no_relatedness_GWAS_viable --geno 0.01 --maf 0.05 --out input_for_prune --make-bed \
--keep ../PHENO_FOR_GWAS_v1_november11.txt
plink --bfile input_for_prune --indep-pairwise 500 10 0.02 --out TEMP_pruning
plink --bfile input_for_prune --extract TEMP_pruning.prune.in --make-bed --out TEMP_pruned_data
flashpca --bfile TEMP_pruned_data --suffix _PD_WGS_PCs_no_relatedness_GWAS_viable.txt --numthreads 19

###
File to move forward with:
pcs_PD_WGS_PCs_no_relatedness_GWAS_viable.txt

# merge with phenotype file...
PHENO_FOR_GWAS_v1_november11.txt

module load R
R
pheno <- read.table("../PHENO_FOR_GWAS_v1_november11.txt",header=T)
pc <- read.table("pcs_PD_WGS_PCs_no_relatedness_GWAS_viable.txt",header=T)
pc$IID <- NULL
MM <- merge(pheno,pc,by="FID")
# add in extra PHENO because PHENO to close to beginning and rvtest cant handle that
MM$PHENO_RV <- MM$PHENO
write.table(MM, file="../PHENO_FOR_GWAS_v1_november11_with_PC.txt", quote=FALSE,row.names=F,sep="\t")

###
File to move forward with:
PHENO_FOR_GWAS_v1_november11_with_PC.txt

6893 samples (3095 females, 3798 males; 6893 founders) remaining after main filters.
```

###  subset data
Subset WGS data using only samples of interest from:
PHENO_FOR_GWAS_v1_november11_with_PC.txt
and including only variants of interest from:
ALL_MISSENSE.txt
ALL_LOF.txt
ALL_CADD_20.txt
ALL_CADD_10.txt
ALL_MISSENSE_and_LOF.txt

```
# make folders for new vcf files

mkdir ALL_MISSENSE
mkdir ALL_LOF
mkdir ALL_CADD_20
mkdir ALL_CADD_10
mkdir ALL_MISSENSE_and_LOF

module load plink/2.0-dev-20191128
module load samtools

# ALL_LOF
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
 	plink2 --pfile pd.june2019.chr"$chnum".freeze9.sqc \
	--extract annotation/ALL_LOF.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--out ALL_LOF/PD_WGS_ALL_LOF_"$chnum" \
	--mac 1 --export vcf bgz id-paste=iid
	tabix -p vcf ALL_LOF/PD_WGS_ALL_LOF_"$chnum".vcf.gz
done

# ALL_CADD_20
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
 	plink2 --pfile pd.june2019.chr"$chnum".freeze9.sqc \
	--extract annotation/ALL_CADD_20.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--out ALL_CADD_20/PD_WGS_ALL_CADD_20_"$chnum" \
	--mac 1 --export vcf bgz id-paste=iid
	tabix -p vcf ALL_CADD_20/PD_WGS_ALL_CADD_20_"$chnum".vcf.gz
done

# ALL_CADD_10
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
 	plink2 --pfile pd.june2019.chr"$chnum".freeze9.sqc \
	--extract annotation/ALL_CADD_10.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--out ALL_CADD_10/PD_WGS_ALL_CADD_10_"$chnum" \
	--mac 1 --export vcf bgz id-paste=iid
	tabix -p vcf ALL_CADD_10/PD_WGS_ALL_CADD_10_"$chnum".vcf.gz
done

# ALL_MISSENSE
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
 	plink2 --pfile pd.june2019.chr"$chnum".freeze9.sqc \
	--extract annotation/all_missense.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--out ALL_MISSENSE/PD_WGS_ALL_MISSENSE_"$chnum" \
	--mac 1 --export vcf bgz id-paste=iid
	tabix -p vcf ALL_MISSENSE/PD_WGS_ALL_MISSENSE_"$chnum".vcf.gz
done

# ALL_MISSENSE_and_LOF
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
 	plink2 --pfile pd.june2019.chr"$chnum".freeze9.sqc \
	--extract annotation/ALL_MISSENSE_and_LOF.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--out ALL_MISSENSE_and_LOF/PD_WGS_ALL_MISSENSE_and_LOF_"$chnum" \
	--mac 1 --export vcf bgz id-paste=iid
	tabix -p vcf ALL_MISSENSE_and_LOF/PD_WGS_ALL_MISSENSE_and_LOF_"$chnum".vcf.gz
done


```


###  burden testing

```
working dir:
cd /data/CARD/PD/WGS/june2019

# output dir:
mkdir BURDEN
cd BURDEN
mkdir ALL_MISSENSE
mkdir ALL_LOF
mkdir ALL_CADD_20
mkdir ALL_CADD_10
mkdir ALL_MISSENSE_and_LOF

# input files:
ALL_MISSENSE_and_LOF/PD_WGS_ALL_MISSENSE_and_LOF_"$chnum".vcf.gz
ALL_MISSENSE/PD_WGS_ALL_MISSENSE_"$chnum".vcf.gz
ALL_CADD_10/PD_WGS_ALL_CADD_10_"$chnum".vcf.gz
ALL_LOF/PD_WGS_ALL_LOF_"$chnum".vcf.gz
ALL_CADD_20/PD_WGS_ALL_CADD_20_"$chnum".vcf.gz

# pheno/covariate:
PHENO_FOR_GWAS_v1_november11_with_PC.txt
names are:
SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 
pheno is:
PHENO

# start burdening...

module load rvtests

# sanity checks...
# GBA
rvtest --noweb --hide-covar --out BURDEN/ALL_MISSENSE/PD_WGS_burden_chr1 --burden cmc  \
--inVcf ALL_MISSENSE_and_LOF/PD_WGS_ALL_MISSENSE_and_LOF_1.vcf.gz \
--pheno PHENO_FOR_GWAS_v1_november11_with_PC.txt --pheno-name PHENO_RV \
--covar PHENO_FOR_GWAS_v1_november11_with_PC.txt --freqUpper 0.05 --imputeCov \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# Loaded 2788 cases, 4105 controls, and 0 missing phenotypes => 6893 samples
P = 4.00621e-13
# LRRK2
rvtest --noweb --hide-covar --out BURDEN/ALL_MISSENSE/PD_WGS_burden_chr12 --burden cmc  \
--inVcf ALL_MISSENSE/PD_WGS_ALL_MISSENSE_12.vcf.gz \
--pheno PHENO_FOR_GWAS_v1_november11_with_PC.txt --pheno-name PHENO_RV \
--covar PHENO_FOR_GWAS_v1_november11_with_PC.txt --freqUpper 0.05 --imputeCov \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr12.txt --gene LRRK2
# Loaded 2788 cases, 4105 controls, and 0 missing phenotypes => 6893 samples
P = 0.101713 ALL_MISSENSE

```

```
# now running all
#!/bin/bash
# sbatch --cpus-per-task=10 --mem=5g --mail-type=END --time=24:00:00 WGS_BURDEN_TESTING_CMC_SKAT_2020.sh variant frequency
# sbatch --cpus-per-task=10 --mem=5g --time=24:00:00 WGS_BURDEN_TESTING_CMC_SKAT_2020.sh ALL_MISSENSE_and_LOF 0.05
# sbatch --cpus-per-task=10 --mem=5g --time=24:00:00 WGS_BURDEN_TESTING_CMC_SKAT_2020.sh ALL_LOF 0.05
VARIANT=$1
FREQUENCY=$2
###
echo "this is"
echo "PD case control" 
echo "whit variant type"
echo $VARIANT
echo "using frequency cut off of"
echo $FREQUENCY
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
rvtest --noweb --hide-covar --out BURDEN/"$VARIANT"/PD_WGS_"$FREQUENCY"_chr$CHNUM --burden cmc --kernel skato \
--inVcf "$VARIANT"/PD_WGS_"$VARIANT"_"$CHNUM".vcf.gz \
--pheno PHENO_FOR_GWAS_v1_november11_with_PC.txt --pheno-name PHENO_RV --imputeCov \
--covar PHENO_FOR_GWAS_v1_november11_with_PC.txt --freqUpper $FREQUENCY \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr$CHNUM.txt
done

```

###  post burden testing file prepping for meta-analyses...


```
Steps:
- merge per chromosome
- sort based on P-value

cd BURDEN
ALL_MISSENSE
ALL_LOF
ALL_CADD_20
ALL_CADD_10
ALL_MISSENSE_and_LOF

rm *_FREQUENCY_*
264 files per folder 

## big loop....

cd /data/CARD/PD/WGS/june2019/BURDEN/

cat variant_types.txt | while read line
do 
	cd $line
	head -1 PD_WGS_0.05_chr8.CMC.assoc > CMC_header.txt
	head -1 PD_WGS_0.05_chr8.SkatO.assoc > SkatO_header.txt
	# SKATO
	cat *_0.05_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../"$line"_0.05.SkatO.assoc
	cat *_0.01_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../"$line"_0.01.SkatO.assoc
	cat *_0.005_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../"$line"_0.005.SkatO.assoc
	cat *_0.001_*.SkatO.assoc | sort -gk 8 | grep -v "Pvalue" > temp.txt
	cat SkatO_header.txt temp.txt > ../"$line"_0.001.SkatO.assoc
	# CMC
	cat *_0.05_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../"$line"_0.05.CMC.assoc
	cat *_0.01_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../"$line"_0.01.CMC.assoc
	cat *_0.005_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../"$line"_0.005.CMC.assoc
	cat *_0.001_*.CMC.assoc | sort -gk 7 | grep -v "Pvalue" > temp.txt
	cat CMC_header.txt temp.txt > ../"$line"_0.001.CMC.assoc
	cd ..
done

# move files to results file
mkdir RESULTS
mv *.assoc RESULTS/

