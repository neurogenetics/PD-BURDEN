## PD WGS data burden testing

```
November 2020
Cornelis

Steps:
- sort out phenotype file
- rename variants
- annotation
- QC (relatedness, ancestry, PC's)
- subset data
- burden testing
- post burden testing file prepping for meta-analyses
- create cumulative frequency list for each group
```

### Phenotype file

```
/data/CARD/PD/GENOMES/august19/genotypes/PHENO_SEPT.txt

merge with psam file from WGS data => pd.june2019.chr*.freeze9.sqc.psam

2810 PD cases
4207 controls
--------------
7017 total

This is excluding pop outliers, genetic registry, other diagnoses (prodomal, MSA, PSP etc)
Note that this does not include relatedness filtering step….

Short sample list saved as:
PHENO_FOR_GWAS_v1_november11.txt
```

### Update variant names in WGS data 

some variants are rs1378555846 instead of 22:10637358:G:A format

```
# First check the size of the header of each pvar file...

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
	grep "#" pd.june2019.chr$chnum.freeze9.sqc.pvar | wc -l
done
# seems to be 31 rows for chr1-22 and 19 rows for chr23

# create stable backup for pvar files:
mkdir PVAR_files
scp pd.june2019.chr*.freeze9.sqc.pvar PVAR_files/
cd PVAR_files
rm pd.june2019.chrX.freeze9.sqc.pvar

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
	grep -v "#" pd.june2019.chr"$chnum".freeze9.sqc.pvar > no_header.txt
	grep "#" pd.june2019.chr"$chnum".freeze9.sqc.pvar > header.txt
	cut -f 1,2,4,5 no_header.txt | sed -e 's/\t/:/g' > variant_names.txt
	cut -f 1,2 no_header.txt > part1
	cut -f 4,5,6,7,8 no_header.txt > part2
	paste part1 variant_names.txt part2 > NEW_no_header.txt
	cat header.txt NEW_no_header.txt > NEW"$chnum".pvar
	wc -l pd.june2019.chr"$chnum".freeze9.sqc.pvar
	wc -l NEW"$chnum".pvar
	wc -l part1
	wc -l part2
	wc -l NEW_no_header.txt
	wc -l variant_names.txt
done
---------
# The above option is OK-ish... cleaner option would be --update-name  in plink... but this doesnt work because multi-allelics in plink mess this up...
similar as above...
# down side here though is that indels have not the same format as UKB.. 22:10687003:C:CTATA vs 22:10687003:C:D4

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
	plink2 --pgen pd.june2019.chr"$chnum".freeze9.sqc.pgen --pvar PVAR_files/NEW"$chnum".pvar --psam pd.june2019.chr"$chnum".freeze9.sqc.psam --freq
	wc -l plink2.afreq >> length.txt
done

# OK done... its possible to use these moving forward using:
plink2 --pgen <filename> --pvar <filename> --psam <filename>

```

### annotation

```
scp hg38.exome_calling_regions.interval_list /data/CARD/PD/WGS/

working folder => /data/CARD/PD/WGS/june2019
```

```
Subset data keeping only coding regions using hg38.exome_calling_regions.interval_list

# This section briefly describes how to filter for exome calling regions only based on gnomAD file
exome_calling_regions.interval_list -> from broad gnomad
wget https://storage.googleapis.com/gnomad-public/intervals/exome_calling_regions.v1.interval_list
convert to hg38...

module load crossmap
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# add chr
awk '{ print "chr",$0 }' exome_calling_regions.interval_list > temp
sed -i 's/chr /chr/g' temp
crossmap bed hg19ToHg38.over.chain.gz temp > temp4.txt
grep -v "Fail" temp4.txt | cut -f 7-11 | grep -v "v1_alt" | grep -v "v2_random" | grep -v "v1_random" > hg38.exome_calling_regions.interval_list
# tested using UK Biobank annotated exome data and >99% of coding variants were retained so works
```


```
cd /data/CARD/PD/WGS/june2019/

module load plink/2.0-dev-20191128

plink2 --pfile pd.june2019.chr22.freeze9.sqc --extract range ../hg38.exome_calling_regions.interval_list --freq

# update chromosome X files
scp pd.june2019.chrX.freeze9.sqc.pgen pd.june2019.chr23.freeze9.sqc.pgen 
scp pd.june2019.chrX.freeze9.sqc.psam pd.june2019.chr23.freeze9.sqc.psam
scp pd.june2019.chrX.freeze9.sqc.pvar pd.june2019.chr23.freeze9.sqc.pvar

# create frequency files
module load plink/2.0-dev-20191128
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
  	plink2 --pgen pd.june2019.chr"$chnum".freeze9.sqc.pgen --pvar PVAR_files/NEW"$chnum".pvar \
	--psam pd.june2019.chr"$chnum".freeze9.sqc.psam --freq --out FREQchr"$chnum" \
    	--extract range ../hg38.exome_calling_regions.interval_list
done


# selecting one test sample to make annotation faster...
echo BF-1101 BF-1101 > textsample.txt

# subset 1 one sample
module load plink/2.0-dev-20191128
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
  	plink2 --pgen pd.june2019.chr"$chnum".freeze9.sqc.pgen --pvar PVAR_files/NEW"$chnum".pvar \
	--psam pd.june2019.chr"$chnum".freeze9.sqc.psam --export vcf id-paste=iid \
    --out anno"$chnum" --keep textsample.txt \
    --extract range ../hg38.exome_calling_regions.interval_list
done
```

```
Running annovar....

# annotate all
module load annovar
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
  	table_annovar.pl anno"$chnum".vcf $ANNOVAR_DATA/hg38 --thread 16 -buildver hg38 \
	-out PD_WGS_anno_chr"$chnum" -remove -protocol refGene,avsnp150,clinvar_20200316,dbnsfp41a \
	-operation g,f,f,f -nastring . -vcfinput
done

```

```

# move data to annotation folder

mv PD_WGS_anno_chr* annotation/
mv FREQchr*afreq annotation/

# clean up folder
cd annotation/
mkdir inputfiles
mv *.avinput inputfiles/
rm *.vcf
# merge data with frequency file...

for chnum in {1..23};
  do
	paste PD_WGS_anno_chr"$chnum".hg38_multianno.txt FREQchr"$chnum".afreq > PD_WGS_anno_chr"$chnum".hg38_multianno.withafreq.txt
done
```

```
### subset "groups" of variants...

# missense
grep exonic PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | grep nonsynonymous | cut -f 73 > ALL_MISSENSE.txt
# n=933810

# LOF (stop, frame)
grep stopgain PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | cut -f 73 > all_stopgain.txt
# n=24960
grep stoploss PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | cut -f 73 > all_stoploss.txt
# n=1193
grep nonframeshift PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | cut -f 73 > all_nonframeshift.txt
# n=16597
grep frame PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | grep -v nonframeshift | grep -v nonsynonymous | cut -f 73 > all_frameshift.txt
# n=26633

# splicing 
grep splicing PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | \
grep -v ncRNA | cut -f 6,73 | grep splicing | cut -f 2 > all_splice_bp.txt
# 13474 splicing

# CADD <10
awk '{ if($42 > 10) { print }}' PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | cut -f 73 > ALL_CADD_10.txt
# n=46372

# CADD <20
awk '{ if($42 > 20) { print }}' PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | cut -f 73 > ALL_CADD_20.txt
# n=38390

#### Prepping final files:

cat all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice_bp.txt > ALL_LOF.txt

cat ALL_MISSENSE.txt all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice_bp.txt > ALL_MISSENSE_and_LOF.txt

### FINAL GROUPS:

ALL_CADD_10.txt
ALL_CADD_20.txt
ALL_MISSENSE_and_LOF.txt
ALL_LOF.txt
ALL_MISSENSE.txt

```

### QC (relatedness, ancestry)

```
cd /data/CARD/PD/WGS/june2019
module load plink/2.0-dev-20191128

# subset samples
mkdir PCA
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
 	plink2 --pfile pd.june2019.chr$chnum.freeze9.sqc \
	--remove qc/pd.june2019.pop.outliers.txt \
	--make-bed --out PCA/geno"$chnum"
done

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


# CONTINUE HERE

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

40 files => 2 tests, 4 variant thresholds, 5 variant groups => 40

### quick check of hits...

cat *.SkatO.assoc | sort -gk 8 > all_skat.txt
cat *.CMC.assoc | sort -gk 7 > all_cmc.txt

# inspect using excel
filter => at least 5 variants "NumVar"

# CMC all genes with E-5:
GBA => good
ATP9B => no link with PD
CDC20 => PARKIN partner
LARS2 => some far fetch role for methylation maybe? and Perrault syndrome?
GALNT18 => no link with PD
BCKDHA => no link with PD

# SKAT-O all genes with E-5:
GBA => good
LRRK2 => good
RPUSD1 => no link with PD
ATP9B => no link with PD
CCDC73 => no strong link wtih PD
CDC20 => PARKIN partner
DCDC2C => restless leg GWAS hit

```


### create cumulative frequency list for each group...

```
Input files needed:
/data/CARD/PD/WGS/june2019/annotation/PD_WGS_chr*.hg38_multianno.withafreq.txt

Variant classes:
all_missense.txt
ALL_LOF.txt
ALL_CADD_20.txt
ALL_CADD_10.txt
ALL_MISSENSE_and_LOF.txt

Case-control files to generate:
module load plink/2.0-dev-20191128

mkdir freq_case_control_per_variant_class
# make plink1 files
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
 	plink2 --pfile pd.june2019.chr"$chnum".freeze9.sqc \
	--extract annotation/ALL_LOF.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--pheno PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--pheno-name PHENO_RV --make-bed \
	--out freq_case_control_per_variant_class/PD_WGS_ALL_LOF_"$chnum"
done
# creating assoc files...
module load plink
cd freq_case_control_per_variant_class
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
plink --bfile PD_WGS_ALL_LOF_"$chnum" \
--assoc --out PD_WGS_ALL_LOF_freq_"$chnum"
done

# prep annotation...
cut -f 1-11,92 UKB_exomes_200K_chr9.hg38_multianno.txt > SHORT.txt

# then merge in R and filter....
module load R
R
require("data.table")
assoc <- read.table("PD_WGS_ALL_LOF_freq_9.assoc",header=T)
assoc = assoc[c("CHR","SNP","F_A","F_U")]
#annotation <- read.table("../annotation/SHORT.txt",header=T)
annotation <- fread("../annotation/SHORT.txt",header=T)
dim(assoc)
dim(annotation)
MM <- merge(assoc,annotation,by.x="SNP",by.y="Otherinfo6")
dim(MM)
# now calculate cumulative case-control frequency
case_feq <- aggregate(x = MM$F_A, by = list(MM$Gene.refGene),FUN = sum)   
control_feq <- aggregate(x = MM$F_U, by = list(MM$Gene.refGene),FUN = sum)   
names(case_feq)[1] <- "GENE"
names(case_feq)[2] <- "case_freq"
names(control_feq)[1] <- "GENE"
names(control_feq)[2] <- "control_freq"
MM2 <- merge(case_feq,control_feq,by="GENE")
write.table(MM2, file="Case_control_cumulative_MAF.txt", quote=FALSE,row.names=F,sep="\t")
## output format:
#    GENE case_freq control_freq
#   ABCA1 0.0005379    0.0006090
#   ABCA2 0.0005379    0.0001218
q()
n

## merge files with result files...

work in progress....


```

