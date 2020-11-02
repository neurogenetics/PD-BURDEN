### PD WGS data

November 2020
Cornelis

Steps:
- annotation
- QC (relatedness, ancestry)
- subset data
- burden testing

#### annotation

scp hg38.exome_calling_regions.interval_list /data/CARD/PD/WGS/

working folder => /data/CARD/PD/WGS/june2019

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
  	plink2 --pfile pd.june2019.chr$chnum.freeze9.sqc --freq --out FREQchr"$chnum" \
    --extract range ../hg38.exome_calling_regions.interval_list
done

# selecting one test sample to make annotation faster...
echo BF-1101 BF-1101 > textsample.txt

# subset 1 one sample
module load plink/2.0-dev-20191128
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
  	plink2 --pfile pd.june2019.chr$chnum.freeze9.sqc --export vcf id-paste=iid \
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
	-out UKB_exomes_200K_chr"$chnum" -remove -protocol refGene,avsnp150,clinvar_20200316,dbnsfp35a \
	-operation g,f,f,f -nastring . -vcfinput
done

# move data to annotation folder

mv UKB_exomes_200K_chr* ../annotation_of_plink_files/
mv FREQchr*afreq ../annotation_of_plink_files/

# clean up folder
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/annotation_of_plink_files/
mkdir inputfiles
mv *.avinput inputfiles/
rm *.vcf 


```

```
### subset "groups" of variants...

# missense
grep exonic UKB_exomes_200K_chr*.hg38_multianno.withafreq.txt | grep nonsynonymous | cut -f 92 > all_missense.txt
# n=4672378

# LOF (stop, frame)
grep stopgain UKB_exomes_200K_chr*.hg38_multianno.withafreq.txt > all_stopgain.txt
# n=159685
grep stoploss UKB_exomes_200K_chr*.hg38_multianno.withafreq.txt > all_stoploss.txt
# n=7085
grep nonframeshift UKB_exomes_200K_chr*.hg38_multianno.withafreq.txt > all_nonframeshift.txt
# n=107106
grep frame UKB_exomes_200K_chr*.hg38_multianno.withafreq.txt | grep -v nonframeshift | grep -v nonsynonymous > all_frameshift.txt
# n=203994

# splicing 
grep splicing UKB_exomes_200K_chr*.hg38_multianno.withafreq.txt | \
grep -v ncRNA | cut -f 6,92 | grep splicing | cut -f 2 > all_splice_15bp.txt
# 24852 exonic;splicing
# 1134001 splicing


# CADD <10
awk '{ if($58 > 10) { print }}' UKB_exomes_200K_chr*.hg38_multianno.withafreq.txt | cut -f 92 > ALL_CADD_10.txt
# n=265920

# CADD <20
awk '{ if($58 > 20) { print }}' UKB_exomes_200K_chr*.hg38_multianno.withafreq.txt | cut -f 92 > ALL_CADD_20.txt
# n=228934


#### Prepping final files:

cat all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice_normal.txt > ALL_LOF.txt

cat all_missense.txt all_frameshift.txt all_stopgain.txt all_stoploss.txt all_splice_normal.txt > ALL_MISSENSE_and_LOF.txt

### FINAL GROUPS:

ALL_CADD_10.txt
ALL_CADD_20.txt
all_frameshift.txt
all_missense.txt
all_nonframeshift.txt


```

#### QC (relatedness, ancestry)

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
# 158317
```



####  subset data




####  burden testing
