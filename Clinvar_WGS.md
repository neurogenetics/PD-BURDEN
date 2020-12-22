# Check for enrichment of ClinVar variants in case-control data...

Cornelis
December 22 2020

General idea:



### Step 1: Find variants and make final lists
```
# search for pathogenic variants....

cd /data/CARD/PD/WGS/june2019/annotation/

grep athogenic PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | cut -f 16 | sort | uniq -c

  18710 Conflicting_interpretations_of_pathogenicity
      1 Conflicting_interpretations_of_pathogenicity,_Affects
     11 Conflicting_interpretations_of_pathogenicity,_other
     11 Conflicting_interpretations_of_pathogenicity,_risk_factor
   1401 Likely_pathogenic
      3 Likely_pathogenic,_risk_factor
   3212 Pathogenic
      3 Pathogenic,_Affects
      1 Pathogenic,_association
      3 Pathogenic,_drug_response
    775 Pathogenic/Likely_pathogenic
      1 Pathogenic/Likely_pathogenic,_drug_response
      1 Pathogenic/Likely_pathogenic,_other
      1 Pathogenic/Likely_pathogenic,_risk_factor
      2 Pathogenic,_other
      4 Pathogenic,_protective
      4 Pathogenic,_risk_factor

# only want:
Conflicting_interpretations_of_pathogenicity
Likely_pathogenic
Pathogenic
Pathogenic/Likely_pathogenic

# header:
cut -f 6,7,9,10,11,12,13,14,15,16,73,85 PD_WGS_anno_chr1.hg38_multianno.withafreq.txt | head -1 > header_for_clinvar.txt

Func.refGene => functional type of variants
Gene.refGene => Gene 
ExonicFunc.refGene => type of mutation
AAChange.refGene => Amino acid change
avsnp150 => rsid if present
CLNALLELEID => clinvar ID
CLNDN => clinvar disease
CLNDISDB => clinvar disease
CLNREVSTAT => Origin of report
CLNSIG => Significant level of gene
Otherinfo6 => variant name
ALT_FREQS => Frequency of allele

# counting lines... 
grep "Conflicting_interpretations_of_pathogenicity" PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | wc -l # 18733
grep -w "Likely_pathogenic" PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | wc -l # 2182
grep -w "Pathogenic" PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | wc -l # 4007
grep -w "Pathogenic/Likely_pathogenic" PD_WGS_anno_chr*.hg38_multianno.withafreq.txt | wc -l # 778

# some manual filtering needed
## Conflicting_interpretations_of_pathogenicity
grep "Conflicting_interpretations_of_pathogenicity" PD_WGS_anno_chr*.hg38_multianno.withafreq.txt \
| grep -v "of_pathogenicity,_" | cut -f 6,7,9,10,11,12,13,14,15,16,73,85 > TEMP.txt

cat header_for_clinvar.txt TEMP.txt > Conflicting_interpretations_of_pathogenicity.txt 
wc -l Conflicting_interpretations_of_pathogenicity.txt # 18711

## Likely_pathogenic
grep "Likely_pathogenic" PD_WGS_anno_chr*.hg38_multianno.withafreq.txt \
| grep -v "Likely_pathogenic,_ri" | grep -v "Pathogenic/Likely" | cut -f 6,7,9,10,11,12,13,14,15,16,73,85 > TEMP.txt

cat header_for_clinvar.txt TEMP.txt > Likely_pathogenic.txt  
wc -l Likely_pathogenic.txt # 1402

## Pathogenic
grep "Pathogenic" PD_WGS_anno_chr*.hg38_multianno.withafreq.txt \
| grep -v "Pathogenic/Li" | grep -v "Pathogenic,_" | cut -f 6,7,9,10,11,12,13,14,15,16,73,85 > TEMP.txt

cat header_for_clinvar.txt TEMP.txt > Pathogenic.txt  
wc -l Pathogenic.txt # 3213

## Pathogenic/Likely_pathogenic
grep "Pathogenic/Likely_pathogenic" PD_WGS_anno_chr*.hg38_multianno.withafreq.txt \
| grep -v "athogenic,_" | cut -f 6,7,9,10,11,12,13,14,15,16,73,85 > TEMP.txt

cat header_for_clinvar.txt TEMP.txt > Pathogenic_Likely_pathogenic.txt  
wc -l Pathogenic_Likely_pathogenic.txt # 776

## Make final variants lists:

1) Conflicting_interpretations_of_pathogenicity.txt # done above
	cut -f 11 Conflicting_interpretations_of_pathogenicity.txt > Conflicting_interpretations_of_pathogenicity_variants.txt
2) Pathogenic.txt # done above
	cut -f 11 Pathogenic.txt > Pathogenic_variants.txt
3) Likely_pathogenic.txt # merging all likely pathogenic ones...
	cat Likely_pathogenic.txt Pathogenic_Likely_pathogenic.txt | cut -f 11 > Likely_pathogenic_variants.txt 
4) Pathogenic.txt and Likely_pathogenic.txt # merging all pathogenic AND likely pathogenic ones...
	cat Pathogenic.txt Likely_pathogenic.txt Pathogenic_Likely_pathogenic.txt | cut -f 11 > ALL_likely_and_pathogenic_variants.txt 
5) All above all_clinvar_of_interest.txt # all four files from above combined...
	cat Likely_pathogenic.txt Pathogenic_Likely_pathogenic.txt Conflicting_interpretations_of_pathogenicity.txt Pathogenic.txt | cut -f 11 > ALL_clinvar_variants.txt 

# counts
  3213 Pathogenic_variants.txt
 24102 ALL_clinvar_variants.txt
  5391 ALL_likely_and_pathogenic_variants.txt
 18711 Conflicting_interpretations_of_pathogenicity_variants.txt
  2178 Likely_pathogenic_variants.txt
```

### Step 2: Prep for burden....

```
cd /data/CARD/PD/WGS/june2019/

mkdir CLINVAR

module load plink/2.0-dev-20191128
module load samtools

# Pathogenic_variants
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  	plink2 --pgen pd.june2019.chr"$chnum".freeze9.sqc.pgen --pvar PVAR_files/NEW"$chnum".pvar \
	--psam pd.june2019.chr"$chnum".freeze9.sqc.psam --extract annotation/Pathogenic_variants.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--out CLINVAR/PD_WGS_Pathogenic_variants_"$chnum" --make-bed
done
# merge files...
cd CLINVAR
ls | grep Pathogenic_variants | grep bim | sed -e 's/.bim//g' > merge_list.txt
module load plink
plink --merge-list merge_list.txt --make-bed --out Pathogenic_variants
module load plink/2.0-dev-20191128
plink2 --bfile Pathogenic_variants --mac 1 --export vcf bgz id-paste=iid --out Pathogenic_variants
tabix -p vcf Pathogenic_variants.vcf.gz
rm PD_WGS_Pathogenic_variants*
# 2718 variants remaining after main filters.

# ALL_clinvar_variants
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  	plink2 --pgen pd.june2019.chr"$chnum".freeze9.sqc.pgen --pvar PVAR_files/NEW"$chnum".pvar \
	--psam pd.june2019.chr"$chnum".freeze9.sqc.psam --extract annotation/ALL_clinvar_variants.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--out CLINVAR/PD_WGS_ALL_clinvar_variants_"$chnum" --make-bed
done
# merge files...
cd CLINVAR
ls | grep ALL_clinvar_variants | grep bim | sed -e 's/.bim//g' > merge_list.txt
module load plink
plink --merge-list merge_list.txt --make-bed --out ALL_clinvar_variants
module load plink/2.0-dev-20191128
plink2 --bfile ALL_clinvar_variants --mac 1 --export vcf bgz id-paste=iid --out ALL_clinvar_variants
tabix -p vcf ALL_clinvar_variants.vcf.gz
rm PD_WGS_ALL_clinvar_variants*
# 19610 variants remaining after main filters.

# ALL_likely_and_pathogenic_variants
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  	plink2 --pgen pd.june2019.chr"$chnum".freeze9.sqc.pgen --pvar PVAR_files/NEW"$chnum".pvar \
	--psam pd.june2019.chr"$chnum".freeze9.sqc.psam --extract annotation/ALL_likely_and_pathogenic_variants.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--out CLINVAR/PD_WGS_ALL_likely_and_pathogenic_variants_"$chnum" --make-bed
done
# merge files...
cd CLINVAR
ls | grep ALL_likely_and_pathogenic_variants | grep bim | sed -e 's/.bim//g' > merge_list.txt
module load plink
plink --merge-list merge_list.txt --make-bed --out ALL_likely_and_pathogenic_variants
module load plink/2.0-dev-20191128
plink2 --bfile ALL_likely_and_pathogenic_variants --mac 1 --export vcf bgz id-paste=iid --out ALL_likely_and_pathogenic_variants
tabix -p vcf ALL_likely_and_pathogenic_variants.vcf.gz
rm PD_WGS_ALL_likely_and_pathogenic_variants*
# 4584 variants remaining after main filters.

# Conflicting_interpretations_of_pathogenicity_variants
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  	plink2 --pgen pd.june2019.chr"$chnum".freeze9.sqc.pgen --pvar PVAR_files/NEW"$chnum".pvar \
	--psam pd.june2019.chr"$chnum".freeze9.sqc.psam --extract annotation/Conflicting_interpretations_of_pathogenicity_variants.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--out CLINVAR/PD_WGS_Conflicting_interpretations_of_pathogenicity_variants_"$chnum" --make-bed
done
# merge files...
cd CLINVAR
ls | grep Conflicting_interpretations_of_pathogenicity_variants | grep bim | sed -e 's/.bim//g' > merge_list.txt
module load plink
plink --merge-list merge_list.txt --make-bed --out Conflicting_interpretations_of_pathogenicity_variants
module load plink/2.0-dev-20191128
plink2 --bfile Conflicting_interpretations_of_pathogenicity_variants --mac 1 --export vcf bgz id-paste=iid --out Conflicting_interpretations_of_pathogenicity_variants
tabix -p vcf Conflicting_interpretations_of_pathogenicity_variants.vcf.gz
rm PD_WGS_Conflicting_interpretations_of_pathogenicity_variants*
# 15026 variants remaining after main filters.

# Likely_pathogenic_variants
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  	plink2 --pgen pd.june2019.chr"$chnum".freeze9.sqc.pgen --pvar PVAR_files/NEW"$chnum".pvar \
	--psam pd.june2019.chr"$chnum".freeze9.sqc.psam --extract annotation/Likely_pathogenic_variants.txt \
	--keep PHENO_FOR_GWAS_v1_november11_with_PC.txt \
	--out CLINVAR/PD_WGS_Likely_pathogenic_variants_"$chnum" --make-bed
done
# merge files...
cd CLINVAR
ls | grep Likely_pathogenic_variants | grep bim | sed -e 's/.bim//g' > merge_list.txt
module load plink
plink --merge-list merge_list.txt --make-bed --out Likely_pathogenic_variants
module load plink/2.0-dev-20191128
plink2 --bfile Likely_pathogenic_variants --mac 1 --export vcf bgz id-paste=iid --out Likely_pathogenic_variants
tabix -p vcf Likely_pathogenic_variants.vcf.gz
rm PD_WGS_Likely_pathogenic_variants*
# 1866 variants remaining after main filters.

-----
Moving forward with two main files:
ALL_clinvar_variants => 19610 variants
ALL_likely_and_pathogenic_variants => 4584 variants


-----
# start testing (sanity)

cd /data/CARD/PD/WGS/june2019

module load rvtests # 2.1.0 
# ALL_clinvar_variants
rvtest --noweb --hide-covar --out BURDEN/CLINVAR/GBA_TEST_ALL_clinvar_variants --burden cmc  \
--inVcf CLINVAR/ALL_clinvar_variants.vcf.gz \
--pheno PHENO_FOR_GWAS_v1_november11_with_PC.txt --pheno-name PHENO_RV \
--covar PHENO_FOR_GWAS_v1_november11_with_PC.txt --freqUpper 0.05 --imputeCov \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# 11 variants => P=6.07555e-07

# ALL_likely_and_pathogenic_variants
rvtest --noweb --hide-covar --out BURDEN/CLINVAR/GBA_TEST_ALL_likely_and_pathogenic_variants --burden cmc  \
--inVcf CLINVAR/ALL_likely_and_pathogenic_variants.vcf.gz \
--pheno PHENO_FOR_GWAS_v1_november11_with_PC.txt --pheno-name PHENO_RV \
--covar PHENO_FOR_GWAS_v1_november11_with_PC.txt --freqUpper 0.05 --imputeCov \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --geneFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/refFlat_HG38_chr1.txt --gene GBA
# 9 variants => P=0.000479568

----
```

### Step 3: creating pathway files...
```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/

mkdir DISGENET
cd DISGENET
# https://www.disgenet.org/downloads
wget https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz

less curated_gene_disease_associations.tsv.gz
gunzip curated_gene_disease_associations.tsv.gz

cut -f 2,6 curated_gene_disease_associations.tsv | sed -e 's/ /_/g' > DISGENET.txt

# reshaping REFFLAT format
# required format...
# set1 1:100-200,1:250-300
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/
gunzip refFlat.txt.gz
cut -f 1 refFlat.txt > TEMP.txt
cut -f 3,5 refFlat.txt | sed -e 's/\t/:/g' > TEMP2.txt
cut -f 6 refFlat.txt > TEMP3.txt
paste TEMP2.txt TEMP3.txt | sed -e 's/\t/-/g' > TEMP4.txt
paste TEMP.txt TEMP4.txt > REFFLAT_reshaped.txt
mv REFFLAT_reshaped.txt DISGENET/

# munging further...
cd DISGENET

module load R
R
digenet <- read.table("DISGENET.txt",header=T)
REFFLAT <- read.table("REFFLAT_reshaped.txt",header=F)
dim(digenet)
#[1] 84038     2
dim(REFFLAT)
#[1] 88819     2
MM <- merge(digenet, REFFLAT, by.x="geneSymbol", by.y="V1")
# remove duplicat lines
MM2 <- unique(MM)
write.table(MM2, file="TEST123.txt", quote=FALSE,row.names=F,sep="\t")
------
# comparing lists...
cut -f 1 DISGENET.txt | sort -u | wc -l
# 9704 => unique genes in input
cut -f 1 TEST123.txt | sort -u | wc -l
# 9563  => unique genes in output so not bad overall...
# = 98.5% will take it
------
# Make final files...
## remove chr from names
cut -f 1,2 TEST123.txt > temp
cut -f 3 TEST123.txt | sed -e 's/chr//g' > temp2
paste temp temp2 > DISGENET_full.txt
# reformat
cut -f 2,3 DISGENET_full.txt > DISGENET_disease_file.txt
# remove GBA and LRRK2
grep -v -w GBA DISGENET_full.txt | grep -v LRRK2 | cut -f 2,3 > DISGENET_disease_file_no_LRRK2_GBA.txt
-------------------
```


### Step 4: # Start testing...
```
cd /data/CARD/PD/WGS/june2019/

# two input .vcf
ALL_clinvar_variants.vcf.gz => 19610 variants
ALL_likely_and_pathogenic_variants.vcf.gz => 4584 variants
# two disease files
DISGENET_disease_file_no_LRRK2_GBA.txt
DISGENET_disease_file.txt

# Loaded 2788 cases, 4105 controls, and 0 missing phenotypes
# Loaded 11116 set to tests

module load rvtests # 2.1.0 
rvtest --noweb --hide-covar --out BURDEN/CLINVAR/ALL_clinvar_variants_full_DISGENET_disease_file --burden cmc \
--inVcf CLINVAR/ALL_clinvar_variants.vcf.gz \
--pheno PHENO_FOR_GWAS_v1_november11_with_PC.txt --pheno-name PHENO_RV \
--covar PHENO_FOR_GWAS_v1_november11_with_PC.txt --freqUpper 0.05 --imputeCov  \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --setFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/DISGENET/DISGENET_disease_file.txt 

rvtest --noweb --hide-covar --out BURDEN/CLINVAR/ALL_clinvar_variants_no_LRRK2_GBA_DISGENET_disease_file --burden cmc \
--inVcf CLINVAR/ALL_clinvar_variants.vcf.gz \
--pheno PHENO_FOR_GWAS_v1_november11_with_PC.txt --pheno-name PHENO_RV \
--covar PHENO_FOR_GWAS_v1_november11_with_PC.txt --freqUpper 0.05 --imputeCov \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --setFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/DISGENET/DISGENET_disease_file_no_LRRK2_GBA.txt 

rvtest --noweb --hide-covar --out BURDEN/CLINVAR/ALL_likely_and_pathogenic_variants_full_DISGENET_disease_file --burden cmc --kernel skato \
--inVcf CLINVAR/ALL_likely_and_pathogenic_variants.vcf.gz \
--pheno PHENO_FOR_GWAS_v1_november11_with_PC.txt --pheno-name PHENO_RV \
--covar PHENO_FOR_GWAS_v1_november11_with_PC.txt --freqUpper 0.05 --imputeCov  \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --setFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/DISGENET/DISGENET_disease_file.txt 

rvtest --noweb --hide-covar --out BURDEN/CLINVAR/ALL_likely_and_pathogenic_variants_no_LRRK2_GBA_DISGENET_disease_file --burden cmc --kernel skato \
--inVcf CLINVAR/ALL_likely_and_pathogenic_variants.vcf.gz \
--pheno PHENO_FOR_GWAS_v1_november11_with_PC.txt --pheno-name PHENO_RV \
--covar PHENO_FOR_GWAS_v1_november11_with_PC.txt --freqUpper 0.05 --imputeCov \
--covar-name SEX,AGE_ANALYSIS,PC1,PC2,PC3,PC4,PC5 --setFile /data/CARD/UKBIOBANK/EXOME_DATA_200K/REFFLAT/DISGENET/DISGENET_disease_file_no_LRRK2_GBA.txt 

# inspect result files...
cut -f 1,3,4,5,6,7 BURDEN/CLINVAR/ALL_clinvar_variants_full_DISGENET_disease_file.CMC.assoc | grep -v nan | sort -gk 6 | head 
cut -f 1,3,4,5,6,7 BURDEN/CLINVAR/ALL_clinvar_variants_no_LRRK2_GBA_DISGENET_disease_file.CMC.assoc | grep -v nan | sort -gk 6 | head 
cut -f 1,3,4,5,6,7 BURDEN/CLINVAR/ALL_likely_and_pathogenic_variants_full_DISGENET_disease_file.CMC.assoc | grep -v nan | sort -gk 6 | head 
cut -f 1,3,4,5,6,7 BURDEN/CLINVAR/ALL_likely_and_pathogenic_variants_no_LRRK2_GBA_DISGENET_disease_file.CMC.assoc | grep -v nan | sort -gk 6 | head 

# create final result files
cut -f 1,3,4,5,6,7 BURDEN/CLINVAR/ALL_clinvar_variants_full_DISGENET_disease_file.CMC.assoc | grep -v nan | sort -gk 6 > BURDEN/CLINVAR/
cut -f 1,3,4,5,6,7 BURDEN/CLINVAR/ALL_clinvar_variants_no_LRRK2_GBA_DISGENET_disease_file.CMC.assoc | grep -v nan | sort -gk 6 > BURDEN/CLINVAR/
cut -f 1,3,4,5,6,7 BURDEN/CLINVAR/ALL_likely_and_pathogenic_variants_full_DISGENET_disease_file.CMC.assoc | grep -v nan | sort -gk 6 > BURDEN/CLINVAR/
cut -f 1,3,4,5,6,7 BURDEN/CLINVAR/ALL_likely_and_pathogenic_variants_no_LRRK2_GBA_DISGENET_disease_file.CMC.assoc | grep -v nan | sort -gk 6 > BURDEN/CLINVAR/
```
