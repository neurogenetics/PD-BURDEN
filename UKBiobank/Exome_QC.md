```
# create variant overview
cat ukb23156_*.stats | grep -e "ukb23156_" -e "records:" > ../STATS/overview_of_vcf_new.gz
cat filtered_*.stats | grep -e "filtered_" -e "records:" > ../STATS/overview_of_filtered_vcf_new.gz

```

```
## make merged plink file
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20210420.zip
unzip plink2_linux_avx2_20210420.zip

ls | grep to_merge_combined_filtered | grep psam | sed -e 's/.psam//g' > MERGE_LIST.txt

OR

ls | grep to_merge_combined_filtered | grep psam | sed -e 's/.psam//g' | grep -v "X" | grep -v "Y" > MERGE_LIST_autosomes.txt

./plink2 --pmerge-list MERGE_LIST.txt --make-pgen --memory 99000 \
--threads 10 --out MERGED_UKB_first_pass

AND

ls | grep to_merge_combined_filtered | grep psam | sed -e 's/.psam//g' | grep "X" > MERGE_LIST_X_chromosome.txt

./plink2 --pmerge-list MERGE_LIST_X_chromosome.txt --make-pgen --memory 99000 \
--threads 10 --out MERGED_UKB_chromosome_X


```

```
# sex check?
./plink2 --pfile MERGED_UKB_chromosome_X --make-bed --chr 23 --from-bp 2781479 --to-bp 156030895 --out first_pass_sex_check
module load plink
plink --bfile first_pass_sex_check --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out gender_check2 

module load R
R
genetic <- read.table("gender_check2.sexcheck",header=T)
reported <- read.table("/data/CARD/UKBIOBANK/PHENOTYPE_DATA/covariates_phenome_to_use.txt",header=T)
MM <- merge(genetic,reported,by.x="IID",by.y="FID")
write.table(MM, file="../COMPARISON_genetic_reported_sex.txt",quote=F,row.names=F,sep="\t")

Comparison results....

SNPSEX => 2
0	100084
1	2476
SNPSEX => 1
0	466
1	50648
SNPSEX => 0
0	10061
1	36895

Some discrapancies but overall looks good...

```

```
#annotation

-----------------

#!/bin/bash
# sbatch --cpus-per-task=10 --mem=10g --time=1:00:00 annotate_UKB_200K_vcf_annotation_July2021.sh ukb23156_c10_b15_v1
# load packages
module load samtools
module load annovar
module load bcftools/1.9
# variables
FILE=$1
echo "this is"
echo $FILE 
# previously done
# ls | grep ukb23156_ | grep _v1.vcf.gz > PVCF_TO_LOOP_OVER.txt
# start here
zcat "$FILE".vcf.gz | cut -f 1-12 | bgzip -c > "$FILE".for.anno.vcf.gz

# for SNPs and 0.20 for indels.
echo "split multiallelics, and get rid of bad variants"
bcftools norm -m- "$FILE".for.anno.vcf.gz > filtered_"$FILE".for.anno.vcf
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' filtered_"$FILE".for.anno.vcf > filtered_"$FILE".for.anno.names.vcf

echo "starting annovar"

table_annovar.pl filtered_"$FILE".for.anno.names.vcf $ANNOVAR_DATA/hg38 --thread 16 -buildver hg38 \
-out VCF_annotation/"$FILE" -remove --otherinfo -polish -protocol refGene,avsnp150,clinvar_20200316 \
-operation g,f,f -nastring . -vcfinput

echo "done"

-----------------

## to start massive loop

cat PVCF_TO_LOOP_OVER.txt  | while read line
do 
   sbatch --cpus-per-task=10 --mem=10g --time=1:00:00 annotate_UKB_200K_vcf_annotation_July2021.sh $line
done
   
```

```
## merge annotation

cd VCF_annotation/
ls | wc -l
2931
# 2931 / 3 = 977 -> makes sense..

cat ukb23156_c*v1.hg38_multianno.txt > MERGED_annotation.txt
wc -l MERGED_annotation.txt
# 17982874 MERGED_annotation.txt

cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/
# autosomes + X
cat to_merge_combined_filtered_c*_b*.pvar | grep -v "#" | grep -v "Y" > MERGED_PLINK.txt
grep -v "##" to_merge_combined_filtered_c19_b45.pvar | grep "#" > header_merged_plink.txt
cat header_merged_plink.txt MERGED_PLINK.txt > MERGED_PLINK_with_header.txt
sed -i -e 's/:/_/g' MERGED_PLINK_with_header.txt
# Y only
cat to_merge_combined_filtered_cY_b*.pvar | grep -v "#" > MERGED_PLINK_Y.txt
grep -v "##" to_merge_combined_filtered_cY_b0.pvar | grep "#" > header_merged_plink_Y.txt
cat header_merged_plink_Y.txt MERGED_PLINK_Y.txt > MERGED_PLINK_Y_with_header.txt
sed -i -e 's/:/_/g' MERGED_PLINK_Y_with_header.txt

cd VCF_annotation/

module load R
R
require(data.table)
anno <- fread("MERGED_annotation.txt",header=T)
anno_Y <- fread("ukb23156_cY_b0_v1.hg38_multianno.txt",header=T)
plink <- fread("../MERGED_PLINK_with_header.txt",header=T)
plink_Y <- fread("../MERGED_PLINK_Y_with_header.txt",header=T)
MM <- merge(plink,anno,by.x="ID",by.y="Otherinfo6")
MM_Y <- merge(plink_Y,anno_Y,by.x="ID",by.y="Otherinfo6")
write.table(MM, file="MERGED_annotation_added_pvar.txt",quote=F,row.names=F,sep="\t")
write.table(MM_Y, file="MERGED_annotation_added_pvar_Y.txt",quote=F,row.names=F,sep="\t")

```

```
Done for now...

```


```
Side-job for Derek...
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/
grep PRKN MERGED_annotation.txt > PRKN_variants_Derek.txt > 
head -1 MERGED_annotation.txt > header.txt
cat header.txt PRKN_variants_Derek.txt > PRKN_variants_Derek_with_header.txt

scp blauwendraatc@biowulf.nih.gov:///data/CARD/UKBIOBANK/EXOME_DATA_200K/PVCF_FILES/VCF_annotation/PRKN_variants_Derek_with_header.txt /Users/blauwendraatc/Desktop/Phase1_INDI/ 
```

