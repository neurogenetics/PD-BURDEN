```
# create variant overview
cat ukb23156_*.stats | grep -e "ukb23156_" -e "records:" > ../STATS/overview_of_vcf.gz
cat filtered_*.stats | grep -e "filtered_" -e "records:" > ../STATS/overview_of_filtered_vcf.gz

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


```

```
#annotation

grep -v "#" MERGED_UKB_first_pass.pvar | cut -f 1,2 > temp1
grep -v "#" MERGED_UKB_first_pass.pvar | cut -f 2 > temp2
grep -v "#" MERGED_UKB_first_pass.pvar | cut -f 3,7 > temp3
grep -v "#" MERGED_UKB_first_pass.pvar | cut -f 4,5 > temp4
paste temp1 temp2 temp4 temp3 > to_annotate_pass1.txt

#!/bin/bash
# sbatch --cpus-per-task=10 --mem=200g --time=24:00:00 annotate_UKB_plink2_format_June2021.sh

module load annovar 
table_annovar.pl to_annotate_pass1.txt $ANNOVAR_DATA/hg38 --thread 16 -buildver hg38 \
-out TESTING -remove --otherinfo -polish -protocol refGene,avsnp150,clinvar_20200316 \
-operation g,f,f -nastring . 


---- below is old...



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

#!/bin/bash
# sbatch --cpus-per-task=10 --mem=200g --time=24:00:00 annotate_UKB_plink2_format.sh $argument

module load annovar
cat to_merge_combined_filtered_c13_b17.pvar | awk '{print $1,$2,$2,$4,$5,$3}\' OFS="\t"| grep -v '##' > to_merge_combined_filtered_c13_b17.avinput

table_annovar.pl to_merge_combined_filtered_c13_b17.avinput $ANNOVAR_DATA/hg38 --thread 16 -buildver hg38 \
-out TESTING -remove --otherinfo -polish -protocol refGene,avsnp150,clinvar_20200316,dbnsfp41a \
-operation g,f,f,f -nastring . 


