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

./plink2 --pmerge-list MERGE_LIST.txt --make-pgen --memory 99000 \
--threads 10 --out MERGED_UKB_first_pass


```
