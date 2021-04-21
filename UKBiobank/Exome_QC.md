```
cat ukb23156_*.stats | grep -e "ukb23156_" -e "records:" > ../STATS/overview_of_vcf.gz
cat filtered_*.stats | grep -e "filtered_" -e "records:" > ../STATS/overview_of_filtered_vcf.gz
```
