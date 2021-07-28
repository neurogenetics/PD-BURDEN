# Other annotation tools

## 1. snpEff/snpSift

Biowulf documentation: https://hpc.nih.gov/apps/snpEff.html

```
cd /data/CARD/PD/AMP_NIH/no_relateds
mkdir test_snpEff
cd test_snpEff

# Copy some example data
cp /data/CARD/PD/AMP_NIH/PD_PHENO/subset_NIH_ONLY_forANNOVAR_chr22.vcf .

module load snpEff
java -Xmx20g -jar $SNPEFF_JAR -v GRCh38.86 subset_NIH_ONLY_forANNOVAR_chr22.vcf > subset_NIH_ONLY_forANNOVAR_chr22.eff.vcf

## Can decide to filter before annotation with snpSift
cat subset_NIH_ONLY_forANNOVAR_chr22.eff.vcf| java -jar $SNPSIFT_JAR filter "((EFF[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE'))" > subset_NIH_ONLY_forANNOVAR_chr22.filtered.vcf

java -jar $SNPSIFT_JAR dbnsfp -v -db /fdb/dbNSFP2/dbNSFP3.2a.txt.gz subset_NIH_ONLY_forANNOVAR_chr22.eff.vcf > subset_NIH_ONLY_forANNOVAR_chr22.annotated.vcf

# See here for all the fields you can extract: https://pcingola.github.io/SnpEff/ss_extractfields/
echo "CHROM
POS
ID
REF
ALT
FILTER
AF
AC
DP
MQ
ANN[*].ALLELE
ANN[*].EFFECT
ANN[*].IMPACT
ANN[*].GENE
ANN[*].GENEID
ANN[*].FEATURE
ANN[*].FEATUREID
ANN[*].BIOTYPE
ANN[*].RANK
ANN[*].HGVS_C
ANN[*].HGVS_P
ANN[*].CDNA_POS
ANN[*].CDNA_LEN
ANN[*].CDS_POS
ANN[*].CDS_LEN
ANN[*].AA_POS
ANN[*].AA_LEN
ANN[*].DISTANCE
ANN[*].ERRORS
LOF[*].GENE
LOF[*].GENEID
LOF[*].NUMTR
LOF[*].PERC
NMD[*].GENE
NMD[*].GENEID
NMD[*].NUMTR
NMD[*].PERC" > snpSift_extract_fields.txt

# Use tr to reformat these fields into a single line with space separation -> can be inputted to snpSift
# Using ANN[*] for example doesn't output the fields in separate columns and is difficult to parse through
java -jar $SNPSIFT_JAR extractFields -s "," -e "." subset_NIH_ONLY_forANNOVAR_chr22.annotated.vcf $(tr '\n' ' ' < snpSift_extract_fields.txt) > subset_NIH_ONLY_forANNOVAR_chr22.annotated.txt

head -2 subset_NIH_ONLY_forANNOVAR_chr22.annotated.txt | cut -f1-10
# CHROM	POS	ID	REF	ALT	FILTER	AF	AC	DP	MQ
# 22	15528133	chr22:15528133:A:G	A	G	PASS	2.643E-5	1	602382	.

head -2 subset_NIH_ONLY_forANNOVAR_chr22.annotated.txt | cut -f11-20
# ANN[*].ALLELE	ANN[*].EFFECT	ANN[*].IMPACT	ANN[*].GENE	ANN[*].GENEID	ANN[*].FEATURE	ANN[*].FEATUREID	ANN[*].BIOTYPE	ANN[*].RANK	ANN[*].HGVS_C
# G,G	upstream_gene_variant,intergenic_region	MODIFIER,MODIFIER	OR11H1,YME1L1P1-OR11H1	ENSG00000130538,ENSG00000236831-ENSG00000130538	transcript,intergenic_region	ENST00000252835.4,ENSG00000236831-ENSG00000130538	protein_coding,.	-1,-1	c.-26A>G,n.15528133A>G

head -2 subset_NIH_ONLY_forANNOVAR_chr22.annotated.txt | cut -f21-37
# ANN[*].HGVS_P	ANN[*].CDNA_POS	ANN[*].CDNA_LEN	ANN[*].CDS_POS	ANN[*].CDS_LEN	ANN[*].AA_POS	ANN[*].AA_LEN	ANN[*].DISTANCE	ANN[*].ERRORS	LOF[*].GENE	LOF[*].GENEID	LOF[*].NUMTR	LOF[*].PERC	NMD[*].GENE	NMD[*].GENEID	NMD[*].NUMTR	NMD[*].PERC
# .,.	-1,-1	-1,-1	-1,-1	-1,-1	-1,-1	-1,-1	25,0	.,.	.	.	.	.	.	.	.	.
```

## 2. Loftee

GitHub: https://github.com/konradjk/loftee

Biowulf documentation: https://hpc.nih.gov/apps/VEP.html

```
ml VEP/101
export PERL5LIB=$PERL5LIB:${VEPCACHEDIR}/Plugins/loftee_GRCh38

vep \
  --offline --cache --dir_cache ${VEPCACHEDIR} \
  --input_file subset_NIH_ONLY_forANNOVAR_chr22.vcf \
  --species human --assembly GRCh38 --fasta ${VEPCACHEDIR}/GRCh38.fa \
  --output_file subset_NIH_ONLY_forANNOVAR_chr22.loftee --stats_file subset_NIH_ONLY_forANNOVAR_chr22_summary.txt \
  --plugin LoF,loftee_path:${VEPCACHEDIR}/Plugins/loftee_GRCh38,\
human_ancestor_fa:${VEPCACHEDIR}/Plugins/loftee_GRCh38/human_ancestor.fa.gz,\
conservation_file:${VEPCACHEDIR}/Plugins/loftee_GRCh38/loftee.sql,\
gerp_bigwig:${VEPCACHEDIR}/Plugins/loftee_GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw

# See the summary html -> /data/CARD/PD/AMP_NIH/no_relateds/test_snpEff/subset_NIH_ONLY_forANNOVAR_chr22_summary.html

# Filter for only 
grep 'LoF=HC' subset_NIH_ONLY_forANNOVAR_chr22.loftee| uniq > refined.LoF.all
grep 'LoF=HC' subset_NIH_ONLY_forANNOVAR_chr22.loftee|cut -f1 | uniq > refined.LoF.rsid

head -2 refined.LoF.all | cut -f1-13
# chr22:15528412:AG:A	22:15528413	-	ENSG00000130538	ENST00000252835	Transcript	frameshift_variant	255	255	85	E/X	gaG/ga	-
# chr22:15528412:AG:A	22:15528413	-	ENSG00000130538	ENST00000643195	Transcript	frameshift_variant	222	222	74	E/X	gaG/ga	-

head -2 refined.LoF.all | cut -f14
# IMPACT=HIGH;STRAND=1;LoF=HC;LoF_flags=SINGLE_EXON,PHYLOCSF_WEAK;LoF_info=PERCENTILE:0.259938837920489,GERP_DIST:1150.95936574936,BP_DIST:726,DIST_FROM_LAST_EXON:-254,50_BP_RULE:FAIL,ANN_ORF:-4470.42,MAX_ORF:-4470.42
# IMPACT=HIGH;STRAND=1;LoF=HC;LoF_flags=SINGLE_EXON,PHYLOCSF_WEAK;LoF_info=PERCENTILE:0.234177215189873,GERP_DIST:1150.95936574936,BP_DIST:726,DIST_FROM_LAST_EXON:-221,50_BP_RULE:FAIL,ANN_ORF:-4255.2,MAX_ORF:-4255.2

head refined.LoF.rsid
# chr22:15528412:AG:A
# chr22:15528422:T:G
# chr22:15528466:A:AG
# chr22:15528966:CTT:C
# chr22:16591038:G:A
# chr22:16591067:T:TC
# chr22:16591178:A:C
# chr22:16591457:C:T
# chr22:16591476:A:AT
# chr22:16591592:C:T
```
