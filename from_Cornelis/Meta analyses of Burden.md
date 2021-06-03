### Meta analyses of Burden

```
November 2020
Cornelis and Mike

Input data:
PD WGS case-control
UK Biobank case-control
UK Biobank proxy-control

Number of input files:
5 variant categories x 4 frequency levels => 20 per dataset x 3 datasets => 60 files total

variant categories:
ALL_MISSENSE
ALL_LOF
ALL_CADD_20
ALL_CADD_10
ALL_MISSENSE_and_LOF
ALL_CADD_20_and_LOF

frequency levels:
0.05
0.01
0.005
0.001

```

### Start prep for Meta analyses

```
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN/
mkdir META_ANAL
cd META_ANAL
mkdir PD

# copy files to working dir...
## UKB
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN/RESULTS/
scp PD_CASE_CONTROL* /data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN/META_ANAL/PD/
scp PD_PARENT_CONTROL* /data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN/META_ANAL/PD/
## WGS
cd /data/CARD/PD/WGS/june2019/BURDEN/RESULTS/
scp *.assoc /data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN/META_ANAL/PD/

# back to working dir 
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN/META_ANAL/PD/
# check inventory
ls | wc -l
# 144
3 datasets x 4 freq levels x 6 variant class x 2 tests => 144 

```

```
# Start prep data for Meta analyses
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/BURDEN/META_ANAL/PD/

ALL_MISSENSE
ALL_LOF
ALL_CADD_20
ALL_CADD_10
ALL_MISSENSE_and_LOF
ALL_CADD_20_and_LOF

# CMC...
cat PD_PARENT_CONTROL_ALL_MISSENSE_0.01.CMC.assoc PD_CASE_CONTROL_ALL_MISSENSE_0.01.CMC.assoc ALL_MISSENSE_0.01.CMC.assoc > TEMP.txt
grep -v -e 'N_INFORMATIVE' -e '-nan' TEMP.txt > INPUT_META_CMC_ALL_MISSENSE_0.01.tab

cat PD_PARENT_CONTROL_ALL_MISSENSE_and_LOF_0.01.CMC.assoc PD_CASE_CONTROL_ALL_MISSENSE_and_LOF_0.01.CMC.assoc ALL_MISSENSE_and_LOF_0.01.CMC.assoc > TEMP.txt
grep -v -e 'N_INFORMATIVE' -e '-nan' TEMP.txt > INPUT_META_CMC_ALL_MISSENSE_and_LOF_0.01.tab

cat PD_PARENT_CONTROL_ALL_MISSENSE_and_LOF_0.05.CMC.assoc PD_CASE_CONTROL_ALL_MISSENSE_and_LOF_0.05.CMC.assoc ALL_MISSENSE_and_LOF_0.05.CMC.assoc > TEMP.txt
grep -v -e 'N_INFORMATIVE' -e '-nan' TEMP.txt > INPUT_META_CMC_ALL_MISSENSE_and_LOF_0.05.tab

# only parent and WGS
cat PD_PARENT_CONTROL_ALL_MISSENSE_and_LOF_0.01.CMC.assoc ALL_MISSENSE_and_LOF_0.01.CMC.assoc > TEMP.txt
grep -v -e 'N_INFORMATIVE' -e '-nan' TEMP.txt > INPUT_META_v2_CMC_ALL_MISSENSE_and_LOF_0.01.tab



# SkatO...
cat PD_PARENT_CONTROL_ALL_MISSENSE_0.01.SkatO.assoc PD_CASE_CONTROL_ALL_MISSENSE_0.01.SkatO.assoc ALL_MISSENSE_0.01.SkatO.assoc > TEMP.txt
grep -v -e 'N_INFORMATIVE' -e '-nan' TEMP.txt > INPUT_META_SKATO_ALL_MISSENSE_0.01.tab

cat PD_PARENT_CONTROL_ALL_MISSENSE_and_LOF_0.01.SkatO.assoc PD_CASE_CONTROL_ALL_MISSENSE_and_LOF_0.01.SkatO.assoc ALL_MISSENSE_and_LOF_0.01.SkatO.assoc > TEMP.txt
grep -v -e 'N_INFORMATIVE' -e '-nan' TEMP.txt > INPUT_META_SKATO_ALL_MISSENSE_and_LOF_0.01.tab

cat PD_PARENT_CONTROL_ALL_MISSENSE_and_LOF_0.05.SkatO.assoc PD_CASE_CONTROL_ALL_MISSENSE_and_LOF_0.05.SkatO.assoc ALL_MISSENSE_and_LOF_0.05.SkatO.assoc > TEMP.txt
grep -v -e 'N_INFORMATIVE' -e '-nan' TEMP.txt > INPUT_META_SKATO_ALL_MISSENSE_and_LOF_0.05.tab


# only parent and WGS
cat PD_PARENT_CONTROL_ALL_MISSENSE_and_LOF_0.01.SkatO.assoc ALL_MISSENSE_and_LOF_0.01.SkatO.assoc > TEMP.txt
grep -v -e 'N_INFORMATIVE' -e '-nan' TEMP.txt > INPUT_META_v2_SKATO_ALL_MISSENSE_and_LOF_0.01.tab

```


### Start Meta analyses

python Z_meta_burdens_CMC.py -i INPUT_META_CMC_ALL_MISSENSE_0.01.tab
python Z_meta_burdens_CMC.py -i INPUT_META_CMC_ALL_MISSENSE_and_LOF_0.01.tab
python Z_meta_burdens_CMC.py -i INPUT_META_CMC_ALL_MISSENSE_and_LOF_0.05.tab
python Z_meta_burdens_CMC.py -i INPUT_META_v2_CMC_ALL_MISSENSE_and_LOF_0.01.tab

python Z_meta_burdens_SKAT.py -i INPUT_META_SKATO_ALL_MISSENSE_0.01.tab
python Z_meta_burdens_SKAT.py -i INPUT_META_SKATO_ALL_MISSENSE_and_LOF_0.01.tab
python Z_meta_burdens_SKAT.py -i INPUT_META_SKATO_ALL_MISSENSE_and_LOF_0.05.tab
python Z_meta_burdens_SKAT.py -i INPUT_META_v2_SKATO_ALL_MISSENSE_and_LOF_0.01.tab


```


# Imports and prep.
import h5py
import numpy as np
import pandas as pd
import math
import sys
import joblib
import subprocess
import statsmodels.api as sm
from scipy import stats

# Read in your combined and filtered results file. This is the only line you need to change to run the routine on a different file. Just pick the set of 3 lines below relevant to your test.
# input_file = "CMC_to_combine_Ps.tab"
# header_text = ['Gene','RANGE','N_INFORMATIVE','NumVar','NumPolyVar','NonRefSite','Pvalue']
# raw_df = pd.read_csv(input_file, delim_whitespace=True, header=0, names=header_text)
# Or...
input_file = "SKAT_to_combine_Ps.tab"
header_text = ['Gene','RANGE','N_INFORMATIVE','NumVar','NumPolyVar','Q','rho','Pvalue']
raw_df = pd.read_csv(input_file, delim_whitespace=True, header=0, names=header_text)

# Make a gene list for the loop.
gene_list = raw_df.Gene.unique()

# Loop it using Fisher to combine p-values.
results = []

for i in range(len(gene_list)):
  this_gene = gene_list[i]
  temp_data = raw_df[raw_df['Gene'] == this_gene]
  combined_gene = stats.combine_pvalues(temp_data['Pvalue'], method='fisher', weights=None)
  test_stat = combined_gene[0]
  p_val = combined_gene[1]
  print(this_gene, test_stat, p_val)
  results.append([this_gene, test_stat, p_val])

output = pd.DataFrame(results, columns=['GENE', 'TEST_STAT', 'P_VAL'])

# Now the export.
output_file = input_file + ".combined_Ps.csv"
output.to_csv(output_file, index=False)

```


### Figure making
```
TBD Manhattan plot

```


