# PD Genetic Burden testing...

```
November 2020
LNG: Cornelis, Andy and Mike
QMUL: Alastair and Ben

Data sources:
PD WGS => Combination of AMP-PD data and internal LNG WGS data (case-control)
~3K cases and ~4K controls
UK Biobank => Combination of PD cases, PD proxies (parent with PD) and controls
~600 cases, ~6K proxies and ~50K controls
```

##### General workflow:
- Individual level data filtering
- Annotation of variants
- Subset variant classes
- Subset genetic data
- Run Burden
- Filter burden for number of variants (min 3 in each gene)
- Meta-analyze burden tests
- Create cumulative case-control frequency files
- Create final results file for interpretation

##### Files in repo
```
UK Biobank Burden.md => describes workflow of UK Biobank data
PD WGS Burden.md => describes workflow of PD WGS data
Meta analyses of Burden.md => describes the meta-analyses of burden results
PD WGS Age at onset Burden.md => describes workflow of PD WGS age at onset burden testing
```

##### Individual level filtering:
```
- European ancestry
- non-related (PIHAT <0.125)
Additional control filtering for UK Biobank => controls have age of recruit >60 and no AD or PD parent and no PD-ism, no dementia diagnosis
```

##### Variant selection after annotation with ANNOVAR and frequency levels:
```
ALL_MISSENSE.txt => 5, 1, 0.5 and 0.1%
ALL_LOF.txt => (splicing, frameshift, stopgain, stoploss) => 5, 1, 0.5 and 0.1%
ALL_CADD_20.txt => 5, 1, 0.5 and 0.1%
ALL_CADD_10.txt => 5, 1, 0.5 and 0.1%
ALL_MISSENSE_and_LOF.txt => All missense + splicing, frameshift, stopgain, stoploss => 5, 1, 0.5 and 0.1%
```



