# Correlation Plots for Annovar vs. SnpEff / LOFTEE

#### The final plots are here: 
`/data/CARD/PD/AMP_NIH/no_relateds/annovar_snpeff_correlation`

## Structure of README:

### [1. Make loop over files](#1-make-loop-over-files-1)
### [2. Make correlation plots](#2-make-correlation-plots-1)
### [3. Check results](#3-check-results-1)

---

### 1. Make loop over files
```
cd /data/CARD/PD/AMP_NIH/no_relateds
mkdir annovar_snpeff_correlation
cd annovar_snpeff_correlation

echo "ALL_CADD_10 ALL_MODERATE_HIGH_IMPACT_SNPEFF
ALL_CADD_20 ALL_HIGH_IMPACT_SNPEFF
ALL_MISSENSE_and_LOF ALL_MISSENSE_and_LOF_SNPEFF
ALL_MISSENSE_and_LOF ALL_MISSENSE_and_ALL_LOF_HC_LOFTEE
ALL_MISSENSE_and_LOF ALL_MISSENSE_and_LOF_SNPEFF_and_HC_LOFTEE
ALL_MISSENSE ALL_MISSENSE_SNPEFF
ALL_LOF ALL_LOF_SNPEFF
ALL_LOF ALL_LOF_HC_LOFTEE
ALL_LOF ALL_LOF_and_HC_LOFTEE
ALL_CADD_20_and_LOF ALL_HIGH_IMPACT_and_LOF_SNPEFF
ALL_CADD_20_and_LOF ALL_HIGH_IMPACT_and_ALL_LOF_HC_LOFTEE
ALL_CADD_20_and_LOF ALL_HIGH_IMPACT_and_LOF_SNPEFF_and_HC_LOFTEE" > variant_categories.txt

echo "AMP_NIH_UKB_ALL_PD_PHENO_META ALL_PD_PHENOTYPES_CONTROL
AMP_NIH_UKB_CASE_CONTROL_META PD_CASE_CONTROL
AMP_NIH_UKB_CASE_PROXIES_META PD_CASE_PROXIES_CONTROL" > pheno_categories.txt
```

### 2. Make correlation plots
```
cat variant_categories.txt | while read line
do
    ANNOVAR=$(echo $line | cut -d" " -f1)
    SNPEFF=$(echo $line | cut -d" " -f2)
    cat pheno_categories.txt | while read PHENO
    do
        PHENO_DIR=$(echo $PHENO | cut -d" " -f1)
        PHENO_FILE=$(echo $PHENO | cut -d" " -f2)
        for NUMVAR in {1..4};
        do
            for MAF in {0.005,0.001,0.05,0.01};
            do
                for TEST in {"SkatO","CMC"};
                do
                    Rscript --vanilla plot_correlation.R \
                    "/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/${PHENO_DIR}/ANNOVAR/minimum_numvar_${NUMVAR}/meta_AMP_NIH_noRelateds_UKB_EXOM_${PHENO_FILE}_2021_${ANNOVAR}_freqUpper${MAF}.combined_Ps.${TEST}.tab" \
                    "/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/${PHENO_DIR}/SNPEFF_LOFTEE/minimum_numvar_${NUMVAR}/meta_AMP_NIH_noRelateds_UKB_EXOM_${PHENO_FILE}_2021_${SNPEFF}_freqUpper${MAF}.combined_Ps.${TEST}.tab"
                done
            done
        done
    done
done
```
```
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# start like this
# Rscript --vanilla plot_correlation.R $ANNOVAR $SNPEFF
ANNOVAR=args[1]
SNPEFF=args[2]

# Load packages
require(dplyr)
require(data.table)
require(stringr)
require(ggpubr)
require(sjmisc)

# Import annovar and SnpEff data
annovar <- fread(ANNOVAR, header=T)
snpeff <- fread(SNPEFF, header=T)

# Calculate the number of genes in each dataset and the number in both
annovar_genes <- length(unique(annovar$GENE))
snpeff_genes <- length(unique(snpeff$GENE))
Mrg <- merge(annovar, snpeff, by="GENE")
common_genes <- length(unique(Mrg$GENE))
print(paste("There are ", annovar_genes, " genes in the annovar dataset, ", snpeff_genes, " genes in the SnpEff/LOFTEE dataset and ", common_genes, " in both datasets"))

# Take the log of the META P-value for plotting
Mrg$META1_logP <- -log(Mrg$META_PVAL.x, base=10)
Mrg$META2_logP <- -log(Mrg$META_PVAL.y, base=10)

# Rename the meta-analysis groups
meta_group <- Mrg$META_ANALYSIS_GROUP.x %>% unique() %>% 
	gsub(pattern="meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_CONTROL_2021", replacement="Cases vs. Controls") %>%
	gsub(pattern="meta_AMP_NIH_noRelateds_UKB_EXOM_ALL_PD_PHENOTYPES_CONTROL_2021", replacement="All PD Phenotypes vs. Controls") %>%
	gsub(pattern="meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021", replacement="Cases and Proxies vs. Controls")

# Reformat the variant group name for the title
variant_group <- unique(Mrg$VARIANT_GROUP.y) %>% gsub(pattern="_", replacement=" ") %>% str_to_title() 

# Make the title with the conditions used in this burden analysis
plt_title <- paste(meta_group, variant_group, sep= " - ") %>% 
	gsub(pattern="All", replacement="") %>% 
	gsub(pattern="Lof Snpeff", replacement="LOF (SnpEff)") %>% 
	gsub(pattern="Lof Hc Loftee", replacement="LOF (LOFTEE)") %>% 
	gsub(pattern="Moderate High Impact Snpeff", replacement="Predicted Moderate or High Impact") %>% 
	gsub(pattern="High Impact Snpeff", replacement="Predicted High Impact") %>% 
	gsub(pattern="Missense Snpeff", replacement="Missense") %>% 
	gsub(pattern="And", replacement="and")
if (str_contains(plt_title, "and Hc Loftee")) {
    plt_title2 <- bquote(bold(.(gsub(pattern="L.*Loftee", replacement="LOF (SnpEff", plt_title))~intersect(LOFTEE)*")"))
} else {
    plt_title2 <- plt_title
}

# Add a second title with the test (Skat-O or CMC), MAF cutoff and minimum number of variants per gene
MAF <- toString(unique(Mrg$MAF.x))
MIN_NUMVAR <- toString(unique(Mrg$MIN_NUMVAR.x))
TEST <- unique(Mrg$TEST.x)
plt_title3 <- bquote(bold(.(TEST)~"Burden Test with MAF" <= .(MAF)~"and"~N[variants~per~gene] >= .(MIN_NUMVAR)))

# Make the name to use for the pdf output, based on the SnpEff annotations
OUTNAME <- SNPEFF %>% gsub(pattern=".*/", replacement="") %>% gsub(pattern=".tab", replacement=paste(".NumVar_", MIN_NUMVAR, sep=""))

# Make a scatterplot with the log10 meta-analysis P-values (x-axis=ANNOVAR, y-axis=SnpEff/LOFTEE)
pdf(paste(OUTNAME, ".pdf", sep=""), onefile=F)
plot.new()
ggscatter(Mrg, x = "META1_logP", y = "META2_logP", add = "reg.line", 
xlab = expression("log"[10]*" Annovar P-value"), ylab = expression("log"[10]*" SnpEff P-value"), 
add.params = list(color = "blue"), color="gray")

# Add titles to plot
mtext(side=3, line = 2.5, "Meta-analysis of AMP-PD x NIH WGS and UK Biobank WES", cex=1, font=2)
mtext(side=3, line = 1, plt_title2, cex=0.8, font=2)
mtext(side=3, line = -0.5, plt_title3, cex=0.8, font=2)

# Add correlations of the meta-analyses
cor1 <- round(cor(Mrg$META1_logP,Mrg$META2_logP, use ="complete.obs"), 2)
mtext(side=3, line = -2.2, bquote(R[Meta] == .(cor1)), cex=1, font=3, col="red")
cor2 <- round(cor(-log(Mrg$AMP_NIH_PVAL.x, base=10), -log(Mrg$AMP_NIH_PVAL.y, base=10), use ="complete.obs"), 2)

# Also add AMPxNIH and UKB specific correlations
if (meta_group == "Cases vs. Controls") {
    cor3 <- round(cor(-log(Mrg$UKB_CASE_PVAL.x, base=10), -log(Mrg$UKB_CASE_PVAL.y, base=10), use ="complete.obs"), 2)
    mtext(side=3, line = -3.5, bquote(R[AMP~x~NIH] == .(cor2)~~~R[UKB] == .(cor3)), cex=0.8, font=1)
} else if (meta_group == "Cases and Proxies vs. Controls") {
    cor3 <- round(cor(-log(Mrg$UKB_CASE_PVAL.x, base=10), -log(Mrg$UKB_CASE_PVAL.y, base=10), use ="complete.obs"), 2)
    cor4 <- round(cor(-log(Mrg$UKB_SIBLING_PROXY_PVAL.x, base=10), -log(Mrg$UKB_SIBLING_PROXY_PVAL.y, base=10), use ="complete.obs"), 2)
    cor5 <- round(cor(-log(Mrg$UKB_PARENT_PROXY_PVAL.x, base=10), -log(Mrg$UKB_PARENT_PROXY_PVAL.y, base=10), use ="complete.obs"), 2)
    mtext(side=3, line = -3.5, bquote(R[AMP~x~NIH] == .(cor2)~~~R[UKB~Case] == .(cor3)), cex=0.8, font=1)
    mtext(side=3, line = -4.5, bquote(R[UKB~Sibling] == .(cor4)~~~R[UKB~Parent] == .(cor5)), cex=0.8, font=1)
} else {
    cor3 <- round(cor(-log(Mrg$UKB_ALL_CASES_PVAL.x, base=10), -log(Mrg$UKB_ALL_CASES_PVAL.y, base=10), use ="complete.obs"), 2)
    mtext(side=3, line = -3.5, bquote(R[AMP~x~NIH] == .(cor2)~~~R[UKB~All~Phenotypes] == .(cor3)), cex=0.8, font=1)
}
dev.off()
```

### 3. Check results
```
ls meta*pdf | wc -l
# 1144

## We expect: 
12 variant categories * 3 phenotypes * 4 MAF cutoffs * 2 tests * 4 variant per gene cutoffs = 1,152 

## Note that eight of the correlation plots failed due to empty files here: 
/data/CARD/PD/AMP_NIH/no_relateds/meta_risk_analysis/AMP_NIH_UKB_CASE_PROXIES_META/SNPEFF_LOFTEE/minimum_numvar_{1..4}/meta_AMP_NIH_noRelateds_UKB_EXOM_PD_CASE_PROXIES_CONTROL_2021_ALL_MODERATE_HIGH_IMPACT_SNPEFF_freqUpper0.001.combined_Ps.{"SkatO","CMC"}.tab
# These files have a very stringent MAF cutoff so this is fine

scp lakejs@biowulf.nih.gov://data/CARD/PD/AMP_NIH/no_relateds/annovar_snpeff_correlation/*pdf /Users/lakejs/Desktop/burden_correlation
```
