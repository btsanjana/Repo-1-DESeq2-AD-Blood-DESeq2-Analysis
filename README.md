# DESeq2 Differential Expression Analysis (GSE270454)

This repository contains DESeq2-based RNA-seq differential gene expression analysis for early Alzheimer’s disease detection using blood transcriptome dataset GSE270454.

## Dataset
- GEO: GSE270454
- Groups: AD, MCI, ASO, ASM

## Method
DESeq2 was used with design:
~ condition

Contrasts performed:
- AD vs MCI
- AD vs ASO
- AD vs ASM

## Output Files
Significant DEGs (padj < 0.05):
- DESeq2_AD_vs_MCI_sig.csv
- DESeq2_AD_vs_ASO_sig.csv
- DESeq2_AD_vs_ASM_sig.csv

All genes:
- DESeq2_AD_vs_MCI_all.csv
- DESeq2_AD_vs_ASO_all.csv
- DESeq2_AD_vs_ASM_all.csv

## Tools
R packages:
- DESeq2
