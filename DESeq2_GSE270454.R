# ============================================================
# Project: AI-based Early Alzheimer's Detection (Blood RNA-seq)
# Dataset: GSE270454
# Step: Differential Expression Analysis using DESeq2
# Author: (Your Name)
# ============================================================

# -------------------------------
# 1. Load Required Libraries
# -------------------------------
library(DESeq2)

# -------------------------------
# 2. Load Count Matrix
# -------------------------------
counts <- read.csv("GSE270454_RNAseq-combined-counts-matrix.csv.gz",
                   row.names = 1,
                   check.names = FALSE)

cat("Counts matrix dimensions:\n")
print(dim(counts))

# -------------------------------
# 3. Create Sample Metadata (coldata)
# -------------------------------
samples <- colnames(counts)

condition <- ifelse(grepl("^AD", samples), "AD",
                    ifelse(grepl("^MCI", samples), "MCI",
                           ifelse(grepl("^ASM", samples), "ASM",
                                  ifelse(grepl("^ASO", samples), "ASO", NA))))

coldata <- data.frame(row.names = samples,
                      condition = factor(condition))

cat("Condition distribution:\n")
print(table(coldata$condition))

# -------------------------------
# 4. Filter Low Count Genes
# -------------------------------
counts <- counts[rowSums(counts) >= 10, ]
cat("Filtered counts matrix dimensions:\n")
print(dim(counts))

# -------------------------------
# 5. Create DESeq2 Dataset Object
# -------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = ~ condition)

# -------------------------------
# 6. Run DESeq2
# -------------------------------
dds <- DESeq(dds)

# -------------------------------
# 7. Differential Expression Results
# -------------------------------

# AD vs MCI
res_AD_MCI <- results(dds, contrast = c("condition", "AD", "MCI"))
summary(res_AD_MCI)

# AD vs ASO
res_AD_ASO <- results(dds, contrast = c("condition", "AD", "ASO"))
summary(res_AD_ASO)

# AD vs ASM
res_AD_ASM <- results(dds, contrast = c("condition", "AD", "ASM"))
summary(res_AD_ASM)

# -------------------------------
# 8. Extract Significant DEGs (padj < 0.05)
# -------------------------------

sig_AD_MCI <- res_AD_MCI[which(res_AD_MCI$padj < 0.05), ]
sig_AD_MCI <- sig_AD_MCI[order(sig_AD_MCI$padj), ]

sig_AD_ASO <- res_AD_ASO[which(res_AD_ASO$padj < 0.05), ]
sig_AD_ASO <- sig_AD_ASO[order(sig_AD_ASO$padj), ]

sig_AD_ASM <- res_AD_ASM[which(res_AD_ASM$padj < 0.05), ]
sig_AD_ASM <- sig_AD_ASM[order(sig_AD_ASM$padj), ]

cat("\nNumber of significant genes (padj < 0.05):\n")
cat("AD vs MCI:", sum(res_AD_MCI$padj < 0.05, na.rm = TRUE), "\n")
cat("AD vs ASO:", sum(res_AD_ASO$padj < 0.05, na.rm = TRUE), "\n")
cat("AD vs ASM:", sum(res_AD_ASM$padj < 0.05, na.rm = TRUE), "\n")

# -------------------------------
# 9. Save Results as CSV
# -------------------------------

write.csv(as.data.frame(res_AD_MCI), "DESeq2_AD_vs_MCI_all.csv")
write.csv(as.data.frame(res_AD_ASO), "DESeq2_AD_vs_ASO_all.csv")
write.csv(as.data.frame(res_AD_ASM), "DESeq2_AD_vs_ASM_all.csv")

write.csv(as.data.frame(sig_AD_MCI), "DESeq2_AD_vs_MCI_sig.csv")
write.csv(as.data.frame(sig_AD_ASO), "DESeq2_AD_vs_ASO_sig.csv")
write.csv(as.data.frame(sig_AD_ASM), "DESeq2_AD_vs_ASM_sig.csv")

cat("\nCSV output files saved successfully.\n")

# -------------------------------
# 10. Display Top 20 DEGs
# -------------------------------

top20_AD_MCI <- head(sig_AD_MCI, 20)
top20_AD_ASO <- head(sig_AD_ASO, 20)
top20_AD_ASM <- head(sig_AD_ASM, 20)

cat("\nTop 20 genes AD vs MCI:\n")
print(top20_AD_MCI)

cat("\nTop 20 genes AD vs ASO:\n")
print(top20_AD_ASO)

cat("\nTop 20 genes AD vs ASM:\n")
print(top20_AD_ASM)

# -------------------------------
# 11. Print Top Gene Names
# -------------------------------
cat("\nTop 10 Gene IDs AD vs MCI:\n")
print(rownames(head(sig_AD_MCI, 10)))

cat("\nTop 10 Gene IDs AD vs ASO:\n")
print(rownames(head(sig_AD_ASO, 10)))

cat("\nTop 10 Gene IDs AD vs ASM:\n")
print(rownames(head(sig_AD_ASM, 10)))

# -------------------------------
# END
# -------------------------------
