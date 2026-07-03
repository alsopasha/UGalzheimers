library(vegan)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

asv_counts <- read.csv("../data/raw/asv_counts.csv", row.names = 1)
metadata <- read.csv("../data/processed/clinical_metadata.csv", row.names = 1)
taxonomy <- read.csv("../data/raw/taxonomy.csv", row.names = 1)

metadata$Cohort <- factor(metadata$Cohort, levels = c("HC", "Early AD"))

rel_abund <- apply(asv_counts, 2, function(x) x / sum(x))
bray_curtis_dist <- vegdist(t(rel_abund), method = "bray")

pcoa_res <- cmdscale(bray_curtis_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(PC1 = pcoa_res$points[, 1], PC2 = pcoa_res$points[, 2], Cohort = metadata$Cohort)

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Cohort)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal()

permanova_res <- adonis2(bray_curtis_dist ~ Cohort, data = metadata, permutations = 9999)

dds <- DESeqDataSetFromMatrix(countData = asv_counts,
                              colData = metadata,
                              design = ~ CAL_mm + Age + Cohort)

dds <- DESeq(dds)
res <- results(dds, contrast = c("Cohort", "Early AD", "HC"))
res_sig <- res[which(res$padj < 0.05), ]

write.csv(as.data.frame(res), "../data/processed/deseq2_results.csv")
write.csv(as.data.frame(permanova_res), "../data/processed/permanova_results.csv")
