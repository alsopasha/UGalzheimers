# ==============================================================================
# Early Alzheimer's Data Analysis
# Author: Yusuf Efe Kivilcim
# Institution: University of Cambridge
#
# Description: Merges clinical, genetic, and 16S microbiome data to evaluate 
# the interaction between salivary 17β-estradiol, P. gingivalis, and APOE ε4.
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(tableone)
library(ggplot2)
library(ggpubr)

# 2. Load and Merge Validation Data
# ------------------------------------------------------------------------------
clinical_data <- read.csv("../data/processed/clinical_genetic_data.csv")
microbiome_data <- read.csv("../data/processed/microbiome_asv_table.csv")

# Merge on Subject_ID
merged_data <- merge(clinical_data, microbiome_data, by = "Subject_ID")

# 3. Data Pre-Processing
# ------------------------------------------------------------------------------
# Convert Group to a factor with specific levels
merged_data$Group <- factor(merged_data$Group, levels = c("HC", "Early_AD"))

# Create a binary APOE e4 Carrier flag (1 if e3/e4 or e4/e4, else 0)
merged_data$APOE_e4_Carrier <- ifelse(grepl("e4", merged_data$APOE_Genotype), 1, 0)
merged_data$APOE_e4_Carrier <- factor(merged_data$APOE_e4_Carrier, levels = c(0, 1), labels = c("No", "Yes"))

# Factorize TREM2 Status
merged_data$TREM2_Status <- factor(merged_data$TREM2_Status, levels = c("Non_Carrier", "Carrier"))

# 4. Table 1: Demographics and Baseline Characteristics
# ------------------------------------------------------------------------------
vars <- c("Age", "BMI", "MoCA", "MMSE", "Estradiol_pg_mL", "PPD_mm", "CAL_mm", 
          "BOP_percent", "P_gingivalis", "APOE_e4_Carrier")
nonnormal_vars <- c("Estradiol_pg_mL", "P_gingivalis", "BOP_percent") 

table1 <- CreateTableOne(vars = vars, strata = "Group", data = merged_data)
print(table1, nonnormal = nonnormal_vars, exact = "APOE_e4_Carrier")

# 5. Non-Parametric Comparison: Salivary Estradiol
# ------------------------------------------------------------------------------
wilcox_result <- wilcox.test(Estradiol_pg_mL ~ Group, data = merged_data)
cat("\nMann-Whitney U Test for Estradiol differences between HC and Early AD:\n")
print(wilcox_result)

# 6. Correlation Analysis: Estradiol vs. Pathogen Abundance
# ------------------------------------------------------------------------------
cor_test <- cor.test(merged_data$Estradiol_pg_mL, merged_data$P_gingivalis, method = "spearman")
cat("\nSpearman Correlation (Estradiol vs P. gingivalis):\n")
print(cor_test)

# Plot the correlation
cor_plot <- ggplot(merged_data, aes(x = Estradiol_pg_mL, y = P_gingivalis, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("HC" = "#00BFC4", "Early_AD" = "#F8766D")) +
  labs(title = "Inverse Correlation: Salivary 17β-Estradiol & P. gingivalis",
       x = "Salivary 17β-Estradiol (pg/mL)",
       y = "Relative Abundance of P. gingivalis") +
  stat_cor(method = "spearman", label.x = 1.5, label.y = 0.25)

# Save plot (uncomment in local environment)
# ggsave("correlation_plot.png", plot = cor_plot, width = 8, height = 6)

# 7. The "Triple Hit" Multivariable Logistic Regression
# ------------------------------------------------------------------------------
# Modeling Early AD risk based on the interaction of hormones, microbes, and genetics
# Adjusted for age and baseline clinical attachment loss (CAL_mm)

model <- glm(Group ~ Estradiol_pg_mL * P_gingivalis * APOE_e4_Carrier + Age + CAL_mm, 
             data = merged_data, family = binomial(link = "logit"))

cat("\nMultivariable Logistic Regression Results:\n")
summary(model)
