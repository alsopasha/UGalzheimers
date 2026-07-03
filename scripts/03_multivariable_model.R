library(dplyr)
library(stats)

clinical_data <- read.csv("../data/processed/clinical_metadata.csv")

spearman_res <- cor.test(clinical_data$Estradiol_pg_mL, 
                         clinical_data$P_gingivalis_Abundance, 
                         method = "spearman")

cat("Spearman Rank-Order Correlation: rho =", spearman_res$estimate, "p-value =", spearman_res$p.value, "\n")

clinical_data$EarlyAD_Status <- ifelse(clinical_data$Cohort == "Early AD", 1, 0)
clinical_data$APOE_e4_Carrier <- as.numeric(clinical_data$APOE_e4_Carrier)

model <- glm(EarlyAD_Status ~ Age + CAL_mm + Estradiol_pg_mL + 
               P_gingivalis_Abundance + APOE_e4_Carrier + 
               P_gingivalis_Abundance:APOE_e4_Carrier,
             data = clinical_data, 
             family = binomial(link = "logit"))

model_summary <- summary(model)
odds_ratios <- exp(coef(model))
conf_intervals <- exp(confint(model))

results_df <- data.frame(
  OddsRatio = odds_ratios,
  CI_Lower = conf_intervals[, 1],
  CI_Upper = conf_intervals[, 2],
  P_Value = model_summary$coefficients[, 4]
)

write.csv(results_df, "../data/processed/logistic_regression_results.csv")
