library(dplyr)
library(ggplot2)
library(ggpubr)

clinical_data <- read.csv("../data/processed/clinical_metadata.csv")

shapiro_results <- clinical_data %>%
  summarise(across(c(Age, BMI, MoCA, MMSE, CAL_mm, PPD_mm, BOP_pct, Estradiol_pg_mL),
                   ~ shapiro.test(.x)$p.value))

baseline_stats <- clinical_data %>%
  group_by(Cohort) %>%
  summarise(across(c(Age, BMI, MoCA, MMSE, CAL_mm, PPD_mm, BOP_pct, Estradiol_pg_mL),
                   list(mean = mean, sd = sd, median = median, iqr = IQR)))

perform_tests <- function(var_name) {
  test_res <- wilcox.test(clinical_data[[var_name]] ~ clinical_data$Cohort)
  return(test_res$p.value)
}

vars_to_test <- c("Age", "BMI", "MoCA", "MMSE", "CAL_mm", "PPD_mm", "BOP_pct", "Estradiol_pg_mL")
p_values <- sapply(vars_to_test, perform_tests)
p_adjusted <- p.adjust(p_values, method = "BH")

results_table <- data.frame(Variable = vars_to_test, P_Value = p_values, P_Adj = p_adjusted)

ggplot(clinical_data, aes(x = Cohort, y = Estradiol_pg_mL, fill = Cohort)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(y = expression("Salivary 17"*beta*"-Estradiol (pg/mL)")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Early AD", "HC")))

ggsave("../data/processed/estradiol_boxplot.png", width = 6, height = 5)
