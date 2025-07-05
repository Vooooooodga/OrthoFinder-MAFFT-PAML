#!/usr/bin/env Rscript

# Objective: Compare RELAX analysis results between two gene sets (test vs. background).
#
# Key Steps:
# 1. Load and clean data:
#    - Remove any genes from the background set that are present in the test set.
#    - Complete missing data in the test set using data from the background set.
# 2. FDR Correction: Apply Benjamini-Hochberg FDR correction to both datasets independently.
# 3. Statistical Summary: Count significant genes with k > 1 (intensified selection) and k < 1 (relaxed selection).
# 4. Hypothesis Testing:
#    - Use a two-sided Fisher's Exact Test to determine if there is a significant difference
#      in the proportion of genes under selection (relaxed, intensified, or any) between the two sets.

# --- Setup ---

if (!requireNamespace("dplyr", quietly = TRUE)) {
  cat("Installing package 'dplyr'...\n")
  install.packages("dplyr")
}
library(dplyr)

cat("R script execution started...\n\n")

# --- Parameters ---

file_background <- "summary_all.csv"
file_test1 <- "summary_test.csv"
file_test2 <- "summary_test2.csv"

alpha <- 0.05

# --- Step 1: Load and Clean Data ---

cat(paste("1. Loading data from", file_background, ",", file_test1, "and", file_test2, "...\n"))

tryCatch({
  df_bg_raw <- read.csv(file_background, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  df_test1_raw <- read.csv(file_test1, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  df_test2_raw <- read.csv(file_test2, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}, error = function(e) {
  cat(paste("Error: Could not read CSV files. Please ensure '", file_background, "', '", file_test1, "', and '", file_test2, "' exist in the current directory.\n", sep=""))
  stop(e)
})

# Merge test files and remove duplicates
df_test_raw <- bind_rows(df_test1_raw, df_test2_raw) %>%
  distinct(Gene, .keep_all = TRUE)

cat("   - Merged and deduplicated test files, resulting in", nrow(df_test_raw), "unique test genes.\n")

# Get the complete list of test genes
test_gene_list <- unique(df_test_raw$Gene)

# Remove all test genes from the background set to ensure its purity
df_bg_cleaned <- df_bg_raw %>%
  filter(!Gene %in% test_gene_list)

cat("   - Removed", nrow(df_bg_raw) - nrow(df_bg_cleaned), "test genes from the background set.\n")

# Complete missing data in the test set
genes_to_complete <- df_test_raw %>%
  filter(is.na(p.value)) %>%
  pull(Gene)

completion_data <- df_bg_raw %>%
  filter(Gene %in% genes_to_complete)

df_test_completed <- df_test_raw %>%
  filter(!is.na(p.value)) %>%
  bind_rows(completion_data) %>%
  distinct(Gene, .keep_all = TRUE)

cat("   - Completed missing data for", nrow(completion_data), "test genes.\n\n")

# --- Step 2: Pre-processing and FDR Correction ---

cat("2. Performing FDR correction...\n")

process_and_correct <- function(df, dataset_name) {
  df <- df %>%
    mutate(
      p.value = as.numeric(as.character(p.value)),
      k = as.numeric(as.character(k))
    ) %>%
    filter(!is.na(p.value) & !is.na(k))

  df$p.adj <- p.adjust(df$p.value, method = "fdr")

  cat(paste("   -", dataset_name, ": Processed", nrow(df), "valid genes.\n"))
  return(df)
}

df_bg_final <- process_and_correct(df_bg_cleaned, "Background set")
df_test_final <- process_and_correct(df_test_completed, "Test set")

cat("\n")

# --- Step 3: Statistical Summary ---

cat("3. Summarizing significant genes (p.adj <", alpha, ")...\n")

summarize_results <- function(df) {
  total_genes <- nrow(df)
  sig_k_gt_1 <- df %>% filter(p.adj < alpha & k > 1) %>% nrow()
  sig_k_lt_1 <- df %>% filter(p.adj < alpha & k < 1) %>% nrow()
  sig_total <- sig_k_gt_1 + sig_k_lt_1
  return(list(
    total = total_genes,
    sig_k_gt_1 = sig_k_gt_1,
    sig_k_lt_1 = sig_k_lt_1,
    sig_total = sig_total
  ))
}

summary_bg <- summarize_results(df_bg_final)
summary_test <- summarize_results(df_test_final)

cat("--------------------------------------------------\n")
cat("Statistical Summary:\n")
cat("--------------------------------------------------\n")
cat(sprintf("Background Set:\n"))
cat(sprintf("  - Total valid genes: %d\n", summary_bg$total))
cat(sprintf("  - Significant & k > 1 (Intensified): %d\n", summary_bg$sig_k_gt_1))
cat(sprintf("  - Significant & k < 1 (Relaxed): %d\n", summary_bg$sig_k_lt_1))
cat(sprintf("  - Total significant genes: %d\n", summary_bg$sig_total))
cat("--------------------------------------------------\n")
cat(sprintf("Test Set:\n"))
cat(sprintf("  - Total valid genes: %d\n", summary_test$total))
cat(sprintf("  - Significant & k > 1 (Intensified): %d\n", summary_test$sig_k_gt_1))
cat(sprintf("  - Significant & k < 1 (Relaxed): %d\n", summary_test$sig_k_lt_1))
cat(sprintf("  - Total significant genes: %d\n", summary_test$sig_total))
cat("--------------------------------------------------\n\n")

# --- Step 4: Hypothesis Testing ---

run_fisher_test <- function(group1_sig, group1_total, group2_sig, group2_total, test_name, hypothesis) {
  cat("--------------------------------------------------\n")
  cat(paste("Test:", test_name, "\n"))
  cat(paste("Hypothesis:", hypothesis, "\n"))
  cat("--------------------------------------------------\n")

  contingency_table <- matrix(c(
    group1_sig, group1_total - group1_sig,
    group2_sig, group2_total - group2_sig
  ), nrow = 2, byrow = TRUE)

  colnames(contingency_table) <- c("Significant", "Not Significant")
  rownames(contingency_table) <- c("Test Set", "Background Set")

  # Perform a two-sided Fisher's Exact Test
  fisher_result <- fisher.test(contingency_table, alternative = "two.sided")

  print(contingency_table)
  cat("\n")
  cat(sprintf("P-value: %g\n", fisher_result$p.value))
  cat(sprintf("Odds Ratio Estimate: %g\n", fisher_result$estimate))
  cat("\n")

  if (fisher_result$p.value < alpha) {
    cat(paste("Conclusion: The result is significant at the", alpha, "level.\n"))
    if (fisher_result$estimate > 1) {
      cat("The proportion of significant genes in the Test Set is significantly HIGHER than in the Background Set.\n")
    } else {
      cat("The proportion of significant genes in the Test Set is significantly LOWER than in the Background Set.\n")
    }
  } else {
    cat(paste("Conclusion: The result is not significant at the", alpha, "level.\n"))
    cat("There is no sufficient evidence to claim a difference in proportions between the Test and Background sets.\n")
  }
  cat("--------------------------------------------------\n\n")
}

# Test 1: Relaxed selection (k < 1)
run_fisher_test(
  summary_test$sig_k_lt_1, summary_test$total,
  summary_bg$sig_k_lt_1, summary_bg$total,
  "Enrichment Analysis for Relaxed Selection (k < 1)",
  "Is the proportion of genes under relaxed selection different in the Test Set compared to the Background Set?"
)

# Test 2: Intensified selection (k > 1)
run_fisher_test(
  summary_test$sig_k_gt_1, summary_test$total,
  summary_bg$sig_k_gt_1, summary_bg$total,
  "Enrichment Analysis for Intensified Selection (k > 1)",
  "Is the proportion of genes under intensified selection different in the Test Set compared to the Background Set?"
)

# Test 3: Any significant selection (k < 1 or k > 1)
run_fisher_test(
  summary_test$sig_total, summary_test$total,
  summary_bg$sig_total, summary_bg$total,
  "Enrichment Analysis for Any Selection",
  "Is the proportion of genes under any selection different in the Test Set compared to the Background Set?"
)

cat("R script execution finished.\n") 