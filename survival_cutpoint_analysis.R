#!/usr/bin/env Rscript
# ==============================================================================
# Biomarker Survival Cutpoint Analysis (Optimal Cutpoint, Tertiles, Quartiles)
# ==============================================================================
# Evaluates survival association (Log-rank P-value & Cox PH Hazard Ratio) for any
# continuous biomarker using maximally selected rank statistics (surv_cutpoint),
# Tertile (T1 vs T3), and Quartile (Q1 vs Q4) stratification.

parse_cli <- function(args) {
  opts <- list(
    input = "data/survival_cohort.csv",
    time_col = "OS_Months",
    event_col = "OS_Status",
    var_col = "PRKDC_mRNA",
    output = "results/survival_cutpoint_summary.csv"
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("-h", "--help")) {
      cat(paste(
        "Usage: Rscript toolkit/survival_cutpoint_analysis.R",
        "--input data/cohort.csv --var-col GENE_EXPR",
        "[--time-col OS_Months] [--event-col OS_Status]",
        "[--output results/survival_cutpoint_summary.csv]\n"
      ))
      quit(status = 0)
    } else if (args[i] == "--input") { opts$input <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--time-col") { opts$time_col <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--event-col") { opts$event_col <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--var-col") { opts$var_col <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--output") { opts$output <- args[i + 1]; i <- i + 2
    } else { i <- i + 1 }
  }
  opts
}

opts <- parse_cli(commandArgs(trailingOnly = TRUE))

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
})

dir.create(dirname(opts$output), recursive = TRUE, showWarnings = FALSE)
d <- read.csv(opts$input, stringsAsFactors = FALSE)
d <- d[!is.na(d[[opts$time_col]]) & !is.na(d[[opts$event_col]]) & !is.na(d[[opts$var_col]]), ]
d$.time <- as.numeric(d[[opts$time_col]])
d$.event <- as.numeric(d[[opts$event_col]])
d$.val <- as.numeric(d[[opts$var_col]])

# 1. Optimal Cutpoint (maxstat)
res_cut <- surv_cutpoint(d, time = ".time", event = ".event", variables = ".val")
cut_val <- as.numeric(summary(res_cut)$cutpoint)
d_cat <- surv_categorize(res_cut)
d_cat$.val <- factor(d_cat$.val, levels = c("low", "high"))
diff_opt <- survdiff(Surv(.time, .event) ~ .val, data = d_cat)
p_opt <- 1 - pchisq(diff_opt$chisq, 1)
cox_opt <- coxph(Surv(.time, .event) ~ .val, data = d_cat)

# 2. Tertiles (T1 vs T3)
q33 <- quantile(d$.val, 1 / 3)
q66 <- quantile(d$.val, 2 / 3)
d$Tertile <- case_when(d$.val <= q33 ~ "Low", d$.val >= q66 ~ "High", TRUE ~ "Mid")
d_t <- subset(d, Tertile %in% c("Low", "High"))
d_t$Tertile <- factor(d_t$Tertile, levels = c("Low", "High"))
diff_t <- survdiff(Surv(.time, .event) ~ Tertile, data = d_t)
p_t <- 1 - pchisq(diff_t$chisq, 1)
cox_t <- coxph(Surv(.time, .event) ~ Tertile, data = d_t)

# 3. Quartiles (Q1 vs Q4)
q25 <- quantile(d$.val, 0.25)
q75 <- quantile(d$.val, 0.75)
d$Quartile <- case_when(d$.val <= q25 ~ "Low", d$.val >= q75 ~ "High", TRUE ~ "Mid")
d_q <- subset(d, Quartile %in% c("Low", "High"))
d_q$Quartile <- factor(d_q$Quartile, levels = c("Low", "High"))
diff_q <- survdiff(Surv(.time, .event) ~ Quartile, data = d_q)
p_q <- 1 - pchisq(diff_q$chisq, 1)
cox_q <- coxph(Surv(.time, .event) ~ Quartile, data = d_q)

summary_df <- data.frame(
  Method = c("Optimal_Cutpoint", "Tertile_High_vs_Low", "Quartile_High_vs_Low"),
  Cutpoint_Info = c(sprintf("%.4f", cut_val), sprintf("<=%.4f vs >=%.4f", q33, q66), sprintf("<=%.4f vs >=%.4f", q25, q75)),
  N_Evaluated = c(nrow(d_cat), nrow(d_t), nrow(d_q)),
  LogRank_P = c(p_opt, p_t, p_q),
  Cox_HR = c(exp(coef(cox_opt)), exp(coef(cox_t)), exp(coef(cox_q))),
  HR_CI95_Low = c(exp(confint(cox_opt))[1], exp(confint(cox_t))[1], exp(confint(cox_q))[1]),
  HR_CI95_High = c(exp(confint(cox_opt))[2], exp(confint(cox_t))[2], exp(confint(cox_q))[2])
)

print(summary_df)
write.csv(summary_df, opts$output, row.names = FALSE)
cat(sprintf("Saved survival cutpoint summary -> %s\n", opts$output))
