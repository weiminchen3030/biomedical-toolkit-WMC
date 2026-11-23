#!/usr/bin/env Rscript
# ==============================================================================
# 3-Tier & 2-Tier Tertile Kaplan-Meier & Cox PH Survival Stratification
# ==============================================================================

parse_cli <- function(args) {
  opts <- list(
    input = "data/survival_cohort.csv",
    time_col = "OS_Months",
    event_col = "OS_Status",
    var_col = "PRKDC_mRNA",
    label = "Biomarker",
    output = "results/survival_tertile_summary.csv"
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("-h", "--help")) {
      cat(paste(
        "Usage: Rscript toolkit/survival_tertile_analysis.R",
        "--input data/cohort.csv --var-col GENE_EXPR",
        "[--time-col OS_Months] [--event-col OS_Status]",
        "[--label Biomarker] [--output results/survival_tertile_summary.csv]\n"
      ))
      quit(status = 0)
    } else if (args[i] == "--input") { opts$input <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--time-col") { opts$time_col <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--event-col") { opts$event_col <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--var-col") { opts$var_col <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--label") { opts$label <- args[i + 1]; i <- i + 2
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

q33 <- quantile(d$.val, 1 / 3)
q66 <- quantile(d$.val, 2 / 3)

low_lbl  <- paste0(opts$label, " Low (<=33%)")
mid_lbl  <- paste0(opts$label, " Mid (33-66%)")
high_lbl <- paste0(opts$label, " High (>=66%)")

d$Tertile3 <- case_when(
  d$.val <= q33 ~ low_lbl,
  d$.val >= q66 ~ high_lbl,
  TRUE ~ mid_lbl
)
d$Tertile3 <- factor(d$Tertile3, levels = c(low_lbl, mid_lbl, high_lbl))

fit3 <- survfit(Surv(.time, .event) ~ Tertile3, data = d)
diff3 <- survdiff(Surv(.time, .event) ~ Tertile3, data = d)
p3 <- 1 - pchisq(diff3$chisq, 2)

d2 <- subset(d, Tertile3 != mid_lbl)
d2$Tertile2 <- factor(d2$Tertile3, levels = c(low_lbl, high_lbl))
diff2 <- survdiff(Surv(.time, .event) ~ Tertile2, data = d2)
p2 <- 1 - pchisq(diff2$chisq, 1)
cox2 <- coxph(Surv(.time, .event) ~ Tertile2, data = d2)

cat(sprintf("=== 3-Tier Tertiles (%s) ===\nOverall Log-rank P = %.4e\n", opts$label, p3))
print(summary(fit3)$table)
cat(sprintf("\n=== 2-Tier Top 33%% vs Bottom 33%% ===\nLog-rank P = %.4e | HR = %.2f (95%% CI: %.2f - %.2f)\n",
            p2, exp(coef(cox2)), exp(confint(cox2))[1], exp(confint(cox2))[2]))

out_df <- data.frame(
  Stratum = rownames(summary(fit3)$table),
  summary(fit3)$table,
  Overall_3Tier_LogRank_P = p3,
  High_vs_Low_LogRank_P = p2,
  High_vs_Low_Cox_HR = exp(coef(cox2)),
  High_vs_Low_CI95_Low = exp(confint(cox2))[1],
  High_vs_Low_CI95_High = exp(confint(cox2))[2],
  check.names = FALSE
)
write.csv(out_df, opts$output, row.names = FALSE)
cat(sprintf("Saved tertile survival summary -> %s\n", opts$output))
