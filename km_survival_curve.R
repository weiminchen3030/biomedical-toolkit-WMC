#!/usr/bin/env Rscript
# ==============================================================================
# Publication-Grade Kaplan-Meier Survival Curve & Cox PH CLI
# ==============================================================================
# Stratifies any continuous biomarker or pre-defined group column by:
#   - optimal   : maximally selected rank statistics (survminer::surv_cutpoint)
#   - median    : 50% median split (Low vs High)
#   - tertile3  : 3-tier tertiles (Low <=33%, Intermediate 33-66%, High >=66%)
#   - tertile2  : extreme tertiles (Low <=33% vs High >=66%)
#   - quartile2 : extreme quartiles (Q1 <=25% vs Q4 >=75%)
#   - group     : existing categorical column in input CSV
#
# Outputs publication-ready 600 DPI PNG + vector PDF with risk table, 95% CI,
# Log-rank P-value, Cox PH Hazard Ratio (95% CI), median survival, and summary CSV.
# ==============================================================================

parse_cli <- function(args) {
  opts <- list(
    input = NULL,
    time_col = "OS_Months",
    event_col = "OS_Status",
    var_col = NULL,
    method = "tertile3",
    label = "Biomarker",
    title = NULL,
    xlab = "Time (Months)",
    ylab = "Overall Survival Probability",
    out_prefix = "results/km_survival",
    dpi = 600
  )
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg %in% c("-h", "--help")) {
      cat(paste(
        "Usage: Rscript km_survival_curve.R --input cohort.csv --var-col GENE_EXPR",
        "[--method optimal|median|tertile3|tertile2|quartile2|group]",
        "[--time-col OS_Months] [--event-col OS_Status]",
        "[--label Biomarker] [--title 'Kaplan-Meier Survival']",
        "[--xlab 'Time (Months)'] [--ylab 'Overall Survival Probability']",
        "[--out-prefix results/km_survival] [--dpi 600]\n"
      ))
      quit(status = 0)
    } else if (arg == "--input") { opts$input <- args[i + 1]; i <- i + 2
    } else if (arg == "--time-col") { opts$time_col <- args[i + 1]; i <- i + 2
    } else if (arg == "--event-col") { opts$event_col <- args[i + 1]; i <- i + 2
    } else if (arg == "--var-col") { opts$var_col <- args[i + 1]; i <- i + 2
    } else if (arg == "--method") { opts$method <- tolower(args[i + 1]); i <- i + 2
    } else if (arg == "--label") { opts$label <- args[i + 1]; i <- i + 2
    } else if (arg == "--title") { opts$title <- args[i + 1]; i <- i + 2
    } else if (arg == "--xlab") { opts$xlab <- args[i + 1]; i <- i + 2
    } else if (arg == "--ylab") { opts$ylab <- args[i + 1]; i <- i + 2
    } else if (arg == "--out-prefix") { opts$out_prefix <- args[i + 1]; i <- i + 2
    } else if (arg == "--dpi") { opts$dpi <- as.numeric(args[i + 1]); i <- i + 2
    } else { i <- i + 1 }
  }
  if (is.null(opts$input) || is.null(opts$var_col)) {
    stop("Both --input and --var-col are required. Run with --help for usage.")
  }
  opts
}

opts <- parse_cli(commandArgs(trailingOnly = TRUE))

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

dir.create(dirname(opts$out_prefix), recursive = TRUE, showWarnings = FALSE)
d <- read.csv(opts$input, stringsAsFactors = FALSE)
d <- d[!is.na(d[[opts$time_col]]) & !is.na(d[[opts$event_col]]) & !is.na(d[[opts$var_col]]), ]
d$.time <- as.numeric(d[[opts$time_col]])
d$.event <- as.numeric(d[[opts$event_col]])

if (opts$method == "group") {
  d$.group <- factor(d[[opts$var_col]])
  palette_cols <- c("#0F4C81", "#E5A000", "#8B2635", "#2E7D32")[seq_along(levels(d$.group))]
} else {
  d$.val <- as.numeric(d[[opts$var_col]])
  if (opts$method == "optimal") {
    res_cut <- surv_cutpoint(d, time = ".time", event = ".event", variables = ".val")
    cut_val <- as.numeric(summary(res_cut)$cutpoint)
    d_cat <- surv_categorize(res_cut)
    d$.group <- factor(ifelse(d_cat$.val == "high", paste0(opts$label, " High"), paste0(opts$label, " Low")),
                       levels = c(paste0(opts$label, " Low"), paste0(opts$label, " High")))
    palette_cols <- c("#0F4C81", "#E5A000")
  } else if (opts$method == "median") {
    med <- median(d$.val, na.rm = TRUE)
    d$.group <- factor(ifelse(d$.val >= med, paste0(opts$label, " High (>=50%)"), paste0(opts$label, " Low (<50%)")),
                       levels = c(paste0(opts$label, " Low (<50%)"), paste0(opts$label, " High (>=50%)")))
    palette_cols <- c("#0F4C81", "#E5A000")
  } else if (opts$method == "tertile3") {
    q33 <- quantile(d$.val, 1 / 3)
    q66 <- quantile(d$.val, 2 / 3)
    lbls <- c(paste0(opts$label, " Low (<=33%)"), paste0(opts$label, " Mid (33-66%)"), paste0(opts$label, " High (>=66%)"))
    d$.group <- factor(case_when(d$.val <= q33 ~ lbls[1], d$.val >= q66 ~ lbls[3], TRUE ~ lbls[2]), levels = lbls)
    palette_cols <- c("#0F4C81", "#7F8C8D", "#E5A000")
  } else if (opts$method == "tertile2") {
    q33 <- quantile(d$.val, 1 / 3)
    q66 <- quantile(d$.val, 2 / 3)
    lbls <- c(paste0(opts$label, " Low (<=33%)"), paste0(opts$label, " High (>=66%)"))
    d <- d[d$.val <= q33 | d$.val >= q66, ]
    d$.group <- factor(ifelse(d$.val >= q66, lbls[2], lbls[1]), levels = lbls)
    palette_cols <- c("#0F4C81", "#E5A000")
  } else if (opts$method == "quartile2") {
    q25 <- quantile(d$.val, 0.25)
    q75 <- quantile(d$.val, 0.75)
    lbls <- c(paste0(opts$label, " Q1 (<=25%)"), paste0(opts$label, " Q4 (>=75%)"))
    d <- d[d$.val <= q25 | d$.val >= q75, ]
    d$.group <- factor(ifelse(d$.val >= q75, lbls[2], lbls[1]), levels = lbls)
    palette_cols <- c("#0F4C81", "#E5A000")
  } else {
    stop(sprintf("Unsupported --method: %s", opts$method))
  }
}

surv_obj <- Surv(time = d$.time, event = d$.event)
fit <- survfit(surv_obj ~ .group, data = d)
diff_res <- survdiff(surv_obj ~ .group, data = d)
p_val <- 1 - pchisq(diff_res$chisq, length(diff_res$n) - 1)
cox_fit <- coxph(surv_obj ~ .group, data = d)
cox_sum <- summary(cox_fit)

hr_last <- exp(coef(cox_fit))[length(coef(cox_fit))]
ci_last <- exp(confint(cox_fit))[length(coef(cox_fit)), ]

annot_text <- sprintf(
  "Log-rank P = %s\nCox HR (High vs Low) = %.2f (95%% CI: %.2f-%.2f)",
  ifelse(p_val < 0.001, sprintf("%.2e", p_val), sprintf("%.4f", p_val)),
  hr_last, ci_last[1], ci_last[2]
)

plot_title <- ifelse(is.null(opts$title), sprintf("%s Survival (%s)", opts$label, opts$method), opts$title)

p_km <- ggsurvplot(
  fit,
  data = d,
  palette = palette_cols,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  tables.theme = theme_cleantable(),
  conf.int = TRUE,
  conf.int.alpha = 0.12,
  censor.shape = 124,
  censor.size = 2.5,
  size = 1.1,
  xlab = opts$xlab,
  ylab = opts$ylab,
  title = plot_title,
  legend.title = "",
  legend.labs = levels(d$.group),
  ggtheme = theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.title = element_text(face = "bold", size = 11),
      legend.position = c(0.75, 0.85),
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA)
    )
)

p_km$plot <- p_km$plot +
  annotate("text", x = max(d$.time, na.rm = TRUE) * 0.04, y = 0.14,
           label = annot_text, hjust = 0, vjust = 0, size = 3.6, fontface = "bold")

png_out <- paste0(opts$out_prefix, ".png")
pdf_out <- paste0(opts$out_prefix, ".pdf")
csv_out <- paste0(opts$out_prefix, "_summary.csv")

png(png_out, width = 7.5, height = 6.2, units = "in", res = opts$dpi)
print(p_km)
dev.off()

pdf(pdf_out, width = 7.5, height = 6.2)
print(p_km)
dev.off()

out_df <- data.frame(
  Stratum = rownames(summary(fit)$table),
  summary(fit)$table,
  Overall_LogRank_P = p_val,
  High_vs_Low_Cox_HR = hr_last,
  High_vs_Low_CI95_Low = ci_last[1],
  High_vs_Low_CI95_High = ci_last[2],
  check.names = FALSE
)
write.csv(out_df, csv_out, row.names = FALSE)
cat(sprintf("Saved KM plots & summary -> %s, %s, %s\n", png_out, pdf_out, csv_out))
