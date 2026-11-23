# Biomedical Toolkit (WMC) — Survival Curve Analysis Tools

Public repository containing reproducible, publication-grade survival curve analysis and plotting tools from Wei-Min Chen (`weiminchen3030`).

## Included Tools

### 1. Kaplan-Meier & Cox PH Survival Curve (`km_survival_curve.R`)
Publication-ready Kaplan-Meier survival curves with risk tables, 95% confidence bands, Log-rank *P*-values, and Cox Proportional Hazards ratios (`HR` + 95% CI).

Supports multiple stratification modes:
- `optimal` (maximally selected rank statistics via `survminer::surv_cutpoint`)
- `median` (50% median split)
- `tertile3` (3-tier tertiles: <=33%, 33–66%, >=66%)
- `tertile2` (extreme tertiles: <=33% vs >=66%)
- `quartile2` (extreme quartiles: Q1 <=25% vs Q4 >=75%)
- `group` (pre-defined categorical column)

**Usage:**
```bash
Rscript km_survival_curve.R \
  --input cohort.csv \
  --time-col OS_Months \
  --event-col OS_Status \
  --var-col PRKDC_mRNA \
  --method tertile3 \
  --label "PRKDC" \
  --out-prefix results/prkdc_km
```

Additional standalone statistical CLI scripts:
- `survival_cutpoint_analysis.R` — Optimal cutpoint, tertile, and quartile Cox/Log-rank summary table.
- `survival_tertile_analysis.R` — 3-tier & 2-tier tertile survival summary table.

### 2. Clonogenic Colony-Formation Survival Curve (`clonogenic_survival_analysis.py`)
End-to-end clonogenic assay curve fitting and plotting CLI:
- Computes Plating Efficiency (`PE`) and Surviving Fraction (`SF`) per replicate (`Mean ± SD / SEM`).
- Fits **Hill** (`SF = 1 / (1 + (D / IC50)^h)`), **LQ-SigmaPlot** (`ln(SF) = -alpha*D - beta*D^2`), and **Constrained LQ** (`alpha >= 0, beta >= 0`) models.
- Extracts `IC50`, `SF10` (`D10`), `IC10`, `AUC`, and Sensitizer Enhancement Ratio (`SER`).
- Exports editable vector `SVG`, TrueType `PDF`, `600 DPI PNG`, multi-sheet Excel, and optional GraphPad Prism `.pzfx`.

**Usage:**
```bash
python clonogenic_survival_analysis.py --input clonogenic_counts.xlsx --out-dir results/clonogenic
```
