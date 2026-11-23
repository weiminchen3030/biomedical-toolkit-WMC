#!/usr/bin/env python3
"""
Clonogenic Survival Assay Analysis, Curve Fitting (IC50 / SF10 / IC10) & Plotting CLI.

Automates the full colony-formation survival workflow from raw colony counts to
SigmaPlot/publication-ready log10 survival curves and multi-model potency metrics:
  1. Parses either lab multi-block Excel/CSV sheets (multiple cell lines/drugs stacked
     with "Concentration", "Plating number", and replicate colony count columns) or
     standard tidy tables.
  2. Computes per-replicate Plating Efficiency (PE = Colonies / Plating_Number) and
     Surviving Fraction (SF = PE / PE_ctrl) along with Avg, SD, and SEM.
  3. Fits survival curves and extracts IC50 (SF=0.50), D10 / SF10 (10% survival,
     classic clonogenic IC10-survival / IC90-kill), IC10_Inhib (SF=0.90), AUC, and
     Sensitizer Enhancement Ratio (SER) across three models:
       - hill          : 2-parameter sigmoidal SF = 1 / (1 + (D / IC50)^h)
       - lq_sigmaplot  : unconstrained log-linear-quadratic ln(SF) = -alpha*D - beta*D^2
       - lq_constrained: monotonic radiobiology LQ with alpha >= 0, beta >= 0
  4. Renders log10 Surviving Fraction curves with clipped safe error bars, clean
     non-scientific tick labels (0.01, 0.05, 0.10, 0.5, 1.00), and exports editable
     SVG, TrueType PDF, 600 DPI PNG, multi-sheet Excel, CSVs, and optional .pzfx.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, lsq_linear

# Enforce editable vector fonts even if standalone
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

# Discover shared plotting & pzfx helpers when available
_REPO_ROOT = Path(__file__).resolve().parent.parent
for _cand in [
    _REPO_ROOT / "skills" / "publication-ready-plotting" / "scripts",
    Path.home() / ".gemini" / "config" / "skills" / "publication-ready-plotting" / "scripts",
    Path(__file__).resolve().parent,
]:
    if _cand.is_dir() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))

try:
    from style_presets import (
        PALETTES as PUB_PALETTES,
        add_panel_labels,
        apply_publication_style,
        mm_to_inches,
        save_publication_figure,
    )
    HAS_STYLE_PRESETS = True
except ImportError:
    HAS_STYLE_PRESETS = False
    PUB_PALETTES = {}

try:
    from export_to_prism_pzfx import dataframe_to_pzfx_bytes
    HAS_PZFX = True
except ImportError:
    HAS_PZFX = False


CLASSIC_COLORS = ["#000000", "#0000FF", "#D24B40", "#009E73", "#E69F00", "#9467BD"]
CLASSIC_MARKERS = ["o", "^", "s", "D", "v", "P"]

CLONOGENIC_PALETTES: Dict[str, List[str]] = {
    "classic": CLASSIC_COLORS,
    "okabe_ito": PUB_PALETTES.get("okabe_ito", ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]),
    "npg": PUB_PALETTES.get("npg", ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4"]),
    "aaas": PUB_PALETTES.get("aaas", ["#3B4992", "#EE0000", "#008B45", "#631879", "#008280", "#BB0021"]),
    "navy_gold": PUB_PALETTES.get("navy_gold", ["#0F4C81", "#F5B700", "#6C757D", "#C0392B", "#2E8B57"]),
}


def _to_float(val: Any) -> Optional[float]:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    s = str(val).strip()
    if not s or s.lower() == "nan":
        return None
    try:
        v = float(s)
        return v if np.isfinite(v) else None
    except ValueError:
        return None


def _read_raw_grid(path: str, sheet: Optional[str] = None) -> pd.DataFrame:
    ext = Path(path).suffix.lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path, sheet_name=sheet or 0, header=None, dtype=object)
    sep = "\t" if ext == ".tsv" else ","
    # Read line-by-line to safely handle irregular trailing columns in block CSVs
    with open(path, "r", encoding="utf-8-sig") as f:
        lines = f.read().splitlines()
    import csv

    rows = list(csv.reader(lines, delimiter=sep))
    max_cols = max((len(r) for r in rows), default=0)
    padded = [r + [None] * (max_cols - len(r)) for r in rows]
    return pd.DataFrame(padded, dtype=object)


def parse_block_format(raw_df: pd.DataFrame) -> pd.DataFrame:
    """
    Parse stacked lab clonogenic tables where each block starts with a header row
    containing 'Concentration' and 'Plating' (with Cell_Line in col 0, Drug in col 1).
    """
    records: List[Dict[str, Any]] = []
    n_rows, n_cols = raw_df.shape

    cell_line = "Sample"
    drug = "Treatment"
    conc_idx: Optional[int] = None
    plate_idx: Optional[int] = None
    current_group: Optional[str] = None

    skip_group_labels = {"cell count", "colony count", "colonies", "count", "counts"}

    for r in range(n_rows):
        row = raw_df.iloc[r].tolist()
        strs = [str(x).strip() if x is not None and not (isinstance(x, float) and np.isnan(x)) else "" for x in row]

        # Detect block header row
        lower_strs = [s.lower() for s in strs]
        c_matches = [i for i, s in enumerate(lower_strs) if s in ("concentration", "conc", "dose")]
        p_matches = [i for i, s in enumerate(lower_strs) if "plating" in s]
        if c_matches and p_matches:
            conc_idx = c_matches[0]
            plate_idx = p_matches[0]
            cell_line = strs[0] if strs[0] else cell_line
            drug = strs[1] if len(strs) > 1 and strs[1] else "Treatment"
            current_group = None
            continue

        if conc_idx is None or plate_idx is None:
            continue

        # Update group name when col 0 holds a real genotype/condition label
        if strs[0] and strs[0].lower() not in skip_group_labels:
            current_group = strs[0]

        conc_val = _to_float(row[conc_idx] if conc_idx < n_cols else None)
        plate_val = _to_float(row[plate_idx] if plate_idx < n_cols else None)
        if conc_val is None or plate_val is None or plate_val <= 0:
            continue

        # Extract the first contiguous block of numeric replicate colony counts after plate_idx
        reps: List[float] = []
        started = False
        for j in range(plate_idx + 1, n_cols):
            num = _to_float(row[j])
            if num is not None:
                started = True
                reps.append(num)
            elif started:
                break

        if not reps:
            continue

        grp_name = current_group or "Group1"
        for rep_idx, col_count in enumerate(reps, start=1):
            records.append(
                {
                    "Cell_Line": cell_line,
                    "Drug": drug,
                    "Group": grp_name,
                    "Concentration": conc_val,
                    "Plating_Number": plate_val,
                    "Replicate": f"Rep{rep_idx}",
                    "Colonies": col_count,
                }
            )

    if not records:
        raise ValueError("Could not locate any valid clonogenic blocks in the input file.")
    return pd.DataFrame(records)


def parse_tidy_format(path: str, sheet: Optional[str] = None) -> pd.DataFrame:
    """Parse tidy CSV/Excel clonogenic tables (long 'Colonies' or wide 'Rep1..RepN')."""
    ext = Path(path).suffix.lower()
    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(path, sheet_name=sheet or 0)
    else:
        df = pd.read_csv(path, sep="\t" if ext == ".tsv" else ",")

    col_map = {c.strip().lower(): c for c in df.columns}
    conc_col = next((col_map[k] for k in ("concentration", "conc", "dose") if k in col_map), None)
    if not conc_col:
        raise ValueError("Tidy format requires a 'Concentration' (or 'Conc' / 'Dose') column.")

    cell_col = next((col_map[k] for k in ("cell_line", "cellline", "cell") if k in col_map), None)
    drug_col = next((col_map[k] for k in ("drug", "treatment", "compound") if k in col_map), None)
    grp_col = next((col_map[k] for k in ("group", "genotype", "condition") if k in col_map), None)
    plate_col = next((col_map[k] for k in ("plating_number", "plating number", "plating", "seeded") if k in col_map), None)
    count_col = next((col_map[k] for k in ("colonies", "colony_count", "count") if k in col_map), None)

    work = df.copy()
    work["Cell_Line"] = work[cell_col].astype(str) if cell_col else "Sample"
    work["Drug"] = work[drug_col].astype(str) if drug_col else "Treatment"
    work["Group"] = work[grp_col].astype(str) if grp_col else "Group1"
    work["Concentration"] = pd.to_numeric(work[conc_col], errors="coerce")
    work["Plating_Number"] = pd.to_numeric(work[plate_col], errors="coerce") if plate_col else 1.0
    work = work.dropna(subset=["Concentration", "Plating_Number"])

    if count_col:
        work["Colonies"] = pd.to_numeric(work[count_col], errors="coerce")
        work = work.dropna(subset=["Colonies"])
        work["Replicate"] = work.groupby(["Cell_Line", "Drug", "Group", "Concentration"]).cumcount() + 1
        work["Replicate"] = work["Replicate"].map(lambda i: f"Rep{i}")
        return work[["Cell_Line", "Drug", "Group", "Concentration", "Plating_Number", "Replicate", "Colonies"]]

    rep_cols = [c for c in df.columns if c.lower().startswith(("rep", "colony", "count")) and c not in (count_col,)]
    if not rep_cols:
        raise ValueError("Tidy format requires either a 'Colonies' column or 'Rep1, Rep2, ...' columns.")
    melted = work.melt(
        id_vars=["Cell_Line", "Drug", "Group", "Concentration", "Plating_Number"],
        value_vars=rep_cols,
        var_name="Replicate",
        value_name="Colonies",
    )
    melted["Colonies"] = pd.to_numeric(melted["Colonies"], errors="coerce")
    return melted.dropna(subset=["Colonies"]).reset_index(drop=True)


def parse_summary_format(path: str, sheet: Optional[str] = None) -> pd.DataFrame:
    """Parse pre-aggregated summary CSV/Excel tables with Avg and SD/SEM columns."""
    ext = Path(path).suffix.lower()
    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(path, sheet_name=sheet or 0)
    else:
        df = pd.read_csv(path, sep="\t" if ext == ".tsv" else ",")

    col_map = {c.strip().lower(): c for c in df.columns}
    conc_col = next((col_map[k] for k in ("concentration", "conc", "dose") if k in col_map), None)
    avg_col = next((col_map[k] for k in ("sf_avg", "sf_mean", "avg", "mean", "sf", "survival") if k in col_map), None)
    std_col = next((col_map[k] for k in ("sf_std", "sf_sd", "std", "sd", "stdev") if k in col_map), None)
    sem_col = next((col_map[k] for k in ("sf_sem", "sem") if k in col_map), None)
    n_col = next((col_map[k] for k in ("n_reps", "n", "reps") if k in col_map), None)
    if not conc_col or not avg_col:
        raise ValueError("Summary format requires 'Concentration' (or 'Dose') and 'SF_Avg' (or 'Avg'/'Mean') columns.")

    cell_col = next((col_map[k] for k in ("cell_line", "cellline", "cell") if k in col_map), None)
    drug_col = next((col_map[k] for k in ("drug", "treatment", "compound") if k in col_map), None)
    grp_col = next((col_map[k] for k in ("group", "genotype", "condition") if k in col_map), None)

    out = pd.DataFrame({
        "Cell_Line": df[cell_col].astype(str) if cell_col else "Sample",
        "Drug": df[drug_col].astype(str) if drug_col else "Treatment",
        "Group": df[grp_col].astype(str) if grp_col else "Group1",
        "Concentration": pd.to_numeric(df[conc_col], errors="coerce"),
        "N_Reps": pd.to_numeric(df[n_col], errors="coerce").fillna(3).astype(int) if n_col else 3,
        "Plating_Number": np.nan,
        "Colonies_Mean": np.nan,
        "PE_0": 1.0,
        "PE_Avg": np.nan,
        "SF_Avg": pd.to_numeric(df[avg_col], errors="coerce"),
    })
    if std_col:
        out["SF_Std"] = pd.to_numeric(df[std_col], errors="coerce").fillna(0.0)
        out["SF_SEM"] = pd.to_numeric(df[sem_col], errors="coerce").fillna(out["SF_Std"] / np.sqrt(out["N_Reps"].clip(lower=1))) if sem_col else (out["SF_Std"] / np.sqrt(out["N_Reps"].clip(lower=1)))
    elif sem_col:
        out["SF_SEM"] = pd.to_numeric(df[sem_col], errors="coerce").fillna(0.0)
        out["SF_Std"] = out["SF_SEM"] * np.sqrt(out["N_Reps"].clip(lower=1))
    else:
        out["SF_Std"] = 0.0
        out["SF_SEM"] = 0.0

    out = out.dropna(subset=["Concentration", "SF_Avg"]).reset_index(drop=True)
    for col_ord in ("Cell_Line", "Drug", "Group"):
        order = list(dict.fromkeys(out[col_ord].tolist()))
        out[col_ord] = pd.Categorical(out[col_ord], categories=order, ordered=True)
    out = out.sort_values(["Cell_Line", "Drug", "Group", "Concentration"], kind="stable").reset_index(drop=True)
    for col_ord in ("Cell_Line", "Drug", "Group"):
        out[col_ord] = out[col_ord].astype(str)
    return out


def load_clonogenic_data(path: str, fmt: str = "auto", sheet: Optional[str] = None) -> pd.DataFrame:
    """Load raw clonogenic colony counts in either multi-block or tidy format."""
    fmt_l = fmt.lower().strip()
    if fmt_l == "block":
        return parse_block_format(_read_raw_grid(path, sheet=sheet))
    if fmt_l == "tidy":
        return parse_tidy_format(path, sheet=sheet)

    raw = _read_raw_grid(path, sheet=sheet)
    try:
        return parse_block_format(raw)
    except ValueError:
        return parse_tidy_format(path, sheet=sheet)


def compute_surviving_fractions(df_rep: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute per-replicate PE & SF and per-dose summary statistics (Avg, SD, SEM)."""
    rep = df_rep.copy()
    rep["PE"] = rep["Colonies"] / rep["Plating_Number"]

    pe0_rows = []
    for (cl, dr, grp), sub in rep.groupby(["Cell_Line", "Drug", "Group"], sort=False):
        ctrl_dose = 0.0 if (sub["Concentration"] == 0).any() else float(sub["Concentration"].min())
        pe0 = float(sub[sub["Concentration"] == ctrl_dose]["PE"].mean())
        pe0_rows.append({"Cell_Line": cl, "Drug": dr, "Group": grp, "Control_Dose": ctrl_dose, "PE_0": pe0})

    df_pe0 = pd.DataFrame(pe0_rows)
    rep = rep.merge(df_pe0, on=["Cell_Line", "Drug", "Group"], how="left")
    rep["SF"] = np.where(rep["PE_0"] > 0, rep["PE"] / rep["PE_0"], np.nan)

    summary = (
        rep.groupby(["Cell_Line", "Drug", "Group", "Concentration"], sort=False, as_index=False)
        .agg(
            N_Reps=("SF", "count"),
            Plating_Number=("Plating_Number", "first"),
            Colonies_Mean=("Colonies", "mean"),
            PE_0=("PE_0", "first"),
            PE_Avg=("PE", "mean"),
            SF_Avg=("SF", "mean"),
            SF_Std=("SF", lambda s: float(np.std(s, ddof=1)) if len(s) > 1 else 0.0),
        )
    )
    summary["SF_SEM"] = summary["SF_Std"] / np.sqrt(summary["N_Reps"].clip(lower=1))
    # Preserve file appearance order for Cell_Line, Drug, and Group while sorting doses ascending
    for col_ord in ("Cell_Line", "Drug", "Group"):
        order = list(dict.fromkeys(summary[col_ord].tolist()))
        summary[col_ord] = pd.Categorical(summary[col_ord], categories=order, ordered=True)
    summary = summary.sort_values(["Cell_Line", "Drug", "Group", "Concentration"], kind="stable").reset_index(drop=True)
    for col_ord in ("Cell_Line", "Drug", "Group"):
        summary[col_ord] = summary[col_ord].astype(str)
    return rep, summary


def _solve_lq_crossing(alpha: float, beta: float, target_sf: float) -> float:
    """
    Solve -alpha*D - beta*D^2 = ln(target_sf) for the descending survival crossing D > 0.
    """
    rhs = float(np.log(target_sf))
    # Equation: beta*D^2 + alpha*D + rhs = 0
    if abs(beta) < 1e-12:
        if alpha > 0:
            d = -rhs / alpha
            return float(d) if d > 0 else np.nan
        return np.nan

    disc = alpha * alpha - 4.0 * beta * rhs
    if disc < 0:
        return np.nan
    sqrt_d = float(np.sqrt(disc))
    roots = [(-alpha - sqrt_d) / (2.0 * beta), (-alpha + sqrt_d) / (2.0 * beta)]
    # Keep positive roots where the survival curve is descending: d/dD(-alpha*D - beta*D^2) < 0 -> alpha + 2*beta*D > 0
    desc_roots = [r for r in roots if r > 0 and (alpha + 2.0 * beta * r) > 0]
    if desc_roots:
        return float(min(desc_roots))
    return np.nan


def fit_hill_survival(doses: np.ndarray, sfs: np.ndarray) -> Dict[str, float]:
    """Fit anchored 2-parameter Hill equation SF(D) = 1 / (1 + (D / IC50)^h)."""
    max_d = float(np.max(doses)) if len(doses) else 1.0
    pos_d = doses[doses > 0]
    p0_ic50 = float(np.median(pos_d)) if len(pos_d) else max_d * 0.5
    below50 = doses[sfs <= 0.5]
    if len(below50):
        p0_ic50 = float(np.min(below50[below50 > 0])) if np.any(below50 > 0) else p0_ic50

    best_popt = None
    best_r2 = -np.inf
    for h_init in (2.5, 1.0, 5.0):
        try:
            popt, _ = curve_fit(
                lambda x, ic50, h: 1.0 / (1.0 + np.power(np.maximum(x, 0.0) / ic50, h)),
                doses,
                sfs,
                p0=[max(p0_ic50, 1e-4), h_init],
                bounds=([1e-8, 0.2], [max(max_d * 100.0, 100.0), 30.0]),
                maxfev=10000,
            )
            pred = 1.0 / (1.0 + np.power(np.maximum(doses, 0.0) / popt[0], popt[1]))
            ss_tot = float(np.sum((sfs - np.mean(sfs)) ** 2))
            r2 = 1.0 - float(np.sum((sfs - pred) ** 2)) / ss_tot if ss_tot > 0 else np.nan
            if np.isnan(r2) or r2 > best_r2:
                best_r2 = r2
                best_popt = popt
        except Exception:
            continue

    if best_popt is None:
        return {
            "Hill_IC50": np.nan,
            "Hill_D10_SF10": np.nan,
            "Hill_IC10_Inhib": np.nan,
            "Hill_Slope_h": np.nan,
            "Hill_R2": np.nan,
        }

    ic50, h = float(best_popt[0]), float(best_popt[1])
    d10_sf10 = ic50 * (9.0 ** (1.0 / h))
    ic10_inhib = ic50 * ((1.0 / 9.0) ** (1.0 / h))
    return {
        "Hill_IC50": ic50,
        "Hill_D10_SF10": d10_sf10,
        "Hill_IC10_Inhib": ic10_inhib,
        "Hill_Slope_h": h,
        "Hill_R2": float(best_r2),
    }


def fit_lq_survival(doses: np.ndarray, sfs: np.ndarray, constrained: bool = False) -> Dict[str, float]:
    """Fit Linear-Quadratic survival model ln(SF) = -alpha*D - beta*D^2 on positive SF points."""
    mask = (sfs > 0) & np.isfinite(sfs) & np.isfinite(doses)
    d_pos = doses[mask]
    sf_pos = sfs[mask]
    if len(d_pos) < 3:
        return {"alpha": np.nan, "beta": np.nan, "sf50": np.nan, "sf10": np.nan, "sf90": np.nan, "r2": np.nan}

    X = np.column_stack([-d_pos, -(d_pos ** 2)])
    y = np.log(sf_pos)
    if constrained:
        res = lsq_linear(X, y, bounds=(0.0, np.inf))
        alpha, beta = float(res.x[0]), float(res.x[1])
    else:
        ab, *_ = np.linalg.lstsq(X, y, rcond=None)
        alpha, beta = float(ab[0]), float(ab[1])

    pred_sf = np.exp(-alpha * doses - beta * (doses ** 2))
    ss_tot = float(np.sum((sfs - np.mean(sfs)) ** 2))
    r2 = 1.0 - float(np.sum((sfs - pred_sf) ** 2)) / ss_tot if ss_tot > 0 else np.nan

    return {
        "alpha": alpha,
        "beta": beta,
        "sf50": _solve_lq_crossing(alpha, beta, 0.50),
        "sf10": _solve_lq_crossing(alpha, beta, 0.10),
        "sf90": _solve_lq_crossing(alpha, beta, 0.90),
        "r2": float(r2),
    }


def fit_all_survival_curves(summary_df: pd.DataFrame, ref_group: Optional[str] = None) -> pd.DataFrame:
    """Fit Hill, SigmaPlot-LQ, and Constrained-LQ models for every (Cell_Line, Drug, Group)."""
    fit_rows: List[Dict[str, Any]] = []
    for (cl, dr, grp), sub in summary_df.groupby(["Cell_Line", "Drug", "Group"], sort=False):
        sub = sub.sort_values("Concentration")
        doses = sub["Concentration"].values.astype(float)
        sfs = sub["SF_Avg"].values.astype(float)
        max_dose = float(np.max(doses))
        min_sf = float(np.min(sfs))

        hill = fit_hill_survival(doses, sfs)
        lq_sp = fit_lq_survival(doses, sfs, constrained=False)
        lq_c = fit_lq_survival(doses, sfs, constrained=True)
        auc = float(np.trapezoid(sfs, doses) if hasattr(np, "trapezoid") else np.trapz(sfs, doses))

        ic50_val = hill["Hill_IC50"]
        sf10_val = hill["Hill_D10_SF10"]

        fit_rows.append(
            {
                "Cell_Line": cl,
                "Drug": dr,
                "Group": grp,
                "PE_Control": float(sub["PE_0"].iloc[0]),
                "Max_Dose_Tested": max_dose,
                "Min_SF_Observed": min_sf,
                "Hill_IC50_SF50": ic50_val,
                "Hill_D10_SF10_IC90": sf10_val,
                "Hill_IC10_Inhib_SF90": hill["Hill_IC10_Inhib"],
                "Hill_Slope_h": hill["Hill_Slope_h"],
                "Hill_R2": hill["Hill_R2"],
                "IC50_Extrapolated": bool(np.isnan(ic50_val) or ic50_val > max_dose or min_sf > 0.5),
                "SF10_Extrapolated": bool(np.isnan(sf10_val) or sf10_val > max_dose or min_sf > 0.1),
                "LQ_SigmaPlot_Alpha": lq_sp["alpha"],
                "LQ_SigmaPlot_Beta": lq_sp["beta"],
                "LQ_SigmaPlot_SF50": lq_sp["sf50"],
                "LQ_SigmaPlot_SF10": lq_sp["sf10"],
                "LQ_SigmaPlot_SF90": lq_sp["sf90"],
                "LQ_SigmaPlot_R2": lq_sp["r2"],
                "LQ_Mono_Alpha": lq_c["alpha"],
                "LQ_Mono_Beta": lq_c["beta"],
                "LQ_Mono_SF50": lq_c["sf50"],
                "LQ_Mono_SF10": lq_c["sf10"],
                "LQ_Mono_R2": lq_c["r2"],
                "Linear_AUC": auc,
            }
        )

    df_fits = pd.DataFrame(fit_rows)

    # Compute Sensitizer Enhancement Ratio (SER = Ref_IC / Group_IC) per (Cell_Line, Drug)
    df_fits["SER_50_vs_Ref"] = np.nan
    df_fits["SER_10_vs_Ref"] = np.nan
    for (cl, dr), idx in df_fits.groupby(["Cell_Line", "Drug"], sort=False).groups.items():
        sub = df_fits.loc[idx]
        groups_in_drug = sub["Group"].tolist()
        if ref_group and ref_group in groups_in_drug:
            r_grp = ref_group
        elif "WT" in groups_in_drug:
            r_grp = "WT"
        else:
            r_grp = groups_in_drug[0]
        ref_row = sub[sub["Group"] == r_grp].iloc[0]
        if np.isfinite(ref_row["Hill_IC50_SF50"]) and ref_row["Hill_IC50_SF50"] > 0:
            df_fits.loc[idx, "SER_50_vs_Ref"] = ref_row["Hill_IC50_SF50"] / df_fits.loc[idx, "Hill_IC50_SF50"]
        if np.isfinite(ref_row["Hill_D10_SF10_IC90"]) and ref_row["Hill_D10_SF10_IC90"] > 0:
            df_fits.loc[idx, "SER_10_vs_Ref"] = ref_row["Hill_D10_SF10_IC90"] / df_fits.loc[idx, "Hill_D10_SF10_IC90"]

    return df_fits


def evaluate_model_curve(fit_row: pd.Series, d_grid: np.ndarray, model: str) -> np.ndarray:
    """Evaluate survival curve y-values on d_grid for the requested fitting model."""
    if model == "lq_sigmaplot":
        a, b = float(fit_row["LQ_SigmaPlot_Alpha"]), float(fit_row["LQ_SigmaPlot_Beta"])
        return np.exp(-a * d_grid - b * (d_grid ** 2))
    if model == "lq_constrained":
        a, b = float(fit_row["LQ_Mono_Alpha"]), float(fit_row["LQ_Mono_Beta"])
        return np.exp(-a * d_grid - b * (d_grid ** 2))
    ic50, h = float(fit_row["Hill_IC50_SF50"]), float(fit_row["Hill_Slope_h"])
    if not np.isfinite(ic50) or not np.isfinite(h) or ic50 <= 0:
        return np.full_like(d_grid, np.nan, dtype=float)
    return 1.0 / (1.0 + np.power(np.maximum(d_grid, 0.0) / ic50, h))


def _format_dose_val(val: float, max_dose: float) -> str:
    if not np.isfinite(val):
        return f">{max_dose:g}"
    if val > max_dose:
        return f"~{val:.2f} (extrap)"
    return f"{val:.2f}"


def plot_clonogenic_panel(
    ax: plt.Axes,
    sub_sum: pd.DataFrame,
    sub_fits: pd.DataFrame,
    drug_label: str,
    unit: str = "nM",
    model: str = "hill",
    palette: str = "classic",
    ymin: float = 0.008,
    show_ref_lines: bool = False,
) -> None:
    """Plot a SigmaPlot / publication-ready log10 clonogenic survival curve onto ax."""
    colors = CLONOGENIC_PALETTES.get(palette.lower(), CLASSIC_COLORS)
    groups = sub_sum["Group"].unique().tolist()
    max_dose = float(sub_sum["Concentration"].max())
    d_grid = np.linspace(0.0, max_dose, 300)

    pos_vals = sub_sum.loc[sub_sum["SF_Avg"] > 0, "SF_Avg"].values
    eff_ymin = min(ymin, float(np.min(pos_vals)) * 0.65) if len(pos_vals) else ymin
    eff_ymax = max(1.8, float(sub_sum["SF_Avg"].max() + sub_sum["SF_Std"].max()) * 1.25)

    for idx, grp in enumerate(groups):
        col = colors[idx % len(colors)]
        mkr = CLASSIC_MARKERS[idx % len(CLASSIC_MARKERS)]
        g_sum = sub_sum[sub_sum["Group"] == grp].sort_values("Concentration")
        fit_row = sub_fits[sub_fits["Group"] == grp].iloc[0]

        # Plot positive points & asymmetric log-safe SD error bars
        vis = g_sum[g_sum["SF_Avg"] > 0]
        if not vis.empty:
            xv = vis["Concentration"].values.astype(float)
            yv = vis["SF_Avg"].values.astype(float)
            sd = vis["SF_Std"].values.astype(float)
            yerr_low = np.minimum(sd, np.maximum(yv - eff_ymin * 1.02, 0.0))
            yerr_high = sd
            ax.errorbar(
                xv,
                yv,
                yerr=[yerr_low, yerr_high],
                fmt=mkr,
                color=col,
                ecolor=col,
                markersize=4.5,
                capsize=2.5,
                elinewidth=0.9,
                capthick=0.9,
                linestyle="none",
                zorder=4,
            )

        # Fitted survival curve
        y_fit = evaluate_model_curve(fit_row, d_grid, model=model)
        mask_fit = np.isfinite(y_fit) & (y_fit >= eff_ymin * 0.9) & (y_fit <= eff_ymax * 1.15)
        if model == "hill":
            ic50_lbl = _format_dose_val(fit_row["Hill_IC50_SF50"], max_dose)
            sf10_lbl = _format_dose_val(fit_row["Hill_D10_SF10_IC90"], max_dose)
        elif model == "lq_sigmaplot":
            ic50_lbl = _format_dose_val(fit_row["LQ_SigmaPlot_SF50"], max_dose)
            sf10_lbl = _format_dose_val(fit_row["LQ_SigmaPlot_SF10"], max_dose)
        else:
            ic50_lbl = _format_dose_val(fit_row["LQ_Mono_SF50"], max_dose)
            sf10_lbl = _format_dose_val(fit_row["LQ_Mono_SF10"], max_dose)

        legend_label = f"{grp} (IC50={ic50_lbl}, SF10={sf10_lbl} {unit})"
        if np.any(mask_fit):
            ax.plot(d_grid[mask_fit], y_fit[mask_fit], color=col, lw=1.3, label=legend_label, zorder=3)
        else:
            ax.plot([], [], color=col, lw=1.3, marker=mkr, label=legend_label)

    if show_ref_lines:
        ax.axhline(0.50, color="#888888", ls="--", lw=0.6, alpha=0.6, zorder=1)
        ax.axhline(0.10, color="#888888", ls=":", lw=0.6, alpha=0.6, zorder=1)

    ax.set_yscale("log")
    ax.set_ylim(eff_ymin, eff_ymax)
    ax.set_xlim(-max_dose * 0.03, max_dose * 1.04)

    ticks = [t for t in (0.001, 0.005, 0.01, 0.05, 0.10, 0.50, 1.00) if eff_ymin * 0.95 <= t <= eff_ymax]
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(mticker.FixedFormatter([f"{t:.2f}" if t != 0.5 else "0.5" for t in ticks]))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())

    ax.set_xlabel(f"{drug_label} ({unit})")
    ax.set_ylabel("Surviving fraction")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="lower left", fontsize=5.8)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze clonogenic colony counts, fit survival curves (IC50 & SF10/IC10), and export plots/Excel/.pzfx."
    )
    parser.add_argument("--input", required=True, help="Input CSV or Excel file with colony counts.")
    parser.add_argument("--sheet", default=None, help="Optional Excel sheet name.")
    parser.add_argument("--format", default="auto", choices=["auto", "block", "tidy", "summary"], help="Input format parser (default: auto).")
    parser.add_argument("--output-dir", default="results/clonogenic_survival", help="Output directory (default: results/clonogenic_survival).")
    parser.add_argument(
        "--model",
        default="hill",
        choices=["hill", "lq_sigmaplot", "lq_constrained"],
        help="Curve model drawn on plots (default: hill). All 3 models are exported in tables.",
    )
    parser.add_argument("--unit", default="nM", help="Concentration unit label for X-axis (default: nM).")
    parser.add_argument(
        "--palette",
        default="classic",
        choices=list(CLONOGENIC_PALETTES.keys()),
        help="Color palette (default: classic SigmaPlot black/blue/red).",
    )
    parser.add_argument("--journal", default="nature", help="Publication journal preset (default: nature).")
    parser.add_argument("--ymin", type=float, default=0.008, help="Lower Y-axis bound on log10 scale (default: 0.008).")
    parser.add_argument("--ref-group", default=None, help="Reference group for SER calculation (default: WT or first group).")
    parser.add_argument("--show-ref-lines", action="store_true", help="Draw dashed horizontal lines at SF=0.50 and SF=0.10.")
    parser.add_argument("--export-pzfx", action="store_true", help="Also export GraphPad Prism (.pzfx) XY tables per assay.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.format == "summary":
        df_sum = parse_summary_format(args.input, sheet=args.sheet)
        df_rep = df_sum[["Cell_Line", "Drug", "Group", "Concentration", "SF_Avg", "SF_Std", "SF_SEM"]].copy()
    else:
        try:
            df_raw = load_clonogenic_data(args.input, fmt=args.format, sheet=args.sheet)
            df_rep, df_sum = compute_surviving_fractions(df_raw)
        except ValueError:
            if args.format == "auto":
                df_sum = parse_summary_format(args.input, sheet=args.sheet)
                df_rep = df_sum[["Cell_Line", "Drug", "Group", "Concentration", "SF_Avg", "SF_Std", "SF_SEM"]].copy()
            else:
                raise
    df_fits = fit_all_survival_curves(df_sum, ref_group=args.ref_group)

    df_rep.to_csv(out_dir / "clonogenic_replicates_sf.csv", index=False)
    df_sum.to_csv(out_dir / "clonogenic_dose_summary.csv", index=False)
    df_fits.to_csv(out_dir / "clonogenic_ic50_ic10_fits.csv", index=False)

    excel_path = out_dir / "clonogenic_survival_summary.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df_fits.to_excel(writer, sheet_name="Curve_Fit_IC50_IC10", index=False)
        df_sum.to_excel(writer, sheet_name="Per_Dose_Summary", index=False)
        df_rep.to_excel(writer, sheet_name="Replicate_Level_SF", index=False)

    print("=== Clonogenic Curve Fits (IC50, D10/SF10, IC10_Inhib) ===")
    disp_cols = [
        "Cell_Line", "Drug", "Group", "Max_Dose_Tested",
        "Hill_IC50_SF50", "Hill_D10_SF10_IC90", "Hill_IC10_Inhib_SF90", "Hill_R2",
        "LQ_SigmaPlot_SF50", "LQ_SigmaPlot_SF10", "SER_50_vs_Ref",
    ]
    print(df_fits[disp_cols].to_string(index=False))

    if HAS_STYLE_PRESETS:
        apply_publication_style(journal=args.journal, palette="okabe_ito")

    # Individual publication plots per (Cell_Line, Drug)
    pairs = list(df_sum.groupby(["Cell_Line", "Drug"], sort=False))
    for (cl, dr), sub_sum in pairs:
        sub_fits = df_fits[(df_fits["Cell_Line"] == cl) & (df_fits["Drug"] == dr)]
        w_in, h_in = mm_to_inches(95.0, 82.0) if HAS_STYLE_PRESETS else (3.74, 3.23)
        fig, ax = plt.subplots(figsize=(w_in, h_in), constrained_layout=True)
        plot_clonogenic_panel(
            ax=ax,
            sub_sum=sub_sum,
            sub_fits=sub_fits,
            drug_label=dr,
            unit=args.unit,
            model=args.model,
            palette=args.palette,
            ymin=args.ymin,
            show_ref_lines=args.show_ref_lines,
        )
        ax.set_title(f"{cl} — {dr}", fontweight="bold")
        stem = f"{cl}_{dr}_clonogenic_survival".replace(" ", "_")
        if HAS_STYLE_PRESETS:
            save_publication_figure(fig, out_dir / stem, formats=("svg", "pdf", "png"), dpi=600, journal=args.journal, check_collisions=False)
        else:
            for ext in ("svg", "pdf", "png"):
                fig.savefig(out_dir / f"{stem}.{ext}", dpi=600 if ext == "png" else 300, bbox_inches="tight")
        plt.close(fig)

        if args.export_pzfx and HAS_PZFX:
            sub_rep = df_rep[(df_rep["Cell_Line"] == cl) & (df_rep["Drug"] == dr)]
            pzfx_bytes = dataframe_to_pzfx_bytes(
                sub_rep,
                mode="xy",
                title=f"{cl}_{dr}_Surviving_Fraction",
                x_col="Concentration",
                y_col="SF",
                group_col="Group",
            )
            (out_dir / f"{stem}.pzfx").write_bytes(pzfx_bytes)

    # Multi-panel combined figure when >1 assay is present
    if len(pairs) > 1:
        n_panels = len(pairs)
        w_in, h_in = mm_to_inches(183.0, 85.0 * ((n_panels + 1) // 2)) if HAS_STYLE_PRESETS else (7.2, 3.35 * ((n_panels + 1) // 2))
        ncols = min(2, n_panels)
        nrows = (n_panels + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(w_in, h_in), constrained_layout=True)
        axes_flat = np.atleast_1d(axes).ravel()
        for idx, ((cl, dr), sub_sum) in enumerate(pairs):
            sub_fits = df_fits[(df_fits["Cell_Line"] == cl) & (df_fits["Drug"] == dr)]
            plot_clonogenic_panel(
                ax=axes_flat[idx],
                sub_sum=sub_sum,
                sub_fits=sub_fits,
                drug_label=dr,
                unit=args.unit,
                model=args.model,
                palette=args.palette,
                ymin=args.ymin,
                show_ref_lines=args.show_ref_lines,
            )
            axes_flat[idx].set_title(f"{cl} — {dr}", fontweight="bold")
        for k in range(n_panels, len(axes_flat)):
            axes_flat[k].set_visible(False)
        if HAS_STYLE_PRESETS:
            add_panel_labels(axes_flat[:n_panels], journal=args.journal)
            save_publication_figure(fig, out_dir / "clonogenic_combined_panels", formats=("svg", "pdf", "png"), dpi=600, journal=args.journal, check_collisions=False)
        else:
            for ext in ("svg", "pdf", "png"):
                fig.savefig(out_dir / f"clonogenic_combined_panels.{ext}", dpi=600 if ext == "png" else 300, bbox_inches="tight")
        plt.close(fig)

    print(f"\nSaved all tables, figures, and Excel summary -> {out_dir}")


if __name__ == "__main__":
    main()
