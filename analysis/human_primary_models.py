#!/usr/bin/env python3
"""Reproduce the three prespecified donor-level human PDR models."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor


REQUIRED = {
    "GSM", "dataset", "ADORA2A_log2cpm", "signature_score",
    "mean_counts", "mean_genes", "pct_mito",
}


def standardize(values: pd.Series) -> np.ndarray:
    values = values.astype(float).to_numpy()
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros(len(values))
    return (values - values.mean()) / sd


def fit_model(data: pd.DataFrame, model: str) -> dict[str, float | int | str]:
    x = standardize(data.ADORA2A_log2cpm)
    y = standardize(data.signature_score)
    columns = [np.ones(len(data)), x]
    names = ["intercept", "ADORA2A"]
    if model in {"M1_DATASET", "M2_TECHNICAL"}:
        levels = sorted(data.dataset.astype(str).unique())
        if len(levels) != 2:
            raise ValueError("The frozen model expects exactly two source datasets")
        columns.append((data.dataset.astype(str) == levels[1]).astype(float).to_numpy())
        names.append("dataset")
    if model == "M2_TECHNICAL":
        technical = np.column_stack([
            standardize(data.mean_counts), standardize(data.mean_genes), standardize(data.pct_mito)
        ])
        _, _, vt = np.linalg.svd(technical, full_matrices=False)
        columns.append(technical @ vt[0])
        names.append("technical_PC1")
    design = np.column_stack(columns)
    rank = int(np.linalg.matrix_rank(design))
    result = sm.OLS(y, design).fit()
    index = names.index("ADORA2A")
    ci = result.conf_int(alpha=0.05)[index]
    vif = float(variance_inflation_factor(design, index))
    return {
        "model": model,
        "n_donors": len(data),
        "beta": float(result.params[index]),
        "ci_low": float(ci[0]),
        "ci_high": float(ci[1]),
        "standard_error": float(result.bse[index]),
        "ols_p": float(result.pvalues[index]),
        "design_rank": rank,
        "design_columns": design.shape[1],
        "residual_df": int(result.df_resid),
        "condition_number": float(np.linalg.cond(design)),
        "exposure_vif": vif,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("workbook", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    data = pd.read_excel(args.workbook, sheet_name="Human_model_input")
    missing = REQUIRED.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    results = pd.DataFrame([fit_model(data, model) for model in
                            ["M0_UNADJUSTED", "M1_DATASET", "M2_TECHNICAL"]])
    expected = np.array([0.48557960664938105, -0.3000491847787868, -1.4016774784030996])
    if not np.allclose(results.beta.to_numpy(), expected, atol=1e-10):
        raise RuntimeError("Model estimates differ from the frozen analysis")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
