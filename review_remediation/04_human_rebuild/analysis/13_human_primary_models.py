#!/usr/bin/env python3
"""Fit donor-level six-gene sensitivity models from fresh pseudobulk inputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor


OUT = Path(__file__).resolve().parents[1]
INPUT = OUT / "human_pseudobulk_donor_audit.tsv"


def z(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    sd = float(values.std(ddof=0))
    if not np.isfinite(sd) or sd == 0:
        raise RuntimeError(f"Cannot standardize {values.name}")
    return (values - float(values.mean())) / sd


def technical_pc1(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    names = ["mean_counts", "mean_genes", "pct_mito"]
    matrix = np.column_stack([z(frame[name]).to_numpy() for name in names])
    _, _, vt = np.linalg.svd(matrix, full_matrices=False)
    loadings = vt[0].copy()
    pivot = int(np.argmax(np.abs(loadings)))
    if loadings[pivot] < 0:
        loadings *= -1
    return matrix @ loadings, loadings


def design(frame: pd.DataFrame, score_column: str, model: str) -> tuple[pd.Series, pd.DataFrame, np.ndarray]:
    y = z(frame[score_column]).rename("score_z")
    x = pd.DataFrame({"intercept": np.ones(len(frame)), "ADORA2A_z": z(frame["ADORA2A_log2cpm_plus_0_5"])}, index=frame.index)
    loadings = np.full(3, np.nan)
    if model in {"M1", "M2"}:
        x["dataset_GSE245561"] = frame["dataset"].astype(str).eq("GSE245561").astype(float)
    if model == "M2":
        x["technical_PC1"], loadings = technical_pc1(frame)
    return y, x, loadings


def fit_model(frame: pd.DataFrame, score_column: str, model: str):
    y, x, loadings = design(frame, score_column, model)
    rank = int(np.linalg.matrix_rank(x.to_numpy()))
    if rank != x.shape[1]:
        raise RuntimeError(f"{model}: rank-deficient design")
    fit = sm.OLS(y.to_numpy(), x.to_numpy()).fit()
    return fit, x, loadings


def model_rows(frame: pd.DataFrame, score_column: str, universe: str) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    regression = []
    influence = []
    lodo = []
    pc_rows = []
    for model in ["M0", "M1", "M2"]:
        fit, x, loadings = fit_model(frame, score_column, model)
        rank = int(np.linalg.matrix_rank(x.to_numpy()))
        condition = float(np.linalg.cond(x.to_numpy()))
        vif = {name: np.nan for name in x.columns}
        for index, name in enumerate(x.columns):
            if name != "intercept":
                vif[name] = float(variance_inflation_factor(x.to_numpy(), index))
        ci = fit.conf_int(alpha=0.05)
        for index, term in enumerate(x.columns):
            regression.append({
                "analysis_universe": universe,
                "model": model,
                "term": term,
                "estimate": float(fit.params[index]),
                "standard_error": float(fit.bse[index]),
                "t_statistic": float(fit.tvalues[index]),
                "residual_df": int(fit.df_resid),
                "exact_ols_p": float(fit.pvalues[index]),
                "ci95_low": float(ci[index, 0]),
                "ci95_high": float(ci[index, 1]),
                "r_squared": float(fit.rsquared),
                "adjusted_r_squared": float(fit.rsquared_adj),
                "design_rank": rank,
                "design_columns": int(x.shape[1]),
                "condition_number": condition,
                "vif": vif[term],
                "interpretive_role": "technical/identification sensitivity" if model == "M2" else "sensitivity association model",
            })

        cooks = fit.get_influence().cooks_distance[0]
        for position, (_, record) in enumerate(frame.iterrows()):
            influence.append({
                "analysis_universe": universe,
                "model": model,
                "dataset": record["dataset"],
                "GSM": record["GSM"],
                "donor_id": record["donor_id"],
                "cooks_distance": float(cooks[position]),
                "leverage": float(fit.get_influence().hat_matrix_diag[position]),
                "externally_studentized_residual": float(fit.get_influence().resid_studentized_external[position]),
            })

        for _, held in frame.iterrows():
            reduced = frame.loc[~frame["GSM"].eq(held["GSM"])].copy()
            try:
                reduced_fit, reduced_x, _ = fit_model(reduced, score_column, model)
                term_index = list(reduced_x.columns).index("ADORA2A_z")
                reduced_ci = reduced_fit.conf_int(alpha=0.05)[term_index]
                record = {
                    "analysis_universe": universe,
                    "model": model,
                    "left_out_dataset": held["dataset"],
                    "left_out_GSM": held["GSM"],
                    "left_out_donor_id": held["donor_id"],
                    "n_remaining": len(reduced),
                    "estimate": float(reduced_fit.params[term_index]),
                    "standard_error": float(reduced_fit.bse[term_index]),
                    "residual_df": int(reduced_fit.df_resid),
                    "exact_ols_p": float(reduced_fit.pvalues[term_index]),
                    "ci95_low": float(reduced_ci[0]),
                    "ci95_high": float(reduced_ci[1]),
                    "design_rank": int(np.linalg.matrix_rank(reduced_x.to_numpy())),
                    "condition_number": float(np.linalg.cond(reduced_x.to_numpy())),
                    "status": "ESTIMABLE",
                }
            except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
                record = {
                    "analysis_universe": universe,
                    "model": model,
                    "left_out_dataset": held["dataset"],
                    "left_out_GSM": held["GSM"],
                    "left_out_donor_id": held["donor_id"],
                    "n_remaining": len(reduced),
                    "estimate": np.nan,
                    "standard_error": np.nan,
                    "residual_df": np.nan,
                    "exact_ols_p": np.nan,
                    "ci95_low": np.nan,
                    "ci95_high": np.nan,
                    "design_rank": np.nan,
                    "condition_number": np.nan,
                    "status": f"NOT_ESTIMABLE: {error}",
                }
            lodo.append(record)

        if model == "M2":
            for covariate, loading in zip(["mean_counts", "mean_genes", "pct_mito"], loadings):
                pc_rows.append({"analysis_universe": universe, "technical_component": "PC1", "covariate": covariate, "loading": float(loading)})
    return regression, influence, lodo, pc_rows


def main() -> None:
    donor = pd.read_csv(INPUT, sep="\t")
    definitions = [
        ("HIGHER_CELL_COUNT_7", donor["donor_quality_status"].eq("HIGHER_CELL_COUNT"), "six_gene_score_primary7"),
        ("EXPANDED_DONOR_9", donor["donor_quality_status"].isin(["HIGHER_CELL_COUNT", "LOW_CELL_SENSITIVITY"]), "six_gene_score_expanded9"),
    ]
    all_regression, all_influence, all_lodo, all_pc = [], [], [], []
    model_inputs = []
    for universe, mask, score_column in definitions:
        frame = donor.loc[mask].copy().reset_index(drop=True)
        frame["analysis_universe"] = universe
        frame["ADORA2A_z"] = z(frame["ADORA2A_log2cpm_plus_0_5"])
        frame["score_z"] = z(frame[score_column])
        frame["technical_PC1"], _ = technical_pc1(frame)
        model_inputs.append(frame)
        regression, influence, lodo, pc = model_rows(frame, score_column, universe)
        all_regression.extend(regression)
        all_influence.extend(influence)
        all_lodo.extend(lodo)
        all_pc.extend(pc)

    regression = pd.DataFrame(all_regression)
    regression.to_csv(OUT / "human_primary_model_full_results.tsv", sep="\t", index=False)
    pd.concat(model_inputs, ignore_index=True).to_csv(OUT / "human_primary_model_input.tsv", sep="\t", index=False)
    pd.DataFrame(all_influence).to_csv(OUT / "human_influence.tsv", sep="\t", index=False)
    lodo = pd.DataFrame(all_lodo)
    lodo.to_csv(OUT / "human_LODO.tsv", sep="\t", index=False)
    pd.DataFrame(all_pc).to_csv(OUT / "human_technical_PC1_loadings.tsv", sep="\t", index=False)

    stability = (
        lodo.loc[lodo["status"].eq("ESTIMABLE")]
        .groupby(["analysis_universe", "model"], as_index=False)
        .agg(
            lodo_n=("estimate", "size"),
            estimate_min=("estimate", "min"),
            estimate_max=("estimate", "max"),
            estimate_median=("estimate", "median"),
            positive_estimates=("estimate", lambda values: int((values > 0).sum())),
            negative_estimates=("estimate", lambda values: int((values < 0).sum())),
            ci_excluding_zero=("ci95_low", lambda low: 0),
        )
    )
    for index, row in stability.iterrows():
        subset = lodo.loc[
            lodo["analysis_universe"].eq(row["analysis_universe"])
            & lodo["model"].eq(row["model"])
            & lodo["status"].eq("ESTIMABLE")
        ]
        stability.loc[index, "ci_excluding_zero"] = int(((subset["ci95_low"] > 0) | (subset["ci95_high"] < 0)).sum())
    stability.to_csv(OUT / "human_coefficient_stability.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
