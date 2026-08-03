#!/usr/bin/env python3
"""Apply the frozen, result-blind human PDR donor eligibility rules."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED = {
    "accession", "GSM", "BioSample", "donor_id", "endothelial_cells",
    "marker_support", "general_qc", "historical_six_gene_coverage",
    "duplicate_status", "metadata_complete_for_biological_unit",
}


def classify(row: pd.Series) -> tuple[str, str]:
    if row.duplicate_status != "NO_EVIDENCE_OF_DUPLICATION":
        return "UNRESOLVED", "donor/specimen duplication status is unresolved"
    if row.metadata_complete_for_biological_unit != "YES":
        return "UNRESOLVED", "independent biological unit cannot be recovered"
    if row.marker_support != "PASS" or row.general_qc != "PASS":
        return "UNRESOLVED", "endothelial marker support or general QC did not pass"
    if row.historical_six_gene_coverage != "6/6 historical score genes available":
        return "INELIGIBLE", "historical primary signature coverage is below 80%"
    if row.endothelial_cells >= 20:
        return "PRIMARY_ELIGIBLE", "independent PDR membrane donor with at least 20 QC-passing endothelial cells"
    if row.endothelial_cells >= 10:
        return "LOW_CELL_SENSITIVITY", "independent PDR membrane donor with 10-19 QC-passing endothelial cells"
    return "DESCRIPTIVE_ONLY", "fewer than 10 QC-passing endothelial cells"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("workbook", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    data = pd.read_excel(args.workbook, sheet_name="Human_eligibility_input")
    missing = REQUIRED.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    forbidden = ("ADORA", "BETA", "P_VALUE", "FOLD_CHANGE", "COMPARATOR_RANK")
    if any(any(token in str(col).upper() for token in forbidden) for col in data.columns):
        raise ValueError("The eligibility input contains a result-bearing field")
    decisions = data.apply(classify, axis=1, result_type="expand")
    decisions.columns = ["eligibility_decision", "rationale"]
    keep = ["accession", "GSM", "BioSample", "donor_id", "endothelial_cells",
            "marker_support", "general_qc", "historical_six_gene_coverage", "duplicate_status"]
    output = pd.concat([data[keep], decisions], axis=1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    counts = output.eligibility_decision.value_counts()
    if counts.get("PRIMARY_ELIGIBLE", 0) != 7 or counts.get("LOW_CELL_SENSITIVITY", 0) != 2:
        raise RuntimeError(f"Unexpected eligibility counts: {counts.to_dict()}")


if __name__ == "__main__":
    main()
