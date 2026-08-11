#!/usr/bin/env python3
"""Restore audited context tables from release-distributed public-data extracts."""

from pathlib import Path
import shutil
import pandas as pd

PROJECT = Path(__file__).resolve().parents[3]
SOURCE = PROJECT / "inputs/context"
OUT = PROJECT / "review_remediation/05_contextual_human"
OUT.mkdir(parents=True, exist_ok=True)

required = [
    "CONTEXTUAL_HUMAN_REGISTRY.tsv", "Context_2026_MNV.tsv",
    "Context_CD31_FVM_effects.tsv", "Context_GSE102485_samples.tsv",
    "Context_GSE146887.tsv", "Context_GSE146887_GSE163090_duplicate_audit.tsv",
    "Context_GSE160306.tsv", "Context_GSE179568_samples.tsv",
    "Context_GSE234047.tsv", "Context_GSE276892.tsv",
    "Context_GSE307925_samples.tsv", "Context_GSE60436_result.tsv",
    "Context_GSE60436_samples.tsv", "Context_GSE94019_samples.tsv",
    "MNV_VERSION_LOCK.tsv",
]
for name in required:
    source = SOURCE / name
    if not source.exists():
        raise FileNotFoundError(source)
    shutil.copy2(source, OUT / name)

registry = pd.read_csv(OUT / "CONTEXTUAL_HUMAN_REGISTRY.tsv", sep="\t")
if not registry["pooling_action"].eq("SEPARATE_CONTEXT_NO_POOLING").all():
    raise RuntimeError("Context release violates the no-pooling boundary")
spatial = pd.read_csv(OUT / "Context_GSE234047.tsv", sep="\t")
if spatial["section"].nunique() != 3 or spatial["donor"].nunique() != 1:
    raise RuntimeError("GSE234047 one-donor/three-section boundary not recovered")
lock = pd.read_csv(OUT / "MNV_VERSION_LOCK.tsv", sep="\t").iloc[0]
if lock["version"] != "v1" or lock["exact_file_name"] != "media-2.xlsx":
    raise RuntimeError("Voigt MNV release input is not locked to v1 media-2.xlsx")
print(f"Restored and validated {len(required)} context tables")
