#!/usr/bin/env bash
set -euo pipefail
release_root="$(cd "$(dirname "$0")" && pwd)"
cd "$release_root"
Rscript review_remediation/02_mouse_unit_audit/recompute_all_mouse_effects.R
Rscript review_remediation/02_mouse_compartment_audit/run_compartment_and_funnel.R
python review_remediation/04_human_rebuild/analysis/13_human_primary_models.py
python review_remediation/05_contextual_human/analysis/rebuild_context_tables_from_release_inputs.py
python review_remediation/07_figures/make_figures.py
python review_remediation/07_figures/make_supplementary_figures.py
python review_remediation/06_source_data/build_source_data.py
python review_remediation/03_systematic_search/build_prisma_checklist.py
python review_remediation/09_release/build_supplementary_tables.py
python validate_outputs.py
python write_checksums.py
