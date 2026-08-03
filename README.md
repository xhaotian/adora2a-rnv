# ADORA2A across mouse OIR and human retinal neovascularization

This repository contains the analysis code supporting the study:

> Cross-dataset transcriptomic assessment of ADORA2A in mouse OIR and human retinal neovascularization

The study synthesizes eight independent public P17 oxygen-induced retinopathy (OIR) studies and evaluates prespecified human retinal-neovascularization evidence. The mouse random-effects estimate was positive, but heterogeneity was substantial and its 95% prediction interval crossed zero. A result-blind, uniformly rebuilt analysis of seven eligible human proliferative diabetic retinopathy donors did not provide reproducible support. The code therefore supports a bounded cross-dataset assessment, not human validation, target designation, or therapeutic efficacy.

## Repository contents

- `analysis/mouse_meta_analysis.R`: Hedges g study effects, REML random-effects synthesis, Hartung–Knapp inference, prediction interval, and leave-one-study-out analysis.
- `analysis/human_primary_models.py`: prespecified donor-level unadjusted, dataset-adjusted, and technically adjusted models, including VIF and residual degrees of freedom.
- `analysis/apply_human_eligibility.py`: deterministic implementation of the result-blind human eligibility rules.
- `docs/DATA_ACCESS.md`: public accession and input-file guidance.
- `docs/INPUT_OUTPUT_MAP.md`: mapping between inputs, scripts, and outputs.
- `code_manifest.tsv`: file-level release manifest.

## Reproduce the principal results

1. Install Python 3.12 and R 4.5 or later.
2. Install the pinned Python packages and the required R packages:

   ```bash
   python -m pip install -r requirements.txt
   Rscript -e 'install.packages(c("metafor", "readxl"))'
   ```

3. Place the article-associated `Source_Data.xlsx` and `Statistical_Results.xlsx` files in `inputs/`.
4. Run:

   ```bash
   Rscript analysis/mouse_meta_analysis.R inputs/Source_Data.xlsx results/mouse
   python analysis/apply_human_eligibility.py inputs/Source_Data.xlsx results/human_eligibility.tsv
   python analysis/human_primary_models.py inputs/Source_Data.xlsx results/human_models.tsv
   ```

The expected frozen estimates are listed in `docs/EXPECTED_RESULTS.md`. Scripts fail with a non-zero exit status if required columns or sheets are absent.

## Data scope

No new sequencing data were generated. Raw public data remain under their source repositories and accession records. Article-associated source data are distributed with the manuscript rather than duplicated in this code repository.

## License and citation

Code is released under the MIT License. Please cite the associated article when using this analysis workflow.
