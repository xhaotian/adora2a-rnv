# Input-output map

| Input | Script | Output |
|---|---|---|
| `Source_Data.xlsx`, sheet `Mouse_meta_input` | `analysis/mouse_meta_analysis.R` | Study effects, pooled REML/Hartung–Knapp result, prediction interval, leave-one-study-out table |
| `Source_Data.xlsx`, sheet `Human_eligibility_input` | `analysis/apply_human_eligibility.py` | Result-blind donor eligibility decisions |
| `Source_Data.xlsx`, sheet `Human_model_input` | `analysis/human_primary_models.py` | Three prespecified human model estimates, confidence intervals, VIF, design rank, and residual df |

The code does not select additional datasets, models, thresholds, or endpoints.
