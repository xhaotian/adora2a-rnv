#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(metafor))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[[1]]))
out_dir <- dirname(script)
project <- normalizePath(file.path(out_dir, "../.."))
effects_path <- file.path(project, "review_remediation/02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_STUDY_EFFECTS.tsv")
registry_path <- file.path(out_dir, "MOUSE_COMPARTMENT_REGISTRY.tsv")

effects <- read.delim(effects_path, check.names = FALSE)
registry <- read.delim(registry_path, check.names = FALSE)
stopifnot(nrow(effects) == 16L, nrow(registry) == 16L,
          setequal(effects$study, registry$accession))
d <- merge(effects, registry, by.x = "study", by.y = "accession", sort = FALSE)
d <- d[match(effects$study, d$study), ]

fit_set <- function(z, label) {
  fit <- rma(yi = hedges_g, vi = sampling_variance, data = z,
             method = "REML", test = "knha")
  pred <- predict(fit)
  data.frame(
    analysis_set = label,
    k = nrow(z),
    positive_effects = sum(z$hedges_g > 0),
    pooled_hedges_g = as.numeric(fit$b),
    standard_error = fit$se,
    residual_df = fit$k - fit$p,
    p = fit$pval,
    ci95_low = fit$ci.lb,
    ci95_high = fit$ci.ub,
    tau2 = fit$tau2,
    I2 = fit$I2,
    prediction_interval_low = pred$pi.lb,
    prediction_interval_high = pred$pi.ub,
    included_accessions = paste(z$study, collapse = ";")
  )
}

sets <- list(
  BROAD_ALL_RETINAL_CONTEXTS = d,
  WHOLE_RETINA_OR_LYSATE = d[d$compartment_category == "WHOLE_RETINA_OR_LYSATE", ],
  ENRICHED_OR_ISOLATED_CELL_COMPARTMENT = d[d$compartment_category == "ENRICHED_OR_ISOLATED_CELL_COMPARTMENT", ]
)
results <- do.call(rbind, Map(fit_set, sets, names(sets)))
write.table(results, file.path(out_dir, "MOUSE_COMPARTMENT_META_RESULTS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Pustejovsky-Rodgers SMD-specific funnel-asymmetry test. This matches
# meta::metabias(method.bias="Pustejovsky"): weighted regression of Hedges g
# on sqrt(1/n_OIR + 1/n_control), inverse sampling-variance weights, with a
# t test for the predictor coefficient on k-2 degrees of freedom.
d$se_star <- sqrt(1 / d$n_disease + 1 / d$n_control)
pr_fit <- lm(hedges_g ~ se_star, data = d, weights = 1 / sampling_variance)
pr_coef <- summary(pr_fit)$coefficients["se_star", ]
pr <- data.frame(
  method = "Pustejovsky-Rodgers weighted regression test for SMD funnel asymmetry",
  k = nrow(d),
  predictor = "sqrt(1/n_OIR + 1/n_control)",
  estimate = unname(pr_coef["Estimate"]),
  standard_error = unname(pr_coef["Std. Error"]),
  test_statistic_t = unname(pr_coef["t value"]),
  degrees_of_freedom = df.residual(pr_fit),
  p = unname(pr_coef["Pr(>|t|)"])
)
write.table(pr, file.path(out_dir, "FUNNEL_ASYMMETRY_DIAGNOSTIC.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

pdf(file.path(out_dir, "Supplementary_Funnel_Plot.pdf"), width = 6.5, height = 6.0,
    family = "sans", useDingbats = FALSE)
par(mar = c(4.6, 4.8, 1.2, 1.0), las = 1)
plot(d$hedges_g, d$se_star, pch = ifelse(d$compartment_category == "ENRICHED_OR_ISOLATED_CELL_COMPARTMENT", 17,
                                         ifelse(d$compartment_category == "OTHER_RETINAL_COMPARTMENT", 15, 16)),
     col = "#0072B2", xlab = "Hedges g (OIR versus normoxia)",
     ylab = expression("Modified SMD standard error  " * sqrt(1/n[OIR] + 1/n[control])),
     ylim = rev(range(c(0, d$se_star)) + c(-0.02, 0.02)))
abline(v = results$pooled_hedges_g[results$analysis_set == "BROAD_ALL_RETINAL_CONTEXTS"],
       lty = 2, col = "#D55E00", lwd = 1.4)
legend("topright", legend = c("Whole retina / lysate", "Enriched / isolated cells", "Other retinal compartment", "Pooled g"),
       pch = c(16, 17, 15, NA), lty = c(NA, NA, NA, 2),
       col = c("#0072B2", "#0072B2", "#0072B2", "#D55E00"), bty = "n", cex = 0.82)
dev.off()

whole <- results[results$analysis_set == "WHOLE_RETINA_OR_LYSATE", ]
enriched <- results[results$analysis_set == "ENRICHED_OR_ISOLATED_CELL_COMPARTMENT", ]
opposite <- sign(whole$pooled_hedges_g) != sign(enriched$pooled_hedges_g)
gate <- data.frame(
  gate = "COMPARTMENT-1",
  verdict = if (opposite) "STOP" else "PASS",
  rationale = if (opposite) {
    "Whole-retina/lysate and enriched/isolated-cell pooled directions are opposite."
  } else {
    "Whole-retina/lysate and enriched/isolated-cell pooled directions are concordant; magnitude and precision differ."
  }
)
write.table(gate, file.path(out_dir, "COMPARTMENT_STOP_GATE.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

print(results)
print(pr)
print(gate)
