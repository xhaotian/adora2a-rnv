#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(metafor))

args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1]]))
out_dir <- dirname(script)
x <- read.delim(file.path(out_dir, "MOUSE_FINAL_SAMPLE_EXPRESSION.tsv"), check.names = FALSE)

effect_one <- function(d) {
  disease <- d$Adora2a[grepl("OIR", d$group, ignore.case = TRUE)]
  control <- d$Adora2a[!grepl("OIR", d$group, ignore.case = TRUE)]
  stopifnot(length(disease) >= 2L, length(control) >= 2L)
  e <- escalc(measure = "SMD", m1i = mean(disease), sd1i = sd(disease), n1i = length(disease),
              m2i = mean(control), sd2i = sd(control), n2i = length(control))
  data.frame(study = unique(d$study), n_disease = length(disease), n_control = length(control),
             mean_disease = mean(disease), mean_control = mean(control),
             sd_disease = sd(disease), sd_control = sd(control),
             hedges_g = as.numeric(e$yi), sampling_variance = as.numeric(e$vi),
             ci95_low = as.numeric(e$yi - qnorm(.975) * sqrt(e$vi)),
             ci95_high = as.numeric(e$yi + qnorm(.975) * sqrt(e$vi)))
}

effects <- do.call(rbind, lapply(split(x, x$study), effect_one))
effects <- effects[order(effects$study), ]
write.table(effects, file.path(out_dir, "MOUSE_FINAL_PRIMARY_STUDY_EFFECTS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

fit_set <- function(d, label) {
  fit <- rma(yi = hedges_g, vi = sampling_variance, data = d, method = "REML", test = "knha")
  pred <- predict(fit)
  data.frame(analysis_set = label, k = nrow(d), included_accessions = paste(d$study, collapse = ";"),
             pooled_hedges_g = as.numeric(fit$b), standard_error = fit$se,
             residual_df = fit$k - fit$p, p = fit$pval, ci95_low = fit$ci.lb, ci95_high = fit$ci.ub,
             tau2 = fit$tau2, I2 = fit$I2, prediction_interval_low = pred$pi.lb,
             prediction_interval_high = pred$pi.ub, positive_studies = sum(d$hedges_g > 0))
}

primary <- fit_set(effects, "FINAL_PRIMARY_DEPOSITED_SAMPLE_SYNTHESIS")
write.table(primary, file.path(out_dir, "MOUSE_FINAL_PRIMARY_META_RESULTS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
loo <- do.call(rbind, lapply(effects$study, function(omitted) {
  z <- fit_set(effects[effects$study != omitted, ], "FINAL_PRIMARY_LOO")
  z$omitted_accession <- omitted
  z
}))
write.table(loo, file.path(out_dir, "MOUSE_FINAL_PRIMARY_LOO.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
strict <- effects[effects$study %in% c("GSE234447", "GSE315511"), ]
strict_result <- fit_set(strict, "STRICT_UNIT_CONFIRMED_SENSITIVITY")
write.table(strict_result, file.path(out_dir, "MOUSE_STRICT_UNIT_CONFIRMED_META_RESULTS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
print(primary)
