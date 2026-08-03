#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(metafor)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: mouse_meta_analysis.R Source_Data.xlsx output_directory")
workbook <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

effects <- as.data.frame(read_excel(workbook, sheet = "Mouse_meta_input"))
required <- c("study", "hedges_g", "vi")
if (!all(required %in% names(effects))) stop("Mouse_meta_input lacks required columns")
if (nrow(effects) != 8) stop("The frozen synthesis expects eight studies")

fit <- rma.uni(yi = effects$hedges_g, vi = effects$vi, method = "REML", test = "knha")
prediction <- predict(fit, level = 95)
summary <- data.frame(
  k = fit$k,
  estimate = as.numeric(fit$b),
  se = fit$se,
  ci_low = fit$ci.lb,
  ci_high = fit$ci.ub,
  p = fit$pval,
  tau2 = fit$tau2,
  I2 = fit$I2,
  prediction_low = prediction$pi.lb,
  prediction_high = prediction$pi.ub,
  method = "REML with Hartung-Knapp inference"
)

leave_one_out <- do.call(rbind, lapply(seq_len(nrow(effects)), function(i) {
  q <- effects[-i, ]
  m <- rma.uni(yi = q$hedges_g, vi = q$vi, method = "REML", test = "knha")
  data.frame(removed = effects$study[[i]], k = m$k, estimate = as.numeric(m$b),
             ci_low = m$ci.lb, ci_high = m$ci.ub, p = m$pval,
             tau2 = m$tau2, I2 = m$I2)
}))

write.table(effects, file.path(out_dir, "study_effects.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(summary, file.path(out_dir, "meta_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(leave_one_out, file.path(out_dir, "leave_one_study_out.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

if (abs(summary$estimate - 1.843) > 0.001 || abs(summary$tau2 - 1.761) > 0.002 ||
    summary$prediction_low >= 0 || any(leave_one_out$estimate <= 0)) {
  stop("Recomputed results differ from the frozen synthesis")
}
