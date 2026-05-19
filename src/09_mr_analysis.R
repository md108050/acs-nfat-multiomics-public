source("src/00_setup.R")

# 出力先
out_dir_figs   <- "results/figures/mr"
out_dir_tables <- "results/tables/mr"

dir.create(out_dir_figs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# 1. 前処理済みデータの読み込み ----
cer_mr_df <- read_csv("data/processed/mr/cer_mr_exposure.csv", show_col_types = FALSE)

# 2. Exposure データのフォーマット整形 ----
exposure_dat <- format_data(
  cer_mr_df,
  type              = "exposure",
  snp_col           = "SNP",
  beta_col          = "beta",
  se_col            = "se",
  eaf_col           = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col  = "other_allele",
  pval_col          = "pval"
)

exposure_dat$exposure    <- "Cer(d18:1/16:0)"
exposure_dat$id.exposure <- "ceramides"

# 3. Outcome データの取得 ----
# Whole body fat mass (MRC-IEU: ukb-b-19393)
outcome_dat <- extract_outcome_data(
  snps     = exposure_dat$SNP,
  outcomes = "ukb-b-19393",
  proxies  = TRUE
)

# 4. データのハーモナイゼーション ----
harmonised_dat <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat  = outcome_dat,
  action       = 2
)

saveRDS(harmonised_dat, "data/processed/mr/harmonised_dat.rds")

# 5. MR解析 ----
mr_results <- mr(harmonised_dat)
write_csv(mr_results, file.path(out_dir_tables, "mr_results.csv"))

# 6. 感度分析 ----
mr_het   <- mr_heterogeneity(harmonised_dat)
mr_pleio <- mr_pleiotropy_test(harmonised_dat)

write_csv(mr_het,                    file.path(out_dir_tables, "mr_heterogeneity.csv"))
write_csv(as.data.frame(mr_pleio),   file.path(out_dir_tables, "mr_pleiotropy.csv"))

res_loo <- mr_leaveoneout(harmonised_dat)
pdf(file.path(out_dir_figs, "mr_leaveoneout.pdf"), width = 8, height = 6)
mr_leaveoneout_plot(res_loo)
dev.off()

# 7. MR-PRESSO ----
mr_presso_df <- data.frame(
  SNP          = harmonised_dat$SNP,
  BetaOutcome  = harmonised_dat$beta.outcome,
  BetaExposure = harmonised_dat$beta.exposure,
  SdOutcome    = harmonised_dat$se.outcome,
  SdExposure   = harmonised_dat$se.exposure
)

set.seed(1234)
mr_presso_result <- mr_presso(
  BetaOutcome     = "BetaOutcome",
  BetaExposure    = "BetaExposure",
  SdOutcome       = "SdOutcome",
  SdExposure      = "SdExposure",
  OUTLIERtest     = TRUE,
  DISTORTIONtest  = TRUE,
  data            = mr_presso_df,
  NbDistribution  = 1000,
  SignifThreshold = 0.05
)

saveRDS(mr_presso_result, "data/processed/mr/mr_presso_result.rds")

# 8. 外れ値除外後の感度分析 ----
outlier_snps        <- harmonised_dat[c(3, 5, 6), "SNP"]
harmonised_filtered <- harmonised_dat[!harmonised_dat$SNP %in% outlier_snps, ]

mr_filtered_results <- mr(harmonised_filtered)
mr_het_filtered     <- mr_heterogeneity(harmonised_filtered)
mr_pleio_filtered   <- mr_pleiotropy_test(harmonised_filtered)

write_csv(mr_filtered_results,              file.path(out_dir_tables, "mr_results_filtered.csv"))
write_csv(mr_het_filtered,                  file.path(out_dir_tables, "mr_heterogeneity_filtered.csv"))
write_csv(as.data.frame(mr_pleio_filtered), file.path(out_dir_tables, "mr_pleiotropy_filtered.csv"))

# 9. Forest plot ----
mr_df <- data.frame(
  method = c("IVW", "Weighted median", "MR Egger"),
  nsnp   = c(3, 3, 3),
  b      = c(0.06714629, 0.04154400, 0.18191118),
  se     = c(0.03384040, 0.01575718, 0.09882273)
)

mr_df$lower  <- mr_df$b - 1.96 * mr_df$se
mr_df$upper  <- mr_df$b + 1.96 * mr_df$se
mr_df$method <- factor(mr_df$method, levels = rev(mr_df$method))

method_colors <- c(
  "IVW"             = "#C72C48",
  "Weighted median" = "#C72C48",
  "MR Egger"        = "#4D4D4D"
)

method_shapes <- c(
  "IVW"             = 15,
  "Weighted median" = 15,
  "MR Egger"        = 15
)

p_forest <- ggplot(mr_df, aes(x = b, y = method)) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper, color = method),
    orientation = "y",
    width     = 0.25,
    linewidth = 1.1
  ) +
  geom_point(aes(shape = method, color = method), size = 4) +
  scale_shape_manual(values = method_shapes) +
  scale_color_manual(values = method_colors) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_x_continuous(
    breaks = seq(-0.05, 0.3, by = 0.1),
    limits = c(-0.05, 0.3),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x        = "Effect per SD increase in genetically-predicted Cer(d18:1/16:0)",
    y        = "",
    subtitle = "Outcome: Whole body fat mass"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y        = element_text(size = 14),
    axis.text.x        = element_text(size = 14),
    axis.title.x       = element_text(hjust = 0.5),
    legend.position    = "none"
  )

ggsave(file.path(out_dir_figs, "mr_forest.pdf"), p_forest, width = 10, height = 3, units = "in")
ggsave(file.path(out_dir_figs, "mr_forest.png"), p_forest, width = 10, height = 3, units = "in", dpi = 300)