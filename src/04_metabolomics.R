source("src/00_setup.R")

# 出力先
out_dir_figs   <- "results/figures/metabolomics"
out_dir_tables <- "results/tables/metabolomics"

dir.create(out_dir_figs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# 1. データの読み込み ----
df_cohort     <- read_parquet("data/meta/cohort_191.parquet")
df_metab_pre  <- read_csv("data/processed/metabolomics/norm_pre.csv",  show_col_types = FALSE)
df_metab_post <- read_csv("data/processed/metabolomics/norm_post.csv", show_col_types = FALSE)
df_metab_men  <- read_csv("data/processed/metabolomics/norm_men.csv",  show_col_types = FALSE)

metab_list <- list(
  "Premenopausal women"  = df_metab_pre,
  "Postmenopausal women" = df_metab_post,
  "Men"                  = df_metab_men
)

# 2. ステロイド・パスウェイの定義 ----
steroids21_list <- c(
  "Preg", "Prog", "Allo_Pregnanolone", "DOC", "DHDOC_5alpha_DHDOC_5beta",
  "B", "DHB_11", "OHB_18", "A", "DOF_11", "F", "E", "OHF_18", "OxoF_18",
  "DHEA_S", "An_G", "OHA4_7alpha", "OHT_7alpha", "OHA4_11beta", "A4_11_keto", "OHT_11beta"
)

pathway_list <- list(
  Mineralocorticoid = c("Preg", "Prog", "Allo_Pregnanolone", "DOC", "DHDOC_5alpha_DHDOC_5beta", "B", "DHB_11", "OHB_18", "A"),
  Glucocorticoid    = c("DOF_11", "F", "E", "OHF_18", "OxoF_18"),
  Androgen          = c("DHEA_S", "An_G", "OHA4_7alpha", "OHT_7alpha", "OHA4_11beta", "A4_11_keto", "OHT_11beta")
)

# 3. 相関分析関数 ----
calculate_correlations <- function(cohort_df, metab_list_input, target_disease) {
  groups  <- c("Premenopausal women", "Postmenopausal women", "Men")
  results <- list()

  for (grp in groups) {
    df_group <- cohort_df |>
      filter(Diagnosis_2groups == target_disease, menopausal_status == grp) |>
      inner_join(metab_list_input[[grp]], by = "Research_ID")

    if (nrow(df_group) > 5) {
      grp_res <- list()

      for (pw in names(pathway_list)) {
        st_list <- pathway_list[[pw]]

        cor_res <- corr.test(
          x      = df_group[, st_list],
          y      = df_group[, "VSR", drop = FALSE],
          method = "spearman",
          adjust = "none"
        )

        fdr <- p.adjust(cor_res$p, method = "BH")

        grp_res[[pw]] <- tibble(
          Pathway = pw,
          Steroid = rownames(cor_res$r),
          r       = as.numeric(cor_res$r),
          P_value = as.numeric(cor_res$p),
          FDR     = as.numeric(fdr),
          Group   = grp
        )
      }
      results[[grp]] <- bind_rows(grp_res)
    }
  }

  bind_rows(results)
}

# 4. ACSにおける相関分析 ----
cor_acs <- calculate_correlations(df_cohort, metab_list, "ACS")
write_csv(cor_acs, file.path(out_dir_tables, "correlation_ACS.csv"))

# 5. Correlation Heatmap (ACS) ----
cor_acs_fig <- cor_acs |>
  mutate(
    Steroid      = factor(Steroid, levels = rev(steroids21_list)),
    Group        = factor(Group, levels = c("Premenopausal women", "Postmenopausal women", "Men")),
    Pathway      = factor(Pathway, levels = c("Mineralocorticoid", "Glucocorticoid", "Androgen")),
    Significance = ifelse(FDR < 0.1, "*", "")
  )

p_heatmap <- ggplot(cor_acs_fig, aes(x = Group, y = Steroid, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Significance), vjust = 0.7, size = 6, color = "black") +
  scale_fill_gradientn(
    colors = c("#053061", "#4393c3", "white", "#f4a582", "#67001f"),
    limits = c(-1, 1),
    name   = "Spearman r"
  ) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 3) +
  scale_x_discrete(position = "top", labels = c("Pre", "Post", "Men")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text  = element_text(size = 14, face = "bold"),
    panel.grid  = element_blank()
  )

ggsave(file.path(out_dir_figs, "heatmap_ACS.pdf"), p_heatmap, width = 12, height = 5)
ggsave(file.path(out_dir_figs, "heatmap_ACS.png"), p_heatmap, width = 12, height = 5, dpi = 300)

# 6. NFATにおける相関分析 ----
cor_nfat <- calculate_correlations(df_cohort, metab_list, "NFAT")
write_csv(cor_nfat, file.path(out_dir_tables, "correlation_NFAT.csv"))

# 7. 散布図 (Premenopausal ACS: VSR vs 11-DOC) ----
df_pre_acs <- df_cohort |>
  filter(Diagnosis_2groups == "ACS", menopausal_status == "Premenopausal women") |>
  inner_join(df_metab_pre, by = "Research_ID")

p_scatter <- ggplot(df_pre_acs, aes(x = DOC, y = VSR)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", col = "red", fill = "gray80") +
  stat_cor(method = "spearman") +
  labs(x = "11-DOC (Normalized Z-score)", y = "VSR") +
  theme_minimal()

ggsave(file.path(out_dir_figs, "scatter_pre_ACS_DOC_VSR.pdf"), p_scatter, width = 5, height = 5)
ggsave(file.path(out_dir_figs, "scatter_pre_ACS_DOC_VSR.png"), p_scatter, width = 5, height = 5, dpi = 300)

# 8. 交絡因子調整重回帰分析 ----
df_post_acs <- df_cohort |>
  filter(Diagnosis_2groups == "ACS", menopausal_status == "Postmenopausal women") |>
  inner_join(df_metab_post, by = "Research_ID")

df_men_acs <- df_cohort |>
  filter(Diagnosis_2groups == "ACS", menopausal_status == "Men") |>
  inner_join(df_metab_men, by = "Research_ID")

results_lm <- list(
  "Premenopausal (11-DOC)"  = tidy(lm(VSR ~ DOC    + cortisol + BMI + smoking + alcohol, data = df_pre_acs))  |> mutate(Group = "Premenopausal women",  Model = "VSR ~ 11-DOC + Covariates"),
  "Postmenopausal (18-OHB)" = tidy(lm(VSR ~ OHB_18 + cortisol + BMI + smoking + alcohol, data = df_post_acs)) |> mutate(Group = "Postmenopausal women", Model = "VSR ~ 18-OHB + Covariates"),
  "Men (11-Deoxycortisol)"  = tidy(lm(VSR ~ DOF_11 + cortisol + BMI + smoking + alcohol, data = df_men_acs))  |> mutate(Group = "Men",                  Model = "VSR ~ 11-Deoxycortisol + Covariates")
)

supp_table_5 <- bind_rows(results_lm) |> relocate(Group, Model)
write_csv(supp_table_5, file.path(out_dir_tables, "multiple_regression_ACS.csv"))
