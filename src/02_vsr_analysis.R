source("src/00_setup.R")

# 出力先
out_dir_figs   <- "results/figures/vsr"
out_dir_tables <- "results/tables/vsr"
dir.create(out_dir_figs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# 1. データの読み込み ----
df_cohort <- read_parquet("data/meta/cohort_191.parquet")

# 2. フィルタリング ----
# slice_omatic == 1: CTで内臓脂肪・皮下脂肪が測定されている症例
df <- df_cohort |> filter(slice_omatic == 1)

# 3. 共変量調整 ----
# 共変量のみのモデルでlog(VSR + 1)を説明
m_cov <- lm(log(VSR + 1) ~ age + sex + BMI + smoking + alcohol, data = df)

# 残差に平均適合値を加えて中心化（logスケール）
y_adj_log <- resid(m_cov) + mean(fitted(m_cov), na.rm = TRUE)

# 原スケールへ逆変換（log1p -> expm1）
df <- df |> mutate(VSR_adj = exp(y_adj_log) - 1)

# 4. 共変量調整後の診断群間比較 ----
fit <- lm(log(VSR + 1) ~ Diagnosis_2groups + age + sex + BMI + smoking + alcohol, data = df)
anova(fit)

# 5. 要約統計量 ----
summ_vsr <- df |>
  summarise(
    n      = sum(!is.na(VSR_adj)),
    median = median(VSR_adj, na.rm = TRUE),
    q1     = quantile(VSR_adj, 0.25, na.rm = TRUE),
    q3     = quantile(VSR_adj, 0.75, na.rm = TRUE),
    .by    = Diagnosis_2groups
  )

write_csv(summ_vsr, file.path(out_dir_tables, "VSR_summary.csv"))

# 6. 可視化 ----
p_box <- ggplot(df, aes(x = Diagnosis_2groups, y = VSR_adj, fill = Diagnosis_2groups)) +
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  labs(x = NULL, y = "Adjusted VSR (residualized)") +
  theme_minimal(base_family = "Helvetica") +
  theme(legend.position = "none")

ggsave(
  filename = file.path(out_dir_figs, "VSR_boxplot.pdf"),
  plot     = p_box,
  width    = 15,
  height   = 15,
  units    = "cm"
)
