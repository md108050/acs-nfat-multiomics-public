source("src/00_setup.R")

# 出力先
out_dir_figs   <- "results/figures/rnaseq/decoupler"
out_dir_tables <- "results/tables/rnaseq"

dir.create(out_dir_figs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# 1. データの読み込み ----
res <- readRDS("data/processed/rnaseq/res_deseq2.rds")

logFC_df <- as.data.frame(res) |>
  rownames_to_column("gene") |>
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange)) |>
  dplyr::select(gene, log2FoldChange) |>
  column_to_rownames("gene")

# 2. ネットワークモデルの取得 ----
progeny_net   <- get_progeny(organism = "human", top = 500)
collectri_net <- get_collectri(organism = "human", split_complexes = FALSE)

# 3. Pathway 活性の推論 (PROGENy) ----
pathway_acts <- run_mlm(
  mat     = logFC_df,
  net     = progeny_net,
  .source = "source",
  .target = "target",
  .mor    = "weight",
  minsize = 5
)

write_csv(pathway_acts, file.path(out_dir_tables, "decoupler_pathway.csv"))

pathway_top7 <- pathway_acts |>
  filter(score > 0) |>
  slice_max(order_by = score, n = 7) |>
  mutate(source = factor(source, levels = source))

p_pathway <- ggplot(pathway_top7, aes(x = reorder(source, score), y = score, color = source)) +
  geom_segment(aes(xend = source, yend = 0), linewidth = 1) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_viridis_d(option = "rocket", begin = 0, end = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.75, colour = "gray21") +
  coord_flip() +
  theme_bw() +
  labs(x = "Pathways", y = "Score") +
  theme(
    text            = element_text(family = "Helvetica", size = 12),
    legend.position = "none",
    axis.title.y    = element_blank()
  )

ggsave(file.path(out_dir_figs, "pathway.pdf"), p_pathway, width = 5, height = 4)
ggsave(file.path(out_dir_figs, "pathway.png"), p_pathway, width = 5, height = 4, dpi = 300)

# 4. Transcription Factor 活性の推論 (CollecTRI) ----
tf_acts <- run_ulm(
  mat     = logFC_df,
  net     = collectri_net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)

write_csv(tf_acts, file.path(out_dir_tables, "decoupler_tf.csv"))

tf_top7 <- tf_acts |>
  filter(score > 0) |>
  slice_max(order_by = score, n = 7) |>
  mutate(source = factor(source, levels = source))

p_tf <- ggplot(tf_top7, aes(x = reorder(source, score), y = score, color = source)) +
  geom_segment(aes(xend = source, yend = 0), linewidth = 1) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_viridis_d(option = "rocket", begin = 0, end = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.75, colour = "gray21") +
  coord_flip() +
  theme_bw() +
  labs(x = "Transcription factor", y = "Score") +
  theme(
    text            = element_text(family = "Helvetica", size = 12),
    legend.position = "none",
    axis.title.y    = element_blank()
  )

ggsave(file.path(out_dir_figs, "tf.pdf"), p_tf, width = 5, height = 4)
ggsave(file.path(out_dir_figs, "tf.png"), p_tf, width = 5, height = 4, dpi = 300)
