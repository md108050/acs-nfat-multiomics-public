source("src/00_setup.R")

# 出力先
out_dir_figs_umap    <- "results/figures/rnaseq/umap"
out_dir_figs_volcano <- "results/figures/rnaseq/volcano"
out_dir_tables       <- "results/tables/rnaseq"

dir.create(out_dir_figs_umap,    recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_figs_volcano, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables,       recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/rnaseq", recursive = TRUE, showWarnings = FALSE)

# 1. データの読み込み ----
meta   <- read.csv("data/meta/rnaseq_meta.csv", row.names = 1)
counts <- read.csv("data/raw/rnaseq_human/rnaseq_count.csv", row.names = 1, check.names = FALSE)

meta$Diagnosis_2groups <- factor(meta$Diagnosis_2groups, levels = c("NFAT", "ACS"))
meta$batch             <- as.factor(meta$batch)

if (!all(rownames(meta) == colnames(counts))) {
  stop("Rownames of metadata and colnames of counts do not match.")
}

# 2. バッチ補正 (ComBat_seq) ----
counts_matrix <- as.matrix(counts)
set.seed(0)
counts_batch_corrected <- ComBat_seq(counts_matrix, batch = meta$batch)

saveRDS(counts_batch_corrected, "data/processed/rnaseq/counts_batch_corrected.rds")

# 3. DESeq2 の実行 ----
dds <- DESeqDataSetFromMatrix(
  countData = counts_batch_corrected,
  colData   = meta,
  design    = ~ Diagnosis_2groups
)

keep <- rowSums(counts(dds)) > 1
dds  <- dds[keep, ]
dds  <- DESeq(dds)
res  <- results(dds, contrast = c("Diagnosis_2groups", "ACS", "NFAT"))

saveRDS(dds, "data/processed/rnaseq/dds.rds")
saveRDS(res, "data/processed/rnaseq/res_deseq2.rds")

# 4. DEG の抽出 (FDR < 0.05, |log2FC| > 1) ----
res_df <- as.data.frame(res) |>
  rownames_to_column("gene_symbol") |>
  mutate(
    padj           = ifelse(is.na(padj), 1, padj),
    log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
    DEG_Status     = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                              ~ "Not_Significant"
    )
  )

write_csv(res_df,                               file.path(out_dir_tables, "DEG_all.csv"))
write_csv(filter(res_df, DEG_Status == "Up"),   file.path(out_dir_tables, "DEG_up.csv"))
write_csv(filter(res_df, DEG_Status == "Down"), file.path(out_dir_tables, "DEG_down.csv"))

# 5. UMAP プロット ----
vsd     <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

set.seed(0)
umap_res <- umap(t(vst_mat))
umap_df  <- data.frame(
  UMAP1     = umap_res$layout[, 1],
  UMAP2     = umap_res$layout[, 2],
  Diagnosis = meta$Diagnosis_2groups
)

umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = Diagnosis)) +
  geom_point(shape = 21, size = 3, color = "black", stroke = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("NFAT" = "darkgray", "ACS" = "lightpink")) +
  theme_bw() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(
    text            = element_text(family = "Helvetica", size = 12),
    legend.position = "top"
  )

ggsave(file.path(out_dir_figs_umap, "umap.pdf"), umap_plot, width = 6, height = 5)
ggsave(file.path(out_dir_figs_umap, "umap.png"), umap_plot, width = 6, height = 5, dpi = 300)

# 6. Volcano Plot ----
highlight_genes <- c("NFKBIA", "MT1M", "PDK4", "PTX3", "SERPINE1", "MT2A", "MERTK", "VEGFA", "GADD45G")

res_df <- res_df |> mutate(log10FDR = -log10(padj))

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = log10FDR)) +
  geom_point(
    data  = filter(res_df, DEG_Status == "Not_Significant"),
    color = "lightgray", alpha = 0.5, size = 1
  ) +
  geom_point(
    data  = filter(res_df, DEG_Status != "Not_Significant" & !gene_symbol %in% highlight_genes),
    color = "black", alpha = 0.6, size = 1.5
  ) +
  geom_point(
    data  = filter(res_df, gene_symbol %in% highlight_genes),
    color = "red", size = 2
  ) +
  geom_text_repel(
    data          = filter(res_df, gene_symbol %in% highlight_genes),
    aes(label     = gene_symbol),
    family        = "Helvetica", size = 3, fontface = "bold",
    box.padding   = 0.5, point.padding = 0.3
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray20", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray20", linewidth = 0.5) +
  theme_bw() +
  labs(
    x = bquote(Log[2](Fold~Change)),
    y = bquote(-Log[10](FDR))
  ) +
  theme(
    text             = element_text(family = "Helvetica", size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir_figs_volcano, "volcano.pdf"), volcano_plot, width = 6, height = 5)
ggsave(file.path(out_dir_figs_volcano, "volcano.png"), volcano_plot, width = 6, height = 5, dpi = 300)
