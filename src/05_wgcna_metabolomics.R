source("src/00_setup.R")

# 出力先
out_dir_figs   <- "results/figures/wgcna"
out_dir_tables <- "results/tables/wgcna"

dir.create(out_dir_figs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# 1. データの読み込み ----
meta                   <- read.csv("data/meta/rnaseq_meta.csv", row.names = 1)
counts_batch_corrected <- readRDS("data/processed/rnaseq/counts_batch_corrected.rds")
df_metab_all           <- read_csv("data/processed/metabolomics/norm_global.csv", show_col_types = FALSE)

# 2. メタボロミクスと共通のサンプルに絞り込む ----
common_samples <- intersect(rownames(meta), df_metab_all$Research_ID)

meta_wgcna   <- meta[common_samples, , drop = FALSE]
counts_wgcna <- counts_batch_corrected[, common_samples]

# 3. 12サンプルのみでDESeq2を実行 ----
meta_wgcna$Diagnosis_2groups <- factor(meta_wgcna$Diagnosis_2groups)

dds_wgcna <- DESeqDataSetFromMatrix(
  countData = counts_wgcna,
  colData   = meta_wgcna,
  design    = ~ Diagnosis_2groups
)

dds_wgcna <- dds_wgcna[rowSums(DESeq2::counts(dds_wgcna)) > 1, ]
dds_wgcna <- DESeq(dds_wgcna)

# 4. VST 変換 ----
vsd     <- vst(dds_wgcna, blind = FALSE)
vst_mat <- assay(vsd)

# 5. MAD によるフィルタリング (上位 500 遺伝子) ----
# MADカットオフ = 500（旧プロジェクトで確認済み）
mad_cutoff <- 500

mad_vals  <- apply(vst_mat, 1, mad)
top_genes <- names(sort(mad_vals, decreasing = TRUE))[1:mad_cutoff]
datExpr   <- t(vst_mat[top_genes, ])

# サンプル・遺伝子の品質確認（旧プロジェクトで全サンプル・遺伝子が良好であることを確認済み）
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 6. Traits データの作成（メタボロミクス） ----
datTraits <- df_metab_all |>
  filter(Research_ID %in% rownames(datExpr)) |>
  dplyr::select(Research_ID, DOC) |>
  column_to_rownames("Research_ID")

datExpr <- datExpr[rownames(datTraits), ]

# 7. ソフト閾値の決定 ----
allowWGCNAThreads()

powers <- c(1:10, seq(12, 30, by = 2))
sft    <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# スケールフリー適合度 >= 0.9 となる最小のソフト閾値を自動選択
softPower <- min(sft$fitIndices$Power[sft$fitIndices$SFT.R.sq >= 0.9])

# ソフト閾値の確認プロット
pdf(file.path(out_dir_figs, "soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     labels = powers, cex = 0.7, col = "red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices$Power, sft$fitIndices$mean.k., labels = powers, cex = 0.7, col = "red")
par(mfrow = c(1, 1))
dev.off()

# 8. モジュール検出 (blockwiseModules) ----
set.seed(0)
net <- blockwiseModules(
  datExpr,
  power             = softPower,
  TOMType           = "unsigned",
  maxModuleSize     = mad_cutoff,
  reassignThreshold = 0,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = TRUE,
  saveTOMFileBase   = "data/processed/rnaseq/TOM",
  verbose           = 3
)

saveRDS(net,       "data/processed/rnaseq/wgcna_net.rds")
saveRDS(datExpr,   "data/processed/rnaseq/wgcna_datExpr.rds")
saveRDS(datTraits, "data/processed/rnaseq/wgcna_datTraits_metabolomics.rds")

# 9. クラスターデンドログラムの出力 ----
pdf(file.path(out_dir_figs, "cluster_dendrogram.pdf"), width = 12, height = 9)
plotDendroAndColors(
  net$dendrograms[[1]],
  labels2colors(net$colors)[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()

# モジュール遺伝子一覧の保存
module_df <- data.frame(
  gene   = names(net$colors),
  colors = labels2colors(net$colors)
)
write_csv(module_df, file.path(out_dir_tables, "module_genes.csv"))

# 各モジュールの遺伝子一覧を保存
for (mod_color in unique(module_df$colors)) {
  module_df |>
    filter(colors == mod_color) |>
    write_csv(file.path(out_dir_tables, paste0("module_", mod_color, ".csv")))
}

# 10. モジュールとメタボロミクスの相関ヒートマップ ----
nSamples <- nrow(datExpr)

moduleColors <- labels2colors(net$colors)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)

set.seed(0)
moduleTraitCor    <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

write_csv(as.data.frame(moduleTraitCor)    |> rownames_to_column("Module"), file.path(out_dir_tables, "module_trait_cor.csv"))
write_csv(as.data.frame(moduleTraitPvalue) |> rownames_to_column("Module"), file.path(out_dir_tables, "module_trait_pvalue.csv"))

textMatrix <- paste(
  round(moduleTraitCor, 2), "\n(p=",
  round(moduleTraitPvalue, 2), ")",
  sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)

pdf(file.path(out_dir_figs, "module_trait_heatmap_metabolomics.pdf"), width = 10, height = 6)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix        = moduleTraitCor,
  xLabels       = names(datTraits),
  yLabels       = names(MEs),
  ySymbols      = names(MEs),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = textMatrix,
  setStdMargins = FALSE,
  cex.text      = 1,
  zlim          = c(-1, 1),
  main          = "Module-Steroid Profile relationships"
)
dev.off()

# 11. モジュールのORA (enrichR) ----
dbs     <- c("GO_Biological_Process_2025", "MSigDB_Hallmark_2020", "NCI-Nature_2016")
modules <- unique(module_df$colors)

enrichr_results <- list()

for (mod_color in modules) {
  genes <- module_df |>
    filter(colors == mod_color) |>
    pull(gene)

  mod_results <- list()
  for (db in dbs) {
    tmp <- enrichr(genes = genes, databases = db)
    mod_results[[db]] <- tmp[[db]] |>
      filter(P.value < 0.05) |>
      as_tibble()
  }

  enrichr_results[[mod_color]] <- bind_rows(mod_results) |>
    mutate(
      Term = str_remove(Term, "\\s*\\(GO:\\d+\\)"),
      Term = str_remove(Term, "\\s*Homo sapiens.*$")
    ) |>
    arrange(P.value)

  enrichr_results[[mod_color]] |>
    dplyr::select(Term, P.value) |>
    write_csv(file.path(out_dir_tables, paste0("ora_enrichr_", mod_color, ".csv")))
}

# 12. ORA point and line プロット (turquoise モジュール) ----
selected_terms_list <- c(
  "Inflammatory Response",
  "KRAS Signaling Up",
  "Allograft Rejection",
  "Microglial Cell Activation",
  "Macrophage Activation",
  "Regulation of ERK1 and ERK2 Cascade",
  "Immune Response-Activating Cell Surface Receptor Signaling Pathway",
  "Complement",
  "Positive Regulation of Reactive Oxygen Species Metabolic Process",
  "Antigen Processing and Presentation of Exogenous Peptide Antigen"
)

plot_module <- "turquoise"

selected_terms <- enrichr_results[[plot_module]] |>
  filter(Term %in% selected_terms_list) |>
  group_by(Term) |>
  slice_min(P.value, n = 1) |>
  ungroup() |>
  mutate(
    logp = -log10(P.value),
    Term = forcats::fct_reorder(Term, logp)
  )

p_ora <- ggplot(selected_terms, aes(x = Term, y = logp, color = logp)) +
  geom_segment(aes(xend = Term, yend = 0), linewidth = 1) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_gradient(low = "#70b3ac", high = "#004d44") +
  geom_hline(yintercept = 1.3, linetype = "dotted", linewidth = 0.5, colour = "gray21") +
  coord_flip() +
  theme_bw() +
  labs(y = bquote(-Log[10](P-value))) +
  theme(
    text            = element_text(family = "Helvetica"),
    legend.position = "none",
    axis.title.y    = element_blank(),
    axis.title.x    = element_text(size = 8),
    axis.text.x     = element_text(size = 7),
    axis.text.y     = element_text(size = 7)
  )

ggsave(file.path(out_dir_figs, "ora_enrichr_turquoise.pdf"), p_ora, width = 18, height = 12, units = "cm")
ggsave(file.path(out_dir_figs, "ora_enrichr_turquoise.png"), p_ora, width = 18, height = 12, units = "cm", dpi = 300)
