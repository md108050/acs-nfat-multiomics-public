source("src/00_setup.R")

# 出力先
out_dir_figs   <- "results/figures/deconvolution"
out_dir_tables <- "results/tables/deconvolution"

dir.create(out_dir_figs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# 1. シングルセルリファレンスデータの読み込み ----
# データ取得元: https://singlecell.broadinstitute.org/single_cell/study/SCP1376
scdata <- readRDS("data/raw/scrnaseq_reference/human_all_lite.rds")

# 2. 細胞型ラベルの再定義 (Tier 3) ----
Idents(scdata) <- scdata$cell_type

new_labels <- case_when(
  Idents(scdata) %in% paste0("hAd",   1:7) ~ "adipocyte",
  Idents(scdata) %in% paste0("hASPC", 1:6) ~ "ASPC",
  Idents(scdata) == "hBcell"               ~ "b_cell",
  Idents(scdata) %in% c("hASDC", "hcDC1", "hcDC2", "hpDC") ~ "dendritic_cell",
  Idents(scdata) %in% c("hEndoA1", "hEndoA2", "hEndoS1", "hEndoS2", "hEndoS3", "hEndoV") ~ "endothelial",
  Idents(scdata) %in% c("hLEC1", "hLEC2")  ~ "LEC",
  Idents(scdata) == "hMac1"                ~ "hMac1",
  Idents(scdata) == "hMac2"                ~ "hMac2",
  Idents(scdata) == "hMac3"                ~ "hMac3",
  Idents(scdata) == "hMast"                ~ "mast_cell",
  Idents(scdata) %in% c("hMes1", "hMes2", "hMes3") ~ "mesothelium",
  Idents(scdata) %in% c("hMono1", "hMono2") ~ "monocyte",
  Idents(scdata) == "hNeu"                 ~ "neutrophil",
  Idents(scdata) == "hNK"                  ~ "nk_cell",
  Idents(scdata) == "hPeri"                ~ "pericyte",
  Idents(scdata) %in% c("hSMC1", "hSMC2") ~ "SMC",
  Idents(scdata) %in% c("hTcell1", "hTcell2", "hTreg") ~ "t_cell",
  TRUE ~ NA_character_
)

# NA（hEndM, hPlasmablast）を除外
keep      <- !is.na(new_labels)
scdata_t3 <- scdata[, keep]
scdata_t3$cell_type <- new_labels[keep]

# 3. VATのダウンサンプリング (最大500細胞/細胞型) ----
cell_meta     <- scdata_t3@meta.data
cell_meta$cell_id <- rownames(cell_meta)

set.seed(123)
downsampled_cells <- cell_meta |>
  filter(depot == "VAT") |>
  group_by(cell_type) |>
  sample_n(size = min(500, n()), replace = FALSE) |>
  pull(cell_id)

scdata_subset <- subset(scdata_t3, cells = downsampled_cells)

# 4. バルクRNA-seqデータの読み込み ----
meta   <- read.csv("data/meta/rnaseq_meta.csv",              row.names = 1)
counts <- read.csv("data/raw/rnaseq_human/rnaseq_count.csv", row.names = 1, check.names = FALSE)

counts <- counts[, rownames(meta)]

# 5. BisqueRNA によるデコンボリューション ----
bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(counts))

Idents(scdata_subset) <- scdata_subset$cell_type
sample.ids <- colnames(scdata_subset)
sc.meta    <- scdata_subset@meta.data

sc.pheno <- data.frame(
  check.names    = FALSE,
  check.rows     = FALSE,
  stringsAsFactors = FALSE,
  row.names      = sample.ids,
  SubjectName    = sc.meta$sample,
  cellType       = sc.meta$cell_type
)

sc.pdata <- new(
  "AnnotatedDataFrame",
  data        = sc.pheno,
  varMetadata = data.frame(
    labelDescription = c("SubjectName", "cellType"),
    row.names        = c("SubjectName", "cellType")
  )
)

sc.eset <- Biobase::ExpressionSet(
  assayData = as.matrix(scdata_subset@assays$RNA@counts),
  phenoData = sc.pdata
)

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers = NULL, use.overlap = FALSE)

ref.based.estimates <- res$bulk.props
res_df  <- data.frame(t(ref.based.estimates))

join_meta <- cbind(meta, res_df)

# 6. 骨髄系細胞スコアの算出 ----
myeloid_cells        <- c("hMac1", "hMac2", "hMac3", "monocyte", "neutrophil", "dendritic_cell", "mast_cell")
join_meta$myeloid_total <- rowSums(join_meta[, myeloid_cells], na.rm = TRUE)

# 7. Table 1 (細胞型割合) ----
listvars <- c(
  "adipocyte", "ASPC", "b_cell", "dendritic_cell", "endothelial",
  "hMac1", "hMac2", "hMac3", "LEC", "mast_cell", "mesothelium",
  "monocyte", "neutrophil", "nk_cell", "pericyte", "SMC", "t_cell",
  "myeloid_total"
)

tab1 <- CreateTableOne(
  vars          = listvars,
  data          = join_meta,
  strata        = "Diagnosis_2groups",
  includeNA     = FALSE,
  addOverall    = FALSE,
  testExact     = fisher.test,
  testNonNormal = kruskal.test
)

res_tab <- print(
  tab1,
  showAllLevels = TRUE,
  nonnormal     = listvars,
  smd           = FALSE,
  explain       = TRUE,
  catDigits     = 3,
  contDigits    = 3,
  pDigits       = 3,
  test          = TRUE,
  format        = "fp"
)

write_csv(as.data.frame(res_tab), file.path(out_dir_tables, "deconvolution_table1.csv"))

# 8. Box plot (hMac1 & Monocyte) ----
plot_cells_main <- c("hMac1", "monocyte")
cell_label_main <- c(
  hMac1    = "Macrophage subcluster hMac1",
  monocyte = "Monocyte"
)

df_box <- join_meta |>
  dplyr::select(Diagnosis_2groups, all_of(plot_cells_main)) |>
  pivot_longer(
    cols      = all_of(plot_cells_main),
    names_to  = "cell_type",
    values_to = "proportion"
  ) |>
  mutate(
    cell_group     = factor(cell_label_main[cell_type], levels = cell_label_main[plot_cells_main]),
    Diagnosis_2groups = factor(Diagnosis_2groups, levels = c("NFAT", "ACS"))
  )

p_box <- ggplot(df_box, aes(x = cell_group, y = proportion, fill = Diagnosis_2groups)) +
  geom_boxplot(
    alpha        = 0.8,
    outlier.shape = NA,
    width        = 0.6,
    position     = position_dodge(0.75)
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
    size     = 1.2,
    alpha    = 0.4
  ) +
  scale_fill_manual(values = c("NFAT" = "#999999", "ACS" = "#E41A1C")) +
  labs(x = "", y = "Estimated proportion", fill = "") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(color = "black", size = 12),
    axis.text.y     = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir_figs, "boxplot_hMac1_monocyte.pdf"), p_box, width = 6, height = 6)
ggsave(file.path(out_dir_figs, "boxplot_hMac1_monocyte.png"), p_box, width = 6, height = 6, dpi = 300)

# 9. Median [IQR] plot (全17細胞型) ----
cell_order <- c(
  "adipocyte", "ASPC", "mesothelium",
  "endothelial", "LEC", "pericyte", "SMC",
  "monocyte", "hMac1", "hMac2", "hMac3",
  "dendritic_cell", "mast_cell", "neutrophil",
  "t_cell", "b_cell", "nk_cell"
)

cell_labels <- c(
  adipocyte      = "Adipocyte",
  ASPC           = "ASPC",
  mesothelium    = "Mesothelium",
  endothelial    = "Endothelial",
  LEC            = "LEC",
  pericyte       = "Pericyte",
  SMC            = "SMC",
  monocyte       = "Monocyte",
  hMac1          = "Macrophage subcluster hMac1",
  hMac2          = "Macrophage subcluster hMac2",
  hMac3          = "Macrophage subcluster hMac3",
  dendritic_cell = "Dendritic cell",
  mast_cell      = "Mast cell",
  neutrophil     = "Neutrophil",
  t_cell         = "T cell",
  b_cell         = "B cell",
  nk_cell        = "NK cell"
)

df_all <- join_meta |>
  dplyr::select(Diagnosis_2groups, all_of(cell_order)) |>
  pivot_longer(
    cols      = all_of(cell_order),
    names_to  = "cell_type",
    values_to = "proportion"
  ) |>
  mutate(
    cell_group    = factor(cell_labels[cell_type], levels = rev(cell_labels[cell_order])),
    Diagnosis_2groups = factor(Diagnosis_2groups, levels = c("NFAT", "ACS"))
  )

summary_df <- df_all |>
  group_by(cell_group, Diagnosis_2groups) |>
  summarise(
    median = median(proportion, na.rm = TRUE),
    q1     = quantile(proportion, 0.25, na.rm = TRUE),
    q3     = quantile(proportion, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p_iqr <- ggplot(summary_df, aes(x = median, y = cell_group, color = Diagnosis_2groups)) +
  geom_errorbar(
    aes(xmin = q1, xmax = q3),
    orientation = "y",
    width    = 0.25,
    position = position_dodge(width = 0.6),
    linewidth = 0.7
  ) +
  geom_point(position = position_dodge(width = 0.6), size = 2.3) +
  scale_color_manual(values = c("NFAT" = "#999999", "ACS" = "#E41A1C")) +
  labs(x = "Estimated proportion, median [IQR]", y = "", color = "") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(color = "black"),
    axis.text.y      = element_text(color = "black", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir_figs, "median_iqr_17celltypes.pdf"), p_iqr, width = 7, height = 6.5)
ggsave(file.path(out_dir_figs, "median_iqr_17celltypes.png"), p_iqr, width = 7, height = 6.5, dpi = 300)
