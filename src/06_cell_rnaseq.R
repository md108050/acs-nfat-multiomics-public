source("src/00_setup.R")

# 出力先
out_dir_figs   <- "results/figures/cell_rnaseq"
out_dir_tables <- "results/tables/cell_rnaseq"

dir.create(out_dir_figs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# 1. データの読み込み ----
count_data <- read.csv("data/raw/rnaseq_cell/rawcounts.csv", header = TRUE, row.names = 1)

meta_data <- data.frame(
  condition = factor(c("CTL", "CTL", "CTL", "DOC", "DOC", "DOC"), levels = c("CTL", "DOC")),
  row.names = colnames(count_data)
)

# 2. DESeq2 による DEG 検出 ----
set.seed(0)

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData   = meta_data,
  design    = ~ condition
)

dds <- dds[rowSums(DESeq2::counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- DESeq2::results(dds)
res <- res[order(res$pvalue), ]

log2FC_threshold <- 1
FDR_threshold    <- 0.05

res_df <- drop_na(as_tibble(res, rownames = "gene"), padj)

deg      <- res_df |> filter(padj < FDR_threshold, abs(log2FoldChange) > log2FC_threshold)
deg_up   <- deg   |> filter(log2FoldChange > 0)
deg_down <- deg   |> filter(log2FoldChange < 0)

deg_up   |> arrange(padj) |> dplyr::select(gene, log2FoldChange, padj) |> write_csv(file.path(out_dir_tables, "deg_up.csv"))
deg_down |> arrange(padj) |> dplyr::select(gene, log2FoldChange, padj) |> write_csv(file.path(out_dir_tables, "deg_down.csv"))

# 3. Volcano Plot ----
gene_of_interest <- c("Fabp4", "Abca6", "Mertk", "C1qa", "Jun", "Dusp1", "Fkbp5", "Tsc22d3")

res_plot <- res_df |>
  mutate(
    negLogP = -log10(padj),
    sig_cat = case_when(
      gene %in% gene_of_interest                                 ~ "highlight",
      log2FoldChange < -log2FC_threshold & padj < FDR_threshold ~ "sig",
      log2FoldChange >  log2FC_threshold & padj < FDR_threshold ~ "sig",
      TRUE                                                       ~ "nonsig"
    )
  )

p_volcano <- ggplot(res_plot, aes(x = log2FoldChange, y = negLogP)) +
  geom_point(aes(color = sig_cat), alpha = 0.6, size = 2.5) +
  scale_color_manual(values = c("highlight" = "red", "sig" = "black", "nonsig" = "grey")) +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed", color = "black") +
  geom_text_repel(
    data          = filter(res_plot, gene %in% gene_of_interest),
    aes(label     = gene),
    size          = 4,
    box.padding   = 0.5,
    point.padding = 0.5,
    max.overlaps  = 50
  ) +
  coord_fixed(ratio = 0.04) +
  labs(x = bquote(Log[2](Fold~Change)), y = bquote(-Log[10](FDR))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text       = element_text(size = 8),
    axis.title      = element_text(size = 10)
  )

ggsave(file.path(out_dir_figs, "volcano.pdf"), p_volcano, width = 15, height = 15, units = "cm")
ggsave(file.path(out_dir_figs, "volcano.png"), p_volcano, width = 15, height = 15, units = "cm", dpi = 300, bg = "white")

# 4. ORA (clusterProfiler) ----
deg_up_entrez   <- bitr(deg_up$gene,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) |> distinct(SYMBOL, .keep_all = TRUE)
deg_down_entrez <- bitr(deg_down$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) |> distinct(SYMBOL, .keep_all = TRUE)

hallmark_mouse <- msigdbr(species = "Mus musculus", category = "H") |>
  dplyr::select(gs_name, entrez_gene)

run_ora <- function(entrez_ids) {
  if (length(entrez_ids) == 0) return(NULL)

  gobp <- enrichGO(
    entrez_ids,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pvalueCutoff  = 0.1,
    pAdjustMethod = "fdr",
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )

  hallmark <- enricher(
    entrez_ids,
    TERM2GENE     = hallmark_mouse,
    pAdjustMethod = "fdr",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.2
  )

  bind_rows(
    if (!is.null(gobp))     as_tibble(gobp)     else tibble(),
    if (!is.null(hallmark)) as_tibble(hallmark) else tibble()
  ) |> arrange(pvalue)
}

ora_up   <- run_ora(deg_up_entrez$ENTREZID)
ora_down <- run_ora(deg_down_entrez$ENTREZID)

if (!is.null(ora_up)   && nrow(ora_up)   > 0) write_csv(ora_up   |> dplyr::select(Description, pvalue), file.path(out_dir_tables, "ora_up.csv"))
if (!is.null(ora_down) && nrow(ora_down) > 0) write_csv(ora_down |> dplyr::select(Description, pvalue), file.path(out_dir_tables, "ora_down.csv"))

# 5. ORA Point and Line Plot ----
selected_terms_list <- c(
  "positive regulation of inflammatory response",
  "microglial cell activation",
  "regulation of myeloid cell differentiation",
  "regulation of p38MAPK cascade",
  "fat cell differentiation",
  "positive regulation of fat cell differentiation",
  "tumor necrosis factor production",
  "complement activation",
  "monocyte chemotaxis",
  "response to steroid hormone"
)

if (!is.null(ora_up) && nrow(ora_up) > 0) {
  selected_terms <- ora_up |>
    filter(Description %in% selected_terms_list) |>
    mutate(
      logp = -log10(pvalue),
      Term = fct_reorder(Description, logp)
    ) |>
    filter(!is.na(Term))

  if (nrow(selected_terms) > 0) {
    p_ora <- ggplot(selected_terms, aes(x = Term, y = logp, color = Term)) +
      geom_segment(aes(xend = Term, yend = 0), linewidth = 1) +
      geom_point(size = 4, alpha = 0.7) +
      scale_color_viridis_d(option = "rocket", begin = 0.7, end = 0) +
      geom_hline(yintercept = 1.3, linetype = "dotted", linewidth = 0.5, colour = "gray21") +
      coord_flip() +
      theme_bw() +
      labs(y = bquote(-Log[10](P-value))) +
      theme(
        legend.position = "none",
        axis.title.y    = element_blank(),
        axis.title.x    = element_text(size = 10),
        axis.text.x     = element_text(size = 8),
        axis.text.y     = element_text(size = 8)
      )

    ggsave(file.path(out_dir_figs, "ora_up.pdf"), p_ora, width = 18, height = 12, units = "cm")
    ggsave(file.path(out_dir_figs, "ora_up.png"), p_ora, width = 18, height = 12, units = "cm", dpi = 300, bg = "white")
  }
}
