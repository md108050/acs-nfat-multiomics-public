source("src/00_setup.R")

# 出力先
out_dir_figs   <- "results/figures/rnaseq/enrichr"
out_dir_tables <- "results/tables/rnaseq"

dir.create(out_dir_figs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# 1. DEGリストの読み込み ----
deg_up   <- read_csv(file.path(out_dir_tables, "DEG_up.csv"),   show_col_types = FALSE)$gene_symbol
deg_down <- read_csv(file.path(out_dir_tables, "DEG_down.csv"), show_col_types = FALSE)$gene_symbol

genes_list <- list(deg_up = deg_up, deg_down = deg_down)

# 2. enrichR の実行 ----
dbs <- c("GO_Biological_Process_2025", "MSigDB_Hallmark_2020", "NCI-Nature_2016")

enrichr_results <- list()

for (deg_name in names(genes_list)) {
  genes <- genes_list[[deg_name]]

  if (length(genes) > 0) {
    tmp_res <- enrichr(genes = genes, databases = dbs)

    combined_res <- bind_rows(
      tmp_res[[dbs[1]]],
      tmp_res[[dbs[2]]],
      tmp_res[[dbs[3]]]
    ) |>
      filter(P.value < 0.05) |>
      arrange(P.value) |>
      mutate(
        Term_Clean = str_remove(Term, "\\s*\\(GO:\\d+\\)"),
        Term_Clean = str_remove(Term_Clean, "\\s*Homo sapiens.*$")
      )

    enrichr_results[[deg_name]] <- combined_res

    write_csv(
      combined_res |> dplyr::select(Term, P.value, Adjusted.P.value, Overlap, Genes),
      file.path(out_dir_tables, paste0("enrichr_", deg_name, ".csv"))
    )
  }
}

# 3. Enrichment Plot ----
selected_terms <- c(
  "TNF-alpha Signaling via NF-kB",
  "UV Response Up",
  "Regulation of p38MAPK Cascade",
  "Hypoxia",
  "Epithelial Mesenchymal Transition",
  "IL-2/STAT5 Signaling",
  "Ceramide signaling pathway",
  "Regulation of Phagocytosis"
)

plot_df <- enrichr_results[["deg_up"]] |>
  filter(Term_Clean %in% selected_terms) |>
  group_by(Term_Clean) |>
  slice_min(P.value, n = 1) |>
  ungroup() |>
  mutate(logp = -log10(P.value)) |>
  arrange(logp) |>
  mutate(Term_Clean = factor(Term_Clean, levels = Term_Clean))

p_enrichr <- ggplot(plot_df, aes(x = logp, y = Term_Clean, color = Term_Clean)) +
  geom_segment(aes(xend = 0, yend = Term_Clean), linewidth = 1) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_viridis_d(option = "rocket", begin = 0.1, end = 0.8, direction = -1) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray20", linewidth = 0.5) +
  theme_bw() +
  labs(
    x = bquote(-Log[10](P-value)),
    y = ""
  ) +
  theme(
    text            = element_text(family = "Helvetica", size = 12),
    legend.position = "none"
  )

ggsave(file.path(out_dir_figs, "enrichr.pdf"), p_enrichr, width = 6, height = 4)
ggsave(file.path(out_dir_figs, "enrichr.png"), p_enrichr, width = 6, height = 4, dpi = 300)
