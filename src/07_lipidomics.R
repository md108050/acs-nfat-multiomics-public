source("src/00_setup.R")

# 出力先
out_dir_figs   <- "results/figures/lipidomics"
out_dir_tables <- "results/tables/lipidomics"

dir.create(file.path(out_dir_figs, "logFC"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir_figs, "vip"),   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables,                   recursive = TRUE, showWarnings = FALSE)

# 1. データの読み込み ----
stats_files <- list.files("data/processed/lipidomics", pattern = "stats_.*\\.csv$", full.names = TRUE)

if (length(stats_files) == 0) stop("No stats files found in data/processed/lipidomics/")

stats_list <- list()
for (f in stats_files) {
  cat_name <- str_replace(basename(f), "stats_(.*)\\.csv", "\\1")
  df       <- read_csv(f, show_col_types = FALSE)

  if (!"FDR" %in% colnames(df) && "p_value" %in% colnames(df)) {
    df <- df |> mutate(FDR = p.adjust(p_value, method = "fdr"))
  }
  if ("SGoF" %in% colnames(df)) {
    df <- df |> mutate(significance = ifelse(SGoF < 0.05, "*", NA_character_))
  }

  if ("category" %in% colnames(df) && "class" %in% colnames(df)) {
    df <- df |> arrange(category, class, Lipid) |> mutate(Lipid = factor(Lipid, levels = unique(Lipid)))
  } else if ("class" %in% colnames(df)) {
    df <- df |> arrange(class, Lipid)           |> mutate(Lipid = factor(Lipid, levels = unique(Lipid)))
  } else if ("category" %in% colnames(df)) {
    df <- df |> arrange(category, Lipid)        |> mutate(Lipid = factor(Lipid, levels = unique(Lipid)))
  }

  stats_list[[cat_name]] <- df
}

# 2. 棒グラフの作成 ----
bar_theme <- theme_minimal() +
  theme(
    text         = element_text(family = "Helvetica"),
    axis.text    = element_text(size = 7),
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title   = element_text(size = 10),
    aspect.ratio = 0.75
  )

# --- Total ---
if ("Total" %in% names(stats_list)) {
  p_total <- ggplot(stats_list$Total, aes(x = Lipid, y = logFC, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = significance, vjust = ifelse(logFC > 0, -0.5, 1.5)), na.rm = TRUE) +
    scale_fill_manual(values = c("#D95F02", "#66A61E", "#1B9E77", "#7570B3", "#E7298A", "#e6ab02", "#a6761d")) +
    bar_theme +
    theme(legend.text = element_text(size = 7), legend.title = element_text(size = 10), plot.margin = margin(10, 10, 10, 10))

  ggsave(file.path(out_dir_figs, "logFC", "total.pdf"), p_total, width = 18, height = 14, units = "cm")
  ggsave(file.path(out_dir_figs, "logFC", "total.png"), p_total, width = 18, height = 14, units = "cm", dpi = 300, bg = "white")
}

# --- FattyAcyls ---
if ("FattyAcyls" %in% names(stats_list)) {
  p_fa <- ggplot(stats_list$FattyAcyls, aes(x = Lipid, y = logFC, fill = class)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = significance, vjust = ifelse(logFC > 0, -0.5, 1.5)), na.rm = TRUE) +
    scale_fill_manual(values = c("#D95F02", "#e67300", "#cc5900")) +
    bar_theme + theme(legend.position = "top")

  ggsave(file.path(out_dir_figs, "logFC", "FattyAcyls.pdf"), p_fa, width = 15, height = 12, units = "cm")
  ggsave(file.path(out_dir_figs, "logFC", "FattyAcyls.png"), p_fa, width = 15, height = 12, units = "cm", dpi = 300, bg = "white")
}

# --- Glycerolipids ---
if ("Glycerolipids" %in% names(stats_list)) {
  p_gl <- ggplot(stats_list$Glycerolipids, aes(x = Lipid, y = logFC, fill = class)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = significance, vjust = ifelse(logFC > 0, -0.5, 1.5)), na.rm = TRUE) +
    scale_fill_manual(values = c("#85d827", "#66A61E", "#477415", "#33550f")) +
    bar_theme + theme(legend.position = "top")

  ggsave(file.path(out_dir_figs, "logFC", "Glycerolipids.pdf"), p_gl, width = 15, height = 12, units = "cm")
  ggsave(file.path(out_dir_figs, "logFC", "Glycerolipids.png"), p_gl, width = 15, height = 12, units = "cm", dpi = 300, bg = "white")
}

# --- Glycerophospholipids ---
if ("Glycerophospholipids" %in% names(stats_list)) {
  format_lipid <- function(x) {
    purrr::map_chr(x, function(s) {
      if (!str_detect(s, "_")) return(s)
      head  <- str_replace(s, "_.*$", "")
      tail  <- str_replace(s, "^[^_]+_", "")
      nums  <- strsplit(tail, "_")[[1]]
      if (length(nums) %% 2 != 0) return(s)
      pairs <- purrr::map_chr(seq(1, length(nums), by = 2), ~ paste0("(", nums[.x], ":", nums[.x + 1], ")"))
      paste0(head, paste0(pairs, collapse = ""))
    })
  }

  p_gp <- ggplot(stats_list$Glycerophospholipids, aes(x = Lipid, y = logFC, fill = class)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = significance, vjust = ifelse(logFC > 0, -0.5, 1.5)), na.rm = TRUE) +
    scale_fill_manual(values = c("#24d09d", "#21c090", "#1eaf84", "#1b9e77", "#188d6a", "#157c5e", "#126c51", "#105b45", "#0d4a38")) +
    scale_x_discrete(labels = format_lipid) +
    bar_theme + theme(legend.position = "top")

  ggsave(file.path(out_dir_figs, "logFC", "Glycerophospholipids.pdf"), p_gp, width = 20, height = 15, units = "cm")
  ggsave(file.path(out_dir_figs, "logFC", "Glycerophospholipids.png"), p_gp, width = 20, height = 15, units = "cm", dpi = 300, bg = "white")
}

# --- Sphingolipids ---
if ("Sphingolipids" %in% names(stats_list)) {
  p_sl <- ggplot(stats_list$Sphingolipids, aes(x = Lipid, y = logFC, fill = class)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = significance, vjust = ifelse(logFC > 0, -0.5, 1.5)), na.rm = TRUE) +
    scale_fill_manual(values = c("#7570B3", "#5B5691", "#403C6E", "#2a2848")) +
    bar_theme + theme(legend.position = "top")

  ggsave(file.path(out_dir_figs, "logFC", "Sphingolipids.pdf"), p_sl, width = 15, height = 12, units = "cm")
  ggsave(file.path(out_dir_figs, "logFC", "Sphingolipids.png"), p_sl, width = 15, height = 12, units = "cm", dpi = 300, bg = "white")
}

# --- SterolLipids ---
if ("SterolLipids" %in% names(stats_list)) {
  p_st <- ggplot(stats_list$SterolLipids, aes(x = Lipid, y = logFC, fill = class)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = significance, vjust = ifelse(logFC > 0, -0.5, 1.5)), na.rm = TRUE) +
    scale_fill_manual(values = c("#E7298A", "#ae1462", "#740d41")) +
    bar_theme + theme(legend.position = "top")

  ggsave(file.path(out_dir_figs, "logFC", "SterolLipids.pdf"), p_st, width = 15, height = 12, units = "cm")
  ggsave(file.path(out_dir_figs, "logFC", "SterolLipids.png"), p_st, width = 15, height = 12, units = "cm", dpi = 300, bg = "white")
}

# 3. Supplementary Table の作成 ----
create_supp_table <- function(df, out_path) {
  if (!"p_value" %in% colnames(df)) return(NULL)

  df |>
    dplyr::select(Lipid, logFC, p_value, FDR) |>
    mutate(
      logFC     = round(logFC, 3),
      `P-value` = case_when(p_value < 0.001 ~ "<0.001", TRUE ~ as.character(round(p_value, 3))),
      FDR       = round(FDR, 3)
    ) |>
    dplyr::select(Lipid, logFC, `P-value`, FDR) |>
    write_csv(out_path)
}

for (cat_name in names(stats_list)) {
  out_path <- file.path(out_dir_tables, sprintf("stats_%s.csv", cat_name))

  if (cat_name == "Total") {
    stats_list[[cat_name]] |>
      arrange(category, class) |>
      dplyr::select(category, Lipid, logFC, p_value, FDR) |>
      mutate(
        logFC     = round(logFC, 3),
        `P-value` = case_when(p_value < 0.001 ~ "<0.001", TRUE ~ as.character(round(p_value, 3))),
        FDR       = round(FDR, 3)
      ) |>
      dplyr::select(Category = category, `Lipid class` = Lipid, LogFC = logFC, `P-value`, FDR) |>
      write_csv(out_path)
  } else if (cat_name != "All") {
    create_supp_table(stats_list[[cat_name]], out_path)
  }
}

# 4. VIP スコアプロット ----
for (cat_name in names(stats_list)) {
  df <- stats_list[[cat_name]]

  if ("VIP_score" %in% colnames(df)) {
    df_vip <- df |>
      filter(!is.na(VIP_score), VIP_score >= 1) |>
      mutate(Lipid = forcats::fct_reorder(Lipid, VIP_score))

    if (nrow(df_vip) > 0) {
      plot_height <- max(8, nrow(df_vip) * 0.8)

      p_vip <- ggplot(df_vip, aes(x = Lipid, y = VIP_score, color = Lipid)) +
        geom_segment(aes(xend = Lipid, yend = 0), linewidth = 1) +
        geom_point(size = 4, alpha = 0.7) +
        scale_color_viridis_d(option = "F", begin = 0.7, end = 0) +
        geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.75, colour = "gray21") +
        coord_flip() +
        theme_bw() +
        labs(y = "VIP score") +
        theme(
          legend.position = "none",
          axis.title.y    = element_blank(),
          axis.title.x    = element_text(size = 14),
          axis.text.x     = element_text(size = 14),
          axis.text.y     = element_text(size = 14)
        )

      ggsave(file.path(out_dir_figs, "vip", paste0(cat_name, ".pdf")), p_vip, width = 15, height = plot_height, units = "cm")
      ggsave(file.path(out_dir_figs, "vip", paste0(cat_name, ".png")), p_vip, width = 15, height = plot_height, units = "cm", dpi = 300, bg = "white")
    }
  }
}
