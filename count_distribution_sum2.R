# miRNA-seq QC plots – per-dataset, free Y-axis scales
# ──────────────────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(patchwork)

# ── 1.  Paths ──────────────────────────────────────────────────────────────
tsv_files <- c(
  breast         = "/on1/MAF/results2/expr_data3/br3.mt_stm_miR.tsv",
  hepatocellular = "/on1/MAF/results2/expr_data3/h3.mt_stm_miR.tsv",
  laryngeal      = "/on1/MAF/results2/expr_data3/lar3.mt_stm_miR.tsv",
  lung           = "/on1/MAF/results2/expr_data3/l3.mt_stm_miR.tsv"
)

meta_files <- c(
  breast         = "/on1/MAF/results2/output_dir_breast/Breast_meta.csv",
  hepatocellular = "/on1/MAF/results2/output_dir_hepato/Hepatocellular_meta.csv",
  laryngeal      = "/on1/MAF/results2/output_dir_laryng/Laryng_meta.csv",
  lung           = "/on1/MAF/results2/output_dir_lung/Lung_meta.csv"
)

OUTPUT_DIR <- "/on1/MAF/results2/count_distribution2"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 2.  Helper: load one dataset ──────────────────────────────────────────
process_dataset <- function(tsv_path, meta_path, dataset_label) {
  
  counts_raw <- read.delim(tsv_path, check.names = FALSE, row.names = 1)
  counts_mat <- as.matrix(counts_raw)
  mode(counts_mat) <- "numeric"
  
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  colnames(meta) <- tolower(trimws(colnames(meta)))
  stopifnot("sample" %in% colnames(meta), "condition" %in% colnames(meta))
  meta$sample    <- trimws(meta$sample)
  meta$condition <- trimws(meta$condition)
  colnames(counts_mat) <- trimws(colnames(counts_mat))
  
  common_samps <- intersect(meta$sample, colnames(counts_mat))
  
  missing_in_counts <- setdiff(meta$sample, colnames(counts_mat))
  if (length(missing_in_counts) > 0)
    warning(dataset_label, ": samples in meta NOT in TSV: ",
            paste(missing_in_counts, collapse = ", "))
  
  counts_sub <- counts_mat[, common_samps, drop = FALSE]
  meta_sub   <- meta[match(common_samps, meta$sample), , drop = FALSE]
  
  data.frame(
    Sample      = common_samps,
    Dataset     = dataset_label,
    LibrarySize = colSums(counts_sub, na.rm = TRUE),
    Detected    = colSums(counts_sub > 0, na.rm = TRUE),
    Condition   = meta_sub$condition,
    stringsAsFactors = FALSE
  )
}

# ── 3.  Helper: smart Y-axis label ────────────────────────────────────────
# Automatically picks "Millions" / "Thousands" / raw counts
# depending on the actual range of the data
lib_axis <- function(lib_sizes) {
  mx <- max(lib_sizes, na.rm = TRUE)
  if (mx >= 1e6)  list(divisor = 1e6, label = "Total Reads (Millions)")
  else if (mx >= 1e3) list(divisor = 1e3, label = "Total Reads (Thousands)")
  else            list(divisor = 1,   label = "Total Reads (counts)")
}

# ── 4.  Helper: build A+B patchwork for one dataset ───────────────────────
plot_one_dataset <- function(df) {
  
  ds_label <- unique(df$Dataset)
  
  # sort: Normal first, then Cancer; preserve that order on x-axis
  df <- df %>%
    mutate(Condition = factor(Condition, levels = c("Normal", "Cancer"))) %>%
    arrange(Condition) %>%
    mutate(Sample = factor(Sample, levels = unique(Sample)))
  
  fill_scale <- scale_fill_manual(
    values = c(Normal = "#43A047", Cancer = "#E53935"),
    drop   = FALSE
  )
  
  base_theme <- theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      strip.background   = element_blank(),
      axis.text.x        = element_text(angle = 90, vjust = 0.5,
                                        hjust = 1, size = 7),
      axis.ticks.x       = element_blank(),
      legend.position    = "top",
      panel.grid.major.x = element_blank()
    )
  
  # ---- Panel A: Library size with smart unit ----
  ax <- lib_axis(df$LibrarySize)
  
  pA <- ggplot(df, aes(x = Sample, y = LibrarySize / ax$divisor,
                       fill = Condition)) +
    geom_col(width = 0.75, colour = "grey20", linewidth = 0.15) +
    fill_scale +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste0("A  Library size — ", ds_label),
      x     = NULL,
      y     = ax$label,
      fill  = "Condition"
    ) +
    base_theme
  
  # ---- Panel B: Detected miRNAs ----
  pB <- ggplot(df, aes(x = Sample, y = Detected, fill = Condition)) +
    geom_col(width = 0.75, colour = "grey20", linewidth = 0.15) +
    fill_scale +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste0("B  Detected miRNAs — ", ds_label),
      x     = NULL,
      y     = "Number of detected miRNAs",
      fill  = "Condition"
    ) +
    base_theme
  
  pA / pB   # patchwork: A on top, B below
}

# ── 5.  Run & save ────────────────────────────────────────────────────────
all_qc <- list()

for (ds in names(tsv_files)) {
  message("Processing ", ds, " ...")
  
  df <- process_dataset(tsv_files[ds], meta_files[ds], ds)
  all_qc[[ds]] <- df
  
  # scale figure width to number of samples (min 8 inches)
  fig_w <- max(8, ncol(read.delim(tsv_files[ds], row.names = 1)) * 0.18)
  
  p <- plot_one_dataset(df)
  
  out_pdf <- file.path(OUTPUT_DIR,
                       paste0("QC_", ds, ".pdf"))
  ggsave(out_pdf, p, width = fig_w, height = 10)
  message("  Saved → ", out_pdf)
}

# ── 6.  Save combined QC table ────────────────────────────────────────────
qc_all <- do.call(rbind, all_qc)
write.csv(qc_all,
          file.path(OUTPUT_DIR, "qc_summary_table.csv"),
          row.names = FALSE)
# ──optional one big pdf for appendix 
all_plots <- lapply(names(all_qc), function(ds) plot_one_dataset(all_qc[[ds]]))

# merge one on another 
combined_appendix <- wrap_plots(all_plots, ncol = 1)

ggsave(
  file.path(OUTPUT_DIR, "QC_all_datasets_appendix.pdf"),
  combined_appendix,
  width  = 16,
  height = 10 * length(all_plots)   
)
message("\n✓ Done. Files saved to: ", normalizePath(OUTPUT_DIR))

