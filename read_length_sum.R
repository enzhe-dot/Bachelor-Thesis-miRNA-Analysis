# ============================================================
# Line plot: Read length distribution (fixed — reads length from tuple)
# ============================================================

file_laryng <- "/on1/MAF/results2/MAFlaryng/01qcAout/multiqc_laryng/multiqc_data/laryng_fastqc_sequence_length_distribution_plot.txt"
file_hepato  <-"/on1/MAF/results2/MAFhepatocellular/01qcAout/multiqc_hepato/multiqc_data/hepato_fastqc_sequence_length_distribution_plot.txt"
outdir <- "/on1/MAF/results2/fastp_summary"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---- Parse: reading length from each tool  ----
parse_robust <- function(filepath, dataset_name) {
  if (!file.exists(filepath)) return(NULL)
  
  lines <- readLines(filepath)
  lines <- lines[-1]
  
  
  all_lengths <- c()
  count_list <- list()
  
  for (i in seq_along(lines)) {
    parts <- strsplit(lines[i], "\t")[[1]]
    n_bins <- length(parts) - 1
    sample_lengths <- numeric(n_bins)
    sample_counts <- numeric(n_bins)
    
    for (k in 1:n_bins) {
      t <- parts[k + 1]
      clean <- gsub("[^0-9.,]", "", t)
      vals <- as.numeric(strsplit(clean, ",")[[1]])
      sample_lengths[k] <- vals[1]
      sample_counts[k] <- vals[2]
    }
    
    # Сохраняем длины из первого сэмпла (они должны быть одинаковы у всех)
    if (i == 1) all_lengths <- sample_lengths
    count_list[[i]] <- sample_counts
  }
  
  cat(dataset_name, ":", length(lines), "samples,", length(all_lengths),
      "bins (", min(all_lengths), "to", max(all_lengths), "bp)\n")
  
  # Матрица процентов
  count_mat <- do.call(rbind, count_list)
  pct_mat <- count_mat
  for (i in 1:nrow(pct_mat)) {
    total <- sum(pct_mat[i, ], na.rm = TRUE)
    if (total > 0) pct_mat[i, ] <- (pct_mat[i, ] / total) * 100
  }
  
  data.frame(length = all_lengths,
             pct = colMeans(pct_mat, na.rm = TRUE),
             dataset = dataset_name,
             stringsAsFactors = FALSE)
}

df_laryng <- parse_robust(file_laryng, "Laryng")
df_hepato <- parse_robust(file_hepato, "Hepatocellular")

# ---- Lung и Breast ----
all_lengths <- sort(unique(c(df_laryng$length, df_hepato$length)))

all_lengths <- all_lengths[all_lengths <= 150]

constant <- data.frame(length = all_lengths, pct = 0, stringsAsFactors = FALSE)
df_lung   <- constant; df_lung$pct[constant$length == 50] <- 100; df_lung$dataset <- "Lung"
df_breast <- constant; df_breast$pct[constant$length == 50] <- 100; df_breast$dataset <- "Breast"

df_all <- rbind(df_laryng, df_hepato, df_lung, df_breast)

# ---- Graph ----
cols <- c("Lung" = "#7B2D8E", "Breast" = "#D7263D",
          "Laryng" = "#2E8B57", "Hepatocellular" = "#1A8FBF")

png(file.path(outdir, "read_length_distribution_lineplot.png"),
    width = 3000, height = 1800, res = 300)

par(mar = c(5, 5, 4, 2))
plot(0, 0, type = "n", xlim = c(0, 150), ylim = c(0, 105),
     xlab = "Read Length (bp)", ylab = "% of Reads",
     main = "Sequence Length Distribution across Datasets", las = 1)

grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

for (ds in unique(df_all$dataset)) {
  sub <- subset(df_all, dataset == ds)
  if (ds %in% c("Lung", "Breast")) {
    segments(50, 0, 50, 100, col = cols[ds], lwd = 4)
  } else {
    lines(sub$length, sub$pct, col = cols[ds], lwd = 2.5)
  }
}

legend("topright", legend = names(cols), col = cols, lwd = 3, bty = "n", cex = 0.8)

dev.off()
cat("Saved to:", file.path(outdir, "read_length_distribution_lineplot.png"), "\n")
