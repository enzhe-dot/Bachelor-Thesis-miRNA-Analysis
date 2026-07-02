library(dplyr)
library(jsonlite)
library(tidyr)
library(ggplot2)

outdir <- "/on1/MAF/results2/fastp_summary"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE) 

base <- "/on1/MAF/results2"
datasets <- c("MAFlung", "MAFbreast", "MAFlaryng", "MAFhepatocellular")

results <- list()

for (ds in datasets) {
  
  pp_path <- file.path(base, ds, "02ppAout")
  star_path <- file.path(base, ds, "03star1out")
  
  samples <- list.dirs(pp_path, full.names = FALSE, recursive = FALSE)
  if (length(samples) == 0) next
  
  raw_lens <- c(); p1_lens <- c(); p2_lens <- c(); star_lens <- c()
  p1_retain <- c(); p2_retain <- c()
  raw_total <- c(); p2_total <- c()
  
  for (s in samples) {
    
    fp1 <- file.path(pp_path, s, "fp1.json")
    if (file.exists(fp1)) {
      d <- fromJSON(fp1)
      bf <- d$summary$before_filtering
      af <- d$summary$after_filtering
      raw_lens <- c(raw_lens, bf$read1_mean_length)
      p1_lens <- c(p1_lens, af$read1_mean_length)
      raw_total <- c(raw_total, bf$total_reads)
      if (bf$total_reads > 0) {
        p1_retain <- c(p1_retain, af$total_reads / bf$total_reads * 100)
      }
    }
    
    fp2 <- file.path(pp_path, s, "fp2.json")
    if (file.exists(fp2)) {
      d <- fromJSON(fp2)
      bf <- d$summary$before_filtering
      af <- d$summary$after_filtering
      p2_lens <- c(p2_lens, af$read1_mean_length)
      p2_total <- c(p2_total, af$total_reads)
      if (bf$total_reads > 0) {
        p2_retain <- c(p2_retain, af$total_reads / bf$total_reads * 100)
      }
    }
    
    slog <- file.path(star_path, s, "Log.final.out")
    if (file.exists(slog)) {
      log_content <- readLines(slog)
      line <- grep("Average input read length", log_content, value = TRUE)
      if (length(line) > 0) {
        star_lens <- c(star_lens, as.numeric(trimws(strsplit(line, "\t")[[1]][2])))
      }
    }
  }
  
  sm4 <- function(x) {
    if (length(x) == 0) return(c(NA, NA, NA, NA))
    c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
  }
  
  raw_stats   <- sm4(raw_lens)
  p1_stats    <- sm4(p1_lens)
  p2_stats    <- sm4(p2_lens)
  star_stats  <- sm4(star_lens)
  p1_ret      <- sm4(p1_retain)
  p2_ret      <- sm4(p2_retain)
  
  cum_retain <- if (length(raw_total) > 0 & length(p2_total) > 0) {
    round(sum(p2_total) / sum(raw_total) * 100, 1)
  } else { NA }
  
  results[[ds]] <- data.frame(
    dataset = ds,
    n_samples = length(samples),
    raw_mean_len  = round(raw_stats[1], 1),
    raw_sd_len    = round(raw_stats[2], 1),
    raw_min_len   = round(raw_stats[3], 0),
    raw_max_len   = round(raw_stats[4], 0),
    p1_mean_len   = round(p1_stats[1], 1),
    p1_sd_len     = round(p1_stats[2], 1),
    p1_retained_pct = round(p1_ret[1], 1),
    p2_mean_len   = round(p2_stats[1], 1),
    p2_sd_len     = round(p2_stats[2], 1),
    p2_retained_pct = round(p2_ret[1], 1),
    cum_retained_pct = cum_retain,
    star_mean_len = round(star_stats[1], 1),
    star_sd_len   = round(star_stats[2], 1)
  )
}

df <- bind_rows(results)
print(df)

write.csv(df, file.path(outdir, "pipeline_summary2.csv"), row.names = FALSE)

# Barplot
df_long <- df %>%
  select(dataset, raw_mean_len, p1_mean_len, p2_mean_len, star_mean_len) %>%
  pivot_longer(-dataset, names_to = "stage", values_to = "mean_length")

df_long$stage <- factor(df_long$stage,
                        levels = c("raw_mean_len", "p1_mean_len", "p2_mean_len", "star_mean_len"),
                        labels = c("Raw", "Pass 1", "Pass 2", "STAR"))

ggplot(df_long, aes(x = dataset, y = mean_length, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Raw" = "#E41A1C", "Pass 1" = "#377EB8",
                               "Pass 2" = "#4DAF4A", "STAR" = "#984EA3")) +
  theme_minimal() +
  labs(x = "Dataset", y = "Mean read length (nt)", fill = "Stage",
       title = "Read length through preprocessing pipeline")

ggsave(file.path(outdir, "read_length_pipeline2.png"),
       width = 10, height = 6)

message("Saved to: ", outdir)
