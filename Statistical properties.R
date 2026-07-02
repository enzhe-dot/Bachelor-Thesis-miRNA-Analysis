library(dplyr)
library(ggplot2)


outdir <- "C:/Users/ПК/OneDrive/Документы/Statistical proterties of dataset"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# === 1. reading 4 datasets  ===
files <- c(
  "C:\\Users\\ПК\\OneDrive\\Документы\\Statistical proterties of dataset\\fastqc_sequence_counts_plot_hepato.tsv" ,
  "C:\\Users\\ПК\\OneDrive\\Документы\\Statistical proterties of dataset\\fastqc_sequence_counts_plot_laryng.tsv",
  "C:\\Users\\ПК\\OneDrive\\Документы\\Statistical proterties of dataset\\fastqc_sequence_counts_plot_lung.tsv" ,
  "C:\\Users\\ПК\\OneDrive\\Документы\\Statistical proterties of dataset\\fastqc_sequence_counts_plot_breast.tsv"
)

names(files) <- c("Hepatocellular","Laryng", "Lung", "Breast" )

df <- bind_rows(lapply(names(files), function(ds) {
  read.table(files[ds], header = TRUE, sep = "\t") %>%
    mutate(dataset = ds)
}), .id = NULL)

# === 2. Total reads rename columns ===
colnames(df) <- gsub("\\.", "_", make.names(colnames(df)))
df$total_reads <- df$Unique_Reads + df$Duplicate_Reads

# === 3. Barplot ab all samples ===
ggplot(df, aes(x = Sample, y = total_reads, fill = dataset)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Sample", y = "Total Reads", fill = "Dataset",
       title = "Total sequencing reads per sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

ggsave(file.path(outdir, "total_reads_per_sample.png"),
       width = 14, height = 6)
# === 4. Boxplot on 4 datasets  ===
ggplot(df, aes(x = dataset, y = total_reads, fill = dataset)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  theme_minimal() +
  labs(x = "Dataset", y = "Total Reads per sample",
       title = "Distribution of total reads across datasets") +
  theme(legend.position = "none")

ggsave(file.path(outdir, "total_reads_by_dataset.png"),
       width = 8, height = 5)

# === 5. Mean for datasets ==
summary_stats <- df %>%
  group_by(dataset) %>%
  summarise(
    n_samples = n(),
    mean_total_reads = mean(total_reads) %>% round(0),
    sd_total_reads = sd(total_reads) %>% round(0),
    median_total_reads = median(total_reads) %>% round(0),
    min_total_reads = min(total_reads),
    max_total_reads = max(total_reads)
  )

print(summary_stats)
write.csv(summary_stats,
          file.path(outdir, "summary_reads_per_dataset.csv"),
          row.names = FALSE)
