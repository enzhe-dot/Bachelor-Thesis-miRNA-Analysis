library(dplyr)
library(tidyr)
library(ggplot2)

outdir <- "C:/Users/ПК/OneDrive/Документы/Statistical proterties of dataset"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# adjust 
files <- c(
  "C:\\Users\\ПК\\Downloads\\fastqc_sequence_duplication_levels_plot_hepato.tsv",
  "C:\\Users\\ПК\\Downloads\\fastqc_sequence_duplication_levels_plot_laryng.tsv",
  "C:\\Users\\ПК\\Downloads\\fastqc_sequence_duplication_levels_plot_lung.tsv",
  "C:\\Users\\ПК\\Downloads\\fastqc_sequence_duplication_levels_plot_breast.tsv"
)

names(files) <- c("Lung", "Breast", "Laryng", "Hepatocellular")


all_data <- data.frame()
for (ds in names(files)) {
  if (!file.exists(files[ds])) { cat("NOT FOUND:", files[ds], "\n"); next }
  df <- read.table(files[ds], header = TRUE, sep = "\t", check.names = FALSE)
  df$dataset <- ds
  all_data <- rbind(all_data, df)
}

# Long format: Sample | duplication_level | percentage | dataset
df_long <- all_data %>%
  pivot_longer(-c(Sample, dataset), names_to = "duplication_level", values_to = "percentage")

# Order on X axes
df_long$duplication_level <- factor(df_long$duplication_level,
                                    levels = c("1", "2", "3", "4", "5", ">5", ">10", ">50"))

# mean across datasets  
df_summary <- df_long %>%
  group_by(dataset, duplication_level) %>%
  summarise(mean_pct = mean(percentage, na.rm = TRUE), .groups = "drop")

# === Stacked barplot ===
ggplot(df_summary, aes(x = dataset, y = mean_pct, fill = duplication_level)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  theme_minimal() +
  labs(x = "Dataset", y = "Mean % of Reads", fill = "Duplication Level",
       title = "Duplication Level Distribution by Dataset")

ggsave(file.path(outdir, "duplication_levels_stacked.png"), width = 8, height = 5)

# === Heatmap ===
ggplot(df_summary, aes(x = duplication_level, y = dataset, fill = mean_pct)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#E41A1C") +
  geom_text(aes(label = round(mean_pct, 1)), size = 3.5) +
  theme_minimal() +
  labs(x = "Duplication Level", y = "Dataset", fill = "Mean %",
       title = "% of Reads per Duplication Level")

ggsave(file.path(outdir, "duplication_levels_heatmap.png"), width = 8, height = 4)

cat("Saved to:", outdir, "\n")
