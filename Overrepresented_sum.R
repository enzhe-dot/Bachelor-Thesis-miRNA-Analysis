# Adjust 
breast_file <- "/on1/MAF/results2/MAFbreast/01qcAout/multiqc_breast/multiqc_data/fastqc_overrepresented_sequences_plot.txt"
lung_file   <- "/on1/MAF/results2/MAFlung/01qcAout/multiqc_lung/multiqc_data/fastqc_overrepresented_sequences_plot.txt"
laryng_file <- "/on1/MAF/results2/MAFlaryng/01qcAout/multiqc_laryng/multiqc_data/fastqc_overrepresented_sequences_plot.txt"
hepato_file <- "/on1/MAF/results2/MAFhepatocellular/01qcAout/multiqc_hepato/multiqc_data/fastqc_overrepresented_sequences_plot.txt"


breast <- read.delim(breast_file, header = TRUE)
lung   <- read.delim(lung_file,   header = TRUE)
laryng <- read.delim(laryng_file, header = TRUE)
hepato <- read.delim(hepato_file, header = TRUE)

# Summary Overrepresented per sample
breast$Total <- breast[,2] + breast[,3]
lung$Total   <- lung[,2]   + lung[,3]
laryng$Total <- laryng[,2] + laryng[,3]
hepato$Total <- hepato[,2] + hepato[,3]

# Table of mean valus
result <- data.frame(
  Dataset = c("Breast", "Lung", "Laryng", "Hepatocellular"),
  Samples = c(nrow(breast), nrow(lung), nrow(laryng), nrow(hepato)),
  Mean_Top_pct         = round(c(mean(breast[,2]), mean(lung[,2]), mean(laryng[,2]), mean(hepato[,2])), 2),
  Mean_Remaining_pct   = round(c(mean(breast[,3]), mean(lung[,3]), mean(laryng[,3]), mean(hepato[,3])), 2),
  Mean_Total_Overrep   = round(c(mean(breast$Total), mean(lung$Total), mean(laryng$Total), mean(hepato$Total)), 2),
  SD_Total_Overrep     = round(c(sd(breast$Total), sd(lung$Total), sd(laryng$Total), sd(hepato$Total)), 2)
)

print(result)


