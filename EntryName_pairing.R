# ============================================================
# Filter DESeq2-результатов with KEGG miRNA
# Format of the KEGG-file: 
# generated manually 
#   00302238          MIR103B1; hsa-miR-103b [KO:K16989]
#   ↑ Entry           ↑ Symbol  ↑ miRNA_name
# ============================================================

library(dplyr)
library(tidyr)


kegg_file <- "C:\\Users\\ПК\\OneDrive\\Документы\\KEGG_05206_207.txt"
deseq_files <- c(
  Breast = "C:\\Users\\ПК\\Downloads\\DESeq2_results_breast (2).csv",
  Larynx = "C:\\Users\\ПК\\Downloads\\DESeq2_results_laryng.csv",
  Lung   = "C:\\Users\\ПК\\Downloads\\DESeq2_results_lung .csv",
  Hepatocellular = "C:\\Users\\ПК\\Downloads\\DESeq2_results_hepato.csv"
)
output_dir <- "C:\\Users\\ПК\\Downloads"
output_names <- c(
  Breast = file.path(output_dir, "DESeq2_breast_05206_miRNA.csv"),
  Larynx = file.path(output_dir, "DESeq2_laryng_05206_miRNA.csv"),
  Lung   = file.path(output_dir, "DESeq2_lung_05206_miRNA.csv"),
  Hepatocellular = file.path(output_dir, "DESeq2_hepato_05206_miRNA.csv")
)

# ============================================================
# Pairing EGG miRNA .txt
# String: 00302238          MIR103B1; hsa-miR-103b [KO:K16989]
# ============================================================
parse_kegg_mirna <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  
  lines <- lines[lines != "" & !grepl("^\\s*$", lines)]
  
  parse_line <- function(line) {
   
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    entry <- parts[1]
    
    
    rest <- paste(parts[-1], collapse = " ")
    
 
    symbol <- sub(";.*", "", rest)
    symbol <- trimws(symbol)
    
   
    mirna <- sub("^[^;]*;\\s*", "", rest)        
    mirna <- sub("\\s*\\[.*\\].*$", "", mirna)   
    mirna <- trimws(mirna)
    
   
    data.frame(
      Entry       = entry,
      Symbol      = symbol,
      miRNA_name  = mirna,
      stringsAsFactors = FALSE
    )
  }
  
  df <- bind_rows(lapply(lines, parse_line))
  df$Entry <- as.character(df$Entry)
  
  cat(sprintf("KEGG file: %d miRNA readings\n", nrow(df)))
  cat("Foe example:\n")
  print(head(df, 10))
  
  
  dups <- df$miRNA_name[duplicated(df$miRNA_name)]
  if (length(dups) > 0) {
    cat("Dublucates of miRNA_name (leave first entry):\n")
    print(unique(dups))
    df <- df[!duplicated(df$miRNA_name), ]
  }
  
  return(df)
}

kegg_df <- parse_kegg_mirna(kegg_file)

# ============================================================
# Filtering DESeq2 with miRNA from  KEGG
# ============================================================
for (cn in names(deseq_files)) {
  path <- deseq_files[cn]
  out  <- output_names[cn]
  
  deseq <- read.csv(path, stringsAsFactors = FALSE)
  cat(sprintf("\n%s: %d rows", cn, nrow(deseq)))
  
  if (!"gene" %in% names(deseq)) {
    stop("In the file ", path, " no such column 'gene'")
  }
  
 
  filtered <- deseq %>%
    inner_join(kegg_df, by = c("gene" = "miRNA_name"))
  
  cat(sprintf("  → %d after filtration (intersect with KEGG)\n", nrow(filtered)))
  
 
  if ("padj" %in% names(filtered)) {
    filtered <- filtered %>% arrange(padj)
  }
  
  write.csv(filtered, out, row.names = FALSE)
  cat(sprintf("  Save: %s\n", out))
}



cat("\n===== Ready  =====\n")
cat("2 files are saved:\n")
for (nm in names(output_names)) cat("  ", output_names[nm], "\n")