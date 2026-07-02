

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(dplyr)


output_dir <- "C:/Users/ПК/OneDrive/Документы/ORA_targets_v3"

cancer_files <- c(
  Breast         = "C:/Users/ПК/OneDrive/Документы/ORA_targets/miRDB_targets_Breast.csv",
  Lung           = "C:/Users/ПК/OneDrive/Документы/ORA_targets/miRDB_targets_Lung.csv",
  Larynx         = "C:/Users/ПК/OneDrive/Документы/ORA_targets/miRDB_targets_Larynx.csv",
  Hepatocellular = "C:/Users/ПК/OneDrive/Документы/ORA_targets/miRDB_targets_Hepatocellular.csv"
)

plots_dir  <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 2. Background genes (universe) ----
all_files_combined <- do.call(rbind, lapply(cancer_files, read.csv, stringsAsFactors = FALSE))
universe_genes <- as.character(unique(all_files_combined$Entrez_ID))
universe_genes <- universe_genes[!is.na(universe_genes) & universe_genes != ""]
cat(sprintf(" Universe genes: %d\n", length(universe_genes)))


# ---- 3. Function ORA for ine Cnacer  ----
run_ORA <- function(cancer_name, file_path, universe) {
  
  cat(sprintf("\n========== %s ==========\n", cancer_name))
  
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  genes <- as.character(unique(df$Entrez_ID))
  genes <- genes[!is.na(genes) & genes != ""]
  cat(sprintf("Генов-таргетов: %d\n", length(genes)))
  
  results <- list()
  
  cat("  GO BP...\n")
  results$GO_BP <- tryCatch(
    enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "BP",
             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500, readable = TRUE),
    error = function(e) { cat("    Error:", e$message, "\n"); NULL }
  )
  
  cat("  GO MF...\n")
  results$GO_MF <- tryCatch(
    enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "MF",
             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500, readable = TRUE),
    error = function(e) { cat("    Error:", e$message, "\n"); NULL }
  )
  
  cat("  GO CC...\n")
  results$GO_CC <- tryCatch(
    enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "CC",
             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
             minGSSize = 10, maxGSSize = 500, readable = TRUE),
    error = function(e) { cat("    Error:", e$message, "\n"); NULL }
  )
  
  cat("  KEGG...\n")
  ekegg <- tryCatch(
    enrichKEGG(gene = genes, organism = "hsa",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
               minGSSize = 10, maxGSSize = 500),
    error = function(e) { cat("    Error:", e$message, "\n"); NULL }
  )
  if (!is.null(ekegg)) {
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  results$KEGG <- ekegg
  
  return(results)
}


# ---- 4. Function saving tables  ----
save_tables <- function(ora_results, cancer_name, tables_dir) {
  for (analysis_name in names(ora_results)) {
    res <- ora_results[[analysis_name]]
    if (is.null(res)) next
    df <- as.data.frame(res)
    if (nrow(df) == 0) next
    df <- df[order(df$p.adjust), ]
    out_path <- file.path(tables_dir, sprintf("%s_%s_ORA.csv", cancer_name, analysis_name))
    write.csv(df, out_path, row.names = FALSE)
    cat(sprintf("  [%s %s] %d pathways → %s\n", cancer_name, analysis_name, nrow(df), basename(out_path)))
  }
}


# ---- 5. visualisation ----
save_plots <- function(ora_results, cancer_name, plots_dir) {
  for (analysis_name in names(ora_results)) {
    res <- ora_results[[analysis_name]]
    if (is.null(res)) next
    df <- as.data.frame(res)
    if (nrow(df) == 0) next
    n_show <- min(20, nrow(df))
    
    tryCatch({
      p <- dotplot(res, showCategory = n_show, title = sprintf("%s | %s", cancer_name, analysis_name)) +
        theme(axis.text.y = element_text(size = 8))
      ggsave(file.path(plots_dir, sprintf("%s_%s_dotplot.png", cancer_name, analysis_name)), p, width = 10, height = 8, dpi = 150)
    }, error = function(e) cat(sprintf("  dotplot %s %s: %s\n", cancer_name, analysis_name, e$message)))
    
    if (nrow(df) >= 2) {
      tryCatch({
        res2 <- pairwise_termsim(res)
        p2 <- emapplot(res2, showCategory = n_show, title = sprintf("%s | %s", cancer_name, analysis_name))
        ggsave(file.path(plots_dir, sprintf("%s_%s_emapplot.png", cancer_name, analysis_name)), p2, width = 10, height = 9, dpi = 150)
      }, error = function(e) cat(sprintf("  emapplot %s %s: %s\n", cancer_name, analysis_name, e$message)))
    }
    
    tryCatch({
      p3 <- cnetplot(res, showCategory = 5, title = sprintf("%s | %s", cancer_name, analysis_name))
      ggsave(file.path(plots_dir, sprintf("%s_%s_cnetplot.png", cancer_name, analysis_name)), p3, width = 12, height = 10, dpi = 150)
    }, error = function(e) cat(sprintf("  cnetplot %s %s: %s\n", cancer_name, analysis_name, e$message)))
    
    tryCatch({
      p4 <- barplot(res, showCategory = n_show, title = sprintf("%s | %s", cancer_name, analysis_name)) +
        theme(axis.text.y = element_text(size = 8))
      ggsave(file.path(plots_dir, sprintf("%s_%s_barplot.png", cancer_name, analysis_name)), p4, width = 10, height = 8, dpi = 150)
    }, error = function(e) cat(sprintf("  barplot %s %s: %s\n", cancer_name, analysis_name, e$message)))
  }
}


# ---- 6. Launching for all cancers  ----
all_results <- list()
for (cn in names(cancer_files)) {
  all_results[[cn]] <- run_ORA(cn, cancer_files[cn], universe_genes)
  save_tables(all_results[[cn]], cn, tables_dir)
  save_plots(all_results[[cn]], cn, plots_dir)
}


# ---- 7. comparative dotplot — all together  GO BP) ----
cat("\n Comparative plot (GO BP)...\n")
tryCatch({
  bp_list <- lapply(names(all_results), function(cn) all_results[[cn]][["GO_BP"]])
  names(bp_list) <- names(all_results)
  bp_list <- bp_list[!sapply(bp_list, is.null)]
  if (length(bp_list) >= 2) {
    merged <- merge_result(bp_list)
    p_compare <- dotplot(merged, showCategory = 20,
                         title = "GO BP — camparison across cancers ") +
      theme(axis.text.y = element_text(size = 8))
    ggsave(file.path(plots_dir, "ALL_cancers_GOBP_comparison_dotplot.png"),
           p_compare, width = 14, height = 10, dpi = 150)
    cat("  Saved comparative dotplot\n")
  }
}, error = function(e) cat("  Comparative plot:", e$message, "\n"))


# ============================================================
# ---- 8. Dotplot on intersect of 4 datasets ----
# ============================================================

cat("\n Intersect (GO BP)...\n")
make_intersection_dotplot <- function(analysis_name, file_prefix) {
  
  res_list <- lapply(names(all_results), function(cn) all_results[[cn]][[analysis_name]])
  names(res_list) <- names(all_results)
  res_list <- res_list[!sapply(res_list, is.null)]
  
  if (length(res_list) < 2) {
    cat(sprintf("  [%s] Not enough results\n", analysis_name))
    return(NULL)
  }
  
  
  pathway_lists <- lapply(res_list, function(x) {
    df <- as.data.frame(x)
    if (nrow(df) == 0) return(character(0))
    df$ID
  })
  
  common_ids <- Reduce(intersect, pathway_lists)
  cat(sprintf("  [%s] пересечений на все 4 рака: %d\n", analysis_name, length(common_ids)))
  
  if (length(common_ids) == 0) {
    cat(sprintf("  [%s] No common pathways - scip \n", analysis_name))
    return(NULL)
  }
  
  #  merge_result, filtering only common ID  
  merged <- merge_result(res_list)
  merged@compareClusterResult <- merged@compareClusterResult[
    merged@compareClusterResult$ID %in% common_ids, ]
  
  if (nrow(merged@compareClusterResult) == 0) return(NULL)
  
  n_show <- min(length(common_ids), 50)
  
  p <- dotplot(merged, showCategory = n_show,
               title = sprintf("Intersection — %s | %d common paths",
                               analysis_name, length(common_ids))) +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(file.path(plots_dir, sprintf("INTERSECT_%s_dotplot.png", file_prefix)),
         p, width = 14, height = max(8, n_show * 0.25), dpi = 150)
  cat(sprintf("  Saved : INTERSECT_%s_dotplot.png (%d путей)\n", file_prefix, n_show))
  
  
  out_tbl <- merged@compareClusterResult %>%
    arrange(Cluster, p.adjust)
  write.csv(out_tbl,
            file.path(tables_dir, sprintf("INTERSECT_%s_detailed.csv", file_prefix)),
            row.names = FALSE)
  cat(sprintf("  Таблица: INTERSECT_%s_detailed.csv\n", file_prefix))
}

# Запускаем для GO BP, GO MF, GO CC, KEGG — пишем отдельные файлы
make_intersection_dotplot("GO_BP", "GO_BP")
make_intersection_dotplot("GO_MF", "GO_MF")
make_intersection_dotplot("GO_CC", "GO_CC")
make_intersection_dotplot("KEGG", "KEGG")


# ---- 9. Ending ----
cat("\n===== Sum up =====\n")
for (cn in names(all_results)) {
  cat(sprintf("\n[%s]\n", cn))
  for (an in names(all_results[[cn]])) {
    res <- all_results[[cn]][[an]]
    if (!is.null(res)) {
      n <- nrow(as.data.frame(res))
      cat(sprintf("  %s: %d significant pathways \n", an, n))
    } else {
      cat(sprintf("  %s: no results \n", an))
    }
  }
}

cat(sprintf("\n Tables: %s\n", tables_dir))
cat(sprintf(" Graphs: %s\n", plots_dir))
