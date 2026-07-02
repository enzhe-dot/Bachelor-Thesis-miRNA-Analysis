# ============================================================
# Venn diagram: ORA pathways across 4 cancers
# GO: all BP+MF+CC combined per cancer → 4-way Venn
# KEGG: one file per cancer → 4-way Venn
# ============================================================

# ── 0. Config ────────────────────────────────────────────────
BASE_DIR <- "C:/Users/ПК/OneDrive/Документы/ORA_targets_v2/tables"
CANCERS <- c("Breast", "Lung", "Larynx", "Hepatocellular")
OUTPUT_DIR <- BASE_DIR # to be adjusted 

library(dplyr)
library(ggVennDiagram)
library(ggplot2)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ── 1. Reading GO (BP+MF+CC) и Assembling + Description ─────

go_paths <- list()
go_dict  <- list()  # ID → description

for (cn in CANCERS) {
  bp <- read.csv(file.path(BASE_DIR, paste0(cn, "_GO_BP_ORA.csv")),
                 stringsAsFactors = FALSE)
  mf <- read.csv(file.path(BASE_DIR, paste0(cn, "_GO_MF_ORA.csv")),
                 stringsAsFactors = FALSE)
  cc <- read.csv(file.path(BASE_DIR, paste0(cn, "_GO_CC_ORA.csv")),
                 stringsAsFactors = FALSE)

  all_go <- rbind(bp, mf, cc)
  
  go_dict[[cn]] <- setNames(all_go[[2]], all_go[[1]])
  go_paths[[cn]] <- unique(all_go[[1]])
  message(sprintf("[GO] %s: %d unique IDs", cn, length(go_paths[[cn]])))
}


go_dict_all <- unlist(go_dict)
go_dict_all <- go_dict_all[!duplicated(names(go_dict_all))]


# ── 2.Reading KEGG ───────────────────────────────────────────

kegg_paths <- list()
kegg_dict  <- list()

for (cn in CANCERS) {
  df <- read.csv(file.path(BASE_DIR, paste0(cn, "_KEGG_ORA.csv")),
                 stringsAsFactors = FALSE)

  
  desc_col <- grep("description|term", names(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(desc_col)) desc_col <- names(df)[2]  # fallback: second column

  kegg_dict[[cn]] <- setNames(df[[desc_col]], df$ID)
  kegg_paths[[cn]] <- unique(df$ID)
  message(sprintf("[KEGG] %s: %d pathways", cn, length(kegg_paths[[cn]])))
}

kegg_dict_all <- unlist(kegg_dict)
kegg_dict_all <- kegg_dict_all[!duplicated(names(kegg_dict_all))]


# ── 3. Intersection function ─────────────────────

make_detailed_table <- function(path_list, dict, prefix) {
  n <- length(path_list)
  cancer_names <- names(path_list)

  # all combiantions 4-way, 3-way, 2-way
  combos <- list()
  combos[[length(combos) + 1]] <- cancer_names                     # 4
  for (i in 1:n) combos[[length(combos) + 1]] <- cancer_names[-i]  # 3
  for (i in 1:(n-1)) {                                              # 2
    for (j in (i+1):n) {
      combos[[length(combos) + 1]] <- cancer_names[c(i, j)]
    }
  }

  res_list <- list()
  for (combo in combos) {
    common <- Reduce(intersect, path_list[combo])
    
    others <- setdiff(cancer_names, combo)
    if (length(others) > 0) {
      for (o in others) common <- setdiff(common, path_list[[o]])
    }
    if (length(common) == 0) next

    intersection_text <- paste(combo, collapse = " & ")

    for (pid in common) {
      desc <- if (pid %in% names(dict)) dict[pid] else "-"
      res_list[[length(res_list) + 1]] <- data.frame(
        Intersection = intersection_text,
        N_cancers    = length(combo),
        Pathway_ID   = pid,
        Description  = desc,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(res_list) == 0) {
    message(sprintf("[%s] No intersections found", prefix))
    return(invisible(NULL))
  }

  res <- bind_rows(res_list)

  # UTF-8 — Excel 
  write.csv(res,
            file.path(OUTPUT_DIR, paste0(prefix, "_intersections.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")

  message(sprintf("[%s] Saved: %s_intersections.csv (%d rows, %d intersection groups)",
                  prefix, prefix, nrow(res),
                  n_distinct(res$Intersection)))
  invisible(res)
}


# ── 4. Venn diagram (GO) ─────────────────────────────────────

p_go <- ggVennDiagram(go_paths, category.names = CANCERS,
                      set_color = c("#FF8F00", "#43A047", "#8E24AA", "#00ACC1"),
                      label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#E53935") +
  labs(title = "GO pathway overlaps (BP + MF + CC)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(OUTPUT_DIR, "Venn_GO.png"), p_go, width = 10, height = 8, dpi = 150)


# ── 5. Venn diagram (KEGG) ───────────────────────────────────

p_kegg <- ggVennDiagram(kegg_paths, category.names = CANCERS,
                        set_color = c("#FF8F00", "#43A047", "#8E24AA", "#00ACC1"),
                        label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#1E88E5") +
  labs(title = "KEGG pathway overlaps") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(OUTPUT_DIR, "Venn_KEGG.png"), p_kegg, width = 10, height = 8, dpi = 150)


# ── 6. Tables of intersection ──────────────────────────────────

make_detailed_table(go_paths,   go_dict_all,   "GO")
make_detailed_table(kegg_paths, kegg_dict_all, "KEGG")

message("\nDone! Files in: ", OUTPUT_DIR)
message("  Venn_GO.png + Venn_KEGG.png")
message("  GO_intersections.csv + KEGG_intersections.csv")
