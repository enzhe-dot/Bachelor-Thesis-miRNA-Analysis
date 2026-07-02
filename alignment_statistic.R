base_dirs <- c(
  Breast = "/on1/MAF/results2/MAFbreast/03star1out",
  Lung   = "/on1/MAF/results2/MAFlung/03star1out",
  Laryng = "/on1/MAF/results2/MAFlaryng/03star1out",
  Hepato = "/on1/MAF/results2/MAFhepatocellular/03star1out"
)

get_star_mean <- function(dir_path) {
  files <- list.files(dir_path, pattern = "Log\\.final\\.out", recursive = TRUE, full.names = TRUE)
  
  vals <- sapply(files, function(f) {
    lines <- readLines(f)
    get_val <- function(p) {
      line <- grep(p, lines, value = TRUE)
      if (length(line) == 0) return(NA)
      val <- gsub(".*\\|\\s*([0-9.]+).*", "\\1", line)
      as.numeric(gsub("%", "", val))
    }
    c(Unique      = get_val("Uniquely mapped reads %"),
      Multi       = get_val("% of reads mapped to multiple loci"),
      TooManyLoci = get_val("% of reads mapped to too many loci"),
      TooShort    = get_val("% of reads unmapped: too short"),
      UnmappedOther = get_val("% of reads unmapped: other ")
      )
  })
  
  round(rowMeans(vals, na.rm = TRUE), 2)
}

result <- data.frame(
  Dataset = names(base_dirs),
  t(sapply(base_dirs, get_star_mean)),
  row.names = NULL
)

print(result)



