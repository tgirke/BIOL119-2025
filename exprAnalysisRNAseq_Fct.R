###########################################
## Helper Functions for RNA-Seq Analysis ##
###########################################
## Last update: 26-Oct-2025
## Author: Thomas Girke
## Includes:
##    GEO download of gene-level read count table 
##    Volcano plot
##    
## =========================================================
## classHelpers.R — compact helpers for GEO RNA-seq labs
## - geo_counts_gq(): download + read GEO counts table (CSV/TSV, gz OK)
## - volcano_simple(): tiny volcano for DESeq2/edgeR
## - heatmap_simple(): tiny heatmap with Euclidean or Pearson distance
## =========================================================

##################################
## Download of Read Count Table ##
##################################
## Function is longer to make it work with many count tables. Problem is
## that they are inconsistently organized and structured. 
geo_counts_gq <- function(acc,
                          download_dir = "GEO_downloads",
                          prefer = c("rnaseq.*counts", "featurecounts", "raw.?counts", "counts"),
                          verbose = TRUE) {
  if (!dir.exists(download_dir)) dir.create(download_dir, recursive = TRUE)
  series_dir <- file.path(download_dir, acc)

  if (verbose) message("Downloading supplementary files for ", acc, " → ", series_dir)
  GEOquery::getGEOSuppFiles(acc, baseDir = download_dir, makeDirectory = TRUE)

  files <- list.files(series_dir, full.names = TRUE, recursive = FALSE)
  if (!length(files)) stop("No supplementary files found for ", acc)

  pick <- character(0)
  for (p in prefer) {
    cand <- grep(p, basename(files), ignore.case = TRUE, value = TRUE)
    if (length(cand)) { pick <- cand[1]; break }
  }
  if (!length(pick)) {
    cand <- grep("\\.(csv|tsv|txt)(\\.gz)?$", basename(files), ignore.case = TRUE, value = TRUE)
    if (!length(cand)) stop("Supplements present but no table looks like counts.")
    pick <- cand[1]
  }
  f <- file.path(series_dir, pick)
  if (verbose) message("Reading: ", basename(f))

  con <- if (grepl("\\.gz$", f, TRUE)) gzfile(f, "rt") else file(f, "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  first <- ""
  repeat {
    first <- readLines(con, n = 1L)
    if (!length(first)) break
    if (!startsWith(first, "#")) break
  }
  if (!length(first)) stop("Empty file: ", f)
  sep <- if (grepl(",", first)) "," else "\t"

  con2 <- if (grepl("\\.gz$", f, TRUE)) gzfile(f, "rt") else f
  tbl <- utils::read.table(con2, header = TRUE, sep = sep, quote = "",
                           comment.char = "#", check.names = FALSE, row.names = NULL)
  names(tbl)[1] <- "geneid"

  gene <- as.character(tbl[["geneid"]])
  if (anyDuplicated(gene)) gene <- make.unique(gene)

  mat <- as.matrix(tbl[, -1, drop = FALSE])
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat <- round(mat)
  storage.mode(mat) <- "integer"
  rownames(mat) <- gene

  rn <- rownames(mat)
  qc_like <- is.na(rn) | rn == "" | grepl("^(__|N_)", rn)
  is_total_row <- FALSE
  if (nrow(mat) > 1) {
    is_total_row <- all(as.numeric(mat[1, ]) == colSums(mat[-1, , drop = FALSE]))
  }
  drop_idx <- which(qc_like | (seq_len(nrow(mat)) == 1 & is_total_row))
  if (length(drop_idx)) {
    if (verbose) message("Dropping ", length(drop_idx), " non-gene QC/summary row(s).")
    mat <- mat[-drop_idx, , drop = FALSE]
  }

  if (verbose) message("Loaded ", nrow(mat), " genes × ", ncol(mat), " samples.")
  invisible(mat)
}

#########################
## Create Volcano Plot ##
#########################
volcano_simple <- function(res,
                           alpha = 0.05,
                           lfc = 1,
                           label_n = 10,
                           device = "screen",
                           out = "plot",
                           title = "Volcano plot") {
  gene <- if ("gene" %in% names(res)) res$gene else rownames(res)
  l <- if ("log2FoldChange" %in% names(res)) res$log2FoldChange else
       if ("logFC" %in% names(res)) res$logFC else stop("No LFC column found.")
  p <- if ("padj" %in% names(res)) res$padj else
       if ("FDR" %in% names(res)) res$FDR else
       if ("adj.P.Val" %in% names(res)) res$adj.P.Val else stop("No adjusted p column found.")
  p[is.na(p)] <- 1
  l[is.na(l)] <- 0
  y <- -log10(p)
  cls <- ifelse(p < alpha & l >=  lfc, "up",
         ifelse(p < alpha & l <= -lfc, "down", "ns"))

  draw <- function() {
    par(mar = c(5, 5, 3, 1))
    plot(l, y, pch = 20, cex = 0.6,
         col = c(down = "#377EB8", ns = "grey70", up = "#E41A1C")[cls],
         xlab = expression(log[2] * " fold change"),
         ylab = expression(-log[10] * " adj p-value"),
         main = title)
    abline(v = c(-lfc, lfc), lty = 2)
    abline(h = -log10(alpha), lty = 2)
    idx <- head(order(p), label_n)
    text(l[idx], y[idx] + 0.05 * ((seq_along(idx) - 1) %% 3),
         labels = gene[idx], pos = ifelse(l[idx] >= 0, 4, 2), cex = 0.65)
    legend("topright", c("Up", "Not Sig.", "Down"),
           col = c("#E41A1C", "grey70", "#377EB8"), pch = 20, bty = "n")
  }

  if (device == "png") {
    png(paste0(out, "_volcano.png"), width = 1400, height = 1000, res = 160); draw(); dev.off()
  } else if (device == "pdf") {
    pdf(paste0(out, "_volcano.pdf"), width = 7, height = 5); draw(); dev.off()
  } else {
    draw()
  }
}

#################
## Annotations ##
#################
## Get additional IDs (gene symbols and Entrez IDs) and gene descriptions
annotate_genes_simple <- function(res_df, gene_col = "gene") {
  stopifnot(gene_col %in% names(res_df) || !is.null(rownames(res_df)))

  ## Extract gene IDs 
  ids_raw <- if (gene_col %in% names(res_df)) res_df[[gene_col]] else rownames(res_df)
  ids <- gsub("\\..*$", "", ids_raw)  # drop version suffixes (e.g., ENSG000001234.7)

  ## Detect keytype
  keytype <- if (any(grepl("^ENSG\\d+$", ids))) {
    "ENSEMBL"
  } else if (any(grepl("^(NM|NR|XM|XR|NP|XP)_", ids))) {
    "REFSEQ"
  } else if (any(grepl("^[0-9]+$", ids))) {
    "ENTREZID"
  } else {
    "SYMBOL"
  }

  ##Primary mapping
  sym    <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids,
                                  column = "SYMBOL", keytype = keytype, multiVals = "first")
  desc   <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids,
                                  column = "GENENAME", keytype = keytype, multiVals = "first")
  entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids,
                                  column = "ENTREZID", keytype = keytype, multiVals = "first")

  ## Fallback for REFSEQ -> ENTREZ -> SYMBOL/GENENAME
  if (keytype == "REFSEQ" && mean(is.na(sym)) > 0.3) {
    entrez_fallback <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids,
                                             column = "ENTREZID", keytype = "REFSEQ", multiVals = "first")
    sym2  <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unname(entrez_fallback),
                                   column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    desc2 <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unname(entrez_fallback),
                                   column = "GENENAME", keytype = "ENTREZID", multiVals = "first")
    ix <- match(ids, names(entrez_fallback))
    fix <- !is.na(ix) & is.na(sym)
    sym[fix]    <- sym2[entrez_fallback[ix[fix]]]
    desc[fix]   <- desc2[entrez_fallback[ix[fix]]]
    entrez[fix] <- entrez_fallback[ix[fix]]
  }

  ## Assemble and return output
  out <- res_df
  if (!("gene" %in% names(out))) out$gene <- ids_raw
  out$gene_clean  <- ids
  out$symbol      <- unname(sym[ids])
  out$entrez_id   <- unname(entrez[ids])
  out$description <- unname(desc[ids])
  out
}

#############################################
## Pretty number printing in DT::datatable ##
#############################################
## Combines rounding and scientific number notation
pretty_num <- function(dt, digits_fixed = 4, digits_sci = 3,
                       sci_lower = 1e-4, sci_upper = 1e4) {
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("data.table package required")

  data.table::setDT(dt)
  dt[, lapply(.SD, function(col) {
    if (!is.numeric(col)) return(col)
    vapply(col, function(val) {
      if (is.na(val)) return(NA_character_)
      a <- abs(val)
      if (a >= sci_upper || (a > 0 && a < sci_lower)) {
        formatC(val, format = "e", digits = digits_sci)       # scientific only for extremes
      } else {
        formatC(round(val, digits_fixed), format = "f", digits = digits_fixed)  # fixed
      }
    }, character(1))
  })]
}

