# Jake Yeung
# Date of Creation: 2019-07-22
# File: ~/projects/scchic_gastru/scripts/Rfunctions/LSI.R
# Functions for LSI


# Main functions ----------------------------------------------------------

require(irlba)
require(hash)
require(igraph)
require(Matrix)

NormalizeMatrix <- function(count.mat, tf.method = "sums", idf.method = "inverselog", remove.empty.features = TRUE, remove.empty.cells = TRUE, use.sweep = TRUE){
  if (remove.empty.cells){
    ncols.orig <- ncol(count.mat)
    cols.keep <- which(Matrix::colSums(count.mat) > 0)
    count.mat <- count.mat[, cols.keep]
    if (ncols.orig - length(cols.keep) > 0){
      print(paste("Removing", ncols.orig - length(cols.keep), "columns. They are empty"))
    }
  }
  if (remove.empty.features){
    # remove empty features
    nrows.orig <- nrow(count.mat)
    rows.keep <- which(Matrix::rowSums(count.mat) > 0)
    count.mat <- count.mat[rows.keep, ]
    if (nrows.orig - length(rows.keep) > 0){
      print(paste("Removing", nrows.orig - length(rows.keep), "rows. They are empty"))
    }
  }
  tf <- TermFreq(count.mat, method = tf.method)
  idf <- InverseDocFreq(count.mat, method = idf.method)
  
  if (use.sweep){
    count.mat.norm <- sweep(count.mat, 
                            MARGIN = 2, 
                            STATS = tf, 
                            FUN = "*")
    # step 2: row transform (inverse document frequency)
    count.mat.norm <- sweep(count.mat.norm, 
                            MARGIN = 1, 
                            # STATS = log10(1 + ncol(count.mat) / Matrix::rowSums(count.mat)), 
                            STATS = idf, 
                            FUN = "*")
  } else {
    count.mat.norm <- count.mat * tf
    count.mat.norm <- count.mat.norm * idf
  }
  return(count.mat.norm)
}


RunLSI <- function(count.mat, tf.method = "sums", idf.method = "inverselog", n.components = 50, .log = FALSE, .center = FALSE, .truncated = TRUE){
  # rows are genes, columns are cells
  
  count.mat.norm <- NormalizeMatrix(count.mat, tf.method, idf.method, remove.empty.features = TRUE)
  
  # # remove empty features
  # rows.keep <- which(Matrix::rowSums(count.mat) > 0)
  # count.mat <- count.mat[rows.keep, ]
# 
  # 
  # tf <- TermFreq(count.mat, method = tf.method)
  # idf <- InverseDocFreq(count.mat, method = idf.method)
  # 
  # count.mat.norm <- sweep(count.mat, 
  #                         MARGIN = 2, 
  #                         STATS = tf, 
  #                         FUN = "*")
  # # step 2: row transform (inverse document frequency)
  # count.mat.norm <- sweep(count.mat.norm, 
  #                         MARGIN = 1, 
  #                         # STATS = log10(1 + ncol(count.mat) / Matrix::rowSums(count.mat)), 
  #                         STATS = idf, 
  #                         FUN = "*")
  if (.log){
    count.mat.norm <- log2(count.mat.norm * 10^6 + 1)
  }
  if (.center){
    count.mat.norm <- sweep(count.mat.norm,
                            MARGIN = 1,
                            STATS = Matrix::rowMeans(count.mat.norm),
                            FUN = "-")
  }
  if (.truncated){
    system.time(
      count.svd <- irlba(A = t(count.mat.norm), nv = n.components, scale = FALSE, center = FALSE)
    )
  } else {
    system.time(
      count.svd <- svd(x = t(count.mat.norm), nu = n.components, nv = n.components)
    )
  }
  rownames(count.svd$u) <- colnames(count.mat.norm)
  rownames(count.svd$v) <- rownames(count.mat.norm)
  return(count.svd)
}

UmapLSIOutput <- function(count.svd, jsettings, cellnames){
  # umap the outputs
  if (missing(jsettings)){
    jsettings <- umap::umap.defaults
  }
  if (missing(cellnames)){
    cellnames <- rownames(count.svd$u)
  }
  umap.out <- umap::umap(count.svd$u, config = jsettings)
  # return as tidy dataframe
  umap.dat <- data.frame(cell = cellnames, umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
  return(umap.dat)
}


# Helper functions --------------------------------------------------------

TermFreq <- function(x, method = "sums"){
  if (method == "sums"){
    jstats <- Matrix::colSums(x)
  }
  return(1 / jstats)  # multiply this to matrix
}

InverseDocFreq <- function(x, method = "inverselog"){
  if (method == "inverselog"){
    jstats <- log10(1 + ncol(x) / Matrix::rowSums(x))
  }
  return(jstats)  # multiply this to matrix 
}


