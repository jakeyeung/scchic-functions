SumAcrossClusters <- function(count.mat, cnames.keep.lst){
  count.mat <- as.matrix(count.mat)
  count.vecs <- lapply(cnames.keep.lst, function(cnames.keep){
    cnames.keep.i <- which(colnames(count.mat) %in% cnames.keep)
    assertthat::assert_that(length(cnames.keep.i) > 0)
    rowSums(count.mat[, cnames.keep.i])
  })
  return(count.vecs)
}
