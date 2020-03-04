CalculateCumSumNegCor <- function(x, y, cname, nsteps = 100, show.plot=FALSE){
  assertthat::assert_that(length(x) == length(y))
  N <- length(x)
  nsteps.vec <- seq(nsteps)
  names(nsteps.vec) <- paste("step", nsteps.vec, sep = "_")
  qvec <- seq(from = 0, to = 1, length.out = nsteps)
  names(qvec) <- paste("q", as.character(signif(qvec, digits = 2)), sep = "_")
  xrange <- seq(from = min(x), to = max(x), length.out = nsteps)
  yrange <- seq(from = min(y), to = max(y), length.out = nsteps)

  cumsum.out <- sapply(nsteps.vec, function(i){
    ythres <- yrange[[i]]
    xthres <- xrange[[i]]
    ykeep <- which(y < ythres)
    xkeep <- which(x < xthres)
    unionkeep <- unique(c(xkeep, ykeep))
    return(length(unionkeep))
  })
  # normalize

  cumsum.norm <- cumsum.out / N

  if (show.plot){
    plot(nsteps.vec, cumsum.norm)
  }

  # make datatable
  df.tmp <- data.frame(cname = cname, counts = cumsum.out, frac = cumsum.norm, i = nsteps.vec, stringsAsFactors = FALSE)
  return(df.tmp)
}

ApplyAcrossClusters <- function(count.mat, cnames.keep.lst, fn){
  count.mat <- as.matrix(count.mat)
  count.vecs <- lapply(cnames.keep.lst, function(cnames.keep){
    cnames.keep.i <- which(colnames(count.mat) %in% cnames.keep)
    assertthat::assert_that(length(cnames.keep.i) > 0)
    apply(count.mat[, cnames.keep.i], 1, fn)
    # rowSums(count.mat[, cnames.keep.i])
  })
  return(count.vecs)
}

SumAcrossClusters <- function(count.mat, cnames.keep.lst){
  count.mat <- as.matrix(count.mat)
  count.vecs <- lapply(cnames.keep.lst, function(cnames.keep){
    cnames.keep.i <- which(colnames(count.mat) %in% cnames.keep)
    assertthat::assert_that(length(cnames.keep.i) > 0)
    rowSums(count.mat[, cnames.keep.i])
  })
  return(count.vecs)
}
