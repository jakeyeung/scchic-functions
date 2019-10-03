# Jake Yeung
# Date of Creation: 2019-08-05
# File: ~/projects/scchic_gastru/scripts/Rfunctions/DoubleStaining.R
# Double staining functions

GetDoubleName <- function(jname.test, cell.lambda.hash, prefix = "Pseudo:"){
  # used for simulation
  return(paste0(prefix, paste(lapply(strsplit(jname.test, ";")[[1]], function(x) cell.lambda.hash[[x]]), collapse = ";")))
}

GetLLMerged <- function(w, cell.count.raw.merged, dat.impute.repress.lst, dat.impute.active, return.mat = FALSE){
  w.vec <- c(w, 1-w)  # active proportion, repress proportion
  if (any(w.vec < 0)){
    print(paste("Error message: w is negative:"))
    print(w.vec)
    w.vec[which.min(w.vec)] <- 0
    w.vec[which.max(w.vec)] <- 1
    print("Fixed w:")
    print(w.vec)
  }
  dat.imputed.merged <- lapply(dat.impute.repress.lst, function(repress.row){
    t(outmat <- apply(dat.impute.active, 1, function(act.row){
      matrixStats::colWeightedMeans(rbind(act.row, repress.row), w = w.vec)  # standard weighted.mean is SUPER slow
    }))
  })
  # Calculate all the likelihoods  ------------------------------------------
  # output will be N active cell by N repress cell
  repress.names <- names(dat.imputed.merged)
  names(repress.names) <- repress.names

  ll.lst <- lapply(repress.names, function(repress.name){
    repress.mat <- dat.imputed.merged[[repress.name]]
    ll <- as.data.frame(apply(repress.mat, 1, function(merged.prob.vec){
      return(dmultinom(x = cell.count.raw.merged, prob = merged.prob.vec, log = TRUE))
    }))
    rownames(ll) <- rownames(repress.mat)
    colnames(ll) <- repress.name
    return(ll)
    # return(repress.mat)
  })
  ll.merged <- do.call(cbind, ll.lst)
  if (return.mat){
    return(ll.merged)
  } else {
    # return negative of maximum likelihood
    # negative because optimizers find minimums
    return(-1 * max(ll.merged))
  }
}

# GetLLMerged <- function(w, cell.count.raw.merged, dat.impute.repress.lst, dat.impute.active, return.mat = FALSE){
#   w.vec <- c(w, 1-w)  # active proportion, repress proportion 
#   dat.imputed.merged <- lapply(dat.impute.repress.lst, function(repress.row){
#     t(outmat <- apply(dat.impute.active, 1, function(act.row){
#       matrixStats::colWeightedMeans(rbind(act.row, repress.row), w = w.vec)  # standard weighted.mean is SUPER slow 
#     }))
#   })
#   # Calculate all the likelihoods  ------------------------------------------
#   # output will be N active cell by N repress cell 
#   repress.names <- names(dat.imputed.merged)
#   names(repress.names) <- repress.names
#   
#   ll.lst <- lapply(repress.names, function(repress.name){
#     repress.mat <- dat.imputed.merged[[repress.name]]
#     ll <- as.data.frame(apply(repress.mat, 1, function(merged.prob.vec){
#       return(dmultinom(x = cell.count.raw.merged, prob = merged.prob.vec, log = TRUE))
#     }))
#     rownames(ll) <- rownames(repress.mat)
#     colnames(ll) <- repress.name
#     return(ll)
#     # return(repress.mat)
#   })
#   ll.merged <- do.call(cbind, ll.lst)
#   if (return.mat){
#     return(ll.merged)
#   } else {
#     # return negative of maximum likelihood 
#     # negative because optimizers find minimums 
#     return(-1 * max(ll.merged))
#   }
# }

FitMixingWeight <- function(cell.count.raw.merged, dat.impute.repress.lst, dat.impute.active, w.init, w.lower, w.upper, jmethod = "Brent"){
  # https://www.psychologie.uni-heidelberg.de/ae/meth/team/mertens/blog/hessian.nb.html 
  # method by Brent is recommended?
  # jmethod <- ""L-BFGS-B""
  optim(par = w.init, 
        fn = GetLLMerged, 
        cell.count.raw.merged = cell.count.raw.merged, 
        dat.impute.repress.lst = dat.impute.repress.lst, 
        dat.impute.active = dat.impute.active, 
        method = jmethod, hessian=TRUE, lower = w.lower, upper = w.upper)
}
