# Jake Yeung
# EntropyFunctions.R
#  
# 2019-12-06


Vectorize(plog2p <- function(p){
  return(ifelse(p == 0, 0, p * log2(p)))
}, vectorize.args = "p") 

CalculateEntropy <- function(p, normalize.p = FALSE){
  if (normalize.p){
    p <- p / sum(p)
  }
  S <- -sum(plog2p(p))
  return(S)
}

GetGeneCounts <- function(x){
  # get number of zeros in x
  return(Matrix::nnzero(x))
}

