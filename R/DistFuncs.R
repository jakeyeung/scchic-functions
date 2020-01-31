# Jake Yeung
# Date of Creation: 2020-01-22
# File: ~/projects/scchicFuncs/R/DistFuncs.R
# Distance Functions

HellingerDistance <- function(x, y){
  # https://en.wikipedia.org/wiki/Hellinger_distance
  d <- 1 / sqrt(2) * sqrt( (sqrt(x) - sqrt(y)) ^ 2 )
  return(d)
}

