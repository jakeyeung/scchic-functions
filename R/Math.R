# Jake Yeung
# Date of Creation: 2019-07-27
# File: ~/projects/scchic_gastru/scripts/Rfunctions/Math.R
# Math functions

SoftMax <- function(x, return.log = TRUE){
  # numericallys table softmax by subtracting the maximum first
  # https://stackoverflow.com/questions/42599498/numercially-stable-softmax
  # x value are in log if .log = TRUE
  numer <- log(exp(x - max(x)))
  denom <- log(sum(exp(x - max(x))))
  plog <- numer - denom
  if (return.log){
    return(plog)
  } else {
    return(exp(plog))
  }
}

