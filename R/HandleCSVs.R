# Jake Yeung
# Date of Creation: 2019-09-05
# File: ~/projects/scchicFuncs/R/HandleCSVs.R
#

CollapseTaggedCountCnames <- function(dat, cnames.indx = 1:3){
  # 2nd row often redundant, merge with first row
  dat <- as.data.frame(dat)
  colnames(dat)[cnames.indx] <- unlist(dat[1, cnames.indx], use.names = FALSE)
  dat <- dat[-1, ]
  return(dat)
}
