# Jake Yeung
# Date of Creation: 2019-09-05
# File: ~/projects/scchicFuncs/R/HandleCSVs.R
#

CollapseTaggedCountCnames <- function(dat, cnames.indx = 1:3){
  # 2nd row often redundant, merge with first row
  if (any(class(dat) == "data.table")){
    colnames(dat)[cnames.indx] <- unlist(dat[1, ..cnames.indx], use.names = FALSE)
  } else {
    colnames(dat)[cnames.indx] <- unlist(dat[1, cnames.indx], use.names = FALSE)
  }
  return(dat)
}
