# Jake Yeung
# Date of Creation: 2019-09-05
# File: ~/projects/scchicFuncs/R/HandleCSVs.R
#

ReadMatTSSFormat <- function(inf, as.sparse = TRUE){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(chromo = sampleName,
                  start = V2,
                  end = V3,
                  coord = V4)
  # remove missing
  dat <- subset(dat, start != "Missing" & end != "Missing")
  # dat$coord <- dat$V4
  # dat$coord <- paste(dat$chromo, dat$startend, sep = ":")
  dat$chromo <- NULL
  dat$start <- NULL
  dat$end <- NULL
  if (as.sparse){
    coords <- dat$coord
    dat$coord <- NULL
    dat <- as.matrix(dat)
    dat[is.na(dat)] <- 0
    rownames(dat) <- coords
    dat <- Matrix::Matrix(dat, sparse = TRUE)
  }
  return(dat)
}

ReadMatSlideWinFormat <- function(inf, as.sparse = TRUE){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(chromo = sampleName,
                  start = V2,
                  end = V3,)
  # remove missing
  # add chr
  dat$chromo <- paste("chr", dat$chromo, sep="")
  dat <- subset(dat, start != "Missing")
  dat$coord <- paste(dat$chromo, paste(dat$start, dat$end, sep = "-"), sep = ":")
  dat$chromo <- NULL
  # dat$startend <- NULL
  dat$start <- NULL
  dat$end <- NULL
  if (as.sparse){
    coords <- dat$coord
    dat$coord <- NULL
    dat <- as.matrix(dat)
    dat[is.na(dat)] <- 0
    rownames(dat) <- coords
    dat <- Matrix::Matrix(dat, sparse = TRUE)
  }
  return(dat)
}

CollapseTaggedCountCnames <- function(dat, cnames.indx = 1:3){
  # 2nd row often redundant, merge with first row
  dat <- as.data.frame(dat)
  colnames(dat)[cnames.indx] <- unlist(dat[1, cnames.indx], use.names = FALSE)
  dat <- dat[-1, ]
  return(dat)
}
