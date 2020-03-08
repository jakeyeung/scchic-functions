# Jake Yeung
# Date of Creation: 2019-09-05
# File: ~/projects/scchicFuncs/R/HandleCSVs.R
#

# GrepAndWriteMat <- function(mat.tmp, jgrp, jgrp.name, outf){
#   cols.i <- grepl(jgrp, colnames(mat.tmp))
#   mat.tmp.filt <- mat.tmp[, cols.i]
#   assertthat::assert_that(ncol(mat.tmp.filt) > 0)
#   print(jgrp.name)
#   print(jgrp)
#   print(dim(mat.tmp.filt))
#   print(1 - Matrix::nnzero(mat.tmp.filt) / length(mat.tmp.filt))
#   # write to output
#   saveRDS(mat.tmp.filt, file = outf)
#   return(mat.tmp.filt)
# }

LoadCellAnnotsEtOH.MorePlates <- function(annots.dir = "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_mflorescu/annotations_EtOH", 
                                          jmarks.annot = c("K27", "K9", "K27_K9")){
  names(jmarks.annot) <- jmarks.annot
  annots <- lapply(jmarks.annot, function(jmark.annot){
    annots1 <- LoadCellAnnotsEtOH(annots.dir = annots.dir, jmarks.annot = jmark.annot) %>%
      ungroup()  # no fltering needed because jmark.annot lilkely will not match two different experis
    if (jmark.annot == "K27"){
      annots2 <- annots1 %>%
        mutate(cell = gsub("BM-EtOH-EndAllAB-K27m3-p9", "BM-EtOH-EndAllAB-K27m3-p8", cell))
      annots3 <- annots1 %>%
        mutate(cell = gsub("BM-EtOH-EndAllAB-K27m3-p9", "BM-EtOH-MNctrl-K27m3-p7", cell))
    } else if (jmark.annot == "K9"){
      annots2 <- annots1 %>%
        mutate(cell = gsub("BM-EtOH-EndAllAB-K9m3-p6", "BM-EtOH-EndAllAB-K9m3-p5", cell))
      annots3 <- annots1 %>%
        mutate(cell = gsub("BM-EtOH-EndAllAB-K9m3-p6", "BM-EtOH-MNctrl-K9m3-p4", cell))
    } else if (jmark.annot == "K27_K9"){
      annots2 <- annots1 %>%
        mutate(cell = gsub("BM-EtOH-EndAllAB-K27m3-K9m3-p3", "BM-EtOH-EndAllAB-K27m3-K9m3-p2", cell))
      annots3 <- annots1 %>%
        mutate(cell = gsub("BM-EtOH-EndAllAB-K27m3-K9m3-p3", "BM-EtOH-MNctrl-K27m3-K9m3-p1", cell))
    } else {
      warning(jmark.annot, "not coded")
      return(NULL)
    }
    annots <- do.call(rbind, list(annots1, annots2, annots3)) %>%
      mutate(cell.orig = cell,
             cell = make.names(cell))
  }) %>%
    bind_rows()
  return(annots)
}

LoadCellAnnotsEtOH <- function(annots.dir = "/Users/yeung/data/dblchic/annotations_EtOH", jmarks.annot = c("K27", "K9", "K27_K9")){
  # get annots
  # annots.dir <- "/Users/yeung/data/dblchic/annotations_EtOH"
  annots.lst <- sapply(jmarks.annot, function(x) list.files(path = file.path(annots.dir), pattern = paste0("*[a-z]_", x, ".txt"), full.names = TRUE)) %>%
    unlist()
  annots.name <- sapply(annots.lst, function(x) strsplit(basename(x), split = "_")[[1]][[1]], USE.NAMES = FALSE)
  annots.dat <- mapply(FUN = LoadCellAnnots, annots.lst, annots.name, SIMPLIFY = FALSE) %>%
    bind_rows()
  annots.dat <- annots.dat[!duplicated(annots.dat), ]
  return(annots.dat)
}


LoadCellAnnots <- function(inf, annot){
  # inf <- "/Users/yeung/data/dblchic/annotations_EtOH"
  dat <- fread(inf, header = FALSE) %>%
    dplyr::rename(cell = V1) %>%
    mutate(celltype = annot) %>%
    as.data.frame()
  return(dat)
}


ReadMatTSSFormat <- function(inf, as.sparse = TRUE, add.coord = FALSE){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(chromo = sampleName,
                  start = V2,
                  end = V3,
                  coord = V4)
  # remove missing
  dat <- subset(dat, start != "Missing" & end != "Missing")
  # dat$coord <- dat$V4
  # dat$coord <- paste(dat$chromo, dat$startend, sep = ":")
  if (add.coord){
      coord2 <- paste(dat$chromo, paste(dat$start, dat$end, sep = "-"), sep = ":")
      dat$coord <- paste(coord2, dat$coord, sep = ";")
  }

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
