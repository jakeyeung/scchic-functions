# Jake Yeung
# Date of Creation: 2020-07-21
# File: 
# description


GetWellPosition <- function(indx, platecols = 24, is.zero.base = TRUE){
  # 0 -> 1,1
  if (is.zero.base){
    indx <- indx + 1
  }
  jcol <- indx %% platecols
  jcol <- ifelse(jcol == 0, platecols, jcol)
  jrow <- ceiling(indx / platecols)
  return(c(jrow, jcol))
}

ReadChrReads <- function(inf, sort.rnames = TRUE, add.chromo = FALSE){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(chromo = sampleName)
  # remove missing
  if (add.chromo){
    # add chr
    dat$chromo <- paste("chr", dat$chromo, sep="")
  }
  if (sort.rnames){
    rnames.i <- gtools::mixedorder(rownames(dat))
    dat <- dat[rnames.i, ]
  }
  return(dat)
}



GetChromoCounts <- function(inf, spikeinchromo = "J02459.1"){
  
  dat <- ReadChrReads(inf, add.chromo = FALSE)
  
  chromos.remove <- grep(pattern = "^KI|^GL", dat$chromo, value = TRUE)
  
  chromos.keep <- dat$chromo[!dat$chromo %in% chromos.remove]
  
  dat.filt <- subset(dat, chromo %in% chromos.keep)
  
  dat.filt.long <- dat.filt %>%
    melt(., id.vars = "chromo", variable.name = "samp", value.name = "counts") %>%
    rowwise() %>%
    mutate(counts = ifelse(is.na(counts), 0, counts), 
           samp = as.character(samp), 
           experi = ClipLast(samp, jsep = "_")) %>%
    # mark = strsplit(experi, "-")[[1]][[3]],
    # conc = strsplit(experi, "-")[[1]][[5]]) %>%
    group_by(samp) %>%
    mutate(totalcounts = sum(counts, na.rm = FALSE))
  dat.spikeincounts <- subset(dat.filt.long, chromo == spikeinchromo) %>%
    dplyr::rename(spikeincounts = counts) %>%
    dplyr::select(samp, spikeincounts)
  dat.chromocounts <- subset(dat.filt.long, chromo %in% jchromos) %>%
    dplyr::rename(chromocounts = counts) %>%
    dplyr::summarise(chromocounts = sum(chromocounts)) %>%
    dplyr::select(samp, chromocounts)
  dat.filt.long <- left_join(dat.filt.long, dat.spikeincounts)
  dat.filt.long <- left_join(dat.filt.long, dat.chromocounts)
}