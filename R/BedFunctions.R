# Jake Yeung
# Date of Creation: 2019-12-07
# File: ~/projects/scchicFuncs/R/BedFunctions.R
#

FilterBinsByBlacklist <- function(rnames.all.common, blfile){
  assertthat::assert_that(file.exists(blfile))
  rnames.gr <- makeGRangesFromDataFrame(data.frame(seqnames = sapply(rnames.all.common, GetChromo, add.chr=FALSE),
                                                   start = sapply(rnames.all.common, GetStart),
                                                   end = sapply(rnames.all.common, GetEnd)))
  bl.gr <- LoadBlacklist(inf = blfile)
  overlaps <- findOverlaps(bl.gr, rnames.gr)
  indx <- seq(length(rnames.gr))
  bl.hits.i <- unique(subjectHits(overlaps))
  bl.hits.l <- !indx %in% bl.hits.i
  rnames.gr.filt <- rnames.gr[bl.hits.l]
  rnames.gr.badbins <- rnames.gr[!bl.hits.l]
  rnames.all.common.blfilt <- names(rnames.gr.filt)
  return(rnames.all.common.blfilt)
}


