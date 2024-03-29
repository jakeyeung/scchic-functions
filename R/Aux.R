

cbind.fill.lst <- function(mats.lst, all.rnames, fill = 0){
  mats.lst.filled <- lapply(mats.lst, function(mat.tmp){
    missing.rnames <- all.rnames[!all.rnames %in% rownames(mat.tmp)]
    mat.tmp.to.fill <- matrix(data = fill, nrow = length(missing.rnames), ncol = ncol(mat.tmp), dimnames = list(missing.rnames, colnames(mat.tmp)))
    mat.tmp.bind <- rbind(mat.tmp, mat.tmp.to.fill)
    mat.tmp.bind <- mat.tmp.bind[all.rnames, ]
    return(mat.tmp.bind)
  })
  return(do.call(cbind, mats.lst.filled))
}



GetCondFromSamp.blood <- function(samp, mark){
  if (mark == "H3K4me1"){
    wt.plates <- c("B6-13W1-BM-H3K4me1", "PZ-ChIC-Bl6-BM-H3K4me1-Index")
    linneg.plates <- c("PZ-ChIC-Bl6-BM-lin-H3K4me1-")
    hsc.plates <- c("-stem-cells-")
  } else if (mark == "H3K4me3"){
    wt.plates <- c("B6-13W1-BM-H3K4me3")
    linneg.plates <- c("-Linneg-")
    hsc.plates <- c("-B6BMSC-")
  } else if (mark == "H3K27me3") {
    wt.plates <- c("B6-13W1-BM-")
    linneg.plates <- c("-Linneg-")
    hsc.plates <- c("-B6BMSC-")
  } else if (mark == "H3K9me3") {
    wt.plates <- c("B6-13W1-BM-")
    linneg.plates <- c("-Linneg-")
    hsc.plates <- c("BMSC-")
  } else {
    print(paste(mark, "not yet coded"))
  }
  blood.plates <- c("-blood-")

  wt.plates.grep <- paste(wt.plates, collapse = "|")
  linneg.plates.grep <- paste(linneg.plates, collapse = "|")
  hsc.plates.grep <- paste(hsc.plates, collapse = "|")
  blood.plates.grep <- paste(blood.plates, collapse = "|")

  is.wt <- grepl(wt.plates.grep, samp)
  is.linneg <- grepl(linneg.plates.grep, samp)
  is.hsc <- grepl(hsc.plates.grep, samp)
  is.blood <- grepl(blood.plates.grep, samp)

  bool.vec <- c(is.wt, is.linneg, is.hsc, is.blood)
  assertthat::assert_that(sum(bool.vec) == 1)

  indx <- which.max(bool.vec)
  if (indx == 1){
    cond <- "Unenriched"
  } else if (indx == 2){
    cond <- "Linneg"
  } else if (indx == 3){
    cond <- "StemCell"
  } else if (indx == 4){
    cond <- "Blood"
  } else {
    stop("must be 1, 2, 3, or 4. Found: ", indx)
  }
}

GetCondFromSamp <- function(samp, mark = "H2K4me1"){
  if (mark == "H3K4me1"){
    wt.plates <- c("B6-13W1-BM-H3K4me1", "PZ-ChIC-Bl6-BM-H3K4me1-Index")
    linneg.plates <- c("PZ-ChIC-Bl6-BM-lin-H3K4me1-")
    hsc.plates <- c("-stem-cells-")
  } else if (mark == "H3K4me3"){
    wt.plates <- c("B6-13W1-BM-H3K4me3")
    linneg.plates <- c("-Linneg-")
    hsc.plates <- c("-B6BMSC-")
  } else if (mark == "H3K27me3") {
    wt.plates <- c("B6-13W1-BM-")
    linneg.plates <- c("-Linneg-")
    hsc.plates <- c("-B6BMSC-")
  } else if (mark == "H3K9me3") {
    wt.plates <- c("B6-13W1-BM-")
    linneg.plates <- c("-Linneg-")
    hsc.plates <- c("BMSC-")
  } else {
    print(paste(mark, "not yet coded"))
  }
  wt.plates.grep <- paste(wt.plates, collapse = "|")
  linneg.plates.grep <- paste(linneg.plates, collapse = "|")
  hsc.plates.grep <- paste(hsc.plates, collapse = "|")
  is.wt <- grepl(wt.plates.grep, samp)
  is.linneg <- grepl(linneg.plates.grep, samp)
  is.hsc <- grepl(hsc.plates.grep, samp)
  bool.vec <- c(is.wt, is.linneg, is.hsc)
  assertthat::assert_that(sum(bool.vec) == 1)

  indx <- which.max(bool.vec)
  if (indx == 1){
    cond <- "Unenriched"
  } else if (indx == 2){
    cond <- "Linneg"
  } else if (indx == 3){
    cond <- "StemCell"
  } else {
    stop("must be 1, 2, or 3. Found: ", indx)
  }
}


SumAcrossChromos <- function(count.mat, jchromos, colfunction = mean){
  names(jchromos) <- jchromos
  rows.i.lst <- lapply(jchromos, function(jchromo){
    terms.keep <- grep(paste0("^", jchromo, ":"), rownames(count.mat))
  })
  mat.sum.by.chromo.merged <- lapply(rows.i.lst, function(terms.keep){
    return(apply(count.mat[terms.keep, ], MARGIN = 2, colfunction))
  })
  mat.sum.by.chromo.merged <- do.call(rbind, mat.sum.by.chromo.merged)
  jcells <- colnames(count.mat)
  names(jcells) <- jcells
  reads.by.chromo <- lapply(jcells, function(jcell){
    reads.by.chromo <- mat.sum.by.chromo.merged[, jcell]
    reads.by.chromo.dat <- data.frame(chromo = names(reads.by.chromo), ncuts = unlist(reads.by.chromo), stringsAsFactors = FALSE) %>%
      mutate(cell = jcell)
    return(reads.by.chromo.dat)
  }) %>%
    bind_rows()
  return(reads.by.chromo)
}

GrepAndWriteMat <- function(mat.tmp, jgrp, jgrp.name, outf, invert=FALSE, save.rds = TRUE){
  cols.i <- grepl(jgrp, colnames(mat.tmp))
  if (invert){
    cols.i <- !cols.i
  }
  mat.tmp.filt <- mat.tmp[, cols.i]
  assertthat::assert_that(ncol(mat.tmp.filt) > 0)
  print(jgrp.name)
  print(jgrp)
  print(dim(mat.tmp.filt))
  # write to output
  if (save.rds){
    saveRDS(mat.tmp.filt, file = outf)
  }
  return(mat.tmp.filt)
}

GetMarkFromStr <- function(x){
  jpat <- "H3K[0-9]*.me[0-9]"
  jmatch <- stringr::str_match(x, jpat)
  assertthat::assert_that(nrow(jmatch) == 1 & ncol(jmatch) == 1)
  return(jmatch[[1]])
}

ClipLast <- function(x, jsep = "-", jsep.out = NULL){
  # B6-13W1-BM-H3K4me3-1_269 -> B6-13W1-BM-H3K4me3
  if (is.null(jsep.out)){
    jsep.out <- jsep
  }
  jsplit <- strsplit(x, jsep)[[1]]
  # remove last one
  N <- length(jsplit) - 1
  return(paste(jsplit[1:N], collapse = jsep.out))
}

KeepLast <- function(x, jsep = "-"){
  # B6-13W1-BM-H3K4me3-1_269 -> 1_269
  jsplit <- strsplit(x, jsep)[[1]]
  # remove last one
  N <- length(jsplit)
  return(jsplit[[N]])
}

GetBins <- function(jchromo, midpt, jname, Q, winsize = 100000L){
  jleft <- midpt - winsize / 2
  jright <- midpt + winsize / 2
  Qsub <- subset(Q, chromo == jchromo & coord > jleft & coord < jright)
  Qsub$name <- jname
  return(as.data.frame(Qsub))
}

Vectorize(AssignHash <- function(x, jhash, null.fill = NA){
  # assign hash key to hash value, handle NULLs
  # null.fill = "original", returns original value x into jhash
  x.mapped <- jhash[[as.character(x)]]
  if (is.null(x.mapped)){
    if (is.na(null.fill)){
      x.mapped <- null.fill
    } else if (as.character(null.fill) == "original"){
      x.mapped <- x
    } else {
      x.mapped <- null.fill
    }
  }
  return(x.mapped)
}, vectorize.args = "x")

GetGOData <- function(out.tb.lst, ontology, jdecreasing=TRUE, order.by="Hyper_Fold_Enrichment"){
  topics <- which(lapply(1:length(out.tb.lst), function(i) !is.null(out.tb.lst[[i]][[ontology]])) == TRUE)
  GOdata <- lapply(topics, function(i) out.tb.lst[[i]][[ontology]])
  GOdata <- lapply(1:length(GOdata), function(i) GOdata[[i]][order(GOdata[[i]][order.by],
                                                                   decreasing = jdecreasing),])
  GOdata <- lapply(1:length(GOdata), function(i) GOdata[[i]][1:top,])

  # annotate each GOdata by topic number
  for (i in topics){
    GOdata[[i]]$topic <- i
  }
  GOdata <- dplyr::bind_rows(GOdata) %>%
    mutate(topic = factor(topic, levels = sort(unique(topic))),
           onto = ontology)
  return(GOdata)
}


GetPeakSize <- function(coord){
  # chr1:3005258-3006803 -> 1545
  jstart <- as.numeric(strsplit(strsplit(coord, ":")[[1]][[2]], "-")[[1]][[1]])
  jend <- as.numeric(strsplit(strsplit(coord, ":")[[1]][[2]], "-")[[1]][[2]])
  return(jend - jstart)
}

CleanCoords <- function(x){
  # "chr11:3,089,532-3,115,355" -> remove commas
  return(gsub(pattern = ",", "", x))
}

BinarizeMatrix <- function(x){
  # https://stackoverflow.com/questions/14526429/turn-a-count-matrix-into-a-binary-existence-matrix
  xbin <- as.numeric(as.matrix(x) > 0)
  xbin <- Matrix::Matrix(xbin, sparse = TRUE, nrow = nrow(x), ncol = ncol(x))
  # get back the column and row names
  rownames(xbin) <- rownames(x)
  colnames(xbin) <- colnames(x)
  return(xbin)
}

ParseCoord <- function(x){
  # chr7:103,796,583-103,857,605 -> chromo, start, end
  out <- list()
  out$chromo <- strsplit(x, split = ":")[[1]][[1]]
  out$start <- strsplit(strsplit(x, split = ":")[[1]][[2]], split = "-")[[1]][[1]]
  out$end <- strsplit(strsplit(x, split = ":")[[1]][[2]], split = "-")[[1]][[2]]
  # remove commas
  out <- lapply(out, function(x) gsub(",", "", x))
  out$start <- as.numeric(gsub(",", "", out$start))
  out$end <- as.numeric(gsub(",", "", out$end))
  return(out)
}

GetChromo <- function(x, add.chr = FALSE){
  # chrY:90799295-90803056 -> chrY
  chromo <- strsplit(x, ":")[[1]][[1]]
  if (add.chr){
    chromo <- paste0("chr", chromo)
  }
  return(chromo)
}

GetStart <- function(x){
  # chrY:90799295-90803056 -> 90799295
  return(strsplit(strsplit(x, ":")[[1]][[2]], "-")[[1]][[1]])
}

GetEnd <- function(x){
  # chrY:90799295-90803056 -> 90799295
  return(strsplit(strsplit(x, ":")[[1]][[2]], "-")[[1]][[2]])
}


