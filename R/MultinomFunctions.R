# Jake Yeung
# MultinomFunctions.R
#  
# 2019-11-14

CollapseRowsByGene <- function(count.mat, as.long = TRUE, track.kept.gene = FALSE){
  # take rowSums and pick gene with most counts
  
  # prepare raw data. Rownames of format  chr1:-17665-32335;rpl24;1
  genes.chic <- sapply(rownames(count.mat), function(x) strsplit(x, ";")[[1]][[2]])
  # get count.filt, sum across same gene
  count.mat.tmp <- data.frame(gene = genes.chic, genefull = rownames(count.mat), data.frame(sum.across.cells = as.matrix(rowSums(count.mat))), stringsAsFactors = FALSE)
  # rownames(count.mat.tmp) <- genes.chic
  count.mat.long <- reshape2::melt(count.mat.tmp, id.vars = c("gene", "genefull"), variable.name = c("cellname"), value.name = "total.cuts") %>%
    group_by(gene) %>%
    summarise(genefull = genefull[which.max(total.cuts)], 
              total.cuts = max(total.cuts))
  
  rnames.keep <- count.mat.long$genefull
  count.mat.filt <- count.mat[rnames.keep, ]
  if (as.long){
    count.mat.filt.long <- reshape2::melt(as.matrix(count.mat.filt))
    colnames(count.mat.filt.long) <- c("genefull", "cellname", "ncuts")
    return(count.mat.filt.long)
  } else {
    return(count.mat.filt)
  }
}


# MultinomFitsToLong <- function(cname, LL.ctype.lst, p.ctype.lst, cell.counts){
MultinomFitsToLong <- function(all.cells, LL.ctype.lst, count.filt){
  names(all.cells) <- all.cells
  dat.long.merged <- lapply(all.cells, MultinomFitsToLong.percell, LL.ctype.lst = LL.ctype.lst, count.filt = count.filt) %>%
    bind_rows()
  return(dat.long.merged)
}

MultinomFitsToLong.percell <- function(cname, LL.ctype.lst, count.filt){
  p.ctype.lst <- lapply(LL.ctype.lst, SoftMax)
  cell.counts <- Matrix::colSums(count.filt)
  LL.vec <- LL.ctype.lst[[cname]]
  p.vec <- p.ctype.lst[[cname]]
  cell.count <- cell.counts[[cname]]
  dat.tmp <- data.frame(cell = cname, LL = LL.vec, p = p.vec, ctype.pred = names(LL.vec), cell.size = cell.count, stringsAsFactors = FALSE)
  return(dat.tmp)
}

SetUpProbs <- function(dat.mat.filt, norm.vec = TRUE){
  # handle zeros
  # dat.mat.filt <- expbase ^ dat.mat.filt
  zero.fill <- min(as.matrix(dat.mat.filt)[which(as.matrix(dat.mat.filt) != 0)])
  dat.mat.filt[which(dat.mat.filt == 0)] <- zero.fill
  
  # make likelihoods
  probs.lst.filt <- as.list(as.data.frame(dat.mat.filt))
  # name the list just to be safe
  probs.lst.filt <- lapply(probs.lst.filt, function(x){
    names(x) <- rownames(dat.mat.filt)
    return(x)
  }) 
  if (norm.vec){
    probs.lst.filt <- lapply(probs.lst.filt, function(x){
	  return(x / sum(x))
    })
  } else {
    return(probs.lst.filt)
  }
}

FitMultinoms <- function(count.filt, all.cells, probs.lst.filt, exppower = 0.5){
  names(all.cells) <- all.cells
  LL.ctype.lst <- lapply(all.cells, function(cell.name){
  cell.vec <- count.filt[, cell.name]
  # cell.vec <- cell.vec[which(cell.vec > 0)]
  LL.vec <- sapply(probs.lst.filt, function(jprob){
    assertthat::assert_that(all(names(cell.vec) == names(jprob)))
    return(dmultinom(x = cell.vec, prob = jprob^(exppower), log = TRUE))
  })
})
}

FitMultinomsFn <- function(count.filt, all.cells, probs.lst.filt, fn){
		# allow jprob to pass through an arbitrary function
  names(all.cells) <- all.cells
  LL.ctype.lst <- lapply(all.cells, function(cell.name){
  cell.vec <- count.filt[, cell.name]
  # cell.vec <- cell.vec[which(cell.vec > 0)]
  LL.vec <- sapply(probs.lst.filt, function(jprob){
    assertthat::assert_that(all(names(cell.vec) == names(jprob)))
    return(dmultinom(x = cell.vec, prob = fn(jprob), log = TRUE))
  })
})
}

SummarizeMultinomFits <- function(LL.ctype.lst, count.filt, all.cells){
  names(all.cells) <- all.cells
  # calculate probability of model given data
  p.ctype.lst <- lapply(LL.ctype.lst, SoftMax)
  cell.counts <- Matrix::colSums(count.filt)
  # summaize
  LL.dat <- lapply(all.cells, function(cname){
    LL.vec <- LL.ctype.lst[[cname]]
    p.vec <- p.ctype.lst[[cname]]
    cell.count = cell.counts[[cname]]
    if (all(is.infinite(LL.vec))){
      LL.max <- NA
      p.max <- NA
      best.ctype <- NA
    } else {
      LL.max <- max(LL.vec)
      p.max <- max(p.vec)
      best.ctype <- names(which.max(LL.vec))
    }
    dat.tmp <- data.frame(cell = cname, LL.max = LL.max, p.max = p.max, ctype.pred = best.ctype, cell.size = cell.count, stringsAsFactors = FALSE)
    return(dat.tmp)
  }) %>%
    bind_rows()
}
