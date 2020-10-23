# Jake Yeung
# Date of Creation: 2020-07-21
# File: 
# description



# Cell cycle functions ----------------------------------------------------

AddCellCycleLabel.bydat <- function(dat){
   # expects samp in colnames
  dat <- dat %>%
    rowwise() %>%
    mutate(rowcoord = AddPlateCoordinates(samp)$rowcoord,
           colcoord = AddPlateCoordinates(samp)$colcoord,
           cellcycle.str = AddCellCycleLabel(colcoord))
  return(dat)
}

AddPlateCoordinates <- function(samp){
  cellid = paste("cell", strsplit(samp, "_")[[1]][[2]], sep = "")
  indx = strsplit(samp, "_")[[1]][[2]]
  rowcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[1]]
  colcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[2]]
  return(list(rowcoord = rowcoord, colcoord = colcoord))
}

AddCellCycleLabel <- function(colcoord){
  # every 8 columns is a new cell cycle
  cellcycle.lst = list("0" = "0_G1", "1" = "1_S", "2" = "2_G2/M")
  cellcycle.hash <- hash::hash(cellcycle.lst)
  
  cellcycle = as.character(floor( ( colcoord - 1 ) / 8))
  cellcycle.str = cellcycle.hash[[cellcycle]]
  return(cellcycle.str)
}



TermFreq2 <- function(x, method = "sums", spikeincounts = NA){
  if (method == "sums"){
    jstats <- Matrix::colSums(x)
  } else if (method == "spikeins"){
    # assertthat::assert_that(!is.na(spikeincounts))
    jstats <- 1 / spikeincounts
  } else {
    print(paste(method, "not yet coded"))
  }
  return(1 / jstats)  # multiply this to matrix
}


NormalizeMatrix2 <- function(count.mat, tf.method = "sums", idf.method = "inverselog", remove.empty.features = TRUE, remove.empty.cells = TRUE, use.sweep = TRUE, spikeincounts = NA){
  if (remove.empty.cells){
    ncols.orig <- ncol(count.mat)
    cols.keep <- which(Matrix::colSums(count.mat) > 0)
    count.mat <- count.mat[, cols.keep]
    if (ncols.orig - length(cols.keep) > 0){
      print(paste("Removing", ncols.orig - length(cols.keep), "columns. They are empty"))
    }
  }
  if (remove.empty.features){
    # remove empty features
    nrows.orig <- nrow(count.mat)
    rows.keep <- which(Matrix::rowSums(count.mat) > 0)
    count.mat <- count.mat[rows.keep, ]
    if (nrows.orig - length(rows.keep) > 0){
      print(paste("Removing", nrows.orig - length(rows.keep), "rows. They are empty"))
    }
  }
  tf <- TermFreq2(count.mat, method = tf.method, spikeincounts = spikeincounts)
  idf <- InverseDocFreq(count.mat, method = idf.method)
  
  if (use.sweep){
    count.mat.norm <- sweep(count.mat, 
                            MARGIN = 2, 
                            STATS = tf, 
                            FUN = "*")
    # step 2: row transform (inverse document frequency)
    count.mat.norm <- sweep(count.mat.norm, 
                            MARGIN = 1, 
                            # STATS = log10(1 + ncol(count.mat) / Matrix::rowSums(count.mat)), 
                            STATS = idf, 
                            FUN = "*")
  } else {
    count.mat.norm <- count.mat * tf
    count.mat.norm <- count.mat.norm * idf
  }
  return(count.mat.norm)
}

RunLSI2 <- function(count.mat, tf.method = "sums", idf.method = "inverselog", n.components = 50, .log = FALSE, .center = FALSE, .truncated = TRUE, spikeincounts = NA){
  # rows are genes, columns are cells
  
  count.mat.norm <- NormalizeMatrix2(count.mat, tf.method, idf.method, remove.empty.features = TRUE, spikeincounts = spikeincounts)
  
  if (.log){
    count.mat.norm <- log2(count.mat.norm * 10^6 + 1)
  }
  if (.center){
    count.mat.norm <- sweep(count.mat.norm,
                            MARGIN = 1,
                            STATS = Matrix::rowMeans(count.mat.norm),
                            FUN = "-")
  }
  if (.truncated){
    system.time(
      count.svd <- irlba(A = t(count.mat.norm), nv = n.components, scale = FALSE, center = FALSE)
    )
  } else {
    system.time(
      count.svd <- svd(x = t(count.mat.norm), nu = n.components, nv = n.components)
    )
  }
  rownames(count.svd$u) <- colnames(count.mat.norm)
  rownames(count.svd$v) <- rownames(count.mat.norm)
  return(count.svd)
}





# Concentration dilution series functions ---------------------------------



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




FitGlmRowSpikeins <- function(input.dat, return.fit.obj = FALSE){
  jfit.glm <- glm(formula = genecounts ~ ncells + offset(log(spikeincounts)), data = input.dat, family = poisson)
  if (return.fit.obj){
    return(jfit.glm)
  }
  jslope.glm.ln <- coefficients(jfit.glm)[["ncells"]]
  jslope.pval.glm <- summary(jfit.glm)$coefficients[, "Pr(>|z|)"][["ncells"]]
  jslope.se.glm <- summary(jfit.glm)$coefficients[, "Std. Error"][["ncells"]]
  jfit.dat <- data.frame(slope.ln = jslope.glm.ln, slope = jslope.glm.ln / log(2), pval = jslope.pval.glm, slope.se.ln = jslope.se.glm, slope.se = jslope.se.glm / log(2), stringsAsFactors = FALSE)
  if (!return.fit.obj){
    return(jfit.dat)
  }
}


FitLmRowSpikeins <- function(input.dat, return.fit.obj = FALSE, pseudocount = 1){
  jfit.lm <- lm(formula = log( (genecounts + pseudocount) / spikeincounts) ~ ncells, data = input.dat)
  if (return.fit.obj){
    return(jfit.lm)
  }
  jslope.lm.ln <- coefficients(jfit.lm)[["ncells"]]
  jslope.pval.lm <- summary(jfit.lm)$coefficients[, "Pr(>|t|)"][["ncells"]]
  jslope.se.lm <- summary(jfit.lm)$coefficients[, "Std. Error"][["ncells"]]
  jfit.dat <- data.frame(slope.ln = jslope.lm.ln, slope = jslope.lm.ln / log(2), pval = jslope.pval.lm, slope.se.ln = jslope.se.lm, slope.se = jslope.se.lm / log(2), stringsAsFactors = FALSE)
  if (!return.fit.obj){
    return(jfit.dat)
  }
}

FitGlmRowChromocounts <- function(input.dat, return.fit.obj = FALSE){
  jfit.glm <- glm(formula = genecounts ~ ncells + offset(log(chromocounts)), data = input.dat, family = poisson)
  if (return.fit.obj){
    return(jfit.glm)
  }
  jslope.glm.ln <- coefficients(jfit.glm)[["ncells"]]
  jslope.pval.glm <- summary(jfit.glm)$coefficients[, "Pr(>|z|)"][["ncells"]]
  jslope.se.glm <- summary(jfit.glm)$coefficients[, "Std. Error"][["ncells"]]
  jfit.dat <- data.frame(slope.ln = jslope.glm.ln, slope = jslope.glm.ln / log(2), pval = jslope.pval.glm, slope.se.ln = jslope.se.glm, slope.se = jslope.se.glm / log(2), stringsAsFactors = FALSE)
  if (!return.fit.obj){
    return(jfit.dat)
  }
}

FitLmRowChromocounts <- function(input.dat, return.fit.obj = FALSE, pseudocount = 1){
  jfit.lm <- lm(formula = log( (genecounts + pseudocount) / chromocounts) ~ ncells, data = input.dat)
  if (return.fit.obj){
    return(jfit.lm)
  }
  jslope.lm.ln <- coefficients(jfit.lm)[["ncells"]]
  jslope.pval.lm <- summary(jfit.lm)$coefficients[, "Pr(>|t|)"][["ncells"]]
  jslope.se.lm <- summary(jfit.lm)$coefficients[, "Std. Error"][["ncells"]]
  jfit.dat <- data.frame(slope.ln = jslope.lm.ln, slope = jslope.lm.ln / log(2), pval = jslope.pval.lm, slope.se.ln = jslope.se.lm, slope.se = jslope.se.lm / log(2), stringsAsFactors = FALSE)
  if (!return.fit.obj){
    return(jfit.dat)
  }
}

FitGlmRowRaw <- function(input.dat, return.fit.obj = FALSE){
  jfit.glm <- glm(formula = genecounts ~ ncells, data = input.dat, family = poisson)  # no offset
  if (return.fit.obj){
    return(jfit.glm)
  }
  jslope.glm.ln <- coefficients(jfit.glm)[["ncells"]]
  jslope.pval.glm <- summary(jfit.glm)$coefficients[, "Pr(>|z|)"][["ncells"]]
  jslope.se.glm <- summary(jfit.glm)$coefficients[, "Std. Error"][["ncells"]]
  jfit.dat <- data.frame(slope.ln = jslope.glm.ln, slope = jslope.glm.ln / log(2), pval = jslope.pval.glm, slope.se.ln = jslope.se.glm, slope.se = jslope.se.glm / log(2), stringsAsFactors = FALSE)
  if (!return.fit.obj){
    return(jfit.dat)
  }
}


FitNormCountsToNcells.lm <- function(input.dat, return.fit.obj = FALSE){
  jfit <- lm(formula = log(chromocounts / spikeinconc) ~ ncells, data = input.dat)
  if (return.fit.obj){
    return(jfit)
  }
  jslope <- coefficients(jfit)[["ncells"]]
  jslope.pval <- summary(jfit)$coefficients[, "Pr(>|t|)"][["ncells"]]
  jslope.se <- summary(jfit)$coefficients[, "Std. Error"][["ncells"]]
  jfit.dat <- data.frame(slope.ln = jslope, slope = jslope / log(2), pval = jslope.pval, slope.se.ln = jslope.se, slope.se = jslope.se / log(2), stringsAsFactors = FALSE)
  if (!return.fit.obj){
    return(jfit.dat)
  }
}

FitNormCountsToNcells.lm <- function(input.dat, return.fit.obj = FALSE){
  # jfit <- lm(formula = log(chromocounts / spikeinconc) ~ ncells, data = input.dat)
  jfit <- lm(formula = log(chromocounts / spikeincounts) ~ ncells, data = input.dat)
  if (return.fit.obj){
    return(jfit)
  }
  jslope <- coefficients(jfit)[["ncells"]]
  jslope.pval <- summary(jfit)$coefficients[, "Pr(>|t|)"][["ncells"]]
  jslope.se <- summary(jfit)$coefficients[, "Std. Error"][["ncells"]]
  jfit.dat <- data.frame(slope.ln = jslope, slope = jslope / log(2), pval = jslope.pval, slope.se.ln = jslope.se, slope.se = jslope.se / log(2), stringsAsFactors = FALSE)
  if (!return.fit.obj){
    return(jfit.dat)
  }
}

FitNormCountsToNcells.lm.naive <- function(input.dat, return.fit.obj = FALSE){
  jfit <- lm(formula = log(chromocounts) ~ ncells, data = input.dat)
  # jfit <- lm(formula = log(chromocounts) ~ ncells, data = input.dat)
  if (return.fit.obj){
    return(jfit)
  }
  jslope <- coefficients(jfit)[["ncells"]]
  jslope.pval <- summary(jfit)$coefficients[, "Pr(>|t|)"][["ncells"]]
  jslope.se <- summary(jfit)$coefficients[, "Std. Error"][["ncells"]]
  jfit.dat <- data.frame(slope.ln = jslope, slope = jslope / log(2), pval = jslope.pval, slope.se.ln = jslope.se, slope.se = jslope.se / log(2), stringsAsFactors = FALSE)
  if (!return.fit.obj){
    return(jfit.dat)
  }
}


FitNormCountsToNcells.glm <- function(input.dat, return.fit.obj = FALSE){
  jfit.glm <- glm(formula = chromocounts ~ ncells + offset(log(spikeincounts)), data = input.dat, family = poisson)
  if (return.fit.obj){
    return(jfit.glm)
  }
  jslope.glm.ln <- coefficients(jfit.glm)[["ncells"]]
  jslope.pval.glm <- summary(jfit.glm)$coefficients[, "Pr(>|z|)"][["ncells"]]
  jslope.se.glm <- summary(jfit.glm)$coefficients[, "Std. Error"][["ncells"]]
  jfit.dat <- data.frame(slope.ln = jslope.glm.ln, slope = jslope.glm.ln / log(2), pval = jslope.pval.glm, slope.se.ln = jslope.se.glm, slope.se = jslope.se.glm / log(2), stringsAsFactors = FALSE)
  if (!return.fit.obj){
    return(jfit.dat)
  }
}



GetChromoCounts <- function(inf, spikeinchromo = "J02459.1", chromos.keep = NA){
  
  dat <- ReadChrReads(inf, add.chromo = FALSE)
  
  
  if (is.na(chromos.keep)){
    chromos.remove <- grep(pattern = "^KI|^GL", dat$chromo, value = TRUE)
    chromos.keep.SpikeAndChromo <- dat$chromo[!dat$chromo %in% chromos.remove]
    chromos.keep <- chromos.keep.SpikeAndChromo[which(chromos.keep.SpikeAndChromo != spikeinchromo)]
  } else {
    chromos.keep.SpikeAndChromo <- c(chromos.keep, spikeinchromo)
  }
  
  dat.filt <- subset(dat, chromo %in% chromos.keep.SpikeAndChromo)
  
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
  dat.chromocounts <- subset(dat.filt.long, chromo %in% chromos.keep) %>%
    dplyr::rename(chromocounts = counts) %>%
    dplyr::summarise(chromocounts = sum(chromocounts)) %>%
    dplyr::select(samp, chromocounts)
  dat.filt.long <- left_join(dat.filt.long, dat.spikeincounts)
  dat.filt.long <- left_join(dat.filt.long, dat.chromocounts)
  return(dat.filt.long)
}