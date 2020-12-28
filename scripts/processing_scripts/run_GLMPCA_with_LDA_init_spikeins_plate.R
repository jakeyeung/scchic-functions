# Jake Yeung
# Date of Creation: 2020-08-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/1-correct_LDA_from_peaks_by_variance.MouseBM.poisson_faster.R
# 

rm(list=ls())
jstart <- Sys.time()

InitGLMPCAfromLDA.plate <- function(count.mat, tm.result, dat.meta, plate.cname = "ncuts.var", bins.keep = 100, do.log = FALSE, svd.on.Yinit = TRUE, sz.cname = "none", bins.keep.manual = "none"){

  ntopics <- ncol(tm.result$topics)
  Y <- count.mat

  # cells.keep <- colnames(Y)
  # dat.meta <- as.data.frame(subset(dat.meta, cell %in% cells.keep))
  # # order dat.meta properly
  # rownames(dat.meta) <- dat.meta$cell
  # dat.meta <- dat.meta[cells.keep, ]

  # print("Size factor:")
  # print(head(size.factor))

  assertthat::assert_that(all(rownames(Y) == colnames(tm.result$terms)))

  if (bins.keep.manual == "none"){
    # keep variable bins
    if (bins.keep > 0){
      print(paste("Keeping only top bins bins.keep:", bins.keep))
      bins.high.i <- as.data.frame(apply(tm.result$terms, MARGIN = 1, function(jcol) order(jcol, decreasing = TRUE)[1:bins.keep])) %>%
        unlist()
      bins.high <- unique(rownames(Y)[bins.high.i])
    } else {
      print(paste("Keeping all bins because bins.keep=", bins.keep))
      bins.high <- colnames(tm.result$terms)
    }
    # bins.high.i <- as.data.frame(apply(tm.result$terms, MARGIN = 1, function(jcol) order(jcol, decreasing = TRUE)[1:bins.keep])) %>%
    #   unlist()
    # bins.high <- unique(rownames(Y)[bins.high.i])
  } else {
      print("Taking bins manually...")
    # take manually
      bins.high.i <- rownames(Y) %in% unique(bins.keep.manual)
      bins.high <- rownames(Y)[bins.high.i]
  }

  print("Dim before")
  print(dim(Y))
  Y.filt <- Y[bins.high, ]
  print("Dim after")
  print(dim(Y.filt))

  # remove empty cells
  print("Removing empty cells dim before")
  print(dim(Y.filt))
  empty.cells <- colSums(Y.filt) == 0
  Y.filt <- Y.filt[, !empty.cells]
  print("Removing empty cells dim after")
  print(dim(Y.filt))

  if (sz.cname == "none"){
    size.factor <- colSums(Y.filt)
  } else {
    size.factor <- dat.meta[[sz.cname]]
  }


  if (args$platecname == "none"){
      X.mat <- NULL
  } else {
      if (length(unique(dat.meta[[plate.cname]])) > 1){
          X.mat <- model.matrix(object = as.formula(paste0("~ ", plate.cname)), data = dat.meta)
      }  else {
          print("Number of factors in 'plate' is only 1, setting to NULL")
          X.mat <- NULL
      }
  }
  # var.vec <- dat.meta[[covar.cname]]
  # names(var.vec) <- dat.meta$cell
  # X <- data.frame(covariate = var.vec, cell = names(var.vec))  # columns of 1s are implicit
  # X.reorder <- X[match(colnames(Y.filt), X$cell), ]
  # X.mat <- matrix(data = X.reorder$covariate, ncol = 1, byrow = TRUE, dimnames = list(X.reorder$cell, covar.cname))

  # factors (U) is c by k
  # loadings (V) is g by k
  if (do.log){
    topics.mat <- log2(tm.result$topics)
    terms.mat <- log2(tm.result$terms[, bins.high])
  } else {
    topics.mat <- tm.result$topics
    terms.mat <- tm.result$terms[, bins.high]
  }

  if (svd.on.Yinit){
    # GLM loglink function for multinom is log( p / (1 - p) )
    # V %*% t(U) on init matrix gives an estimate of p
    # after estimate p / (1 - p), THEN remove mean and finally do SVD to get factors and loadings estimate
    p <- t(terms.mat) %*% t(topics.mat)
    logodds <- log(p / (1 - p))
    # remove mean and SVD
    logodds.centered <- t(scale(t(logodds), center = TRUE, scale = FALSE))
    # logodds.centered.check <- sweep(logodds, MARGIN = 1, STATS = rowMeans(logodds), FUN = "-")
    logodds.pca <- prcomp(t(logodds.centered), center = FALSE, scale. = FALSE, rank. = ntopics)
    U.init <- logodds.pca$x  # cells by k
    V.init <- logodds.pca$rotation  # genes by k, no need to transpose
  } else {
    U.init <- scale(topics.mat, center = TRUE, scale = FALSE)  # c by k
    V.init <- scale(terms.mat, center = TRUE, scale = FALSE)  # k by g
    V.init <- t(V.init)  # now its g by k
  }
  return(list(Y.filt = Y.filt, U.init = U.init, V.init = V.init, X.mat = X.mat, size.factor = size.factor, ntopics = ntopics))
}

suppressPackageStartupMessages(library("argparse"))
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(scchicFuncs)
library(hash)
library(igraph)
library(umap)
# library(DescTools)
library(glmpca)
# library(devtools)
# dev_mode(on = TRUE)
# devtools::install_github("willtownes/glmpca")
library(glmpca)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infpeaks', metavar='INFILE', help='LDA output')
parser$add_argument('-infmeta', metavar='INFILE', help='metadata containing plate colname and sizefactor colname')
parser$add_argument('-outbase', metavar='OUTBASE', help='GLMPCA base output, will add .pdf and .RData to this outbase to make outf')
parser$add_argument('-platecname', metavar='COLNAME', default='none', help='Colname used to set up X matrix. Set none to not include plate')
parser$add_argument('-sizefactorcname', metavar='cname', default = 'none', help='Colname from metadata to take sizefactors, default none, which takes size factors from raw count mat')
parser$add_argument('-bincutoff', metavar='BINCUTOFF', type = 'integer', help='Bin cutoff', default = 0)
parser$add_argument('-niter', metavar='N', type = 'integer', help='Number of iterations', default = 1000)
parser$add_argument('-tol', metavar='Float', type = 'double', help='Tolerance', default = 1e-8)
parser$add_argument('-penalty', metavar='Float', type = 'double', help='Penalty of the fit', default = 1.5)
parser$add_argument('-binskeep', metavar='N', type = 'integer', help='Nuber of bins to keep for GLMPCA. Set to 0 to keep all bins', default = 250)
parser$add_argument('-genesfile', metavar='INFILE', help='Input file with gene as column name for rows to keep. binskeep automatically set to 0', default = "none")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

if (endsWith(x = args$infmeta, suffix = ".RData")){
  load(args$infmeta, verbose=TRUE)  # dat.spikeins.mat
} else if (endsWith(x = args$infmeta, suffix = ".rds")){
  dat.spikeins.mat <- readRDS(args$infmeta)
  print("Loading from rds rownames shoudl be samp names")
  print(head(rownames(dat.spikeins.mat)))
} else if (endsWith(x = args$infmeta, suffix = ".txt")){
  dat.spikeins.mat <- as.data.frame(fread(args$infmeta))
  rownames(dat.spikeins.mat) <- dat.spikeins.mat$cell
  print("Peeking in dat spikeins mat")
  print(head(dat.spikeins.mat))
} else {
  stop("Not RData, rds of txt")
}
dat.meta <- dat.spikeins.mat
 
if (args$genesfile == "none"){
  jbins.keep.manual <- "none"
} else {
  jbins.keep.dat <- fread(args$genesfile)  # expect gene column name
  jbins.keep.manual <- unique(jbins.keep.dat$gene)
}

print(paste("Found", length(jbins.keep.manual), "genes to filter"))
print(head(jbins.keep.manual))

assertthat::assert_that(length(jbins.keep.manual) > 0)

# Constants ---------------------------------------------------------------

# hubprefix <- "/hpc/hub_oudenaarden"

# jcovar.cname <- "cell.var.within.sum.norm.log2.CenteredAndScaled"
jcovar.cname <- args$platecname
jsize.cname <- args$sizefactorcname

# bincutoff <- 10000
bincutoff <- args$bincutoff

# niter <- 1000
niter <- args$niter
# jtol <- 1e-8
jtol <- args$tol
jbins.keep <- args$binskeep
# calculating var raw
# jpenalty <- 1.5
jpenalty <- args$penalty


# ldadate <- "2020-02-11"
# jdate <- "2020-08-26"
# jsuffix <- "KeepBestPlates2.LDAfromPeaks"

# outdir <- args$outdir
# dir.create(outdir)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

  # Setup output paths ------------------------------------------------------
  
  # outbase <- paste0("PZ_", jmark, ".", jcond, ".", jsuffix, ".GLMPCA_var_correction.binskeep_", jbins.keep, ".covar_", jcovar.cname, ".niter_", niter, ".tol_", jtol, ".", jdate)
  outbase <- args$outbase  # nclude full direcotry path
  outdircheck <- dirname(outbase)
  assertthat::assert_that(dir.exists(outdircheck))

  outname <- paste0(outbase, ".RData")
  outname.pdf <- paste0(outbase, ".pdf")
  outf <- outname
  outf.pdf <- outname.pdf
  
  # print(file.exists(outf))
  # assertthat::assert_that(!file.exists(outf))
  # assertthat::assert_that(!file.exists(outf.pdf))
  
  # if (file.exists(outf)){
  #   return(NULL)
  # }
  
  # Load data  --------------------------------------------------------------
  
  inf.peaks <- args$infpeaks
  assertthat::assert_that(file.exists(inf.peaks))
  
  
  # load raw counts from peaks 
  load(inf.peaks, v=T)
  count.mat.peaks <- as.matrix(count.mat)
  # remove bins that have zero variance
  genevars <- apply(count.mat.peaks, 1, var)
  geneskeep <- rownames(count.mat.peaks)[which(genevars > 0)]
  count.mat.peaks <- count.mat.peaks[geneskeep, ]
  tm.result.peaks <- AddTopicToTmResult(posterior(out.lda), jsep = "_")
  tm.result.peaks$terms <- tm.result.peaks$terms[, geneskeep]

  # rearrange meta dat
  dat.meta <- subset(dat.meta, cell %in% colnames(count.mat.peaks))
  dat.meta <- dat.meta[colnames(count.mat.peaks), ]
  
  # get distances
  peaks.all <- colnames(tm.result.peaks$terms)
  coords.all <- sapply(peaks.all, function(p) strsplit(p, ";")[[1]][[1]])
  peaks.dist <- data.frame(rname = peaks.all, 
                           coord = coords.all, 
                           chromo = sapply(coords.all, GetChromo),
                           jstart = as.numeric(sapply(coords.all, GetStart)),
                           jend = as.numeric(sapply(coords.all, GetEnd)), 
                           stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(jdist = jend - jstart)
  
  peaks.keep <- subset(peaks.dist, jdist > bincutoff)$rname
  print(paste("Filtering peaks of length less than", bincutoff))
  print("NPeaks before...")
  print(length(peaks.all))
  print("NPeaks after...")
  print(length(peaks.keep))
  
  tm.result.peaks$terms <- tm.result.peaks$terms[, peaks.keep]
  
  print("dim count.mat.peaks before...")
  print(dim(count.mat.peaks))
  
  count.mat.peaks <- count.mat.peaks[peaks.keep, ]
  
  print("dim count.mat.peaks after...")
  print(dim(count.mat.peaks))
  
  topics.mat.peaks <- tm.result.peaks$topics
  
  
  # start -----------------
  pdf(outf.pdf, useDingbats = FALSE)
  
  # print(jmark)
  print("Current time elapsed:")
  print(Sys.time() - jstart)
  
  # Plot it all -------------------------------------------------------------
  
  umap.out <- umap(topics.mat.peaks, config = jsettings)
  dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
  dat.umap.long <- DoLouvain(topics.mat.peaks, jsettings, dat.umap.long)
  # cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#41f0d9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115")
  ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  
  # Plot variance -----------------------------------------------------------
  
  # m.plates <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = plate)) + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #   geom_point()
  # print(m.plates)
  
  dev.off()
  
  # set up GLMPCA
  
  print("Running GLMPCA")


  
  glm.inits <- InitGLMPCAfromLDA.plate(count.mat = count.mat.peaks, tm.result = tm.result.peaks, dat.meta = dat.meta, plate.cname = jcovar.cname, bins.keep = jbins.keep, do.log = FALSE, svd.on.Yinit = TRUE, sz.cname = jsize.cname, bins.keep.manual = jbins.keep.manual)

  # set up X independently 
  
  system.time(
    glm.out <- glmpca(Y = glm.inits$Y.filt, L = glm.inits$ntopics, 
                      fam = "poi", 
                      ctl = list(maxIters = niter, tol = jtol, penalty = jpenalty, verbose=TRUE), 
                      init = list(factors = glm.inits$U.init, loadings = glm.inits$V.init), 
                      # minibatch = "none",
                      # optimizer = "fisher", 
                      minibatch = "stochastic",
                      optimizer = "avagrad", 
                      X = glm.inits$X.mat, Z = NULL, sz = glm.inits$size.factor)
  )
  print(traceback())
  save(glm.out, glm.inits, dat.meta, file = outf)
  
  print(Sys.time() - jstart)
  
