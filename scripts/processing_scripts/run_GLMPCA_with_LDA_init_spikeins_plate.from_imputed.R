# Jake Yeung
# Date of Creation: 2020-08-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/1-correct_LDA_from_peaks_by_variance.MouseBM.poisson_faster.R
# 

rm(list=ls())
jstart <- Sys.time()

InitGLMPCAfromImputed.plate <- function(count.mat, mat.imputed, dat.meta, plate.cname = "ncuts.var", sz.cname = "none", ntopics=30){

  # ntopics <- ncol(tm.result$topics)
  ntopics <- 30

  # set rownames same as mat.imputed
  # Y <- count.mat[rnames.keep, ]
  Y <- count.mat

  bins.high <- rownames(mat.imputed)
  cells.keep <- colnames(mat.imputed)

  print("Filter bins and cells")

  print("Dim before")
  print(dim(Y))
  Y.filt <- Y[bins.high, cells.keep]
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
  # factors (U) is c by k
  # loadings (V) is g by k

  # GLM loglink function for multinom is log( p / (1 - p) )
  # V %*% t(U) on init matrix gives an estimate of p
  # after estimate p / (1 - p), THEN remove mean and finally do SVD to get factors and loadings estimate
  # p <- mat.imputed
  # logodds <- log(p / (1 - p))
  logodds <- mat.imputed  # alreay in -Inf Inf
  # remove mean and SVD
  logodds.centered <- t(scale(t(logodds), center = TRUE, scale = FALSE))
  # logodds.centered.check <- sweep(logodds, MARGIN = 1, STATS = rowMeans(logodds), FUN = "-")
  logodds.pca <- prcomp(t(logodds.centered), center = FALSE, scale. = FALSE, rank. = ntopics)
  U.init <- logodds.pca$x  # cells by k
  V.init <- logodds.pca$rotation  # genes by k, no need to transpose
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

parser$add_argument('-infimputed', metavar='INFILE', help='RData containing mat.adj or rds')
parser$add_argument('-infcounts', metavar='INFILE', help='RData containing count.mat or rds')
parser$add_argument('-infmeta', metavar='INFILE', help='metadata containing plate colname and sizefactor colname')
parser$add_argument('-outbase', metavar='OUTBASE', help='GLMPCA base output, will add .pdf and .RData to this outbase to make outf')
parser$add_argument('-platecname', metavar='COLNAME', default='none', help='Colname used to set up X matrix. Set none to not include plate')
parser$add_argument('-sizefactorcname', metavar='cname', default = 'none', help='Colname from metadata to take sizefactors, default none, which takes size factors from raw count mat')
parser$add_argument('-niter', metavar='N', type = 'integer', help='Number of iterations', default = 1000)
parser$add_argument('-jntopics', metavar='Topics', type = 'integer', help='Number of topics from matadj', default = 30)
parser$add_argument('-tol', metavar='Float', type = 'double', help='Tolerance', default = 1e-8)
parser$add_argument('-penalty', metavar='Float', type = 'double', help='Penalty of the fit', default = 1.5)
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
 

# Constants ---------------------------------------------------------------


jcovar.cname <- args$platecname
jsize.cname <- args$sizefactorcname


niter <- args$niter
jtol <- args$tol
# calculating var raw
jpenalty <- args$penalty



  # Setup output paths ------------------------------------------------------
  
  # outbase <- paste0("PZ_", jmark, ".", jcond, ".", jsuffix, ".GLMPCA_var_correction.binskeep_", jbins.keep, ".covar_", jcovar.cname, ".niter_", niter, ".tol_", jtol, ".", jdate)
  outbase <- args$outbase  # nclude full direcotry path
  outdircheck <- dirname(outbase)
  assertthat::assert_that(dir.exists(outdircheck))

  outname <- paste0(outbase, ".RData")
  outf <- outname
  
  # Load data  --------------------------------------------------------------
  
  inf.counts <- args$infcounts
  assertthat::assert_that(file.exists(inf.counts))

  inf.imputed <- args$infimputed
  assertthat::assert_that(file.exists(inf.imputed))
  
  # load raw counts from peaks 
  if (endsWith(x = inf.counts, suffix = ".RData")){
    load(inf.counts, v=T)
    count.mat.peaks <- as.matrix(count.mat)
  } else {
    count.mat.peaks <- as.matrix(readRDS(inf.counts))
  }

  if (endsWith(x = inf.imputed, suffix = ".RData")){
    load(inf.imputed, v=T)
    mat.imputed <- mat.adj
  } else {
    mat.imputed <- as.matrix(readRDS(inf.imputed))
  }
  # # load mat.adj
  # load(args$infimputed, v=T)  # mat.adj
  # mat.imputed <- mat.adj

  # set up GLMPCA
  
  print("Running GLMPCA")


  
  glm.inits <- InitGLMPCAfromImputed.plate(count.mat = count.mat.peaks, mat.imputed = mat.imputed, dat.meta = dat.meta, plate.cname = jcovar.cname, sz.cname = jsize.cname, ntopics = args$jntopics)

  print("Glm initialized")
  print(lapply(glm.inits, dim))

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
  
