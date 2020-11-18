# Jake Yeung
# Date of Creation: 2020-08-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/1-correct_LDA_from_peaks_by_variance.MouseBM.poisson_faster.R
# 

rm(list=ls())
jstart <- Sys.time()

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
parser$add_argument('-outbase', metavar='OUTBASE', help='GLMPCA base output, will add .pdf and .RData to this outbase to make outf')
parser$add_argument('-cname', metavar='COLNAME', default='cell.var.within.sum.norm.log2.CenteredAndScaled', help='Colname used to regress out')
parser$add_argument('-bincutoff', metavar='BINCUTOFF', type = 'integer', help='Bin cutoff', default = 0)
parser$add_argument('-niter', metavar='N', type = 'integer', help='Number of iterations', default = 1000)
parser$add_argument('-tol', metavar='Float', type = 'double', help='Tolerance', default = 1e-8)
parser$add_argument('-penalty', metavar='Float', type = 'double', help='Penalty of the fit', default = 1.5)
parser$add_argument('-binskeep', metavar='N', type = 'integer', help='Nuber of bins to keep for GLMPCA. Set to 0 to keep all bins', default = 250)
parser$add_argument('-chromos', metavar='SpaceSepString', help='Space separated string of chromosomes for calculating intrachr var', nargs='+')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")
parser$add_argument("--winsorize", action="store_true",
                        help="Winsorize yes or not")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}


# Constants ---------------------------------------------------------------

hubprefix <- "/hpc/hub_oudenaarden"

# jcovar.cname <- "cell.var.within.sum.norm.log2.CenteredAndScaled"
jcovar.cname <- args$cname
winsorize <- args$winsorize

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

# jchromos <- paste("chr", c(seq(19)), sep = "")
jchromos <- args$chromos

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
  
  dat.impute.log <- log2(t(tm.result.peaks$topics %*% tm.result.peaks$terms))
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos) %>%
    rowwise() %>%
    mutate(plate = ClipLast(as.character(cell), jsep = "_"))
  
  dat.var.merge <- left_join(dat.umap.long, dat.var)
  
  m.plates <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = plate)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point()
  
  m.var <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c(direction = -1)
  
  m.var.plates <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c(direction = -1) + facet_wrap(~plate)
  
  # PCA on topics, show variance
  pca.out <- prcomp(topics.mat.peaks, center = TRUE, scale. = TRUE)
  dat.pca <- data.frame(cell = rownames(pca.out$x), PC1 = pca.out$x[, 1], PC2 = pca.out$x[, 2], stringsAsFactors = FALSE) %>%
    left_join(dat.var.merge)
  
  m.pca.var <- ggplot(dat.pca, aes(x = PC1, y = PC2, color = cell.var.within.sum.norm)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c(direction = -1)
  
  m.pca.var.plates <- ggplot(dat.pca, aes(x = PC1, y = PC2, color = cell.var.within.sum.norm)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c(direction = -1) + facet_wrap(~plate)
  
  print(m.plates)
  print(m.var)
  print(m.var.plates)
  print(m.pca.var)
  print(m.pca.var.plates)
  
  
  # Calculate raw varaince and compare with imputed variance  ---------------
  
  
  # dat.var.raw <- CalculateVarRaw(count.mat, merge.size = mergesize, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE)
  # # center and scale ncuts.var
  # dat.var.raw$ncuts.var.log2 <- log2(dat.var.raw$ncuts.var)
  # dat.var.raw$ncuts.var.CenteredAndScaled <- (dat.var.raw$ncuts.var - mean(dat.var.raw$ncuts.var)) / sd(dat.var.raw$ncuts.var)
  # dat.var.raw$ncuts.var.log2.CenteredAndScaled <- (dat.var.raw$ncuts.var.log2 - mean(dat.var.raw$ncuts.var.log2)) / sd(dat.var.raw$ncuts.var.log2)
  
  # dat.merge2 <- left_join(dat.var.merge, dat.var.raw)
  dat.merge2 <- dat.var.merge
  dat.merge2$cell.var.within.sum.norm.CenteredAndScaled <- (dat.merge2$cell.var.within.sum.norm - mean(dat.merge2$cell.var.within.sum.norm)) / sd(dat.merge2$cell.var.within.sum.norm)
  dat.merge2$cell.var.within.sum.norm.log2 <- log2(dat.merge2$cell.var.within.sum)
  
  dat.merge2$cell.var.within.sum.norm.log2.CenteredAndScaled <- (dat.merge2$cell.var.within.sum.norm.log2 - mean(dat.merge2$cell.var.within.sum.norm.log2)) / sd(dat.merge2$cell.var.within.sum.norm.log2)
  
  # winsorize?
  if (winsorize){
    dat.merge2[[jcovar.cname]] <- DescTools::Winsorize(dat.merge2[[jcovar.cname]], probs = c(0.01, 0.99))
  }
  
  # m1 <- ggplot(dat.merge2, aes(x = ncuts.var, y = cell.var.within.sum.norm)) + geom_point() + 
  #   scale_x_log10() + scale_y_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # m2 <- ggplot(dat.merge2, aes(x = ncuts, y = ncuts.var)) + geom_point() + 
  #   scale_x_log10() + scale_y_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m3 <- ggplot(dat.merge2, aes_string(x = jcovar.cname, y = "cell.var.within.sum.norm")) + geom_point() + scale_y_log10() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # print(m1)
  # print(m2)
  print(m3)
  
  dev.off()
  
  # if (file.exists(outf)){
  #   print(paste("Outf exists, skipping...", outf))
  #   return(NULL)
  # }
  
  # set up GLMPCA
  
  print("Running GLMPCA")
  
  glm.inits <- InitGLMPCAfromLDA(count.mat.peaks, tm.result.peaks, dat.merge2, covar.cname = jcovar.cname, bins.keep = jbins.keep, do.log = FALSE, svd.on.Yinit = TRUE, use.orig.sz = TRUE)
  
  
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
  
  save(glm.out, glm.inits, dat.merge2, file = outf)
  
  print(Sys.time() - jstart)
  
