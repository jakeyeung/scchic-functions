# Jake Yeung
# make_merged_exprs_mat_for_MARA.args.R
# 2020-08-25
# DESCRIPTION
# 
#     prepare exprs matrix for MARA from inputed LDA 
# 
# FOR HELP
# 
#     Rscript make_merged_exprs_mat_for_MARA.args.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-08-25
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

jstart <- Sys.time()  

suppressPackageStartupMessages(library("argparse"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(topicmodels)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Iput LDA object (.Robj or .RData)')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output matrix txt of imputed values')
parser$add_argument('-keepNbins', metavar='N bins', type = 'integer', default=0, 
                                            help='Number of bins to keep for each topic')
parser$add_argument('--AddChr', action="store_true", default=FALSE,
                                            help='Add chr to front of rownames')
parser$add_argument("--NoNormalize", action="store_true", 
                        help="Do not normalize matrix. If you do this you have to add an intercept in your model")
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

# Load data  --------------------------------------------------------------

# jmark <- "H3K4me1"
# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

inf <- args$infile
load(inf, v=T)  # contains count.mat and out.lda
tm.result <- posterior(out.lda)
dat.impute <- t(tm.result$topics %*% tm.result$terms)
# now we have
# dat.impute = dat.impute, count.mat = count.mat, tm.result = tm.result

# dat.imputes <- lapply(dat.out, function(x) x$dat.impute)
# dat.mats <- lapply(dat.out, function(x) x$dat.impute)
# 
# dat.imputes.merge <- do.call(cbind, dat.imputes)
# dat.mats.merge <- do.call(cbind, dat.mats)

tm.result <- AddTopicToTmResult(tm.result)
ntopics <- nrow(tm.result$topics)
ntopics.vec <- rownames(tm.result$terms); names(ntopics.vec) <- ntopics.vec
if (args$keepNbins > 0){
  # for each topic, filter out top N binsA
  terms.keep.lst <- lapply(ntopics.vec, function(jtop){
    loadings <- sort(tm.result$terms[jtop, ], decreasing=TRUE)
    loadings.keep.tmp <- loadings[1:args$keepNbins]
    terms.keep.tmp <- names(loadings.keep.tmp)
    return(terms.keep.tmp)
  })
  terms.keep <- unique(unlist(terms.keep.lst))
} else {
  terms.keep <- colnames(tm.result$terms)
}

# Center and rorate -------------------------------------------------------

# GLM loglink function for multinom is log( p / (1 - p) )
# V %*% t(U) on init matrix gives an estimate of p
# after estimate p / (1 - p), THEN remove mean and finally do SVD to get factors and loadings estimate
# p <- dat.imputes.merge

print("Dim before:")
print(dim(dat.impute))
p <- dat.impute[terms.keep, ]
print("Dim after filtering terms:")
print(dim(p))

logodds <- log(p / (1 - p))
# remove mean and SVD
logodds.centered <- t(scale(t(logodds), center = TRUE, scale = FALSE))
# # logodds.centered.check <- sweep(logodds, MARGIN = 1, STATS = rowMeans(logodds), FUN = "-")
# logodds.pca <- prcomp(t(logodds.centered), center = FALSE, scale. = FALSE, rank. = ntopics)
# U.init <- logodds.pca$x  # cells by k
# V.init <- logodds.pca$rotation  # genes by k, no need to transpose

# make nice
geneids <- sapply(rownames(logodds), function(x) strsplit(x, ";")[[1]][[1]])

if (args$AddChr){
    print("Adding chr:")
    print("geneids before")
    print(head(geneids))
    geneids <- paste("chr", geneids, sep = "")
    print("geneids after")
    print(head(geneids))
}

logodds.out <- data.frame(Gene.ID = geneids, logodds, stringsAsFactors = FALSE)
logodds.centered.out <- data.frame(Gene.ID = geneids, logodds.centered, stringsAsFactors = FALSE)

if (!args$NoNormalize){
    fwrite(logodds.centered.out, file = args$outfile, sep = "\t", col.names = TRUE, row.names = FALSE)
} else {
    fwrite(logodds.out, file = args$outfile, sep = "\t", col.names = TRUE, row.names = FALSE)
}

print(Sys.time() - jstart)

# # write logodds to output
# outname <- paste0("ldaOut.marks_merged.BM_AllMerged.merged_by_clusters_no_NAs.K-30.txt")
# outdir <- "/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/count_mats_peaks_norm_merged"
# outdir2 <- "/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/count_mats_peaks_unnorm_LDA_merged"
# outtxt <- file.path(outdir, outname)
# outtxt2 <- file.path(outdir2, outname)
# 
# fwrite(logodds.centered.out, file = outtxt, sep = "\t", col.names = TRUE, row.names = FALSE)
# fwrite(logodds.out, file = outtxt2, sep = "\t", col.names = TRUE, row.names = FALSE)
