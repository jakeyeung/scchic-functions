# Jake Yeung
# run_LDA_model2.R
# 2019-11-10
# DESCRIPTION
# 
#     Run LDA clean up params
# 
# FOR HELP
# 
#     Rscript run_LDA_model2.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-11-10
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

jstart <- Sys.time() 
suppressPackageStartupMessages(library("argparse"))
library(topicmodels)
library(dplyr)
library(ggplot2)
library(Matrix)
library(here)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument('inpath', metavar='INFILE',
                                            help='.RData where count mat is in count.dat$counts or .rds object to count mat')
parser$add_argument('outdir', metavar='OUTDIR',
                                            help='Out directory')
parser$add_argument("-t", "--topics", metavar='Comma sep string', required=TRUE,
                                            help='CSV of topics to iterate')
parser$add_argument("-b", "--binarizemat", action="store_true", default=FALSE,
                        help="Binarize matrix")
parser$add_argument("-n", "--projname", metavar='Name of project', default="MyProj",
                        help="Name of project for naming pdf and Robj output")
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

setwd(here())
print(paste("Work directory: ", getwd()))
source("scripts/Rfunctions/ParseStrings.R")
source("scripts/Rfunctions/Aux.R")

# args <- commandArgs(trailingOnly=TRUE)

print("input args:")
print(args)

inpath <- args$inpath
outdir <- args$outdir
topic.vec <- as.numeric(StrToVector(args$topics, delim = ","))
binarizemat <- args$binarizemat
projname <- args$projname  # helps with writing pdf and Robj output

print(paste("Will iterate through", length(topic.vec), "Ks"))
print(topic.vec)

plotpath <- file.path(outdir, paste0("plots.", projname, ".pdf"))
outpath <- file.path(outdir, paste0("ldaOut.", projname, ".Robj"))

# Load counts -------------------------------------------------------------

if (endsWith(inpath, ".rds")){
  count.mat <- readRDS(inpath)
} else {
  # assume .RData object with count.dat$counts
  load(inpath, v=T)
  count.mat <- count.dat$counts
}

print(dim(count.mat))

# Plot mean and variance --------------------------------------------------

dat.meanvar <- data.frame(Sum = Matrix::rowSums(count.mat), 
                          Mean = Matrix::rowMeans(count.mat),
                          Var = apply(count.mat, 1, var),
                          peak = rownames(count.mat),
                          stringsAsFactors=FALSE)
dat.meanvar <- dat.meanvar %>%
  rowwise() %>%
  mutate(CV = sqrt(Var) / Mean,
         peaksize = GetPeakSize(peak))

p1 <- ggplot(dat.meanvar, aes(x = peaksize)) + geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


p2 <- ggplot(dat.meanvar, aes(x = log10(Mean), y = log10(CV), size = peaksize)) + geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = -0.5)

p3 <- ggplot(dat.meanvar, aes(x = peaksize, y = Sum)) + geom_point(alpha = 0.1) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_log10()

# save plots
pdf(plotpath, useDingbats = FALSE)
  print(p1)
  print(p2)
  print(p3)
dev.off()

# Run LDA on count matrix -------------------------------------------------

print("Running LDA")

# binarize matrix
if (binarizemat){
  count.mat.orig <- count.mat
  count.mat <- BinarizeMatrix(count.mat)
  print(paste('Max count after binarizing', max(count.mat)))
} else {
  count.mat.orig <- NA
}

if (length(topic.vec) > 1){
  print("Running multicore LDA for topics:")
  print(topic.vec)
  out.lda <- parallel::mclapply(topic.vec, function(nc){
          LDA(x = t(count.mat), k = nc, method = "Gibbs", control=list(seed=0))            
  }, mc.cores = length(topic.vec))
} else {
  print("Running single LDA for topics:")
  print(topic.vec)
  out.lda <- LDA(x = t(count.mat), k = topic.vec, method = "Gibbs", control=list(seed=0))
}

# save output
print("Saving LDA")
save(out.lda, count.mat, count.mat.orig, file = outpath)
print("Time elapsed after LDA")
print(Sys.time() - jstart)


print("Time elapsed after tuning")
print(Sys.time() - jstart)

