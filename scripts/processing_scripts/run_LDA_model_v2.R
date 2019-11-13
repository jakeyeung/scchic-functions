# Jake Yeung
# run_LDA_model2.R
# Clean up LDA after several months of usage we do things differently

jstart <- Sys.time() 

library(topicmodels)
library(dplyr)
library(ggplot2)
# library(ldatuning)  # dont do LDA tuning anymore, hard to find easy consensus
library(Matrix)
library(here)

setwd(here())

print(paste("Work directory: ", getwd()))

source("scripts/Rfunctions/ParseStrings.R")
source("scripts/Rfunctions/Aux.R")

args <- commandArgs(trailingOnly=TRUE)

print("input args:")
print(args)

inpath <- args[[1]]
outdir <- args[[2]]
# nclst <- StrToNumeric(args[[3]])
topic.vec <- as.numeric(StrToVector(args[[3]], delim = ","))
binarizemat <- StrToBool(args[[4]])
projname <- args[[5]]  # helps with writing pdf and Robj output
# tunemodels <- StrToBool(args[[5]])  # input a topic vec directly
# remove bad cells and peaks BEFORE hand 
# meanmax <- StrToNumeric(args[[6]])  # remove suspicious peaks 
# cellmin <- StrToNumeric(args[[7]])  # remove cells with low counts
# cellmax <- StrToNumeric(args[[8]])  # remove suspiciious cells 

print(paste("Will iterate through", length(topic.vec), "Ks"))
print(topic.vec)

plotpath <- file.path(outdir, paste0("plots.", projname, ".pdf"))
outpath <- file.path(outdir, paste0("ldaOut.", projname, ".Robj"))

# print args
print("Args:")
print(paste(inpath, outdir, nclst, topic.vec, tunemodels, cellmin, cellmax, binarizemat))

# Load counts -------------------------------------------------------------

load(inpath, v=T)

count.mat <- count.dat$counts

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

# suspicious peaks
bad.peaks <- dat.meanvar %>%
  filter(Mean > meanmax)
print(head(bad.peaks$peak))
print(paste("There are", nrow(bad.peaks), "peaks with more than", meanmax, "counts. Removing them..."))

# filter them out before running count.mat
print("Dimensions before filtering peaks...")
print(dim(count.mat))

count.mat <- count.mat[which(!rownames(count.mat) %in% bad.peaks$peak), ]
print("Dimensions after filtering peaks...")
print(dim(count.mat))

# remove peaks in chrM
# M.peaks <- grep("chrM|chrY|chrX", dat.meanvar$peak, value=TRUE)
M.peaks <- grep("chrM", dat.meanvar$peak, value=TRUE)

count.mat <- count.mat[which(!rownames(count.mat) %in% M.peaks), ]
print("Dimensions after filtering peaks M chromos")
print(dim(count.mat))


# Remove cells with zero entries
print("Dimensions before filtering cells")
print(dim(count.mat))
count.mat <- count.mat[, which(Matrix::colSums(count.mat) > cellmin)]
count.mat <- count.mat[, which(Matrix::colSums(count.mat) < cellmax)]
print("Dimensions after filtering cells...")
print(dim(count.mat))

# Run LDA on count matrix -------------------------------------------------

# nclst <- 10
print("Running LDA")

# binarize matrix
if (binarizemat){
  count.mat.orig <- count.mat
  count.mat <- BinarizeMatrix(count.mat)
  print(paste('Max count after binarizing', max(count.mat)))
} else {
  count.mat.orig <- NA
}

# print(head(count.mat[1:5, 1:5]))
if (!tunemodels){
  print("Running single LDA for topics:")
  print(nclst)
  out.lda <- LDA(x = t(count.mat), k = nclst, method = "Gibbs", control=list(seed=0))
} else {
  print("Running multicore LDA for topics:")
  print(topic.vec)
  out.lda <- parallel::mclapply(topic.vec, function(nc){
          LDA(x = t(count.mat), k = nc, method = "Gibbs", control=list(seed=0))            
  }, mc.cores = length(topic.vec))
}

# save output
print("Saving LDA")
save(out.lda, count.mat, count.mat.orig, file = outpath)
print("Time elapsed after LDA")
print(Sys.time() - jstart)

# tune LDA 

# if (tunemodels){
#     # topic.vec <- c(4, 9, 11, 14, 16, 18)
#     print("Running tuning")
#     optimal.topics <- FindTopicsNumber(t(as.matrix(count.mat)), topics=topic.vec, mc.cores = length(topic.vec), method="Gibbs", metrics=c("Arun2010", "CaoJuan2009", "Griffiths2004", "Deveaud2014"), control = list(seed=0))
#     FindTopicsNumber_plot(optimal.topics)
# 
#     save(optimal.topics, file = tunepath)
#     print("Saving optimal topics")
# }

print("Time elapsed after tuning")
print(Sys.time() - jstart)

