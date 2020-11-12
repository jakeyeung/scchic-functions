# Jake Yeung
# make_clusters_from_topics_mixturemodel.R
# 2020-04-18
# DESCRIPTION
# 
#     Get clusters from mixture model
# 
# FOR HELP
# 
#     Rscript make_clusters_from_topics_mixturemodel.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-04-18
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(mixtools)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(igraph)
library(umap)
library(topicmodels)
library(scchicFuncs)

library(ggrepel)
library(forcats)

# # for annotating bin to gene
# if (args$species == "mouse"){
#   library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#   library(org.Mm.eg.db)
# } else {
#   print(paste("Species:", args$species, "not yet coded"))
# }
library(ChIPseeker)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infile', metavar='INFILE',
                                            help='lda output .RData (expects out.lda)')
parser$add_argument('-outprefix', metavar='OUTPREFIX include dir',
                                            help='Output prefix for output .RData and .pdf objects')
parser$add_argument('-annotfile', metavar='CSVFILE', 
                                            help='Tab-delimited file with columns cell and cluster which defines cells into each cluster')
# parser$add_argument('-topicskeep', metavar='List of numbers', nargs = "+",
#                                             help='list of topics to keep (prefixes add later? e.g. -topicskeep 1 2 3 4)')
parser$add_argument('-mark', metavar='HistoneMarkName', 
                                            help='Name of histone mark')
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

# Contants ----------------------------------------------------------------

jsep <- ""  # separation between topic and number

# Get clusters from topics
dat.annot <- fread(args$annotfile)

topics.keep <- unique(dat.annot$cluster)
topics.keep <- topics.keep[!is.na(topics.keep)]
names(topics.keep) <- topics.keep

# topics.keep <- paste("topic", args$topicskeep, sep = jsep)
# topics.names <- topics.keep
# names(topics.keep) <- topics.names

# these settings alter the clustering and UMAP
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")  # let's hope you don't have more clusters than cbPalettes... otherwise need to recycle these 

# Load LDA outputs for single stains --------------------------------------------------------

jmark <- args$mark

inf.lda <- args$infile
assertthat::assert_that(file.exists(inf.lda))

bname <- ClipLast(basename(inf.lda), jsep = "\\.", jsep.out = ".")
outprefix <- args$outprefix

print(paste("Writing output files using this prefix:", outprefix))

# Load data ---------------------------------------------------------------

outf <- paste0(outprefix, ".RData")
outpdf <- paste0(outprefix, ".pdf")

load(inf.lda, v=T)
tm.result <- posterior(out.lda)

# dont separate with _
colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = jsep)
rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = jsep)

topics.mat <- tm.result$topics

dat.umap.long <- DoUmapAndLouvain(topics.mat, jsettings)
  # left_join(., annots.dat)  # in general case dont add annots

dat.umap.long <- dat.umap.long %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

# add var
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
(jchromos <- paste("chr", seq(19), sep = ""))
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

print("Check datumaplong and datvar")
print(head(dat.umap.long))
print(head(dat.var))
dat.umap.long <- left_join(dat.umap.long, dat.var)

print("Check datumaplong again")
print(head(dat.umap.long))
dat.umap.long.merge.tmp <- left_join(dat.umap.long, dat.annot %>% dplyr::select(c(cell, cluster)))

dat.umap.long$mark <- jmark

m.louv <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + facet_wrap(~plate) + ggtitle(jmark)

m.cluster <- ggplot(dat.umap.long.merge.tmp, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + facet_wrap(~plate) + ggtitle(jmark)

m.var <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~plate) + ggtitle(jmark)


# Define topics  ----------------------------------------------------------


# plot and make objects

pdf(outpdf, useDingbats = FALSE)

print(m.cluster)
print(m.louv)
print(m.var)

jthres <- 0.5
# label cells by topic
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

mm.celltype.lst <- lapply(topics.keep, function(jtopic){
  print(jtopic)
  cells.keep <- subset(dat.annot, cluster == jtopic)$cell
  # tvec.raw <- sort(tm.result$topics[, jtopic])
  # # transform
  # tvec <- log(tvec.raw / (1 - tvec.raw))
  # # xline <- quantile(tvec, )
  # mm <- normalmixEM(x = tvec, lambda = c(0.9, 0.1), mu = c(-5, -1), sigma = c(2, 1), k = 2)
  # (xline <- min(mm$x[which(mm$posterior[, 1] < jthres)]))
  # xcells <- names(tvec)[which(tvec > xline)]
  # print(paste(length(xcells), "/", length(tvec), "assigned to", jtopic))
  # # get topic value for assigned cells
  # tvec.raw.filt <- tvec.raw[xcells]

  # (xline <- min(mm$x[which(mm$posterior[, 1] < jthres)]))
  # # xline <- max(mm$x[indx.btwn][which(post.filt[, 2] > 0.5)])
  # plot(density(tvec), main = paste(jtopic, jthres), xlab = "Log Odds [log(p / (1 - p))]")
  # abline(v = xline, col = 'blue')
  # plot.mixEM(mm, whichplots = 2, xlab2 = "Log Odds [log(p / (1 - p))]", main2 = paste(jtopic, jthres))
  # abline(v = xline, col = 'blue')

  # cells.keep <- xcells
  m.check <- PlotXYWithColor(dat.umap.long %>% mutate(is.celltype = cell %in% cells.keep), xvar = "umap1", yvar = "umap2", cname = "is.celltype", jtitle = paste(jtopic, jthres), cont.color = FALSE, col.palette = cbPalette, use.ggrastr=FALSE)
  print(m.check)
  return(list(topic = jtopic, cells.keep = cells.keep))
  # return(list(topic = jtopic, topic.weight = tvec.raw.filt, celltype = xcells, mm = mm, threshold = xline))
})

celltypes <- lapply(names(mm.celltype.lst), function(ctype){
  data.frame(cluster = ctype, cell = mm.celltype.lst[[ctype]]$cells.keep,
             topic = mm.celltype.lst[[ctype]]$topic,
             stringsAsFactors = FALSE)
}) %>%
  bind_rows()

# # handle replicates
# celltypes.dedup <- celltypes %>%
#   group_by(cell) %>%
#   filter(topic.weight == max(topic.weight))
# 
# assertthat::assert_that(!all(duplicated(celltypes.dedup)))
celltypes.dedup <- celltypes

print("Check celltypes datumap celltypes")
print(head(celltypes))
print(head(dat.umap.long))
print(head(celltypes.dedup))
dat.merge <- left_join(dat.umap.long, celltypes.dedup)

dat.merge$plate <- sapply(dat.merge$cell, function(x) ClipLast(x, jsep = "_"))

# handle duplicates

# show final
m.final <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = topic)) + geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = cbPalette, na.value = "grey85")
print(m.final)


dev.off()

# save annotated dat to output
save(dat.merge, count.mat, mm.celltype.lst, out.lda, file = outf)


