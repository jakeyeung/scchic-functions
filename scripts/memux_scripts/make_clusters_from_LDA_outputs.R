# Jake Yeung
# make_clusters_from_LDA_outputs.R
# 2020-04-04
# DESCRIPTION
# 
#     After LDA looks good, cluster single stains for training data
# 
# FOR HELP
# 
#     Rscript make_clusters_from_LDA_outputs.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-04-04
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infile', metavar='INFILE',
                                            help='lda output .RData (expects out.lda)')
parser$add_argument('-outprefix', metavar='OUTPREFIX include dir',
                                            help='Output prefix for output .RData and .pdf objects')
parser$add_argument('-mark', metavar='HistoneMarkName', 
                                            help='Name of histone mark')
parser$add_argument('-topn', metavar='NaturalNmber', type="integer", default=150,
                                            help='Number of bins to plot for each topic. Default 150 is good usually')
parser$add_argument("--WriteTopicLoadings", action="store_true", default=FALSE,
                        help="Write topic loadings to pdf (requires inftss, only available for certain species (Mouse)")
parser$add_argument('-inftss', metavar='Path to TSS gene annotations', required = FALSE, 
                                            help='Path to TSS annotations for annotating bin to gene, only relevant if WriteTopicLoadings is true. Example path: /hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed')
parser$add_argument('-species', metavar='Species', required = FALSE, default = "mouse",
                                            help='For now mouse')
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

# for annotating bin to gene
if (args$species == "mouse"){
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
} else {
  print(paste("Species:", args$species, "not yet coded"))
}
library(ChIPseeker)

# Contants ----------------------------------------------------------------

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
# outprefix <- file.path(outdir, paste0("UmapAnnots", bname))

print(paste("Writing output files using this prefix:", outprefix))

# Load data ---------------------------------------------------------------

outf <- paste0(outprefix, ".RData")
outpdf <- paste0(outprefix, ".pdf")

load(inf.lda, v=T)
tm.result <- posterior(out.lda)

# dont separate with _
colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "")
rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "")

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

dat.umap.long <- left_join(dat.umap.long, dat.var)

dat.umap.long$mark <- jmark

m.louv <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + facet_wrap(~plate) + ggtitle(jmark)

# m.celltype <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = celltype)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_manual(values = cbPalette) + facet_wrap(~plate)

m.var <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~plate) + ggtitle(jmark)



# Define topics  ----------------------------------------------------------


# plot and make objects

pdf(outpdf, width = 1240/72, height = 715/72, useDingbats = FALSE)

# m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = celltype)) + geom_point(size = 3) + theme_minimal() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_manual(values = cbPalette, na.value = "grey90") + ggtitle(jmark)
# print(m)
# 
# m.plate <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = celltype)) + geom_point(size = 3) + theme_minimal() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_manual(values = cbPalette, na.value = "grey90") + ggtitle(jmark) + facet_wrap(~plate)
# print(m.plate)

print(m.louv)
print(m.var)


par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
dat.merge <- dat.umap.long
dat.merge$cluster <- paste("louvain", dat.merge$louvain, sep = "")
dat.merge$topic <- dat.merge$cluster
dat.merge$topic.weight <- NA  # only relevant if we use LDA weights
dat.merge$plate <- sapply(dat.merge$cell, function(x) ClipLast(x, jsep = "_"))

# handle duplicates

# show final 
m.final <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = topic)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  facet_wrap(~plate) + ggtitle(jmark)
print(m.final)

# write topic loadings

if (args$WriteTopicLoadings){
    # jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed"
    jinf.tss <- args$inftss
    topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)
    if (args$species == "mouse"){
      annots.out <- AnnotateBins2(tm.result$terms, top.thres = 0.995, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene)
    } else {
      print(paste("Species:", args$species, "not yet coded. Try mouse"))
      stop(paste(args$species, "not yet implemented. Try mouse"))
    }
    annots.out$terms.annot$gene <- sapply(annots.out$terms.annot$termgene, function(x) ifelse(is.na(x), NA, strsplit(x, ";")[[1]][[2]]))

    dat.topics <- data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE)
    dat.umap.long.merge <- left_join(dat.merge, dat.topics)
      
    topn <- args$topn

    for (jtopic in topics.sum$topic){
      print(jtopic)
      m.umap <- PlotXYWithColor(dat.umap.long.merge, xvar = "umap1", yvar = "umap2", cname = jtopic, use.ggrastr = FALSE, jsize=4)  # can't install ggrastr on hpc easily
      terms.sub <- subset(annots.out$terms.annot, topic == jtopic)
      top.genes <- terms.sub$gene[1:topn]
      # integrating with public data not yet implemented
      # dat.sum.sub <- subset(dat.sum.long, gene %in% top.genes)
      # m.exprs <- ggplot(dat.sum.sub,
      #                   aes(x = forcats::fct_reorder(celltype, zscore, .desc=TRUE), y = zscore)) +
      #   geom_boxplot(outlier.shape = NA) +
      #   geom_jitter(width = 0.1, size = 0.5) +
      #   theme_classic() +
      #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      #   ggtitle(paste(jtopic, "Top:", topn, "N Unique Genes", length(top.genes))) + 
      #   xlab("")
      # print(m.exprs)
      
      # plot top 150 genes?
      jsub.terms <- subset(terms.sub, topic == jtopic & rnk <= topn) %>%
        ungroup() %>%
        mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
      m.top <- jsub.terms %>%
        ggplot(aes(x = term, y = log10(weight), label = gene)) +
        geom_point(size = 0.25) +
        theme_bw(8) +
        geom_text_repel(size = 2, segment.size = 0.1, segment.alpha = 0.25) +
        theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
        xlab("") + ylab("Log10 Bin Weight") +
        ggtitle(paste("Top peak weights for:", jtopic))
      print(m.umap)
      print(m.top)
    }
}


dev.off()

# save annotated dat to output
save(dat.merge, count.mat, out.lda, file = outf)

