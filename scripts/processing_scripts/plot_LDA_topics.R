# Jake Yeung
# plot_LDA_topics.R
# 2020-04-17
# DESCRIPTION
# 
#     Plot downstream topics from ldaoutput
# 
# FOR HELP
# 
#     Rscript plot_LDA_topics.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-04-17
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(topicmodels)
library(dplyr)
library(ggplot2)
library(Matrix)
library(here)
library(scchicFuncs)

library(tidyr)

library(hash)
library(igraph)
library(umap)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(data.table)

library(ggrepel)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input .RData or .RObj from run_lda_model script')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output pdf')
parser$add_argument("-keeptop", metavar="Integer", type="integer", default = 150,
                        help="Keep top n ")
parser$add_argument("-inftss", metavar="Path", default="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed",
                        help="Path to TSS for annotating peak to gene")
parser$add_argument("--TopicsOnly", action="store_true", default=FALSE,
                        help="Plot topic weights only (does not care about finding/matching to genes)")
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

keeptop <- args$keeptop
inf.tss <- args$inftss

plotpath <- args$outfile
load(args$infile, v=T)  # out.lda, count.mat, count.mat.orig

# write plots to output
# save plots
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

tm.result <- posterior(out.lda)

tm.result <- AddTopicToTmResult(tm.result, jsep="_")

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- sort(unique(sapply(colnames(tm.result$terms), function(x) strsplit(x, ":")[[1]][[1]])))
pdf(plotpath, width = 1240/72, height = 815/72, useDingbats = FALSE)
    # do UMAP, plot imputed intrachromosomal variance 
    dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings) %>%
      rowwise() %>%
      mutate(plate = ClipLast(as.character(cell), jsep = "_"))
    dat.var <- CalculateVarAll(dat.impute.log, jchromos)
    dat.merge <- left_join(dat.umap, dat.var)
    m.louv <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      facet_wrap(~plate)
    m.intrachrom <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      facet_wrap(~plate) + 
      scale_color_viridis_c(direction = -1)
    print(m.louv)
    print(m.intrachrom)

    # plot topics
    topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)


    dat.merged.topics <- left_join(dat.umap, data.frame(cell = rownames(tm.result$topics), tm.result$topics))

    terms.mat <- tm.result$terms

    if (!args$TopicsOnly){
      annot.out <- AnnotateBins(terms.mat = terms.mat, inf.tss = inf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")
      terms.filt <- annot.out$terms.filt
    }

    for (jtop in topics.sum$topic){

      print(jtop)
      m.umap <- PlotXYWithColor(dat.merged.topics, xvar = "umap1", yvar = "umap2", cname = jtop, use.ggrastr=FALSE)



      if (!args$TopicsOnly){
        top.genes <- subset(terms.filt, topic == jtop & rnk <= keeptop)$gene
        assertthat::assert_that(length(top.genes) > 0)
        jsub.terms <- subset(terms.filt, topic == jtop & rnk < keeptop) %>%
          ungroup() %>%
          mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
        m.top <- jsub.terms %>%
          # mutate(term = forcats::fct_reorder(term, dplyr::desc(weight))) %>%
          ggplot(aes(x = term, y = log10(weight), label = gene)) +
          geom_point(size = 0.25) +
          theme_bw(8) +
          geom_text_repel(size = 2, segment.size = 0.1, segment.alpha = 0.25) +
          theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
          xlab("") + ylab("Log10 Bin Weight") +
          ggtitle(paste("Top peak weights for:", jtop))
        # plot everything
      }

      print(m.umap)
      if (!args$TopicsOnly){
        print(m.top)
      }
    }
dev.off()
