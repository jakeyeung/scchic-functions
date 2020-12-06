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
library(glmpca)

library(hash)
library(igraph)
library(umap)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infile', metavar='INFILE',
                                            help='Iput LDA object (.Robj or .RData)')
parser$add_argument('-outfile', metavar='OUTFILE',
                                            help='Output matrix txt of imputed values')
parser$add_argument('-outpdf', metavar='OUTFILE',
                                            help='Output pdfs of the GLM output')
parser$add_argument('--AddChr', action="store_true", default=FALSE,
                                            help='Add chr to front of rownames')
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


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# Load data  --------------------------------------------------------------


inf <- args$infile
load(inf, v=T)  # contains count.mat and out.lda

# Loading objects:
#   glm.out
#   glm.inits
#   dat.merge2

# make plots

pdf(args$outpdf, useDingbats=FALSE)

dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings)

# dat.merge2$louvain <- NULL
# dat.merge2$umap1 <- NULL
# dat.merge2$umap2 <- NULL
# 
# # dat.merge <- left_join(dat.umap, subset(dat.merge2, select = -c(umap1, umap2, louvain))
# dat.merge <- left_join(dat.umap, dat.merge2)
# 
m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)

# if ("cell.var.within.sum.norm" %in% colnames(dat.merge)){
#     m.var <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
#       geom_point() + 
#       scale_color_viridis_c(direction = -1) + 
#       theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#     print(m.var)
# }
# 

dev.off()

dat.impute <- as.matrix(glm.out$loadings) %*% as.matrix(t(glm.out$factors))

# remove rows with no variance
genevars <- apply(dat.impute, 1, var)
geneskeep <- which(genevars > 0)
dat.impute <- dat.impute[geneskeep, ]

exprs.mat <- data.frame(Gene.ID = rownames(dat.impute), dat.impute, stringsAsFactors = FALSE)

exprs.mat$Gene.ID <- sapply(exprs.mat$Gene.ID, function(x) strsplit(x, ";")[[1]][[1]])
# remove gene names 

if (args$AddChr){
    geneids <- exprs.mat$Gene.ID
    print("Adding chr:")
    print("geneids before")
    print(head(geneids))
    geneids <- paste("chr", geneids, sep = "")
    print("geneids after")
    print(head(geneids))
    exprs.mat$Gene.ID <- geneids
}


fwrite(exprs.mat, file = args$outfile, sep = "\t", col.names = TRUE, row.names = FALSE)
# gzip file?

print(Sys.time() - jstart)

