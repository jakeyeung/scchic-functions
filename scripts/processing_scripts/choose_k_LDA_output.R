# Jake Yeung
# 3-convert_LDA_to_bins.R
# 2019-09-30
# DESCRIPTION
# 
#     Description
# 
# FOR HELP
# 
#     Rscript 3-convert_LDA_to_bins.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-09-30
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(scchicFuncs)
library(topicmodels)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input LDA Output. .RData object named out.lda containing a list of LDA outputs')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output LDA bins. .rds output. Chosen K is added to file name')
parser$add_argument("-k", "--kchoose", type="integer", help="K to filter LDA of interest. Set nothing (NULL) for choosing best k based on likelihood", default = NULL)
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

load(args$infile, v=T)  # out.lda

if (is.null(args$kchoose)){
    out.lda <- ChooseBestLDA(out.lda)
} else {
    kvec <- sapply(out.lda, function(x) x@k)
    i <- which(kvec == args$kchoose)
    assertthat::assert_that(length(i) == 1)
    out.lda <- out.lda[[i]]
}
print(paste("K chosen:", out.lda@k))

outfile.noext <- tools::file_path_sans_ext(args$outfile)
# add K and then .rds
jsuffix <- paste0(".K_", out.lda@k, ".rds")
outfile <- paste0(outfile.noext, jsuffix)

saveRDS(out.lda, file = outfile)
