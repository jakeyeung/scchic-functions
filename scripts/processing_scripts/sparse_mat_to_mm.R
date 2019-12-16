# Jake Yeung
# sparse_mat_to_mm.R
# 2019-07-11
# DESCRIPTION
# 
#     Read sparse mat as RData, write to mm (expect count.dat$counts)
# 
# FOR HELP
# 
#     Rscript sparse_mat_to_mm.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-07-11
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(Matrix)

BinarizeMatrix <- function(x){
  # https://stackoverflow.com/questions/14526429/turn-a-count-matrix-into-a-binary-existen
  xbin <- as.numeric(as.matrix(x) > 0)
  xbin <- Matrix::Matrix(xbin, sparse = TRUE, nrow = nrow(x), ncol = ncol(x))
  # get back the column and row names
  rownames(xbin) <- rownames(x)
  colnames(xbin) <- colnames(x)
  return(xbin)
}

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input RData object with count.dat$counts')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='File output to mm')
parser$add_argument("-b", "--binarize", action="store_true", default=FALSE,
                        help="Binarize the matrix")
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

if (endsWith(args$infile, ".rds")){
  count.dat <- list()
  print(paste("reading rds", args$infile))
  count.dat$counts <- readRDS(args$infile)
} else {
  print(paste("Assuming", args$infile, "is an .RData file"))
  load(args$infile, v=T)  # count.dat$counts
}

if ( args$binarize ) {
  print("Binarizing matrix...")
  count.dat$counts <- BinarizeMatrix(count.dat$counts)   
}

writeMM(count.dat$counts, sparse = TRUE, file = args$outfile)

# write rownames and column names
write.table(rownames(count.dat$counts), file = paste0(args$outfile, ".rownames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(colnames(count.dat$counts), file = paste0(args$outfile, ".colnames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
