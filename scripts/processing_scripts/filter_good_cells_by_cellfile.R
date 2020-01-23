# Jake Yeung
# filter_good_cells.R
# 2019-06-24
# DESCRIPTION
# 
#     Take .RData and filter out good cells, may need to modify sample names?
# 
# FOR HELP
# 
#     Rscript filter_good_cells.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-06-24
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='RData')
parser$add_argument('cellfile', metavar='CELLFILE',
                                            help='File list of cells')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='RData filtered out')
parser$add_argument('--invert', action="store_true", default=FALSE,
                                            help='Use cellfile as bad cells')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")
parser$add_argument("-p", "--preprocess", action="store_true", default=FALSE,
                        help="Print extra output [default]")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

good.cells <- unlist(read.csv2(args$cellfile, header=FALSE, stringsAsFactors=FALSE), use.names = FALSE)

if (endsWith(args$infile, suffix=".rds")){
    count.dat <- list()
    count.dat$counts <- readRDS(args$infile)
} else {
    load(args$infile, v=T)  # count.dat$counts
}

# preprocess sample names
print(head(colnames(count.dat$counts)))

if (args$preprocess){
  source("scripts/Rfunctions/QCFunctionsGastru.R")  # PreprocessSamp
  print("Preprocessing sample (warning, may not always be wise)")
  colnames(count.dat$counts) <- sapply(colnames(count.dat$counts), PreprocessSamp)
  print(head(colnames(count.dat$counts)))
}

print(paste("Number of cells before filtering...", ncol(count.dat$counts)))
if (!args$invert){
    # default
  samps.i <- which(colnames(count.dat$counts) %in% good.cells)
} else {
    # invert
  samps.i <- which(!colnames(count.dat$counts) %in% good.cells)
}

count.dat$counts <- count.dat$counts[, samps.i]
print(paste("Number of cells after filtering...", ncol(count.dat$counts)))

save(count.dat, file = args$outfile)
