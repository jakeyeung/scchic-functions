# Jake Yeung
# csv_to_rds_TSS.R
# 2019-06-24
# DESCRIPTION
# 
#     CSV from bamToCountTable (TSS) to RDS
# 
# FOR HELP
# 
#     Rscript csv_to_rds_TSS.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-06-24
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(data.table)
library(Matrix)
library(tidyr)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='CSV file from bamToCountTable using bedfile option')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='RDS output file')
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

# two rows of headers!
dat <- fread(args$infile, header=TRUE)[-1, ] %>%
  dplyr::rename(chromo = sampleName, start=V2, end=V3, coord=V4)
# # merge the two headers
# cnames.i <- which(!is.na(dat[1, ]))
# colnames(dat)[cnames.i] <- dat[1, cnames.i]
# dat <- dat[-1, ]
coord <- dat$coord
dat$chromo <- NULL
dat$start <- NULL
dat$end <- NULL
dat$coord <- NULL
dat <- as.matrix(dat)
dat[is.na(dat)] <- 0
rownames(dat) <- coord
dat <- Matrix::Matrix(dat, sparse = TRUE)
count.dat <- list()
count.dat$counts <- dat
save(count.dat, file = args$outfile)
# saveRDS(dat, file = args$outfile)


