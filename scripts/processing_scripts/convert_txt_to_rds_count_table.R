#!/bin/sh
# Jake Yeung
# convert_txt_to_rds_count_table.sh
#  
# 2020-02-15
# 
#     Description
# 
# FOR HELP
# 
#     Rscript convert_txt_to_rds_count_table.sh --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-02-15
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(scchicFuncs)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='TXT', help='Output of bam to count table using bed or bin strategy')
parser$add_argument('outfile', metavar='OUTFILE', help='Output count table in .rds output')
parser$add_argument('-format', metavar='bin or bed', help='Should be bin or bed, default bed', default = "bed")
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

if (args$format == "bed"){
  mat <- ReadMatTSSFormat(args$infile, as.sparse=TRUE, add.coord = TRUE)
} else if (args$format == "bin"){
  mat <- ReadMatSlideWinFormat(args$infile, as.sparse=TRUE)
} else {
  stop("Format must be bed or bin, found: ", args$format)
}

saveRDS(mat, file = args$outfile)
