# Jake Yeung
# filter_good_cells_good_bins.R
# 2019-12-24
# DESCRIPTION
# 
#     Filter good cells and good bins, optionally allow blacklist filter
# 
# FOR HELP
# 
#     Rscript filter_good_cells_good_bins.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-12-24
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infilerz', metavar='INFILESRZ', nargs = "+",
                                            help='infile containing TA frac')
parser$add_argument('-infilecounts', metavar='INFILESCOUNTMAT', nargs = "+",
                                            help='infile containing mat')
parser$add_argument('-names', metavar="space delim strings", nargs = "+",
                                            help='Label for each infile')
parser$add_argument('-outdir', metavar='OUTDIR',
                                            help='outdir')
parser$add_argument('-countcutoff', metavar='INTEGER', type = 'integer', 
                                            help='Minimum counts')
parser$add_argument('-TAcutoff', metavar='TAcutoff', type = 'character',
                                            help='Minimum TA fraction')
parser$add_argument('-blfile', metavar='BEDFILE',
                                            help='List of blacklist regons to exclude')
parser$add_argument('--overwrite', action="store_true", default=FALSE, help="Force overwrite")
parser$add_argument("--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")

                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

print(args)

# check files
assertthat::assert_that(length(args$infilerz) == length(args$infilecounts))
assertthat::assert_that(length(args$names) == length(args$infilecounts))

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

jstart <- Sys.time()


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(JFuncs)
library(scchicFuncs)

library(GenomicRanges)



# Constants set to arguments  ---------------------------------------------


# blfile <- "/home/jyeung/hpc/databases/blacklists/mm10.blacklist.copy.bed.gz"
blfile <- args$blfile
assertthat::assert_that(file.exists(blfile))

# outdir <- paste0("/home/jyeung/data/from_rstudioserver/dblchic/quality_control_dbl.mixed_with_unenriched_again.", Sys.Date())
outdir <- args$outdir
dir.create(outdir)

cutoff.counts <- args$countcutoff
cutoff.TA <- args$TAcutoff
overwrite <- args$overwrite

infs.rz <- as.list(args$infilerz)
infs.mat <- as.list(args$infilecounts)

# indir <- "/home/jyeung/hpc/dblchic/from_rstudio/bonemarrow/countTables_RZ_wthunenriched_again"
# infs.rz <- c(file.path(indir, "BM-B6-H3K4me1_with_dblexperi.2019-12-23.LHcounts.csv.gz"),
#              file.path(indir, "BM-B6-H3K27me3_with_dblexperi.2019-12-23.LHcounts.csv.gz"),
#              file.path(indir, "BM-B6-K4m1-K27m3-191008-merged.LHcounts.csv.gz"))
# infs.mat <- c(file.path(indir, "BM-B6-H3K4me1_with_dblexperi.2019-12-23.mq_40.bsize_100000.step_20000.csv.gz"),
#              file.path(indir, "BM-B6-H3K27me3_with_dblexperi.2019-12-23.mq_40.bsize_100000.step_20000.csv.gz"),
#              file.path(indir, "BM-B6-K4m1-K27m3-191008-merged.mq_40.bsize_100000.step_20000.csv.gz"))

lapply(infs.rz, function(x) assertthat::assert_that(file.exists(x)))
lapply(infs.mat, function(x) assertthat::assert_that(file.exists(x)))

# jnames <- c("H3K4me1", "H3K27me3", "H3K4me1_H3K27me3")
jnames <- args$names
names(jnames) <- jnames

outdir <- args$outdir

outpaths <- lapply(jnames, function(jname){
  jout <- file.path(outdir, paste0("count_mat.",jname, ".countcutoff_", cutoff.counts, ".TAcutoff_", cutoff.TA, ".rds"))
  if (!overwrite){
    assertthat::assert_that(!file.exists(jout), msg = paste("Outfile exists, not overwriting for safety:", jout))
  }
  return(jout)
})
pdfout <- file.path(outdir, paste0("qc_plots.", paste(jnames, collapse = "-"), ".pdf"))
if (!overwrite){
  assertthat::assert_that(!file.exists(pdfout), msg = paste("Pdfout exists, not overwriting for safety:", pdfout))
}

names(infs.mat) <- jnames
names(infs.rz) <- jnames

print(infs.rz)

dat.rz <- ReadLH.SummarizeTA(infs.rz, remove.nones = FALSE, na.to.zero = TRUE, bind.rows = FALSE)

# add jname to dat.rz
dat.rz.filt <- lapply(jnames, function(jmark){
  jtmp <- dat.rz[[jmark]] %>%
    rowwise() %>%
    mutate(mark = jmark)
}) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(experi = ClipLast(samp, jsep = "_"),
         cellindx = paste("cell", strsplit(samp, "_")[[1]][[2]], sep = ""))



# Get good cells  ---------------------------------------------------------

empty.wells <- GetEmptyWells()

dat.rz.filt <- dat.rz.filt %>%
  rowwise() %>%
  mutate(empty.well = cellindx %in% empty.wells, 
         good.cell = total.count > cutoff.counts & TA.frac > cutoff.TA & !empty.well)


m.density.bymark <- ggplot(dat.rz.filt, aes(x = total.count, fill = mark)) + geom_density(alpha = 0.3)  + 
  scale_x_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_viridis_d()

m.scatter.bymark.col <- ggplot(dat.rz.filt, aes(x = total.count, y = TA.frac, color = good.cell, size = empty.well, shape = empty.well)) + geom_point()  + 
  scale_x_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, ncol = 1)

# for each mark, plot each plate 


# how many good cells?
dat.rz.filt %>%
  group_by(experi) %>%
  filter(good.cell) %>%
  summarise(ncells = length(total.count))

dat.rz.filt %>%
  group_by(experi) %>%
  summarise(ncells = length(total.count))

cells.keep <- subset(dat.rz.filt, good.cell)$samp
print(paste("Number of cells before:", nrow(dat.rz.filt)))
print(paste("Number of cells after:", length(cells.keep)))


# Plot pdf outputs --------------------------------------------------------

pdf(pdfout, useDingbats = FALSE)
  print(m.scatter.bymark.col)
  print(m.density.bymark)
dev.off()


# Load mats ---------------------------------------------------------------

# filter good cells
mats <- lapply(infs.mat, function(inf){
  mat <- ReadMatSlideWinFormat(inf)
  cols.i <- colnames(mat) %in% cells.keep
  mat.filt <- mat[, cols.i]
})

# remove blacklist
print("Before removing blacklist")
print(lapply(mats, dim))
mats <- lapply(mats, function(mat){
  rnames.filt <- FilterBinsByBlacklist(rownames(mat), blfile = blfile)
  rows.i <- rownames(mat) %in% rnames.filt
  mat.filt <- mat[rows.i, ]
})
print("After removing blacklist")
print(lapply(mats, dim))

# save each mat separately 
for (jname in jnames){
  outpath <- outpaths[[jname]]
  jtmp <- mats[[jname]]
  assertthat::assert_that(nrow(jtmp) > 0 & ncol(jtmp) > 0)
  saveRDS(jtmp, file = outpath)
}

