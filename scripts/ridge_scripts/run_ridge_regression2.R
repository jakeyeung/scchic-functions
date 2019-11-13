# source('/home/yeung/projects/ridge-regression/ridgeInR.R')
# source("/home/hub_oudenaarden/jyeung/projects/from_PhD/ridge-regression/ridgeInR.R")
source("/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/ridge_scripts/ridgeInR.R")
# source("ridgeInR.R")

# exp: matrix of expression, row centered. First column is Gene.ID with gene names. Column names are sample names

args <- commandArgs(trailingOnly = TRUE)
exp <- args[1]
site <- args[2]
lambda <- args.to.numeric(args[3], default=FALSE)
if (is.numeric(lambda)){
  print(lambda)
  print(paste("Using global lambda:"))
} else {
  print(lambda)
  print(paste("Using lambda from CV"))
}

# E = read.table(exp, header=T, row.names = 1)
# N = read.table(site, header=T, row.names = 1)

E <- read.table.handlerows(exp)
N <- read.table.handlerows(site)

print(head(E))
print(head(N))
# if (ncol(E) > 5) print(E[1:5, 1:5])
# if (ncol(N) > 5) print(N[1:5, 1:5])

# E = log2(as.matrix(E) + 1.0)  # convert to log2 BEFOREhand
E = as.matrix(E)
N = as.matrix(N)

# Es = center.cols(E)
# Ns = center.cols(N)

# take common rows
common.genes <- intersect(rownames(E), rownames(N))
max.genes <- max(nrow(E), nrow(N))
sprintf('Using %s/%s genes because they are in common.', length(common.genes), max.genes)
E <- E[common.genes, ]
N <- N[common.genes, ]

if (lambda == FALSE){
  opt = optimize.lambda(N, E)
} else {
  print("Setting global optimum")
  opt = list(lambda.opt = lambda)
}
r = ridge.regression(N, E, opt$lambda.opt)
top20 = sort(r$combined.Zscore, decreasing=TRUE)[1:30]
# print(dim(r$Ahat))
# print(rownames(r$Ahat))

# if (ncol(E) == 24){
#   x = seq(18,64,by=2)      
# } else if (ncol(E) == 8){
#   x = seq(22,64,by=6)
# } else if (ncol(E) == 32){
#   x = c(seq(18,64,by=2), seq(22,64,by=6))
# } else {
#   print(paste("N columns in E:", ncol(E)))
#   warning("Expected 24, 6, or 32 columns in E")
#   x = rep("", ncol(E))
# }
# i = 1
# for(n in names(top20)) {
#       pdf(file=paste(i, '_', n, '.pdf', sep=''), height=5, width=7)
#       activities = r$Ahat[n, ]
#       par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#       plot(x, activities, type='l', xlab='CT', ylab='Activities', lwd=4, lty=2, col='blue', xaxt='n', main=n)
#       legend('topright', inset=c(-0.27,0), legend=toString(round(r$combined.Zscore[n], 4)), title='Z-score')
#       points(x, activities, cex=1.5, col='red', pch=19)
#       axis(1, at=x, las=3)
#       dev.off()
#       i = i + 1
# }

print("Printing colnames:")
print(colnames(E))
for (n in colnames(E)){
      write(n, "Colnames", append=TRUE)
}

for(n in names(r$combined.Zscore)){
      write(paste(n, r$combined.Zscore[n], sep='\t'), "Zscores", append=TRUE)
}

for(n in names(r$combined.Zscore)){
      x = paste( r$Ahat[n, ], collapse="\t")
      write(paste(n, x, sep='\t'), "Activities", append=TRUE)
}

# print(r$AhatSE)
# print(str(r$AhatSE))
for(n in names(r$combined.Zscore)){
      x = paste( r$AhatSE[n, ], collapse="\t")
      write(paste(n, x, sep='\t'), "StandardError", append=TRUE)
}

write(toString(mean(r$fov)), 'FOV')
write(toString(opt$lambda.opt), "Lambda")

save(r, E, N, exp, site, file = "r.Robj")

