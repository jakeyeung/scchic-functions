FitMixtureModelLabelCells <-function(topics.mat, topics.keep, jthres = 0.5, show.plots = TRUE, dat.umap.long = NULL, quantile.if.fail = 0.99, fail.if.nfrac = 0.75){
  # topics.mat: rows are samples, columns are latent variables / topics
  # dat.umap.long requires columns umap1, umap2 and cell name for plotting mixturre fitting output on umap
  assertthat::assert_that(all(topics.keep %in% colnames(topics.mat)))
  
  if (is.null(names(topics.keep))){
    names(topics.keep) <- topics.keep
  }
  # label cells by topic 
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  mm.celltype.lst <- lapply(topics.keep, function(jtopic){
    print(jtopic)
    tvec.raw <- sort(topics.mat[, jtopic])
    # transform
    tvec <- log(tvec.raw / (1 - tvec.raw))
    # xline <- quantile(tvec, )
    mm <- mixtools::normalmixEM(x = tvec, lambda = c(0.9, 0.1), mu = c(-5, -1), sigma = c(2, 1), k = 2)
    (xline <- min(mm$x[which(mm$posterior[, 1] < jthres)]))
    # xline <- -2
    xcells <- names(tvec)[which(tvec > xline)]
    print(paste(length(xcells), "/", length(tvec), "assigned to", jtopic))
    if (length(xcells) / length(tvec) >= fail.if.nfrac){
      print(paste0("MM failed: too many cells assigned in mixture model, resorting to taking top: ", quantile.if.fail))
      xcells <- names(tvec)[which(tvec > quantile(tvec, quantile.if.fail))]
      print(paste("Manual threshold:", length(xcells), "/", length(tvec), "assigned to", jtopic))
    }
    if (length(xcells) == 0){
      print(paste0("No cells assigned, resorting to taking top: ", quantile.if.fail))
      xcells <- names(tvec)[which(tvec > quantile(tvec, quantile.if.fail))]
      print(paste("Manual threshold:", length(xcells), "/", length(tvec), "assigned to", jtopic))
    }
    # get topic value for assigned cells
    tvec.raw.filt <- tvec.raw[xcells]
    
    # xline <- max(mm$x[indx.btwn][which(post.filt[, 2] > 0.5)])
    if (show.plots){
      plot(density(tvec), main = paste(jtopic, jthres), xlab = "Log Odds [log(p / (1 - p))]")
      abline(v = xline, col = 'blue')
      plot.mixEM(mm, whichplots = 2, xlab2 = "Log Odds [log(p / (1 - p))]", main2 = paste(jtopic, jthres))
      abline(v = xline, col = 'blue')
    }
    
    cells.keep <- xcells
    if (show.plots){
      m.check <- PlotXYWithColor(dat.umap.long %>% mutate(is.celltype = cell %in% cells.keep), xvar = "umap1", yvar = "umap2", cname = "is.celltype", jtitle = paste(jtopic, jthres), cont.color = FALSE, col.palette = cbPalette)
      print(m.check)
    }
    
    return(list(topic = jtopic, topic.weight = tvec.raw.filt, celltype = xcells, mm = mm, threshold = xline))
  })
}

TidyMixtureModelOutputs <- function(mm.celltype.lst, dedup = TRUE){
  celltypes <- lapply(names(mm.celltype.lst), function(ctype){
    data.frame(cluster = ctype, cell = mm.celltype.lst[[ctype]]$celltype, 
               topic = mm.celltype.lst[[ctype]]$topic, 
               topic.weight = mm.celltype.lst[[ctype]]$topic.weight, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  
  # handle replicates
  if (dedup){
    celltypes <- celltypes %>%
      group_by(cell) %>%
      filter(topic.weight == max(topic.weight))
  }
  return(celltypes)
}