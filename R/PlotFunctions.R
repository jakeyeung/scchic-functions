# Jake Yeung
# Date of Creation: 2019-06-16
# File: ~/projects/scchic_gastru/functions/PlotFunctions.R
# 

GetGeneAnnotsHash <- function(inf.annot){
  dat.annot <- data.table::fread(inf.annot, col.names = c("chromo", "start", "end", "bname"))
  # add chr
  dat.annot$chromo <- paste("chr", dat.annot$chromo, sep = "")
  rnames.old <- paste(dat.annot$chromo, paste(dat.annot$start, dat.annot$end, sep = "-"), sep = ":")
  rnames.new <- dat.annot$bname
  annots.hash <- hash::hash(rnames.old, rnames.new)
}

AddGeneNameToRows <- function(mat, annots.hash){
  # mat rownmaes got stripped of gene names, add them back
  rnames.old <- rownames(mat)
  rnames.new <- sapply(rnames.old, function(x) annots.hash[[x]])
  rownames(mat) <- rnames.new
  return(mat)
}

PlotDecreasingWeights <- function(jsub.terms, jtitle = "", order.term.by.weight = TRUE, textsize = 3, themesize = 5){
  if (order.term.by.weight){
    jsub.terms <-  jsub.terms %>%
      mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
  }
  m.top <- jsub.terms %>%
    ggplot(aes(x = term, y = log10(weight), label = gene)) +
    geom_point(size = 0.25) +
    theme_bw(themesize) +
    geom_text_repel(size = textsize, segment.size = 0.1, segment.alpha = 0.25) +
    theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
    xlab("") + ylab("Log10 Bin Weight") +
    ggtitle(jtitle)
  return(m.top)
}

PlotPseudobulkZscore <- function(dat.bulk.sub, order.celltype.by.zscore = TRUE, xlabsize = 8, themesize = 12){
  if (order.celltype.by.zscore){
    dat.bulk.sub <- dat.bulk.sub %>%
      mutate(celltype = forcats::fct_reorder(celltype, dplyr::desc(zscore), .fun = median))
  }
  m.exprs <- ggplot(dat.bulk.sub,
                    aes(x = celltype , y = zscore)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.5) +
    theme_classic(themesize) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = xlabsize)) +
    ggtitle(paste("topic", jtop, "Top:", keeptop, "N Unique Genes", length(top.genes)))
  return(m.exprs)

}


PlotSpatialGene <- function(dat.adj, jfits.nlm.all, jgene){
  jgene.short <- strsplit(jgene, "_")[[1]][[2]]
  jsub <- dat.adj %>% filter(gene == jgene)
  jfit.sub <- subset(jfits.nlm.all %>% filter(gene == jgene))
  best.model <- jfit.sub$model
  jfit <- GetFit(jfit.sub, colnames(jfit.sub), model = best.model)[[1]]
  desmat <- model.matrix(~1 + plate, jsub)
  assertthat::assert_that(ncol(desmat) == 2)
  if (best.model == "gauss"){
    jsub$exprs.pred <- GetGauss(jfit$par, jsub$x, desmat[, 2])
  } else if (best.model == "flat") {
    jsub$exprs.pred <- GetFlat(jfit$par, jsub$x, desmat[, 2])
  } else {
    warning("each model must be coded separately")
  }
  jsub <- jsub %>%
    rowwise() %>%
    mutate(exprs.pred.adj = ifelse(plate == "p02", exprs.pred - batch.effect, exprs.pred))
  m1 <- ggplot(jsub, aes(x = x, y = exprs, color = plate)) + geom_point() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    ggtitle(paste(jgene.short, "")) + 
    ylab("Log2 Expression") + xlab("Position [A -> P]") + 
    geom_line(mapping = aes(x = x, y = exprs.pred), data = jsub, inherit.aes = FALSE, linetype = "dotted") 
  m2 <- ggplot(jsub, aes(x = x, y = exprs.adj, color = plate)) + geom_point() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(paste(jgene.short, "Batch-Corrected")) + 
    ylab("Adjusted Log2 Expression") + xlab("Position [A -> P]") + 
    geom_line(mapping = aes(x = x, y = exprs.pred.adj), data = jsub, inherit.aes = FALSE, linetype = "dotted") 
  return(list(m1, m2))
  # multiplot(m1, m2, cols = 2)
}

PlotXYNoColor <- function(jsub, xvar, yvar, jcol = "gray80", jsize = 1){
  m <- ggplot(jsub, aes_string(x = xvar, y = yvar)) +
    ggrastr::geom_point_rast(size = jsize, color = jcol) +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank()) +
    xlab("") + ylab("")
  return(m)
}
PlotXYWithColor <- function(jsub, xvar = "X1", yvar = "X2", cname = "activity", jcol = scales::muted("darkblue"), jtitle = "", jcol.low = "gray85", jcol.mid = "gray50", jsize = 1, leg.name = NULL, jjrange = "auto",
                            cont.color = TRUE, col.palette = NA, strip.ticks = FALSE, manual.mid = NA, remove.axis.info = TRUE){
  if (is.null(leg.name)){
    leg.name <- cname
  }
  cname.str <- cname
  if (strip.ticks){
    cname <- StripTicks(cname)
  }
  jsub <- RankOrder(jsub, cname = cname, out.cname = "orderrank")
  m1 <- ggplot(jsub, aes_string(x = xvar, y = yvar, col = cname.str, order = "orderrank")) +
    ggrastr::geom_point_rast(size = jsize) +
    theme_bw() +
    xlab("") + ylab("") + ggtitle(jtitle)
  if (remove.axis.info){
    m1 <- m1 + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                         axis.ticks=element_blank(),
                         axis.text.x=element_blank(),
                         axis.text.y=element_blank(),
                         panel.border=element_blank())
  }  else {
    m1 <- m1 + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  }
  if (cont.color){
    jrange <- range(jsub[[cname]])
    if (is.na(manual.mid)){
      jmid <- min(jsub[[cname]]) + diff(range(jsub[[cname]])) / 2
    } else {
      jmid <- manual.mid
    }
    # print(jmid)
    if (jjrange != "auto"){
      jrange <- jjrange
    }
    m1 <- m1 +
      scale_color_gradient2(low = jcol.low, mid = jcol.mid, high = jcol, midpoint = jmid, limit = jrange, name = leg.name)
  }  else {
    m1 <- m1 + scale_color_manual(values = col.palette)
  }
  return(m1)
}

RankOrder <- function(dat.tmp, cname, out.cname = "orderrank"){
  dat.tmp[[out.cname]] <- rank(dat.tmp[[as.character(cname)]], ties.method = "first")
  dat.tmp <- dat.tmp %>%
    arrange_(.dots = out.cname)
  return(dat.tmp)
}
