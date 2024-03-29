# Jake Yeung
# Date of Creation: 2019-08-25
# File: ~/projects/scchicFuncs/R/AuxLDA.R
#

AddTopicToTmResult <- function(tm.result, jsep = ""){
  colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = jsep)
  rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = jsep)
  return(tm.result)
}

AnnotateCoordsFromList.GeneWise <- function(coords.vec,
                                            inf.tss="/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed",
                                            txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                                            annodb = "org.Mm.eg.db",
                                            chromos.keep=c(paste("chr", seq(19), sep = ""), "chrX", "chrY")){
  library(ChIPseeker)
  regions <- data.frame(seqnames = sapply(coords.vec, GetChromo),
                        start = sapply(coords.vec, GetStart),
                        end = sapply(coords.vec, GetEnd),
                        stringsAsFactors = FALSE)
  rownames(regions) <- coords.vec
  # regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))
  print("Chromos to keep")
  print(chromos.keep)
  regions <- subset(regions, seqnames %in% chromos.keep)
  
  regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
  regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                  # TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                  TxDb=txdb,
                                                  # annoDb='org.Mm.eg.db'))
                                                  annoDb=annodb))
  regions.annotated$region_coord <- names(regions.range)
  
  if (is.character(inf.tss)){
    assertthat::assert_that(file.exists(inf.tss))
    tss.dat <- fread(inf.tss, col.names = c("seqnames", "start", "end", "tssname"))
  } else {
    print(paste("inf.tss not a filepath, interpreting as tss.dat has been directly added"))
    tss.dat <- inf.tss
  }
  tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, ";")[[1]][[2]])
  
  annots.biomart <- regions.annotated %>%
    mutate(midpt = start + (end - start) / 2)
  
  annots.gr <- makeGRangesFromDataFrame(annots.biomart %>% dplyr::select(seqnames, start, end, SYMBOL, region_coord), keep.extra.columns = TRUE)
  annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)
  
  out2 <- findOverlaps(annots.tss.gr, annots.gr, type = "any")
  # out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")
  
  out2.df = data.frame(annots.gr[subjectHits(out2),], annots.tss.gr[queryHits(out2),]) %>%
    mutate(midpt = start + round(width / 2),
           midpt.1 = start.1 + round(width.1 / 2),
           dist.to.tss = midpt.1 - midpt)
  
  # # filter closest
  # out2.df.closest <- out2.df %>%
  #   group_by(region_coord) %>%
  #   filter(abs(dist.to.tss) == min(abs(dist.to.tss)))
  
  # terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, sep = ";")
  # terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)
  
  return(list('regions.annotated' = regions.annotated, 'out2.df' = out2.df, 'out2' = out2, "annots.tss.gr" = annots.tss.gr, 'annots.gr' = annots.gr))
}


AnnotateBins2 <- function(terms.mat, top.thres=0.995, inf.tss="/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep=c(paste("chr", seq(19), sep = ""), "chrX", "chrY"), skip.split = FALSE){
  assertthat::assert_that(file.exists(inf.tss))
  assertthat::assert_that(class(terms.mat) == "matrix")
  regions <- data.frame(seqnames = sapply(colnames(terms.mat), GetChromo),
                        start = sapply(colnames(terms.mat), GetStart),
                        end = sapply(colnames(terms.mat), GetEnd),
                        stringsAsFactors = FALSE)
  rownames(regions) <- colnames(terms.mat)
  # regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))
  print("Chromos to keep")
  print(chromos.keep)
  regions <- subset(regions, seqnames %in% chromos.keep)
  
  regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
  regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                  # TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                  TxDb=txdb,
                                                  # annoDb='org.Mm.eg.db'))
                                                  annoDb=annodb))
  regions.annotated$region_coord <- names(regions.range)
  
  topic.regions <- lapply(seq(nrow(terms.mat)), function(clst){
    return(SelectTopRegions(terms.mat[clst, ], colnames(terms.mat), method = "thres", method.val = top.thres))
  })
  
  print(paste("Using TSS definitions from:", inf.tss))
  terms.long <- data.frame(term = colnames(terms.mat), as.data.frame(t(terms.mat)), stringsAsFactors = FALSE) %>%
    gather(key = "topic", value = "weight", -term) %>%
    mutate(topic = gsub("X", "", topic)) %>%
    group_by(topic) %>%
    arrange(desc(weight)) %>%
    mutate(rnk = seq(length(weight))) %>%
    rowwise()
  terms.filt.top <- terms.long %>%
    # filter(rnk < 1000) %>%  # DO GENOME WIDE
    rowwise()
  tss.dat <- fread(inf.tss, col.names = c("seqnames", "start", "end", "tssname"))
  tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, ";")[[1]][[2]])
  
  annots.biomart <- regions.annotated %>%
    mutate(midpt = start + (end - start) / 2) %>%
    filter(region_coord %in% terms.filt.top$term)
  
  annots.gr <- makeGRangesFromDataFrame(annots.biomart %>% dplyr::select(seqnames, start, end, SYMBOL, region_coord, ENSEMBL), keep.extra.columns = TRUE)
  annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)
  
  out <- findOverlaps(annots.tss.gr, annots.gr, type = "within")
  out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")
  
  out2.df = data.frame(annots.gr[queryHits(out2),], annots.tss.gr[subjectHits(out2),]) %>%
    mutate(midpt = start + round(width / 2),
           midpt.1 = start.1 + round(width.1 / 2),
           dist.to.tss = midpt.1 - midpt)
  
  # filter closest
  out2.df.closest <- out2.df %>%
    group_by(region_coord) %>%
    filter(abs(dist.to.tss) == min(abs(dist.to.tss)))
  
  terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, sep = ";")
  terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)
  
  terms.annot <- terms.filt.top %>%
    mutate(termgene = ifelse(!is.null(terms.hash[[term]]), terms.hash[[term]], NA))
  
  terms.filt <- terms.filt.top %>%
    mutate(termgene = ifelse(!is.null(terms.hash[[term]]), terms.hash[[term]], NA)) %>%
    filter(!is.na(termgene))
  if (!skip.split){
    terms.filt <- terms.filt %>%
      mutate(gene = sapply(termgene, function(x) strsplit(x, ";")[[1]][[2]])) %>%
      group_by(gene)  # dont filter use
  }
  return(list('topic.regions' = topic.regions, 'regions.annotated' = regions.annotated, 'terms.annot' = terms.annot, 'out2.df.closest' = out2.df.closest, 'terms.filt' = terms.filt))
}



AnnotateCoordsFromList <- function(coords.vec,
                                   inf.tss="/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed",
                                   txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                                   annodb = "org.Mm.eg.db",
                                   chromos.keep=c(paste("chr", seq(19), sep = ""), "chrX", "chrY")){
  # library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  # library(org.Mm.eg.db)
  library(ChIPseeker)
  assertthat::assert_that(file.exists(inf.tss))
  regions <- data.frame(seqnames = sapply(coords.vec, GetChromo),
                        start = sapply(coords.vec, GetStart),
                        end = sapply(coords.vec, GetEnd),
                        stringsAsFactors = FALSE)
  rownames(regions) <- coords.vec
  # regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))
  print("Chromos to keep")
  print(chromos.keep)
  regions <- subset(regions, seqnames %in% chromos.keep)

  regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
  regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                  # TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                  TxDb=txdb,
                                                  # annoDb='org.Mm.eg.db'))
                                                  annoDb=annodb))
  regions.annotated$region_coord <- names(regions.range)

  tss.dat <- fread(inf.tss, col.names = c("seqnames", "start", "end", "tssname"))
  tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, ";")[[1]][[2]])

  annots.biomart <- regions.annotated %>%
    mutate(midpt = start + (end - start) / 2)

  annots.gr <- makeGRangesFromDataFrame(annots.biomart %>% dplyr::select(seqnames, start, end, SYMBOL, region_coord), keep.extra.columns = TRUE)
  annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)

  out <- findOverlaps(annots.tss.gr, annots.gr, type = "within")
  out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")

  out2.df = data.frame(annots.gr[queryHits(out2),], annots.tss.gr[subjectHits(out2),]) %>%
    mutate(midpt = start + round(width / 2),
           midpt.1 = start.1 + round(width.1 / 2),
           dist.to.tss = midpt.1 - midpt)

  # filter closest
  out2.df.closest <- out2.df %>%
    group_by(region_coord) %>%
    filter(abs(dist.to.tss) == min(abs(dist.to.tss)))

  terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, sep = ";")
  terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)

  return(list('regions.annotated' = regions.annotated, 'out2.df.closest' = out2.df.closest))
}


AnnotateBins2 <- function(terms.mat, top.thres=0.995, inf.tss="/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep=c(paste("chr", seq(19), sep = ""), "chrX", "chrY"), skip.split = FALSE){
  assertthat::assert_that(file.exists(inf.tss))
  assertthat::assert_that(class(terms.mat) == "matrix")
  regions <- data.frame(seqnames = sapply(colnames(terms.mat), GetChromo),
                        start = sapply(colnames(terms.mat), GetStart),
                        end = sapply(colnames(terms.mat), GetEnd),
                        stringsAsFactors = FALSE)
  rownames(regions) <- colnames(terms.mat)
  # regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))
  print("Chromos to keep")
  print(chromos.keep)
  regions <- subset(regions, seqnames %in% chromos.keep)

  regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
  regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                  # TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                  TxDb=txdb,
                                                  # annoDb='org.Mm.eg.db'))
                                                  annoDb=annodb))
  regions.annotated$region_coord <- names(regions.range)

  topic.regions <- lapply(seq(nrow(terms.mat)), function(clst){
    return(SelectTopRegions(terms.mat[clst, ], colnames(terms.mat), method = "thres", method.val = top.thres))
  })

  print(paste("Using TSS definitions from:", inf.tss))
  terms.long <- data.frame(term = colnames(terms.mat), as.data.frame(t(terms.mat)), stringsAsFactors = FALSE) %>%
    gather(key = "topic", value = "weight", -term) %>%
    mutate(topic = gsub("X", "", topic)) %>%
    group_by(topic) %>%
    arrange(desc(weight)) %>%
    mutate(rnk = seq(length(weight))) %>%
    rowwise()
  terms.filt.top <- terms.long %>%
    # filter(rnk < 1000) %>%  # DO GENOME WIDE
    rowwise()
  tss.dat <- fread(inf.tss, col.names = c("seqnames", "start", "end", "tssname"))
  tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, ";")[[1]][[2]])

  annots.biomart <- regions.annotated %>%
    mutate(midpt = start + (end - start) / 2) %>%
    filter(region_coord %in% terms.filt.top$term)

  annots.gr <- makeGRangesFromDataFrame(annots.biomart %>% dplyr::select(seqnames, start, end, SYMBOL, region_coord, ENSEMBL), keep.extra.columns = TRUE)
  annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)

  out <- findOverlaps(annots.tss.gr, annots.gr, type = "within")
  out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")

  out2.df = data.frame(annots.gr[queryHits(out2),], annots.tss.gr[subjectHits(out2),]) %>%
    mutate(midpt = start + round(width / 2),
           midpt.1 = start.1 + round(width.1 / 2),
           dist.to.tss = midpt.1 - midpt)

  # filter closest
  out2.df.closest <- out2.df %>%
    group_by(region_coord) %>%
    filter(abs(dist.to.tss) == min(abs(dist.to.tss)))

  terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, sep = ";")
  terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)

  terms.annot <- terms.filt.top %>%
    mutate(termgene = ifelse(!is.null(terms.hash[[term]]), terms.hash[[term]], NA))

  terms.filt <- terms.filt.top %>%
    mutate(termgene = ifelse(!is.null(terms.hash[[term]]), terms.hash[[term]], NA)) %>%
    filter(!is.na(termgene))
  if (!skip.split){
    terms.filt <- terms.filt %>%
      mutate(gene = sapply(termgene, function(x) strsplit(x, ";")[[1]][[2]])) %>%
      group_by(gene)  # dont filter use
  }
  return(list('topic.regions' = topic.regions, 'regions.annotated' = regions.annotated, 'terms.annot' = terms.annot, 'out2.df.closest' = out2.df.closest, 'terms.filt' = terms.filt))
}

DoUmapAndLouvain <- function(topics.mat, jsettings){
  umap.out <- umap(topics.mat, config = jsettings)
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
  dat.umap.long <- DoLouvain(topics.mat = topics.mat, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long)
  return(dat.umap.long)
}

CollapseMatByLouvains <- function(count.mat, dat.umap.longs){
  jmat <- left_join(melt(as.matrix(count.mat)), dat.umap.longs %>% dplyr::select(c(cell, louvmark)), by = c("Var2" = "cell")) %>%
    group_by(louvmark, Var1) %>%
    summarise(count = sum(value)) %>%
    dcast(data = ., Var1 ~ louvmark) %>%
    as.data.frame()
  rownames(jmat) <- jmat$Var1
  jmat$Var1 <- NULL
  return(jmat)
}

GetTmResultFromGensim <- function(inf.topics, inf.terms, inf.cellnames, inf.binnames){
  cellnames <- fread(inf.cellnames, header=FALSE)$V1
  binnames <- fread(inf.binnames, header=FALSE)$V1
  topics.mat.vi <- as.data.frame(fread(inf.topics, header=FALSE))
  terms.mat.vi <- as.data.frame(fread(inf.terms, header=FALSE))
  colnames(topics.mat.vi) <- gsub("^V", "Topic_", colnames(topics.mat.vi))
  rownames(topics.mat.vi) <- cellnames
  colnames(terms.mat.vi) <- binnames
  rownames(terms.mat.vi) <- paste("Topic_", rownames(terms.mat.vi), sep = "")

  tm.result <- list(topics = as.matrix(topics.mat.vi), terms = as.matrix(terms.mat.vi))
  return(tm.result)
}

AnnotateBins <- function(terms.mat, top.thres=0.995, inf.tss="/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knowngene, annodb = "org.Mm.eg.db"){
  # assertthat::assert_that(is.list(tm.result))  # expect terms
  # kchoose <- out.lda@k
  # tm.result <- posterior(out.lda)
  assertthat::assert_that(file.exists(inf.tss))
  assertthat::assert_that(class(terms.mat) == "matrix")
  regions <- data.frame(seqnames = sapply(colnames(terms.mat), GetChromo),
                        start = sapply(colnames(terms.mat), GetStart),
                        end = sapply(colnames(terms.mat), GetEnd),
                        stringsAsFactors = FALSE)
  rownames(regions) <- colnames(terms.mat)
  regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))

  regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
  regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                  # TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                  TxDb=txdb,
                                                  # annoDb='org.Mm.eg.db'))
                                                  annoDb=annodb))
  regions.annotated$region_coord <- names(regions.range)

  topic.regions <- lapply(seq(nrow(terms.mat)), function(clst){
    return(SelectTopRegions(terms.mat[clst, ], colnames(terms.mat), method = "thres", method.val = top.thres))
  })

  print(paste("Using TSS definitions from:", inf.tss))
  terms.long <- data.frame(term = colnames(terms.mat), as.data.frame(t(terms.mat)), stringsAsFactors = FALSE) %>%
    gather(key = "topic", value = "weight", -term) %>%
    mutate(topic = gsub("X", "", topic)) %>%
    group_by(topic) %>%
    arrange(desc(weight)) %>%
    mutate(rnk = seq(length(weight))) %>%
    rowwise()
  terms.filt.top <- terms.long %>%
    # filter(rnk < 1000) %>%  # DO GENOME WIDE
    rowwise()
  tss.dat <- fread(inf.tss, col.names = c("seqnames", "start", "end", "tssname"))
  tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, ";")[[1]][[2]])

  annots.biomart <- regions.annotated %>%
    mutate(midpt = start + (end - start) / 2) %>%
    filter(region_coord %in% terms.filt.top$term)

  annots.gr <- makeGRangesFromDataFrame(annots.biomart %>% dplyr::select(seqnames, start, end, SYMBOL, region_coord), keep.extra.columns = TRUE)
  annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)

  out <- findOverlaps(annots.tss.gr, annots.gr, type = "within")
  out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")

  out2.df = data.frame(annots.gr[queryHits(out2),], annots.tss.gr[subjectHits(out2),]) %>%
    mutate(midpt = start + round(width / 2),
           midpt.1 = start.1 + round(width.1 / 2),
           dist.to.tss = midpt.1 - midpt)

  # filter closest
  out2.df.closest <- out2.df %>%
    group_by(region_coord) %>%
    filter(abs(dist.to.tss) == min(abs(dist.to.tss)))

  terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, sep = ";")
  terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)

  terms.filt <- terms.filt.top %>%
    mutate(termgene = ifelse(!is.null(terms.hash[[term]]), terms.hash[[term]], NA)) %>%
    filter(!is.na(termgene)) %>%
    mutate(gene = sapply(termgene, function(x) strsplit(x, ";")[[1]][[2]])) %>%
    group_by(gene)  # dont filter use
  # return(terms.filt)
  # annotate terms to nearest gene
  return(list('topic.regions' = topic.regions, 'regions.annotated' = regions.annotated, 'terms.filt' = terms.filt))
}


AnnotateTermsToNearestGene <- function(terms.mat, regions.annotated, inf.tss="/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed"){
  print(paste("Using TSS definitions from:", inf.tss))
  assertthat::assert_that(file.exists(inf.tss))
  terms.long <- data.frame(term = colnames(terms.mat), as.data.frame(t(terms.mat)), stringsAsFactors = FALSE) %>%
    gather(key = "topic", value = "weight", -term) %>%
    mutate(topic = gsub("X", "", topic)) %>%
    group_by(topic) %>%
    arrange(desc(weight)) %>%
    mutate(rnk = seq(length(weight))) %>%
    rowwise()
  terms.filt.top <- terms.long %>%
    # filter(rnk < 1000) %>%  # DO GENOME WIDE
    rowwise()
  tss.dat <- fread(inf.tss, col.names = c("seqnames", "start", "end", "tssname"))
  tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, ";")[[1]][[2]])
  annots.biomart <- regions.annotated %>%
    mutate(midpt = start + (end - start) / 2) %>%
    filter(region_coord %in% terms.filt.top$term)

  annots.gr <- makeGRangesFromDataFrame(annots.biomart %>% dplyr::select(seqnames, start, end, SYMBOL, region_coord), keep.extra.columns = TRUE)
  annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)

  out <- findOverlaps(annots.tss.gr, annots.gr, type = "within")
  out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")

  out2.df = data.frame(annots.gr[queryHits(out2),], annots.tss.gr[subjectHits(out2),]) %>%
    mutate(midpt = start + round(width / 2),
           midpt.1 = start.1 + round(width.1 / 2),
           dist.to.tss = midpt.1 - midpt)

  # filter closest
  out2.df.closest <- out2.df %>%
    group_by(region_coord) %>%
    filter(abs(dist.to.tss) == min(abs(dist.to.tss)))

  terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, sep = ";")
  terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)

  terms.filt <- terms.filt.top %>%
    mutate(termgene = ifelse(!is.null(terms.hash[[term]]), terms.hash[[term]], NA)) %>%
    filter(!is.na(termgene)) %>%
    mutate(gene = sapply(termgene, function(x) strsplit(x, ";")[[1]][[2]])) %>%
    group_by(gene)  # dont filter use
  return(terms.filt)
}

OrderTopicsByEntropy <- function(tm.result, jquantile = 0.99){
  # analyze topic matrix across cells
  topics.long <- data.frame(cell = rownames(tm.result$topics), as.data.frame(tm.result$topics)) %>%
    gather(key = "topic", value = "weight", -cell) %>%
    rowwise() %>%
    group_by(topic) %>%
    mutate(zscore = scale(weight, center = TRUE, scale = TRUE))
  topics.sum <- topics.long %>%
    group_by(topic) %>% # do entropy on 1 to 99% of cells
    filter(zscore < quantile(zscore, jquantile)) %>%
    mutate(zscore.prob = exp(zscore) / sum(exp(zscore))) %>%
    summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
    arrange(entropy)
  return(topics.sum)
}

GetPeaksFromGene <- function(jgene, regions.annot, dist = 50000){
  jsub <- subset(regions.annot, grepl(jgene, SYMBOL)) %>% arrange(abs(distanceToTSS)) %>% filter(distanceToTSS <= dist)
  jpeaks <- jsub$region_coord
  return(list(regions.sub = jsub, peaks = jpeaks))
}

SelectBestPeak <- function(jpeaks, regions.annot = NULL, tm.result){
  if (length(jpeaks) == 1){
    warning("Only one peak entered, returning only option")
    # no need to select
    return(jpeaks)
  }
  terms.sub <- tm.result$terms[, jpeaks]
  jmax <- apply(terms.sub, 2, max)
  jmax.i <- apply(terms.sub, 2, which.max)  # which topic max occurs. Often they agree.

  print(unique(jmax.i))

  jpeak <- jpeaks[which(jmax == max(jmax))]
  return(jpeak)
}

PlotAllMarks <- function(jgene, jpeak, jmarks, out.objs, custom.settings){
  out <- lapply(jmarks, function(jmark){
    print(jmark)
    m <- PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmark, show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene)
  })
  return(out)
}


LoadLDABins <- function(jmark, jbin=TRUE, top.thres=0.995, inf = NULL, convert.chr20.21.to.X.Y = TRUE, add.chr.prefix = FALSE, choose.k = "auto"){
  # jbin <- "TRUE"
  # top.thres <- 0.995
  if (is.null(inf)){
    if (jbin){
      inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin, "/lda_out_meanfilt.BM-", jmark, ".CountThres0.K-5_10_15_20_25.Robj")
    } else {
      inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin, "/lda_out_meanfilt.BM-", jmark, ".CountThres0.K-5_15_25.Robj")
    }
  }
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  if (choose.k == "auto"){
    out.lda <- ChooseBestLDA(out.lda)
  } else {
    kvec <- lapply(out.lda, function(x) x@k)
    i <- which(kvec == as.numeric(choose.k))
    out.lda <- out.lda[[i]]
  }


  if (add.chr.prefix){
    out.lda@terms <- paste0("chr", out.lda@terms)
  }

  kchoose <- out.lda@k
  tm.result <- posterior(out.lda)

  if (convert.chr20.21.to.X.Y){
    colnames(tm.result$terms) <- gsub("chr20", "chrX", colnames(tm.result$terms))
    colnames(tm.result$terms) <- gsub("chr21", "chrY", colnames(tm.result$terms))
  }

  regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
                        start = sapply(colnames(tm.result$terms), GetStart),
                        end = sapply(colnames(tm.result$terms), GetEnd),
                        stringsAsFactors = FALSE)
  rownames(regions) <- colnames(tm.result$terms)
  regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))

  regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
  regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                  TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                  annoDb='org.Mm.eg.db'))
  regions.annotated$region_coord <- names(regions.range)

  topic.regions <- lapply(seq(kchoose), function(clst){
    return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
  })

  return(list('out.lda' = out.lda, 'tm.result' = tm.result, 'topic.regions' = topic.regions, 'regions.annotated' = regions.annotated, 'count.mat' = count.mat))
}

GetVar <- function(tm.result, regions.annotated){
  mat.norm <- t(tm.result$topics %*% tm.result$terms)
  # find peaks with largest range
  maxmin <- sort(apply(mat.norm, 1, function(jcol) mad(jcol)), decreasing = TRUE)
  # what are the sizes?
  peak.size <- sapply(names(maxmin), function(x) ParseCoord(x)$end - ParseCoord(x)$start)
  # normalize maxmin by peaksize
  maxmin.norm <- maxmin / peak.size
  dat.var <- data.frame(Var = maxmin, peak.size = peak.size, coord = names(maxmin))
  regions.annotated <- dplyr::left_join(regions.annotated, dat.var, by = c("region_coord"="coord"))
  return(regions.annotated)
}

ChooseBestLDA <- function(out.lda){
  # pick best k
  if (length(out.lda) > 1){
    out.lda.lst <- out.lda
    # we did multicore, so which one to choose?
    Kvec <- sapply(out.lda.lst, function(x) x@k)
    best.K <- Kvec[which.max(sapply(out.lda.lst, function(x) x@loglikelihood))]
    # plot loglikelihood
    par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    plot(Kvec, sapply(out.lda.lst, function(x) x@loglikelihood), 'o')
    kchoose <- best.K
    out.lda <- out.lda.lst[[which(Kvec == kchoose)]]
    print(paste("Likelihood: ", out.lda@loglikelihood))
  } else {
    kchoose <- out.lda@k
  }
  return(out.lda)
}


GetUmapSettings <- function(nn, jmetric, jmindist, seed=123){
  # nn <- 5
  # jmetric <- 'euclidean'
  # # jmetric <- 'cosine'
  # jmindist <- 0.1
  custom.settings <- umap.defaults
  custom.settings$n_neighbors <- nn
  custom.settings$metric <- jmetric
  custom.settings$min_dist <- jmindist
  custom.settings$random_state <- seed
  return(custom.settings)
}

# color by loadings on Kvec
ColorsByGamma <- function(topic, topics.mat, cols.vec = c("pink", "red", "darkred")){
  # jcol <- out.lda@gamma[, topic]
  # jcol <- tmResult$topics[, topic]
  jcol <- topics.mat[, topic]
  colorPal <- grDevices::colorRampPalette(cols.vec)
  jcol.rgb <- colorPal(200)[as.numeric(cut(jcol,breaks = 200))]
  return(jcol.rgb)
}

ColorsByCounts <- function(counts.vec, nbreaks=100, colvec = c("pink", "red", "darkred")){
  colorPal <- grDevices::colorRampPalette(colvec)
  jcol.rgb <- colorPal(nbreaks)[as.numeric(cut(counts.vec,breaks = nbreaks))]
  return(jcol.rgb)
}

SelectTopRegions <- function(beta.row, regions, method = "thres", method.val = 0.01){
  # take betas (in log scale) and select top X fraction
  # or do simple cutoff, set method = "cutoff"
  if (method == "cutoff"){
    return(regions[which(beta.row > method.val)])
  } else if (method == "thres"){
    return(regions[which(beta.row > quantile(beta.row, method.val))])
  } else {
    stop(paste("Method", method, "not yet implemented"))
  }
}

GetCountMatFromLDA <- function(out.lda){
  # https://stackoverflow.com/questions/20004493/convert-simple-triplet-matrixslam-to-sparse-matrixmatrix-in-r
  count.mat <- Matrix::sparseMatrix(i=out.lda@wordassignments$i,
                                  j=out.lda@wordassignments$j,
                                  x=out.lda@wordassignments$v,
                                  dims=c(out.lda@wordassignments$nrow, out.lda@wordassignments$ncol))
  return(count.mat)
}


.modelMatSelection <- function(
  # from cisTopics package
  object,
  target,
  method,
  all.regions=FALSE
){
  # Check info
  if (length(object@selected.model) < 1){
    stop('Please, run selectModel() first.')
  }

  if (target == 'cell'){
    if (method == 'Z-score'){
      modelMat <- scale(object@selected.model$document_expects, center=TRUE, scale=TRUE)
    }

    else if (method == 'Probability'){
      alpha <- object@calc.params[['runModels']]$alpha/length(object@selected.model$topic_sums)
      modelMat <- apply(object@selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
    }
    else{
      stop('Incorrect method selected. Chose method between "Z-score" and "Probability".')
    }
    colnames(modelMat) <- object@cell.names
    rownames(modelMat) <- paste0('Topic', 1:nrow(modelMat))
  }

  else if (target == 'region'){
    if (!all.regions){
      if (length(object@binarized.cisTopics) < 1){
        stop('Please, use binarizecisTopics() first for defining the high confidence regions for dimensionality reduction!')
      }
      else {
        regions <- unique(unlist(lapply(object@binarized.cisTopics, rownames)))
      }
    }

    topic.mat <- object@selected.model$topics

    if (method == 'NormTop'){
      normalizedTopics <- topic.mat/(rowSums(topic.mat) + 1e-05)
      modelMat <- apply(normalizedTopics, 2, function(x) x * (log(x + 1e-05) - sum(log(x + 1e-05))/length(x)))
    }

    else if (method == 'Z-score'){
      modelMat <- scale(object@selected.model$topics, center=TRUE, scale=TRUE)
    }

    else if (method == 'Probability'){
      beta <- object@calc.params[['runModels']]$beta
      topic.mat <- object@selected.model$topics
      modelMat <-  (topic.mat + beta)/rowSums(topic.mat + beta)
    }

    else{
      stop('Incorrect method selected. Chose "NormTop", "Z-score" and "Probability".')
    }

    colnames(modelMat) <- object@region.names
    rownames(modelMat) <- paste0('Topic', 1:nrow(modelMat))

    if (!all.regions){
      modelMat <- modelMat[,regions]
    }
  }

  else{
    stop('Please, provide target="cell" or "region".')
  }

  return(modelMat)
}
