DoLouvain <- function(topics.mat, custom.settings.louv, dat.umap.long = NULL, clstr.cname = "louvain"){
  # Do Louvain for clustering
  dat.umap.louv <- umap(topics.mat, config = custom.settings.louv)
  dat.umap.louv.long <- data.frame(umap1 = dat.umap.louv$layout[, 1], umap2 = dat.umap.louv$layout[, 2], cell = rownames(dat.umap.louv$layout),
                                   stringsAsFactors = FALSE)
  cell.indx <- hash(rownames(dat.umap.louv$knn$indexes), dat.umap.louv$knn$indexes[, 1])
  cell.indx.rev <- hash(dat.umap.louv$knn$indexes[, 1], rownames(dat.umap.louv$knn$indexes))
  nr <- nrow(dat.umap.louv$knn$indexes)
  nc <- ncol(dat.umap.louv$knn$indexes)
  edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
  colnames(edgelist) <- c("from", "to")
  for (vertex.i in seq(nr)){
    istart <- nc*(vertex.i - 1)+1
    iend <- nc*vertex.i
    edgelist[istart : iend, 1] <- cell.indx.rev[[as.character(vertex.i)]]
    edgelist[istart : iend, 2] <- sapply(dat.umap.louv$knn$indexes[vertex.i, 1:nc], function(x) cell.indx.rev[[as.character(x)]])
    # edgelist[istart : iend, 3] <- 1 / (dat.umap$knn$distances[vertex.i, 1:nc] + 0.1)
  }
  g <- graph_from_data_frame(edgelist, directed=FALSE)
  g.out <- cluster_louvain(g, weights = NULL)
  V(g)$color <- g.out$membership
  clstr <- hash(g.out$names, g.out$membership)
  if (is.data.frame(dat.umap.long)){
    # dat.umap.long$louvain <- as.character(sapply(dat.umap.long$cell, function(x) clstr[[as.character(x)]]))
    if (clstr.cname == "louvain"){
      dat.umap.long[[clstr.cname]] <- as.character(sapply(dat.umap.long$cell, function(x) clstr[[as.character(x)]]))
    } else {
      #add louvain prefix
      dat.umap.long[[clstr.cname]] <- as.character(sapply(dat.umap.long$cell, function(x) paste0("louvain_", clstr[[as.character(x)]])))
    }
  } else {
    dat.umap.long <- clstr
  }
  return(dat.umap.long)
}

DoLeiden <- function(topics.mat, K=30, res_param=10^seq(-5, 0), dat.umap.long = NULL, random_seed = 123, weight=FALSE, verbose = TRUE){
  # requires jclustering, leidenbase
  pdat <- data.frame(cell = rownames(topics.mat))
  rownames(pdat) <- pdat$cell
  rownames(pdat) <- pdat$cell
  l.out <- leiden_clustering(topics.mat, pdat, random_seed = random_seed, weight = weight, verbose=verbose, resolution_parameter = res_param, k = K)
  dat.leiden <- data.frame(cell = names(l.out$optim_res$membership), cluster = paste("leiden", as.character(l.out$optim_res$membership), sep = "_"))
  if (!is.null(dat.umap.long)){
    dat.leiden <- left_join(dat.umap.long, dat.leiden)
  }
  return(dat.leiden)
}

