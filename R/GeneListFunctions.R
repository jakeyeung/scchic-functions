# Jake Yeung
# Date of Creation: 2020-06-09
# File: ~/projects/scchic-functions/R/GeneListFunctions.R
# Hardcoded gene names for easy access


# Check other genes -------------------------------------------------------

GetWKMgenes <- function(inf.kob = "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Kobayashi_2019_zf_genelists/hspc_genes_list.txt"){
  
  genesets.zf <- list()
  # genesets.ens.zf <- list()
  
  kobayashi.genes <- data.table::fread(inf.kob, header = FALSE)$V1
  
  # HSPCs
  jgenes.choose <- c("meis1b", "pmp22b", "ahnak", "krt8", "anxa5b", "mrc1b", "pa2g4a", "npm1a", "ahcy", "adh5", "fabp3", "myb")
  # add kobayashi?
  jgenes.choose <- c(jgenes.choose, kobayashi.genes)
  genesets.zf[["HSPC"]] <- jgenes.choose
  # jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)
  
  # lymphs
  jgenes.choose <- c("pax5", "cd79a", "bhlhe40", "cd83", "cxcr4a", "cd74b", "cd74a", "CD37", "zfp36l1a")  # lymphs
  genesets.zf[["Blymph"]] <- jgenes.choose
  # jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)
  
  # monocytes
  jgenes.choose1 <- c("adam8a", "odc1", "lta4h", "thy1", "scpp8", "illr4", "timp2b", "mmp9", "mmp13a", "scinlb")  # monocytes
  jgenes.choose2 <- c("cpa5", "lyz", "lect2l", "npsn", "sms", "abcb9", "ch25hl2", "papss2b", "hsd3b7", "cfd")  # neutros
  jgenes.choose <- c(jgenes.choose1, jgenes.choose2)
  genesets.zf[["granuMono"]] <- jgenes.choose
  # jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)
  
  # eryths
  jgenes.choose <- c("rhag", "prdx2", "epor", "gata1a", "tspo", "slc4a1a", "sptb", "cahz", "ba1", "alas2", "epb41b", "nt5c2l1")
  genesets.zf[["eryth"]] <- jgenes.choose
  
  genesets.ens.zf <- lapply(genesets.zf, function(gset){
    JFuncs::Gene2Ensembl.ZF(gset, return.original=TRUE, species="drerio")
  })
  return(list(gset.zf = genesets.zf, gset.ens.zf = genesets.ens.zf))
}

