# Jake Yeung
# Date of Creation: 2020-12-16
# File: ~/projects/scchic-functions/R/GenesetFunctions.R
# 


GetGeneSets <- function(inf.genes.k4me1 = "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1/geneset_H3K4me1_TSS_topics2.filt.colorcoded2.txt", 
                        inf.genes.k4me3 = "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/celltype_specific_genes_defined_by_K4me3_TSS.txt"){
  
  # inf.genes.k4me1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1/geneset_H3K4me1_TSS_topics2.filt.colorcoded2.txt"
  dat.genes.k4me1 <- fread(inf.genes.k4me1)
  
  # inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/celltype_specific_genes_defined_by_K4me3_TSS.txt")
  dat.genes <- fread(inf.genes.k4me3)
  
  dat.genes.all <- bind_rows(dat.genes, subset(dat.genes.k4me1, select = c(-topic)))
  
  dat.genes.sub <- subset(dat.genes, !(jset == "Basophils" & rnk > 100))
  dat.genes.sub <- subset(dat.genes.sub, !(jset == "pDCs" & rnk > 100))
  
  dat.genes.sub <- subset(dat.genes, rnk < 150)
  
  # swap basophils
  dat.genes.sub <- subset(dat.genes.sub, jset != "Basophils" & jset != "pDCs")
  
  basos.k4me3 <- subset(dat.genes, jset %in% c("Basophils"))$gene
  basos.k4me1 <- subset(dat.genes.k4me1, jset %in% c("Basophils"))$gene
  
  pdcs.k4me3 <- subset(dat.genes, jset %in% c("pDCs"))$gene
  pdcs.k4me1 <- subset(dat.genes.k4me1, jset %in% c("pDCs"))$gene
  
  basos.manual <- c("Il4", "Il6", "Cpa3", "Il1r1")
  basos.rname <- subset(dat.genes.all, symbol %in% basos.manual)$gene
  basos.both <- c(intersect(basos.k4me1, basos.k4me3), basos.rname)
  pdcs.both <- intersect(pdcs.k4me1, pdcs.k4me3)
  
  dat.to.add <- subset(dat.genes.all, gene %in% c(basos.both, pdcs.both))
  
  dat.genes.sub.join <- bind_rows(dat.genes.sub, dat.to.add)
  
  dat.genes.sub.join$jset <- factor(dat.genes.sub.join$jset, levels = c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs"))
  
  dat.genes.sub.join <- dat.genes.sub.join %>%
    arrange(jset, rnk)
  
  return(dat.genes.sub.join)
  
}



# H3K9me3 bin function ----------------------------------------------------


GetParamsWideFormat <- function(params.long.filt, jvalue.var = "estimate"){
  # params.dat.wide <- data.table::dcast(subset(params.lst$H3K9me3, bin %in% k9.bins.names), formula = bin ~ param, value.var = "estimate") %>%
  params.dat.wide <- data.table::dcast(params.long.filt, formula = bin ~ param, value.var = jvalue.var) %>%
    rowwise() %>%
    mutate(ClusterBcells.Estimate = ClusterBcells.Estimate / log(2),
           ClusterGranulocytes.Estimate = ClusterGranulocytes.Estimate / log(2),
           ClusterEryths.Estimate = ClusterEryths.Estimate / log(2),
           ClusterBcells.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
           ClusterEryths.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
           ClusterBcells.Estimate_ClusterEryths.Estimate =  mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
           Bcells.effect = ClusterBcells.Estimate - ClusterEryths.Estimate_ClusterGranulocytes.Estimate,
           Eryths.effect = ClusterEryths.Estimate - ClusterBcells.Estimate_ClusterGranulocytes.Estimate,
           Granulocytes.effect = ClusterGranulocytes.Estimate - ClusterBcells.Estimate_ClusterEryths.Estimate,
           HSPCs.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))
  return(params.dat.wide)
}

GetK9CelltypeBins <- function(params.dat.wide, low.in.k9 = TRUE, keeptop = 150){
  if (low.in.k9){
    jsort.hspcs <- params.dat.wide %>%
      group_by(bin) %>%
      # arrange(HSPCs.effect)
      arrange(desc(HSPCs.effect))
    jbins.hspcs <- jsort.hspcs$bin[1:keeptop]
    
    jsort.bcell <- params.dat.wide %>%
      group_by(bin) %>%
      # arrange(desc(Bcells.effect)) 
      arrange(Bcells.effect)
    jbins.bcell <- jsort.bcell$bin[1:keeptop]
    
    jsort.granu <- params.dat.wide %>%
      group_by(bin) %>%
      # arrange(desc(Granulocytes.effect))
      arrange(Granulocytes.effect)
    jbins.granu <- jsort.granu$bin[1:keeptop]
    
    jsort.eryth <- params.dat.wide %>%
      group_by(bin) %>%
      # arrange(descEryths.effect)) 
      arrange(Eryths.effect)
    jbins.eryth <- jsort.eryth$bin[1:keeptop]
  } else {
    jsort.hspcs <- params.dat.wide %>%
      group_by(bin) %>%
      arrange(HSPCs.effect)
    jbins.hspcs <- jsort.hspcs$bin[1:keeptop]
    
    jsort.bcell <- params.dat.wide %>%
      group_by(bin) %>%
      arrange(desc(Bcells.effect))
    jbins.bcell <- jsort.bcell$bin[1:keeptop]
    
    jsort.granu <- params.dat.wide %>%
      group_by(bin) %>%
      arrange(desc(Granulocytes.effect))
    jbins.granu <- jsort.granu$bin[1:keeptop]
    
    jsort.eryth <- params.dat.wide %>%
      group_by(bin) %>%
      arrange(desc(Eryths.effect))
    jbins.eryth <- jsort.eryth$bin[1:keeptop]
  }
  
  
  
  bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)
  
  bins.keep.lst <- list("Eryths" = jbins.eryth,
                        "Bcells" = jbins.bcell,
                        "Granulocytes" = jbins.granu,
                        "HSPCs" = jbins.hspcs)
  return(bins.keep.lst)
}


