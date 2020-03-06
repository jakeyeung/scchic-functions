# Jake Yeung
# Date of Creation: 2020-03-06
# File: ~/projects/scchic-functions/R/PublicDataFuncs.R
# 


MarkerToCelltype <- function(){
  # get correspondance between marker to celltype
  x <- list("Car1" = "Erythroblast",
            "core" = "HSCs",
            "Vpreb1" = "Bcell",
            "Siglech" = "pDendritic",
            "Prg2" = "Eosinophil",
            "Gstm1" = "Neutrophil",
            "Ly86" = "Monocyte",
            "Ccl5" = "NKcell",
            "Prss34" = "Basophil",
            "Cd74" = "cDendritic",
            "Pf4" = "Megakaryocyte",
            "Fcrla" = "Bcell",
            "Fcnb" = "Neutrophil",
            "Hba.a2" = "Erythroblast",
            "Ltf" = "Neutrophil")
  return(x)
}

