
GetPlateCoord2 <- function(cell, platecols = 24, is.zero.base = FALSE){
  # cell0 -> 1,1
  indx <- as.numeric(strsplit(cell, "cell")[[1]][[2]])
  if (is.zero.base){
    indx <- indx + 1
  }
  jcol <- indx %% platecols
  jcol <- ifelse(jcol == 0, platecols, jcol)
  jrow <- ceiling(indx / platecols)
  return(c(jrow, jcol))
}

BatchColumn2Ctype <- function(Batch, Column, Row = NA){
  
  if (Batch == "SL1"){
    if (Column >= 1 & Column <= 4){
      ctype <- "Tcells"
    } else if (Column >= 5 & Column <= 8){
      ctype <- "Bcells"
    } else if (Column >= 9 & Column <= 12){
      ctype <- "NKs"
    } else if (Column >= 13 & Column <= 16){
      ctype <- "Eryths"
    } else if (Column >= 17 & Column <= 24)
      ctype <- "AllCells"
  }
  
  if (Batch == "SL2"){
    if (Column >= 1 & Column <= 4){
      ctype <- "Granulocytes"
    } else if (Column >= 5 & Column <= 8){
      ctype <- "Monocytes"
    } else if (Column >= 9 & Column <= 12){
      ctype <- "DCs"
    } else if (Column >= 13 & Column <= 24){
      ctype <- "AllCells"
    }
  }
  
  if (Batch == "SL3"){
    if (Column == 1){
      ctype <- "HSCs"
    } else if (Column == 2){
      ctype <- "MPPs"
    } else if (Column == 3){
      ctype <- "LT"
    } else if (Column == 4){
      ctype <- "ST"
    } else if (Column >= 5 & Column <= 24){
      ctype <- "LSK"
    } else {
      warning("Error: unknown column:", Column)
    }
  }
  
  if (Batch == "SL4"){
    if (Column == 1){
      ctype <- "GMP"
    } else if (Column == 2){
      ctype <- "CMP"
    } else if (Column == 3){
      ctype <- "MEP"
    } else if (Column >= 4 & Column <= 24){
      ctype <- "LSK"
    }
  }
  
  if (Batch == "SL5"){
    if (Column == 1 & Row >= 1 & Row <= 8){
      ctype <- "pDCs"
    } else if (Column == 1 & Row >= 9 & Row <= 16){
      ctype <- "IL7RLinNeg"
    } else if (Column >= 2 & Column <= 12){
      ctype <- "LinNeg"
    } else if (Column >= 13 & Column <= 24){
      ctype <- "LSK"
    }
  }
  
  return(ctype)
}
