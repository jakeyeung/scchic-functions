# Jake Yeung
# Date of Creation: 2020-09-13
# File: ~/projects/scchic-functions/R/PlateLayoutFunctions.R
# 

AnnotateSortFromLayout.dat <- function(jdat){
  # expects cell in colname
  jdat <- bind_rows(jdat) %>%
    rowwise() %>%
    mutate(experi = ClipLast(cell, jsep = "_"), 
           plate = as.numeric(strsplit(experi, "-")[[1]][[6]]),
           rowcoord = AddPlateCoordinates(cell)$rowcoord,
           colcoord = AddPlateCoordinates(cell)$colcoord,
           stype = AnnotateSortFromLayout(plate, rowcoord, colcoord))
  return(jdat)
}


AnnotateSortFromLayout <- function(plate, rowcoord, colcoord){
  assertthat::assert_that(is.numeric(plate))
  assertthat::assert_that(is.numeric(rowcoord))
  assertthat::assert_that(is.numeric(colcoord))
  if (rowcoord >= 1 & rowcoord <= 8 & colcoord == 1){
    print("Empty cell")
    warning("Empty cell")
    ctype <- "Empty"
    return(ctype)
  }
  if (plate >= 1 & plate <= 7){
    if (colcoord >= 1 & colcoord <= 11){
      ctype <- "Unenriched"
    } else if (colcoord >= 12 & colcoord <= 18){
      ctype <- "LinNeg"
    } else if (colcoord >= 19 & colcoord <= 24){
      ctype <- "LSK"
    }
    
  } else if (plate >= 8 & plate <= 13){
    if (colcoord >= 1 & colcoord <= 12){
      ctype <- "Unenriched"
    } else if (colcoord >= 13 & colcoord <= 24){
      ctype <- "LinNeg"
    }
  } else {
    print(paste("Unknown plate:", plate))
    ctype <- "Unknown"
  }
  return(ctype)
}


