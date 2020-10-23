# Jake Yeung
# Date of Creation: 2020-09-13
# File: ~/projects/scchic-functions/R/PlateLayoutFunctions.R
# 

SplitGetLast <- function(x, jsplit = "-"){
  # get last element after splitting
  xsplit <- strsplit(x, jsplit)[[1]]
  xnew <- xsplit[[length(xsplit)]]
  return(xnew)
}

GetRepBM <- function(experiname){
  if (grepl("rep3", experiname)){
    return("rep3")
  } else if (grepl("rep2", experiname)){
    return("rep2")
  } else {
    print("Neither rep2 nor rep3 in name: Returning rep2")
    return("rep2")
  }
}


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



AnnotateSortFromLayoutBMall <- function(plate, rowcoord, colcoord, jrep, jmark){
  assertthat::assert_that(is.numeric(plate))
  assertthat::assert_that(is.numeric(rowcoord))
  assertthat::assert_that(is.numeric(colcoord))
  # assertthat::assert_that(is.numeric(jrep))
  
  if (jrep == "rep3"){
    ctype <- AnnotateSortFromLayout(plate = plate, rowcoord = rowcoord, colcoord = colcoord)
  } else if (jrep == "rep2"){
    ctype <- AnnotateSortFromLayout.rep2(plate = plate, rowcoord = rowcoord, colcoord = colcoord, jmark = jmark)
  }
  return(ctype)
}



AnnotateSortFromLayout.rep2 <- function(plate, rowcoord, colcoord, jmark){
  assertthat::assert_that(is.numeric(plate))
  assertthat::assert_that(is.numeric(rowcoord))
  assertthat::assert_that(is.numeric(colcoord))
  if (rowcoord >= 1 & rowcoord <= 8 & colcoord == 1){
    print("Empty cell")
    warning("Empty cell")
    ctype <- "Empty"
    return(ctype)
  }
  if ( (plate == 1 | plate == 2) & jmark == "H3K4me1" ){
    if (colcoord >= 1 & colcoord <= 11){
      ctype <- "Unenriched"
    } else if (colcoord >= 12 & colcoord <= 18){
      ctype <- "LinNeg"
    } else if (colcoord >= 19 & colcoord <= 24){
      ctype <- "LSK"
    }
  } else {
    if (colcoord >= 1 & colcoord <= 11){
        ctype <- "Unenriched"
      } else if (colcoord >= 12 & colcoord <= 19){
        ctype <- "LinNeg"
      } else if (colcoord >= 20 & colcoord <= 24){
        ctype <- "LSK"
      }
  }
  # } else if (plate >= 3 & plate <= 6) {
  #   if (colcoord >= 1 & colcoord <= 11){
  #     ctype <- "Unenriched"
  #   } else if (colcoord >= 12 & colcoord <= 19){
  #     ctype <- "LinNeg"
  #   } else if (colcoord >= 20 & colcoord <= 24){
  #     ctype <- "LSK"
  #   }
  # } else {
  #   print(paste("Unknown plate:", plate))
  #   ctype <- "Unknown"
  # }
  return(ctype)
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


