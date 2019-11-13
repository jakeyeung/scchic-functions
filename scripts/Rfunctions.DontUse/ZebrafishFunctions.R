GetTimeHpf <- function(cell.indx, start.indx = 0, nb.cols = 24){
  # columns 1 to 12 are hpf 6, 13 to 24 are hpf 8, continues for every row up to 384
  cellname <- paste0("cell", cell.indx)
  
  rowcol <- GetPlateCoord(cellname, platecols = nb.cols, is.zero.base = TRUE)
  jtime <- ifelse(rowcol[[2]] <= nb.cols / 2, "6hpf", "8hpf")
  return(jtime)
}
