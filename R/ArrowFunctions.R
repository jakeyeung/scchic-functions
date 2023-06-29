GetArrows <- function(dat.pca.merge.wide.filt, grid.n = 40, jfactor = 1){
  # grid.n  <- 25
  # jfactor <- 0.5
  pos <- as.data.frame(subset(dat.pca.merge.wide.filt, select = c(pc1.x, pc2.x)))
  ppos <- as.data.frame(subset(dat.pca.merge.wide.filt, select = c(pc1.y, pc2.y)))
  
  # arrow estimates for each cell
  ars <- data.frame(pos[,1],pos[,2],ppos[,1],ppos[,2])
  colnames(ars) <- c('x0','y0','x1','y1')
  arsd <- data.frame(xd=ars$x1-ars$x0,yd=ars$y1-ars$y0)
  rownames(ars) <- rownames(arsd) <- rownames(pos);
  
  rownames(pos) <- dat.pca.merge.wide.filt$cell
  
  rx <- range(c(range(ars$x0),range(ars$x1)))
  ry <- range(c(range(ars$y0),range(ars$y1)))
  gx <- seq(rx[1],rx[2],length.out=grid.n)
  gy <- seq(ry[1],ry[2],length.out=grid.n)
  
  
  grid.sd <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)/2
  min.arrow.size <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)*1e-2
  
  min.grid.cell.mass <- 1
  garrows <- do.call(rbind,lapply(gx,function(x) {
    # cell distances (rows:cells, columns: grid points)
    cd <- sqrt(outer(pos[,2],-gy,'+')^2 + (x-pos[,1])^2)
    cw <- dnorm(cd,sd=grid.sd)
    # calculate x and y delta expectations
    gw <- Matrix::colSums(cw)
    cws <- pmax(1,Matrix::colSums(cw));
    gxd <- jfactor * Matrix::colSums(cw*arsd$xd)/cws
    gyd <- jfactor * Matrix::colSums(cw*arsd$yd)/cws
    
    al <- sqrt(gxd^2+gyd^2);
    vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
    
    out <- cbind(rep(x, sum(vg)), gy[vg], x+gxd[vg], gy[vg]+gyd[vg])
  }))
  colnames(garrows) <- c('pc1.x','pc2.x','pc1.y','pc2.y')
  
  print(dim(garrows))
  dat.garrows <- data.frame(garrows)
  return(dat.garrows)
}
