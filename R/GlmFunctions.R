# Jake Yeung
# Date of Creation: 2020-01-15
# File: ~/projects/scchic-functions/R/GlmFunctions.R
# Mostly stolen from Will Townes 
# https://github.com/willtownes/scrna2019/blob/master/util/functions.R

Norm<-function(v){sqrt(sum(v^2))}

ColNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  apply(x,2,Norm)
}


l2norms.frac <- function(factors){
  l2norms <- ColNorms(factors)
  l2norms.frac <- l2norms / sum(l2norms)
  return(l2norms.frac)
}

# l2norms.lst <- funlapply(factors, function(jmark){
#   jfactors <-glmpca.factors[[jmark]]
#   l2norms <- colNorms(jfactors)
#   l2norms.frac <- l2norms / sum(l2norms)
#   return(signif(l2norms.frac, digits = 2) * 100)
# })



InitGLMPCAfromLDA <- function(count.mat, tm.result, dat.var.merge, covar.cname = "ncuts.var", bins.keep = 100, do.log = FALSE, svd.on.Yinit = TRUE, use.orig.sz = TRUE){
  
  ntopics <- ncol(tm.result$topics)
  Y <- count.mat
  
  if (use.orig.sz){
    size.factor <- colSums(Y)
  } else {
    size.factor <- NULL
  }
  
  # print("Size factor:")
  # print(head(size.factor))
  
  assertthat::assert_that(all(rownames(Y) == colnames(tm.result$terms)))
  
  # keep variable bins
  if (bins.keep > 0){
    print(paste("Keeping only top bins bins.keep:", bins.keep))
    bins.high.i <- as.data.frame(apply(tm.result$terms, MARGIN = 1, function(jcol) order(jcol, decreasing = TRUE)[1:bins.keep])) %>%
      unlist()
    bins.high <- unique(rownames(Y)[bins.high.i])
  } else {
    print(paste("Keeping all bins because bins.keep=", bins.keep))
    bins.high <- colnames(tm.result$terms)
  }
  Y.filt <- Y[bins.high, ]
  
  var.vec <- dat.var.merge[[covar.cname]]
  names(var.vec) <- dat.var.merge$cell
  
  X <- data.frame(covariate = var.vec, cell = names(var.vec))  # columns of 1s are implicit
  X.reorder <- X[match(colnames(Y.filt), X$cell), ]
  X.mat <- matrix(data = X.reorder$covariate, ncol = 1, byrow = TRUE, dimnames = list(X.reorder$cell, covar.cname))
  
  # factors (U) is c by k
  # loadings (V) is g by k
  if (do.log){
    topics.mat <- log2(tm.result$topics)
    terms.mat <- log2(tm.result$terms[, bins.high])
  } else {
    topics.mat <- tm.result$topics
    terms.mat <- tm.result$terms[, bins.high]
  }
  
  if (svd.on.Yinit){
    # GLM loglink function for multinom is log( p / (1 - p) )
    # V %*% t(U) on init matrix gives an estimate of p
    # after estimate p / (1 - p), THEN remove mean and finally do SVD to get factors and loadings estimate
    p <- t(terms.mat) %*% t(topics.mat)
    logodds <- log(p / (1 - p))
    # remove mean and SVD
    logodds.centered <- t(scale(t(logodds), center = TRUE, scale = FALSE))
    # logodds.centered.check <- sweep(logodds, MARGIN = 1, STATS = rowMeans(logodds), FUN = "-")
    logodds.pca <- prcomp(t(logodds.centered), center = FALSE, scale. = FALSE, rank. = ntopics)
    U.init <- logodds.pca$x  # cells by k
    V.init <- logodds.pca$rotation  # genes by k, no need to transpose
  } else {
    U.init <- scale(topics.mat, center = TRUE, scale = FALSE)  # c by k
    V.init <- scale(terms.mat, center = TRUE, scale = FALSE)  # k by g
    V.init <- t(V.init)  # now its g by k
  }
  return(list(Y.filt = Y.filt, U.init = U.init, V.init = V.init, X.mat = X.mat, size.factor = size.factor, ntopics = ntopics))
}

##### Deviance functions #####

poisson_deviance<-function(x,mu,sz){
  #assumes log link and size factor sz on the same scale as x (not logged)
  #stopifnot(all(x>=0 & sz>0))
  2*sum(x*log(x/(sz*mu)),na.rm=TRUE)-2*sum(x-sz*mu)
}

multinomial_deviance<-function(x,p){
  -2*sum(x*log(p))
}

binomial_deviance<-function(x,p,n){
  term1<-sum(x*log(x/(n*p)), na.rm=TRUE)
  nx<-n-x
  term2<-sum(nx*log(nx/(n*(1-p))), na.rm=TRUE)
  2*(term1+term2)
}


gof_func<-function(x,sz,mod=c("binomial","multinomial","poisson","geometric")){
  #Let n=colSums(original matrix where x is a row)
  #if binomial, assumes sz=n, required! So sz>0 for whole vector
  #if poisson, assumes sz=n/geometric_mean(n), so again all of sz>0
  #if geometric, assumes sz=log(n/geometric_mean(n)) which helps numerical stability. Here sz can be <>0
  #note sum(x)/sum(sz) is the (scalar) MLE for "mu" in Poisson and "p" in Binomial
  mod<-match.arg(mod)
  fit<-list(deviance=0,df.residual=length(x)-1,converged=TRUE)
  if(mod=="multinomial"){
    fit$deviance<-multinomial_deviance(x,sum(x)/sum(sz))
  } else if(mod=="binomial"){
    fit$deviance<-binomial_deviance(x,sum(x)/sum(sz),sz)
  } else if(mod=="poisson"){
    fit$deviance<-poisson_deviance(x,sum(x)/sum(sz),sz)
  } else if(mod=="geometric"){
    if(any(x>0)) {
      fit<-glm(x~offset(sz),family=MASS::negative.binomial(theta=1))
    }
  } else { stop("invalid model") }
  if(fit$converged){
    dev<-fit$deviance
    df<-fit$df.residual #length(x)-1
    pval<-pchisq(dev,df,lower.tail=FALSE)
    res<-c(dev,pval)
  } else {
    res<-rep(NA,2)
  }
  names(res)<-c("deviance","pval")
  res
}



compute_size_factors<-function(m,mod=c("binomial","multinomial","poisson","geometric")){
  #given matrix m with samples in the columns
  #compute size factors suitable for the discrete model in 'mod'
  mod<-match.arg(mod)
  sz<-Matrix::colSums(m) #base case, multinomial or binomial
  if(mod %in% c("multinomial","binomial")){ return(sz) }
  sz<-log(sz)
  sz<-sz - mean(sz) #make geometric mean of sz be 1 for poisson, geometric
  if(mod=="poisson"){ return(exp(sz)) }
  sz #geometric, use log scale size factors
}



##### Null Residuals functions #####

poisson_deviance_residuals<-function(x,xhat){
  #x,xhat assumed to be same dimension
  #sz<-exp(offsets)
  #xhat<-mu*sz
  term1<-x*log(x/xhat)
  term1[is.nan(term1)]<-0 #0*log(0)=0
  s2<-2*(term1-(x-xhat))
  sign(x-xhat)*sqrt(abs(s2))
}

binomial_deviance_residuals<-function(X,p,n){
  #X a matrix, n is vector of length ncol(X)
  #if p is matrix, must have same dims as X
  #if p is vector, its length must match nrow(X)
  if(length(p)==nrow(X)){
    mu<-outer(p,n)
  } else if(!is.null(dim(p)) && dim(p)==dim(X)){
    mu<-t(t(p)*n)
  } else { stop("dimensions of p and X must match!") }
  term1<-X*log(X/mu)
  # term1[is.nan(term1)]<-0 #0*log(0)=0
  term1[is.na(term1)]<-0 #0*log(0)=0
  nx<- t(n-t(X))
  term2<-nx*log(nx/outer(1-p,n))
  term2<-nx*log(nx/outer(1-p,n))
  # term2[is.nan(term2)]<-0
  sign(X-mu)*sqrt(2*(term1+term2))
}

# multinomial_deviance_residuals<-function(X,p,n){
#   #not clear if this actually makes sense. Don't use this function!
#   #X a matrix, n is vector of length ncol(X)
#   #if p is matrix, must have same dims as X
#   #if p is vector, its length must match nrow(X)
#   if(length(p)==nrow(X)){
#     mu<-outer(p,n)
#   } else if(!is.null(dim(p)) && dim(p)==dim(X)){
#     mu<-t(t(p)*n)
#   } else { stop("dimensions of p and X must match!") }
#   sign(X-mu)*sqrt(-2*X*log(p))
# }

null_residuals<-function(m,mod=c("binomial","multinomial","poisson","geometric"),type=c("deviance","pearson")){
  mod<-match.arg(mod)
  type<-match.arg(type)
  sz<-compute_size_factors(m,mod)
  if(mod %in% c("multinomial","binomial")) {
    phat<-Matrix::rowSums(m)/sum(sz)
    if(type=="pearson"){ #pearson resids same for multinomial as binomial
      mhat<-outer(phat,sz)
      return((m-mhat)/sqrt(mhat*(1-phat)))
    } else { #deviance residuals
      if(mod=="multinomial"){
        stop("multinomial deviance residuals not implemented yet")
        #return(multinomial_deviance_residuals(m,phat,sz))
      } else { #binomial
        return(binomial_deviance_residuals(m,phat,sz))
      }
    }
  } else if(mod=="poisson"){
    mhat<-outer(Matrix::rowSums(m)/sum(sz), sz) #first argument is "lambda hat" (MLE)
    if(type=="deviance"){ 
      return(poisson_deviance_residuals(m,mhat))
    } else { #pearson residuals
      return((m-mhat)/sqrt(mhat))
    }
  } else { #geometric
    gfunc<-function(g){
      fit<-glm(m[g,]~offset(sz),family=MASS::negative.binomial(theta=1))
      residuals(fit,type=type)
    }
    return(t(vapply(1:nrow(m),gfunc,rep(0.0,ncol(m)))))
  }
}
