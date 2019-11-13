
# Functions ---------------------------------------------------------------

GetBIC.lm <- function(jfit){
  L <- as.numeric(logLik(jfit))
  k <- attributes(logLik(jfit))$df
  N <- nobs(jfit)
  BIC <- -2*L + k * log(N)
  return(BIC)
}

SelectBestModel <- function(fit.row, cnames){
  # Find best model by grepping "BIC." and
  # returning the model that has the lowest bic
  col.i <- grepl("^bic.", cnames)
  fit.rowsub <- fit.row[col.i]
  cnames.sub <- cnames[col.i]
  # print(fit.rowsub)
  best.i <- which.min(fit.rowsub)
  models.sub <- sapply(cnames.sub, function(cname) strsplit(cname, split = "bic.")[[1]][[2]])
  model.best <- models.sub[best.i]
  return(model.best)
}

FitLM <- function(jsub, jmodel){
  jfit <- lm(formula = as.formula(jmodel), data = jsub)
  return(jfit)
}


GetFit <- function(fit.row, cnames, model){
  fit.cname <- paste0("fit.", model)
  indx <- which(cnames == fit.cname)
  assertthat::assert_that(length(indx) > 0)
  return(fit.row[[indx]])
}

GetBatchEffectFromFit <- function(jfit, varname = "platep02"){
  be <- coefficients(jfit)[[varname]]
}

GetBEFromBestModel <- function(fit.row, cnames, varname = "platep02"){
  indx <- which(cnames == "model")
  assertthat::assert_that(length(indx) > 0)
  best.model <- fit.row[[indx]]
  jfit <- GetFit(fit.row, cnames, best.model)
  # print(jfit)
  if (class(jfit) == "lm"){
    be <- GetBatchEffectFromFit(jfit, varname = varname)
  } else {
    be <- jfit$par[["be"]]
    assertthat::assert_that(!is.null(be))
  }
  return(be)
}



GetBIC.nlm <- function(optim.out, N){
  LogL <- -1 * optim.out$value  # optim minimizes -logL
  k <- length(optim.out$par)
  BIC <- -2 * LogL + k * log(N)
}

GetLRTStat <- function(optim.out){
  LogL <- -1 * optim.out$value
  return(-2 * LogL)
}

GetChiSqrDof <- function(optim.out2, optim.out1){
  # Dgrees of freedom for LRT test comparing two models 
  return(abs(length(optim.out2$par) - length(optim.out1$par)))
}

# analyze the parameters
GetGaussPars <- function(jsub){
  jfit <- jsub$fit.gauss
  out <- data.frame(int = jfit$par[[1]], plate = jfit$par[[2]], scale = jfit$par[[3]], mu = jfit$par[[4]], sig = jfit$par[[5]])
  # out <- as.data.frame(t(data.frame(jfit$par)))
  return(out)
}



# Fit Gaussian Peak -------------------------------------------------------

FitGaussPeak <- function(jsub, params.init, lims.lower, lims.upper){
  # int <- params.init[[1]]
  # plate <- params.init[[2]]
  # scale <- params.init[[3]]
  # gauss center <- params.init[[4]]
  # gauss width <- params.init[[5]]
  desmat <- model.matrix(~ 1 + plate, data = jsub)
  assertthat::assert_that(ncol(desmat) == 2)  # only 1 plate effect modeled, separate coding for 2 plate effects
  plate.x <- desmat[, 2]
  optim(par = params.init, fn = Gauss.logL, x = jsub$x, plate.x = plate.x, y = jsub$exprs, method = "L-BFGS-B", hessian=TRUE, lower = lims.lower, upper = lims.upper)
}

Gauss.logL <-  function(p, x, plate.x, y){
  # int <- params.init[[1]]
  # plate <- params.init[[2]]
  # scale <- params.init[[3]]
  # gauss center <- params.init[[4]]
  # gauss width <- params.init[[5]]
  mu <- GetGauss(p, x, plate.x)
  # mu <- p[1] + p[2] * plate.x + p[3]*dnorm(x,mean=p[4],sd=p[5])
  sig2.dat <- sum((y - mu) ^ 2) / length(y)
  loss <- -sum(dnorm(y, mean=mu, sd=sqrt(sig2.dat), log=T))
}

GetGauss <- function(p, x, plate.x){
  mu = p[1] + p[2] * plate.x + p[3]*dnorm(x,mean=p[4],sd=p[5])
  return(mu)
}


# Fit flat model ----------------------------------------------------------

FitFlatNLM <- function(jsub, params.init, lims.lower, lims.upper){
  # int <- params.init[[1]]
  # plate <- params.init[[2]]
  desmat <- model.matrix(~ 1 + plate, data = jsub)
  assertthat::assert_that(ncol(desmat) == 2)  # only 1 plate effect modeled, separate coding for 2 plate effects
  plate.x <- desmat[, 2]
  optim(par = params.init, fn = Flat.logL, x = jsub$x, plate.x = plate.x, y = jsub$exprs, method = "L-BFGS-B", hessian=TRUE, lower = lims.lower, upper = lims.upper)
}

Flat.logL <- function(p, x, plate.x, y){
  # int <- params.init[[1]]
  # plate <- params.init[[2]]
  mu <- GetFlat(p, x, plate.x)
  sig2.dat <- sum((y - mu) ^ 2) / length(y)
  loss <- -sum(dnorm(y, mean=mu, sd=sqrt(sig2.dat), log=T))
}

GetFlat <- function(p, x, plate.x){
  mu <- p[1] + p[2] * plate.x
  return(mu)
}


# Fit line model  -------------------------------------------------------

FitLineNLM <- function(jsub, params.init, lims.lower, lims.upper){
  # int <- params.init[[1]]
  # plate <- params.init[[2]]
  # x effect <- params.init[[3]]
  desmat <- model.matrix(~ 1 + plate, data = jsub)
  assertthat::assert_that(ncol(desmat) == 2)  # only 1 plate effect modeled, separate coding for 2 plate effects
  plate.x <- desmat[, 2]
  optim(par = params.init, fn = Line.logL, x = jsub$x, plate.x = plate.x, y = jsub$exprs, method = "L-BFGS-B", hessian=TRUE, lower = lims.lower, upper = lims.upper)
}

Line.logL <- function(p, x, plate.x, y){
  # int <- params.init[[1]]
  # plate <- params.init[[2]]
  # x effect <- params.init[[3]]
  mu <- GetLine(p, x, plate.x)
  sig2.dat <- sum((y - mu) ^ 2) / length(y)
  loss <- -sum(dnorm(y, mean=mu, sd=sqrt(sig2.dat), log=T))
}

GetLine <- function(p, x, plate.x){
  mu <- p[1] + p[2] * plate.x + p[[3]] * x
  return(mu)
}

# Fit expo model  -------------------------------------------------------

FitExpoNLM <- function(jsub, params.init, lims.lower, lims.upper){
  # int <- params.init[[1]]
  # plate <- params.init[[2]]
  # x effect <- params.init[[3]]
  desmat <- model.matrix(~ 1 + plate, data = jsub)
  assertthat::assert_that(ncol(desmat) == 2)  # only 1 plate effect modeled, separate coding for 2 plate effects
  plate.x <- desmat[, 2]
  optim(par = params.init, fn = Expo.logL, x = jsub$x, plate.x = plate.x, y = jsub$exprs, method = "L-BFGS-B", hessian=TRUE, lower = lims.lower, upper = lims.upper)
}

Expo.logL <- function(p, x, plate.x, y){
  # int <- params.init[[1]]
  # plate <- params.init[[2]]
  # x effect <- params.init[[3]]
  mu <- GetExpo(p, x, plate.x)
  sig2.dat <- sum((y - mu) ^ 2) / length(y)
  loss <- -sum(dnorm(y, mean=mu, sd=sqrt(sig2.dat), log=T))
}

GetExpo <- function(p, x, plate.x){
  mu <- p[1] + p[2] * plate.x + p[[3]] * exp(x)
  return(mu)
}


# 
# # Fit flat model  ---------------------------------------------------------
# 
# FitGaussPeak <- function(jsub, params.init, lims.lower, lims.upper){
#   # int <- params.init[[1]]
#   # plate <- params.init[[2]]
#   # scale <- params.init[[3]]
#   # gauss center <- params.init[[4]]
#   # gauss width <- params.init[[5]]
#   x <- jsub$x
#   plate.x <- model.matrix(~ 1 + plate, data = jsub)[, 2]
#   y <- jsub$exprs
#   optim(par = params.init, fn = Gauss.logL, x = x, plate.x = plate.x, method = "L-BFGS-B", hessian=TRUE, lower = lims.lower, upper = lims.upper)
# }
# 
# Gauss.logL <-  function(p, x, plate.x){
#   mu <- GetGauss(p, x, plate.x)
#   mu = p[1] + p[2] * plate.x + p[3]*dnorm(x,mean=p[4],sd=p[5])
#   sig2.dat <- sum((y - mu) ^ 2) / length(y)
#   loss <- -sum(dnorm(y, mean=mu, sd=sqrt(sig2.dat), log=T))
# }
# 
# GetGauss <- function(p, x, plate.x){
#   mu = p[1] + p[2] * plate.x + p[3]*dnorm(x,mean=p[4],sd=p[5])
#   return(mu)
# }
# 
