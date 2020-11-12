# Jake Yeung
# Date of Creation: 2020-06-06
# File: ~/projects/scchic-functions/R/FitFunctions.R
# 

FitGlmRowClustersPlate <- function(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = NULL, returnobj=FALSE){
  # use Offset by size of library
  # https://stats.stackexchange.com/questions/66791/where-does-the-offset-go-in-poisson-negative-binomial-regression
  # fit GLM for a row of a sparse matrix, should save some space?
  
  # pvalue by deviance goodness of fit: https://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
  # offset is in log because the model says the log counts is equal to RHS
  
  if (!is.null(nrow(jrow))){
    # probably a matrix of many rows, sum them up
    print(paste("Merging", nrow(jrow), "rows"))
    row <- Matrix::colSums(jrow)
  }
  dat <- data.frame(cell = cnames, ncuts = jrow, stringsAsFactors = FALSE) %>%
    left_join(., dat.annots.filt.mark, by = "cell") %>%
    left_join(., ncuts.cells.mark, by = "cell")
  
  # m1.pois <- glm(ncuts ~ 1 + Cluster + offset(ncuts.total), data = dat, family = "poisson")
  m1.pois <- glm(ncuts ~ 1 + Plate + Cluster + offset(log(ncuts.total)), data = dat, family = "poisson")
  mnull.pois <- glm(ncuts ~ 1 + Plate + offset(log(ncuts.total)), data = dat, family = "poisson")
  
  if (!returnobj){
    jsum <- anova(mnull.pois, m1.pois)
    pval <- pchisq(jsum$Deviance[[2]], df = jsum$Df[[2]], lower.tail = FALSE)
    out.dat <- data.frame(pval = pval, 
                          dev.diff = jsum$Deviance[[2]],
                          df.diff = jsum$Df[[2]],
                          t(as.data.frame(coefficients(m1.pois))), 
                          stringsAsFactors = FALSE)
    if (!is.null(jbin)){
      out.dat$bin <- jbin
      rownames(out.dat) <- jbin
    }
    return(out.dat)
  } else {
    return(list(fit.full = m1.pois, fit.null = mnull.pois, dat.input = dat))
  }
}




RefitPoissonForPlot <- function(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, return.means = TRUE){
  jfit <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, returnobj = TRUE)
  # https://stackoverflow.com/questions/27669101/strange-error-in-glm
  # confint gives error sometimes
  ci <- confint.default(jfit$fit.full)
  
  params.fc.dat <- data.frame(param = rownames(ci), logLambdaLower = ci[, 1], logLambdaUpper = ci[, 2]) %>%
    filter(param != "(Intercept)") %>%
    rowwise() %>%
    mutate(Cluster = ifelse(param == "(Intercept)", "aHSPCs", gsub("Cluster", "", param)),
           logLambda = mean((logLambdaLower + logLambdaUpper) / 2))
  
  params.int.dat <- data.frame(param = rownames(ci)[[1]], logLambdaLower = ci[1, 1], logLambdaUpper = ci[1, 2]) %>%
    rowwise() %>%
    mutate(Cluster = ifelse(param == "(Intercept)", "aHSPCs", gsub("Cluster", "", param)),
           logLambda = mean((logLambdaLower + logLambdaUpper) / 2))
  
  if (return.means){
    # use for plotting against raw data 
    params.mean.dat <- params.fc.dat %>%
      mutate(logLambda = params.int.dat$logLambda + logLambda,
             logLambdaLower = params.int.dat$logLambdaLower + logLambdaLower,
             logLambdaUpper = params.int.dat$logLambdaUpper + logLambdaUpper) %>%
      bind_rows(., params.int.dat)
    
    # plot raw 
    input.dat <- jfit$dat.input %>% rowwise() %>% mutate(logLambda = log(ncuts) - log(ncuts.total))
    return(list(input.dat = input.dat, params.mean.dat = params.mean.dat, params.int.dat = params.int.dat))
  } else {
    # useful for just plotting CI of FC parameters
    return(list(params.fc.dat = params.fc.dat, params.int.dat = params.int.dat))
  }
}


# RefitPoissonForPlot.old <- function(jrow, cnames, datannots.filt.mark, ncuts.cells.mark){
#   jfit <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, returnobj = TRUE)
#   ci <- confint(jfit$fit.full)
#   
#   params.fc.dat <- data.frame(param = rownames(ci), logLambdaLower = ci[, 1], logLambdaUpper = ci[, 2]) %>%
#     filter(param != "(Intercept)") %>%
#     rowwise() %>%
#     mutate(Cluster = ifelse(param == "(Intercept)", "aHSPCs", gsub("Cluster", "", param)),
#            logLambda = mean((logLambdaLower + logLambdaUpper) / 2))
#   
#   params.int.dat <- data.frame(param = rownames(ci)[[1]], logLambdaLower = ci[1, 1], logLambdaUpper = ci[1, 2]) %>%
#     rowwise() %>%
#     mutate(Cluster = ifelse(param == "(Intercept)", "aHSPCs", gsub("Cluster", "", param)),
#            logLambda = mean((logLambdaLower + logLambdaUpper) / 2))
#   
#   params.mean.dat <- params.fc.dat %>%
#     mutate(logLambda = params.int.dat$logLambda + logLambda,
#            logLambdaLower = params.int.dat$logLambdaLower + logLambdaLower,
#            logLambdaUpper = params.int.dat$logLambdaUpper + logLambdaUpper) %>%
#     bind_rows(., params.int.dat)
#   
#   # plot raw 
#   input.dat <- jfit$dat.input %>% rowwise() %>% mutate(logLambda = log(ncuts) - log(ncuts.total))
#   return(list(input.dat = input.dat, params.mean.dat = params.mean.dat, params.int.dat = params.int.dat))
# }

FitGlmRowClusters <- function(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = NULL, returnobj=FALSE){
  # use Offset by size of library
  # https://stats.stackexchange.com/questions/66791/where-does-the-offset-go-in-poisson-negative-binomial-regression
  # fit GLM for a row of a sparse matrix, should save some space?
  
  # pvalue by deviance goodness of fit: https://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
  # offset is in log because the model says the log counts is equal to RHS
  
  if (!is.null(nrow(jrow))){
    # probably a matrix of many rows, sum them up
    print(paste("Merging", nrow(jrow), "rows"))
    row <- Matrix::colSums(jrow)
  }
  dat <- data.frame(cell = cnames, ncuts = jrow, stringsAsFactors = FALSE) %>%
    left_join(., dat.annots.filt.mark, by = "cell") %>%
    left_join(., ncuts.cells.mark, by = "cell")
  
  # m1.pois <- glm(ncuts ~ 1 + Cluster + offset(ncuts.total), data = dat, family = "poisson")
  m1.pois <- glm(ncuts ~ 1 + Cluster + offset(log(ncuts.total)), data = dat, family = "poisson")
  mnull.pois <- glm(ncuts ~ 1 + offset(log(ncuts.total)), data = dat, family = "poisson")
  
  if (!returnobj){
    jsum <- anova(mnull.pois, m1.pois)
    pval <- pchisq(jsum$Deviance[[2]], df = jsum$Df[[2]], lower.tail = FALSE)
    out.dat <- data.frame(pval = pval, 
                          dev.diff = jsum$Deviance[[2]],
                          df.diff = jsum$Df[[2]],
                          t(as.data.frame(coefficients(m1.pois))), 
                          stringsAsFactors = FALSE)
    if (!is.null(jbin)){
      out.dat$bin <- jbin
      rownames(out.dat) <- jbin
    }
    return(out.dat)
  } else {
    return(list(fit.full = m1.pois, fit.null = mnull.pois, dat.input = dat))
  }
}

