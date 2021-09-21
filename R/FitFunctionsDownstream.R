# Jake Yeung
# Date of Creation: 2020-12-27
# File: ~/projects/scchic-functions/R/FitFunctionsDownstream.R
# 


SummarizeParamsPvalues <- function(jfits.lst, jmark, paramname = "Cluster"){
  # jmark <- "H3K4me1"
  jnames <- names(jfits.lst); names(jnames) <- jnames
  # https://stats.stackexchange.com/questions/315311/how-to-find-p-value-using-estimate-and-standard-error
  params.dat.all.withpval <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xkeep <- grepl(paste0("^", paramname, ".*.Estimate$"), x = names(x))
    jparams <- x[xkeep]
    xkeep.se <- grepl(paste0("^", paramname, ".*.StdError$"), x = names(x))
    jparams.se <- x[xkeep.se]
    jout <- data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), se = unlist(jparams.se), mark = jmark, stringsAsFactors = FALSE) %>%
      rowwise() %>%
      mutate(z = estimate / se, 
             pval.param = exp(-0.717 * z - 0.416 * z ^ 2))
  }) %>%
    bind_rows() 
  return(params.dat.all.withpval)
}


SummarizeMeanValue <- function(jfits.lst, jmark){
  # jmark <- "H3K4me1"
  jnames <- names(jfits.lst); names(jnames) <- jnames
  # https://stats.stackexchange.com/questions/315311/how-to-find-p-value-using-estimate-and-standard-error
  params.dat.all.withpval <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xkeep <- grepl(paste0("X.Intercept..Estimate"), x = names(x))
    jparams <- x[xkeep]
    jout <- data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE) %>%
      rowwise()
  }) %>%
    bind_rows() 
  return(params.dat.all.withpval)
}

