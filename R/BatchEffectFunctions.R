AssignEffect <- function(batch, cluster, effects, jtype = c("plate", "cluster"), plateprefix = "jrep2", clusterprefix = "cluster"){
  if (jtype == "plate"){
    jname <- paste0(plateprefix, batch)
    if (jname %in% names(effects)){
      # print("Adjusting plate nonzero")
      indx <- names(effects) %in% jname
      adj <- effects[indx]
      # adj <-  effects[[effects[jname]]]
    }  else {
      # print("Adjusting plate  zero")
      adj <- 0
    }
  } else if (jtype == "cluster"){
    jname <- paste0(plateprefix, batch, ":", clusterprefix,  cluster)
    if (jname %in% names(effects)){
      # print(jname)
      # print(names(effects))
      indx <- names(effects) %in% jname
      adj <- effects[indx]
      # print("Adjusting cluster nonzero")
    } else {
      # print("Adjusting cluster zero")
      adj <- 0
    }
  } else {
    print(paste("jtype must be plate or cluster", jtype))
  }
  return(adj)
}

AdjustBatchEffect <- function(jdat, plateprefix = "jrep2", clusterprefix = "cluster", platesuffix = "old"){
  jout <- lm(data = jdat, formula = log2exprs ~ 1 + jrep2 + cluster + jrep2:cluster)
  plateeffect.i <- grep(paste0(plateprefix, ".*.", platesuffix, "$"), names(jout$coefficients), value = FALSE)
  plateeffect <- jout$coefficients[plateeffect.i]
  # interactioneffect.i <- grep(plateprefix, ".*.:", names(jout$coefficients), value = FALSE)
  interactioneffect.i <- grep(paste0(plateprefix, ".*.:", clusterprefix), names(jout$coefficients), value = FALSE)
  interactioneffect <- jout$coefficients[interactioneffect.i]
  
  jdat <- jdat %>%
    rowwise() %>%
    mutate(plateadj2 = AssignEffect(jrep2, cluster, plateeffect, jtype = "plate"),
           clstradj2 = AssignEffect(jrep2, cluster, interactioneffect, jtype = "cluster"),
           log2exprsadj = log2exprs - plateadj2 - clstradj2)
  return(jdat)
}
