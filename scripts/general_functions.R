# general functions for working with the RH-RE DMS data

# make factors with different levels for ordering of REs in plotting
REs <- list()
REs[[1]] <- factor(c("SRE1 (AA)", "SRE2 (GA)", "ERE (GT)", "AC", 
                     "AG", "AT", "CA", "CC", 
                     "CG", "CT", "GC", "GG", 
                     "TA", "TC", "TG", "TT"), 
                   levels=c("ERE (GT)", "SRE1 (AA)", "SRE2 (GA)", "AC", 
                            "AG", "AT", "CA", "CC", 
                            "CG", "CT","GC", "GG", 
                            "TA", "TC", "TG", "TT"))
REs[[2]] <- factor(REs[[1]], 
                   levels=c("SRE1 (AA)", "SRE2 (GA)", "AT", "ERE (GT)",
                            "CA", "TA", "CT", "TT",
                            "AC", "GC", "AG", "GG",
                            "CC", "TC", "CG", "TG"))
REs[[3]] <- factor(REs[[1]], 
                   levels=c("SRE1 (AA)", "CA", "AC", "CC",
                            "SRE2 (GA)", "TA", "GC", "TC",
                            "AT", "CT", "AG", "CG",
                            "ERE (GT)", "TT", "GG", "TG"))
REs[[4]] <- factor(REs[[1]], 
                   levels=c("SRE1 (AA)", "CA", "AC", "CC",
                            "TC", "GC", "TA", "SRE2 (GA)",
                            "AT", "CT", "AG", "CG",
                            "TG", "GG", "TT", "ERE (GT)"))
REs[[5]] <- factor(REs[[1]], 
                   levels=c("SRE1 (AA)", "SRE2 (GA)", "TA", "CA",
                            "CC", "TC", "GC", "AC",
                            "ERE (GT)", "TT", "CT", "AT",
                            "TG", "GG", "AG", "CG"))

# convert REBC number to RE variant
REBC_to_RE <- function(REBC, levels=1) {
  if(is.character(REBC)) {
    REBC <- sub("AncSR._REBC", "", REBC)
    REBC <- sub("REBC", "", REBC)
    REBC <- as.integer(REBC)
  }
  RE <- REs[[levels]][REBC]
  return(RE)
}

# calculate coefficient of determination for a model fit to data
calc_R2 <- function(pred, y) {
  ssr <- sum((y - pred)^2)
  y.mu <- mean(y)
  sst <- sum((y - y.mu)^2)
  return(1 - ssr/sst)
}