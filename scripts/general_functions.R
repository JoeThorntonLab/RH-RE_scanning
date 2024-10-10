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
REs[[6]] <- factor(c("SRE (AA)", "GA", "ERE (GT)", "AC", 
                     "AG", "AT", "CA", "CC", 
                     "CG", "CT", "GC", "GG", 
                     "TA", "TC", "TG", "TT"), 
                   levels=c("ERE (GT)", "SRE (AA)", "AC", "AG", 
                            "AT", "CA", "CC", "CG", 
                            "CT", "GA", "GC", "GG", 
                            "TA", "TC", "TG", "TT"))

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

# breaks for log10plus1 transformation
log_breaksplus1 <- function (n = 5, base = 10) {
  # force_all(n, base)
  n_default <- n
  function(x, n = n_default) {
    raw_rng <- suppressWarnings(range(x + 1, na.rm = TRUE))
    if (any(!is.finite(raw_rng))) {
      return(numeric())
    }
    rng <- log(raw_rng, base = base)
    min <- floor(rng[1])
    max <- ceiling(rng[2])
    if (max == min) {
      return(base^min)
    }
    by <- floor((max - min)/n) + 1
    breaks <- base^seq(min, max, by = by)
    relevant_breaks <- base^rng[1] <= breaks & breaks <= 
      base^rng[2]
    if (sum(relevant_breaks) >= (n - 2)) {
      return(breaks)
    }
    while (by > 1) {
      by <- by - 1
      breaks <- base^seq(min, max, by = by)
      relevant_breaks <- base^rng[1] <= breaks & breaks <= 
        base^rng[2]
      if (sum(relevant_breaks) >= (n - 2)) {
        return(breaks)
      }
    }
    log_sub_breaks(rng, n = n, base = base)
  }
}

# log10 of count + 1 transformation
log10plus1 <- scales::trans_new(name = "log10plus1",
                                transform = function(x) log10(x+1),
                                inverse = function(x) 10^x-1,
                                breaks = log_breaksplus1(base = 10),
                                domain = c(1e-100, Inf))


# get colors for ancestral backgrounds
# bg: character vector of ancestral backgrounds
bg_color <- function(bg = c("AncSR1", "AncSR2")) {
  colors <- list(AncSR1 = "#5E0DEC", AncSR2 = "#9bec0d")
  return(sapply(bg, function(x) colors[[x]]))
}

levels=c("ERE (GT)", "SRE (AA)", "AC", "AG", 
         "AT", "CA", "CC", "CG", 
         "CT", "GA", "GC", "GG", 
         "TA", "TC", "TG", "TT")

# colorblind friendly palette for REs
# from https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
RE_color <- c("#490092", "#24ff24", "#ffb6db", "#b6dbff",
                       "#009292", "#006ddb", "#000000", "#ff6db6",
                       "#004949", "#ffff6d", "#b66dff", "#924900",
                       "#920000", "#000000", "#db6d00", "#6db6ff")
names(RE_color) <- levels(REs[[6]])
