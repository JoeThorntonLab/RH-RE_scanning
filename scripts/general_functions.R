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


# get colors for ancestral backgrounds
# bg: character vector of ancestral backgrounds
bg_color <- function(bg = c("AncSR1", "AncSR2")) {
  colors <- list(AncSR1 = "#5E0DEC", AncSR2 = "#A2CB5B")
  return(sapply(bg, function(x) colors[[x]]))
}


# get colors palette for REs
RE_color <- function(RE = REs[[1]]) {
  rotate <- function(x) t(apply(x, 2, rev))
  n <- 4
  mm <- tcrossprod(seq(1,0,length.out = n))
  
  tmp1 <- sapply(col2rgb("#9bec0d")/255, function(x) 1-mm*(1-x))
  tmp2 <- sapply(col2rgb("#ec0d2b")/255, function(x) 1-rotate(mm)*(1-x))
  tmp3 <- sapply(col2rgb("#0deccd")/255, function(x) 1-rotate(rotate(mm))*(1-x))
  tmp4 <- sapply(col2rgb("#5e0dec")/255, function(x) 1-rotate(rotate(rotate(mm)))*(1-x))
  
  tmp <- (tmp1*tmp2*tmp3*tmp4)
  RE_palette <- matrix(rgb(tmp), nrow = n)
  RE_palette_legend <- matrix(c("SRE1 (AA)", "CA", "AC", "CC", 
                                "SRE2 (GA)", "TA", "GC", "TC", 
                                "AT", "CT","AG", "CG", 
                                "ERE (GT)", "TT", "GG", "TG"), 
                              nrow = n, byrow = TRUE)
  RE_palette <- as.vector(RE_palette)
  names(RE_palette) <- as.vector(RE_palette_legend)
  return(RE_palette[as.character(RE)])
}


# determine whether two amino acid sequences can be connected by a single 
# nucleotide change given the genetic code
connected <- function(seqs) {
  seq1 <- unlist(strsplit(seqs[1], ""))
  seq2 <- unlist(strsplit(seqs[2], ""))
  
  if(sum(seq1 != seq2) == 1) {
    i <- which(seq1 != seq2)
    codons1 <- names(GENETIC_CODE[GENETIC_CODE==seq1[i]])
    codons2 <- names(GENETIC_CODE[GENETIC_CODE==seq2[i]])
    pairs <- expand.grid(list(codons1, codons2), stringsAsFactors=FALSE)
    connected <- apply(pairs, 1, function(x) 
      sum(do.call("!=", strsplit(as.character(x), ""))) == 1)
    return(sum(connected) > 0)
  } else return(FALSE)
}
