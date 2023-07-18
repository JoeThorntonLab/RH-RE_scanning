library(Biostrings)
library(MASS)
library(Matrix)
library(dplyr)

# global variables
# GENETIC_CODE_NO_STOP <- GENETIC_CODE[GENETIC_CODE != "*"]
# AA <- sort(unique(GENETIC_CODE_NO_STOP))
AA <- sort(unique(GENETIC_CODE))

# returns TRUE if seq1 and seq2 have Hamming distance 1
# seq1 can be a vector of multiple character strings; seq2 can only have 1
# if return.diff = TRUE, returns the differences for seqs with Hamming
# distance 1
adjacent <- function(seq1, seq2, return.diff = FALSE) {
  seq1 <- strsplit(seq1, "")
  seq2 <- strsplit(seq2, "")[[1]]
  
  distance <- sapply(seq1, function(x) sum(x != seq2))
  hd1 <- distance == 1
  # distance <- sum(seq1 != seq2)
  if(return.diff) {
    diff <- sapply(seq1[hd1], function(x) {
      i <- which(x != seq2)
      return(c(x[i], seq2[i]))
    })
    return(list(adj = hd1, diff = diff))
  } else return(hd1)
}

# return all possible nucleotide sequences that code for a given amino acid 
# sequence
aa_to_codon <- function(aaseq) {
  aaseq <- strsplit(aaseq, "")[[1]]
  codons <- lapply(aaseq, function(x) names(GENETIC_CODE)[GENETIC_CODE == x])
  ntseqs <- expand.grid(codons)
  ntseqs <- apply(ntseqs, 1, function(x) do.call(paste0, as.list(x)))
  return(ntseqs)
}

# Calculate number of possible mutational paths in nucleotide space between two
# amino acid variants that are either identical or one step apart in aa space.
nt_paths <- function(aaseq1, aaseq2, aa_wadj_mat) {
  aaseq1 <- strsplit(aaseq1, "")[[1]]
  aaseq2 <- strsplit(aaseq2, "")[[1]]
  idiff <- which(aaseq1 != aaseq2)  # indices at which aa differences occur
  hd <- length(idiff)  # amino acid Hamming distance
  if(hd > 1) {
    stop("Sequences must have Hamming distance <= 1")
  }
  else if(hd == 1) {
    # get number of nucleotide paths for the single amino acid difference
    aadiff <- c(aaseq1[idiff], aaseq2[idiff])
    npathsdiff <- aa_wadj_mat[match(aadiff[1], rownames(aa_wadj_mat)),
                              match(aadiff[2], colnames(aa_wadj_mat))]
    # multiply by the number of backgrounds the difference occurs in to get
    # total number of possible nucleotide paths
    bgs <- lapply(aaseq1[-idiff], 
                  function(x) names(GENETIC_CODE)[GENETIC_CODE == x])
    nbg <- do.call(prod, lapply(bgs, length))
    npaths <- npathsdiff * nbg
    return(npaths)
  }
  else if(hd == 0) {
    # get number of synonymous mutations for each amino acid
    npathsaa <- diag(aa_wadj_mat[match(aaseq1, rownames(aa_wadj_mat)),
                                 match(aaseq1, rownames(aa_wadj_mat))])
    # multiply each by the number of backgrounds it can occur in to get
    # total number of possible nucleotide paths
    bgs <- lapply(aaseq1, function(x) names(GENETIC_CODE)[GENETIC_CODE == x])
    nbgaa <- lapply(bgs, length)
    npaths <- npathsaa * sapply(1:length(aaseq1), 
                                function(x) do.call(prod, nbgaa[-x]))
    npaths <- sum(npaths)
    return(npaths)
  }
}

# calculate total number of possible paths in nucleotide space starting from a
# given amino acid variant
nt_paths_total <- function(aaseq) {
  aaseq <- strsplit(aaseq, "")[[1]]
  # get number of synonymous nt variants for aaseq
  codons <- lapply(aaseq, function(x) names(GENETIC_CODE)[GENETIC_CODE == x])
  nntvars <- do.call(prod, lapply(codons, length))
  # each variant can access nt length * 3 other nt variants
  npaths <- nntvars * length(aaseq) * 3 * 3
  return(npaths)
}

# matrix exponentiation, capable of handling sparse matrices
mat_exp <- function(mat, k) {
  out <- mat
  for(i in 2:k) {
    out <- out %*% mat
  }
  return(out)
}


### build codon adjacency matrix, including STOP codons (64*64)
codon_adj_mat <- matrix(0, length(GENETIC_CODE), length(GENETIC_CODE))
colnames(codon_adj_mat) <- names(GENETIC_CODE)
rownames(codon_adj_mat) <- names(GENETIC_CODE)

for(i in seq(nrow(codon_adj_mat))) {
  adj <- sapply(colnames(codon_adj_mat), adjacent, 
                rownames(codon_adj_mat)[i])
  codon_adj_mat[i, adj] <- 1
}

# convert to sparse matrix
codon_adj_mat <- as(codon_adj_mat, "sparseMatrix")


### build weighted amino acid adjacency matrix
# weights represent the number of synonymous mutations that can convert amino 
# acid i into amino acid j, including STOP (21*21)
aa_wadj_mat <- matrix(0, length(AA), length(AA))
colnames(aa_wadj_mat) <- AA
rownames(aa_wadj_mat) <- AA

for(i in AA) {
  # get codons for AA i
  codons1 <- names(GENETIC_CODE)[GENETIC_CODE == i]
  for(j in AA) {
    # get codons for AA j
    codons2 <- names(GENETIC_CODE)[GENETIC_CODE == j]
    
    codon_submat <- codon_adj_mat[rownames(codon_adj_mat) %in% codons1, 
                                  colnames(codon_adj_mat) %in% codons2]
    # sum number of edges between codons for AA i and AA j
    aa_wadj_mat[rownames(aa_wadj_mat) == i, colnames(aa_wadj_mat) == j] <- 
      sum(codon_submat)
  }
}

# convert to sparse matrix
aa_wadj_mat <- as(aa_wadj_mat, "sparseMatrix")


### build weighted amino acid library adjacency matrices where entries represent 
# the number of nucleotide mutations that can convert between library variants

# read in amino acid variants to be included in the network
if(!file.exists(file.path("..", "results", "mutation_effects_model", 
                          "meanF_data.rda"))) {
  meanF_data <- read.csv(file.path("..", "results", "cleaned_data", 
                                   "meanF_data_corrected_NovaSeq.csv.gz"),
                         row.names = 1, stringsAsFactors = TRUE)
  save(meanF_data, file = file.path("..", "results", "mutation_effects_model", 
                                    "meanF_data.rda"))
} else load(file.path("..", "results", "mutation_effects_model", 
                      "meanF_data.rda"))
aa_lib_AncSR1 <- meanF_data %>% 
  filter(type == "exp", sig == "significant", bg == "AncSR1") %>% 
  pull(AA_var) %>% as.character() %>% unique() %>% sort()
aa_lib_AncSR2 <- meanF_data %>% 
  filter(type == "exp", sig == "significant", bg == "AncSR2") %>% 
  pull(AA_var) %>% as.character() %>% unique() %>% sort()

## AncSR1
# first compute number of nt mutations for variants with aa Hamming distance 1
i_AncSR1 <- numeric()
j_AncSR1 <- numeric()
x_AncSR1 <- numeric()
for(k in 1:length(aa_lib_AncSR1)) {
  aa1 <- aa_lib_AncSR1[k]
  adj <- adjacent(aa_lib_AncSR1, aa1)
  iadj <- which(adj)
  nadj <- sum(adj)
  aa2 <- aa_lib_AncSR1[iadj]
  
  i_AncSR1 <- append(i_AncSR1, iadj)
  j_AncSR1 <- append(j_AncSR1, rep(k, nadj))
  x_AncSR1 <- append(x_AncSR1, 
                     as.numeric(sapply(aa2, nt_paths, aa1, aa_wadj_mat)))
}

## construct matrices
# unweighted matrix (no synonymous mutations)
aa_lib_AncSR1_adj_mat <- sparseMatrix(i=i_AncSR1, j=j_AncSR1, 
                                      x=replace(x_AncSR1, x_AncSR1 > 0, 1), 
                                      dims = rep(length(aa_lib_AncSR1), 2))
colnames(aa_lib_AncSR1_adj_mat) <- aa_lib_AncSR1
rownames(aa_lib_AncSR1_adj_mat) <- aa_lib_AncSR1
# codon bias-weighted matrix (no synonymous mutations)
aa_lib_AncSR1_wadj_mat <- sparseMatrix(i=i_AncSR1, j=j_AncSR1, x=x_AncSR1, 
                                       dims = rep(length(aa_lib_AncSR1), 2))
colnames(aa_lib_AncSR1_wadj_mat) <- aa_lib_AncSR1
rownames(aa_lib_AncSR1_wadj_mat) <- aa_lib_AncSR1
# codon bias-weighted matrix (with synonymous mutations)
aa_lib_AncSR1_swadj_mat <- aa_lib_AncSR1_wadj_mat
diag(aa_lib_AncSR1_swadj_mat) <- sapply(1:length(aa_lib_AncSR1), 
                                        function(x) nt_paths(aa_lib_AncSR1[x], 
                                                             aa_lib_AncSR1[x], 
                                                             aa_wadj_mat))

## add "nonfunctional" category
# calculate number of accessible nonfunctional nucleotide variants for each aa 
# variant
naccessible_AncSR1 <- sapply(aa_lib_AncSR1, nt_paths_total)
naccessible_fxnal_AncSR1 <- rowSums(aa_lib_AncSR1_swadj_mat)
naccessible_nfxnal_AncSR1 <- naccessible_AncSR1 - naccessible_fxnal_AncSR1

# codon bias-weighted matrix with synonymous and nonfunctional mutations:
# these are immediately purged, so the path remains on the same node (same as
# synonymous substitutions)
aa_lib_AncSR1_nswadj_mat <- 
  aa_lib_AncSR1_swadj_mat + diag(naccessible_nfxnal_AncSR1)

## AncSR2
# first compute number of nt mutations for variants with aa Hamming distance 1
i_AncSR2 <- numeric()
j_AncSR2 <- numeric()
x_AncSR2 <- numeric()
for(k in 1:length(aa_lib_AncSR2)) {
  aa1 <- aa_lib_AncSR2[k]
  adj <- adjacent(aa_lib_AncSR2, aa1)
  iadj <- which(adj)
  nadj <- sum(adj)
  aa2 <- aa_lib_AncSR2[iadj]
  
  i_AncSR2 <- append(i_AncSR2, iadj)
  j_AncSR2 <- append(j_AncSR2, rep(k, nadj))
  x_AncSR2 <- append(x_AncSR2, 
                     as.numeric(sapply(aa2, nt_paths, aa1, aa_wadj_mat)))
}

## construct matrices
# unweighted matrix (no synonymous mutations)
aa_lib_AncSR2_adj_mat <- sparseMatrix(i=i_AncSR2, j=j_AncSR2, 
                                      x=replace(x_AncSR2, x_AncSR2 > 0, 1), 
                                      dims = rep(length(aa_lib_AncSR2), 2))
colnames(aa_lib_AncSR2_adj_mat) <- aa_lib_AncSR2
rownames(aa_lib_AncSR2_adj_mat) <- aa_lib_AncSR2
# codon bias-weighted matrix (no synonymous mutations)
aa_lib_AncSR2_wadj_mat <- sparseMatrix(i=i_AncSR2, j=j_AncSR2, x=x_AncSR2, 
                                       dims = rep(length(aa_lib_AncSR2), 2))
colnames(aa_lib_AncSR2_wadj_mat) <- aa_lib_AncSR2
rownames(aa_lib_AncSR2_wadj_mat) <- aa_lib_AncSR2
# codon bias-weighted matrix (with synonymous mutations)
aa_lib_AncSR2_swadj_mat <- aa_lib_AncSR2_wadj_mat
diag(aa_lib_AncSR2_swadj_mat) <- sapply(1:length(aa_lib_AncSR2), 
                                        function(x) nt_paths(aa_lib_AncSR2[x], 
                                                             aa_lib_AncSR2[x], 
                                                             aa_wadj_mat))

## add "nonfunctional" category
# calculate number of accessible nonfunctional nucleotide variants for each aa 
# variant
naccessible_AncSR2 <- sapply(aa_lib_AncSR2, nt_paths_total)
naccessible_fxnal_AncSR2 <- rowSums(aa_lib_AncSR2_swadj_mat)
naccessible_nfxnal_AncSR2 <- naccessible_AncSR2 - naccessible_fxnal_AncSR2

# codon bias-weighted matrix with synonymous and nonfunctional mutations:
# these are immediately purged, so the path remains on the same node (same as
# synonymous substitutions)
aa_lib_AncSR2_nswadj_mat <- 
  aa_lib_AncSR2_swadj_mat + diag(naccessible_nfxnal_AncSR2)


# convert adjacency matrices to transition probability matrices
aa_lib_AncSR1_trans_mat <- t(apply(aa_lib_AncSR1_adj_mat, 1, function(x) x / sum(x)))
aa_lib_AncSR1_trans_mat <- as(aa_lib_AncSR1_trans_mat, "sparseMatrix")
aa_lib_AncSR1_wtrans_mat <- t(apply(aa_lib_AncSR1_wadj_mat, 1, function(x) x / sum(x)))
aa_lib_AncSR1_wtrans_mat <- as(aa_lib_AncSR1_wtrans_mat, "sparseMatrix")
aa_lib_AncSR1_swtrans_mat <- t(apply(aa_lib_AncSR1_swadj_mat, 1, function(x) x / sum(x)))
aa_lib_AncSR1_swtrans_mat <- as(aa_lib_AncSR1_swtrans_mat, "sparseMatrix")
aa_lib_AncSR1_nswtrans_mat <- t(apply(aa_lib_AncSR1_nswadj_mat, 1, function(x) x / sum(x)))
aa_lib_AncSR1_nswtrans_mat <- as(aa_lib_AncSR1_nswtrans_mat, "sparseMatrix")
aa_lib_AncSR2_trans_mat <- t(apply(aa_lib_AncSR2_adj_mat, 1, function(x) x / sum(x)))
aa_lib_AncSR2_trans_mat <- as(aa_lib_AncSR2_trans_mat, "sparseMatrix")
aa_lib_AncSR2_wtrans_mat <- t(apply(aa_lib_AncSR2_wadj_mat, 1, function(x) x / sum(x)))
aa_lib_AncSR2_wtrans_mat <- as(aa_lib_AncSR2_wtrans_mat, "sparseMatrix")
aa_lib_AncSR2_swtrans_mat <- t(apply(aa_lib_AncSR2_swadj_mat, 1, function(x) x / sum(x)))
aa_lib_AncSR2_swtrans_mat <- as(aa_lib_AncSR2_swtrans_mat, "sparseMatrix")
aa_lib_AncSR2_nswtrans_mat <- t(apply(aa_lib_AncSR2_nswadj_mat, 1, function(x) x / sum(x)))
aa_lib_AncSR2_nswtrans_mat <- as(aa_lib_AncSR2_nswtrans_mat, "sparseMatrix")

# remove rows/columns for unconnected nodes
unconnected_AncSR1 <- rowSums(aa_lib_AncSR1_adj_mat) == 0
aa_lib_AncSR1_trans_mat <- aa_lib_AncSR1_trans_mat[!unconnected_AncSR1, !unconnected_AncSR1]
aa_lib_AncSR1_wtrans_mat <- aa_lib_AncSR1_wtrans_mat[!unconnected_AncSR1, !unconnected_AncSR1]
aa_lib_AncSR1_swtrans_mat <- aa_lib_AncSR1_swtrans_mat[!unconnected_AncSR1, !unconnected_AncSR1]
aa_lib_AncSR1_nswtrans_mat <- aa_lib_AncSR1_nswtrans_mat[!unconnected_AncSR1, !unconnected_AncSR1]

unconnected_AncSR2 <- rowSums(aa_lib_AncSR2_adj_mat) == 0
aa_lib_AncSR2_trans_mat <- aa_lib_AncSR2_trans_mat[!unconnected_AncSR2, !unconnected_AncSR2]
aa_lib_AncSR2_wtrans_mat <- aa_lib_AncSR2_wtrans_mat[!unconnected_AncSR2, !unconnected_AncSR2]
aa_lib_AncSR2_swtrans_mat <- aa_lib_AncSR2_swtrans_mat[!unconnected_AncSR2, !unconnected_AncSR2]
aa_lib_AncSR2_nswtrans_mat <- aa_lib_AncSR2_nswtrans_mat[!unconnected_AncSR2, !unconnected_AncSR2]

# export transition matrices
save(aa_lib_AncSR1_trans_mat, file = file.path("Matrix2_AncSR1.rda"))
save(aa_lib_AncSR1_wtrans_mat, file = file.path("Matrix3_AncSR1.rda"))
save(aa_lib_AncSR1_swtrans_mat, file = file.path("Matrix4_AncSR1.rda"))
save(aa_lib_AncSR1_nswtrans_mat, file = file.path("Matrix5_AncSR1.rda"))


## how well does the weighted matrix capture the effects of a nucleotide matrix?
codon_waa_cor <- numeric(length=100)

for(t in 1:100) {
  codon_trans_mat <- t(apply(codon_adj_mat, 1, function(x) x / sum(x)))
  aa_wtrans_mat <- t(apply(aa_wadj_mat, 1, function(x) x / sum(x)))
  codon_tstep_mat <- matrix(nrow = length(AA), ncol = length(AA))
  colnames(codon_tstep_mat) <- AA
  rownames(codon_tstep_mat) <- AA
  
  
  for(i in 1:length(AA)) {
    codons <- names(GENETIC_CODE)[GENETIC_CODE == AA[i]]
    v0 <- as.numeric(rownames(codon_trans_mat) %in% codons)
    v0 <- v0 / sum(v0)
    v1 <- v0 %*% mat_exp(codon_trans_mat, t)
    codon_tstep_mat[i,] <- sapply(AA, function(x) 
      sum(v1[colnames(v1) %in% names(GENETIC_CODE)[GENETIC_CODE == x]]))
  }
  
  codon_waa_cor[t] <- cor(as.vector(codon_tstep_mat),
                          as.vector(mat_exp(aa_wtrans_mat, t)))
}

plot(1:50, codon_waa_cor[1:50], type = 'l', xlab = "steps", ylab = "correlation", 
     main = "nt transition matrix vs. weighted aa transition matrix")


# how do the stationary distributions of each matrix compare?
# AncSR1
aa_lib_AncSR1_stat_dist <- (mat_exp(aa_lib_AncSR1_trans_mat, 100))[1,]
range(aa_lib_AncSR1_stat_dist - aa_lib_AncSR1_stat_dist %*% aa_lib_AncSR1_trans_mat)
names(aa_lib_AncSR1_stat_dist) <- colnames(aa_lib_AncSR1_trans_mat)

aa_lib_AncSR1_wstat_dist <- (mat_exp(aa_lib_AncSR1_wtrans_mat, 100))[1,]
range(aa_lib_AncSR1_wstat_dist - aa_lib_AncSR1_wstat_dist %*% aa_lib_AncSR1_wtrans_mat)
names(aa_lib_AncSR1_wstat_dist) <- colnames(aa_lib_AncSR1_wtrans_mat)

aa_lib_AncSR1_swstat_dist <- (mat_exp(aa_lib_AncSR1_swtrans_mat, 100))[1,]
range(aa_lib_AncSR1_swstat_dist - aa_lib_AncSR1_swstat_dist %*% aa_lib_AncSR1_swtrans_mat)
names(aa_lib_AncSR1_swstat_dist) <- colnames(aa_lib_AncSR1_swtrans_mat)

aa_lib_AncSR1_nswstat_dist <- (mat_exp(aa_lib_AncSR1_nswtrans_mat, 100))[1,]
range(aa_lib_AncSR1_nswstat_dist - aa_lib_AncSR1_nswstat_dist %*% aa_lib_AncSR1_nswtrans_mat)
names(aa_lib_AncSR1_nswstat_dist) <- colnames(aa_lib_AncSR1_nswtrans_mat)

plot(aa_lib_AncSR1_stat_dist, aa_lib_AncSR1_wstat_dist, 
     main = "AncSR1 stationary", xlab = "Matrix 2 (unweighted)", 
     ylab = "Matrix 3 (codon bias)")
abline(a=0, b=1, col="red")
mtext(paste("Pearson's R =", 
            round(cor(aa_lib_AncSR1_stat_dist, aa_lib_AncSR1_wstat_dist), 2)), 
      side = 3, adj = 0)
plot(aa_lib_AncSR1_wstat_dist, aa_lib_AncSR1_swstat_dist, 
     main = "AncSR1 stationary", xlab = "Matrix 3 (codon bias)", 
     ylab = "Matrix 4 (codon bias + synonymous mutations)")
abline(a=0, b=1, col="red")
mtext(paste("Pearson's R =", 
            round(cor(aa_lib_AncSR1_wstat_dist, aa_lib_AncSR1_swstat_dist), 2)), 
      side = 3, adj = 0)
plot(aa_lib_AncSR1_swstat_dist, aa_lib_AncSR1_nswstat_dist, 
     main = "AncSR1 stationary", xlab = "Matrix 4 (codon bias + synonymous mutations)", 
     ylab = "Matrix 5 (codon bias + synonymous + nonfunctional mutations)")
abline(a=0, b=1, col="red")
mtext(paste("Pearson's R =", 
            round(cor(aa_lib_AncSR1_swstat_dist, aa_lib_AncSR1_nswstat_dist), 2)), 
      side = 3, adj = 0)

# AncSR2
aa_lib_AncSR2_stat_dist <- (mat_exp(aa_lib_AncSR2_trans_mat, 100))[1,]
range(aa_lib_AncSR2_stat_dist - aa_lib_AncSR2_stat_dist %*% aa_lib_AncSR2_trans_mat)
names(aa_lib_AncSR2_stat_dist) <- colnames(aa_lib_AncSR2_trans_mat)

aa_lib_AncSR2_wstat_dist <- (mat_exp(aa_lib_AncSR2_wtrans_mat, 100))[1,]
range(aa_lib_AncSR2_wstat_dist - aa_lib_AncSR2_wstat_dist %*% aa_lib_AncSR2_wtrans_mat)
names(aa_lib_AncSR2_wstat_dist) <- colnames(aa_lib_AncSR2_wtrans_mat)

aa_lib_AncSR2_swstat_dist <- (mat_exp(aa_lib_AncSR2_swtrans_mat, 100))[1,]
range(aa_lib_AncSR2_swstat_dist - aa_lib_AncSR2_swstat_dist %*% aa_lib_AncSR2_swtrans_mat)
names(aa_lib_AncSR2_swstat_dist) <- colnames(aa_lib_AncSR2_swtrans_mat)

aa_lib_AncSR2_nswstat_dist <- (mat_exp(aa_lib_AncSR2_nswtrans_mat, 100))[1,]
range(aa_lib_AncSR2_nswstat_dist - aa_lib_AncSR2_nswstat_dist %*% aa_lib_AncSR2_nswtrans_mat)
names(aa_lib_AncSR2_nswstat_dist) <- colnames(aa_lib_AncSR2_nswtrans_mat)

plot(aa_lib_AncSR2_stat_dist, aa_lib_AncSR2_wstat_dist, 
     main = "AncSR2 stationary", xlab = "Matrix 2 (unweighted)", 
     ylab = "Matrix 3 (codon bias)")
abline(a=0, b=1, col="red")
mtext(paste("Pearson's R =", 
            round(cor(aa_lib_AncSR2_stat_dist, aa_lib_AncSR2_wstat_dist), 2)), 
      side = 3, adj = 0)
plot(aa_lib_AncSR2_wstat_dist, aa_lib_AncSR2_swstat_dist, 
     main = "AncSR2 stationary", xlab = "Matrix 3 (codon bias)", 
     ylab = "Matrix 4 (codon bias + synonymous mutations)")
abline(a=0, b=1, col="red")
mtext(paste("Pearson's R =", 
            round(cor(aa_lib_AncSR2_wstat_dist, aa_lib_AncSR2_swstat_dist), 2)), 
      side = 3, adj = 0)
plot(aa_lib_AncSR2_swstat_dist, aa_lib_AncSR2_nswstat_dist, 
     main = "AncSR2 stationary", xlab = "Matrix 4 (codon bias + synonymous mutations)", 
     ylab = "Matrix 5 (codon bias + synonymous + nonfunctional mutations)")
abline(a=0, b=1, col="red")
mtext(paste("Pearson's R =", 
            round(cor(aa_lib_AncSR2_swstat_dist, aa_lib_AncSR2_nswstat_dist), 2)), 
      side = 3, adj = 0)

# simulating evolution
# AncSR1
EGKA_AncSR1i <- which(rownames(aa_lib_AncSR1_trans_mat) == "EGKA")
EGKA_AncSR1_start <- rep(0, nrow(aa_lib_AncSR1_trans_mat))
EGKA_AncSR1_start[EGKA_AncSR1i] <- 1


# simulate evolution for 5 steps
AncSR1_5step <- as.vector(EGKA_AncSR1_start %*% (mat_exp(aa_lib_AncSR1_trans_mat, 5)))
AncSR1_w_5step <- as.vector(EGKA_AncSR1_start %*% (mat_exp(aa_lib_AncSR1_wtrans_mat, 5)))
AncSR1_sw_5step <- as.vector(EGKA_AncSR1_start %*% (mat_exp(aa_lib_AncSR1_swtrans_mat, 5)))
AncSR1_nsw_5step <- as.vector(EGKA_AncSR1_start %*% (mat_exp(aa_lib_AncSR1_nswtrans_mat, 5)))

plot(AncSR1_5step, AncSR1_w_5step, 
     main = "AncSR1 EGKA 5 steps", xlab = "Matrix 2 (unweighted)", 
     ylab = "Matrix 3 (codon bias)")
abline(a=0, b=1, col="red")
points(AncSR1_5step[EGKA_AncSR1i], AncSR1_w_5step[EGKA_AncSR1i], col = "magenta", pch=19)
mtext(paste("Pearson's R =", 
            round(cor(AncSR1_5step, AncSR1_w_5step), 2)), 
      side = 3, adj = 0)
plot(AncSR1_w_5step, AncSR1_sw_5step, 
     main = "AncSR1 EGKA 5 steps", xlab = "Matrix 3 (codon bias)", 
     ylab = "Matrix 4 (codon bias + synonymous mutations)")
abline(a=0, b=1, col="red")
points(AncSR1_w_5step[EGKA_AncSR1i], AncSR1_sw_5step[EGKA_AncSR1i], col = "magenta", pch=19)
mtext(paste("Pearson's R =", 
            round(cor(AncSR1_w_5step, AncSR1_sw_5step), 2)), 
      side = 3, adj = 0)
plot(AncSR1_sw_5step, AncSR1_nsw_5step, 
     main = "AncSR1 EGKA 5 steps", xlab = "Matrix 4 (codon bias + synonymous mutations)", 
     ylab = "Matrix 5 (codon bias + synonymous + nonfunctional mutations)")
abline(a=0, b=1, col="red")
points(AncSR1_sw_5step[EGKA_AncSR1i], AncSR1_nsw_5step[EGKA_AncSR1i], col = "magenta", pch=19)
mtext(paste("Pearson's R =", 
            round(cor(AncSR1_sw_5step, AncSR1_nsw_5step), 2)), 
      side = 3, adj = 0)


# AncSR2
EGKA_AncSR2i <- which(rownames(aa_lib_AncSR2_trans_mat) == "EGKA")
EGKA_AncSR2_start <- rep(0, nrow(aa_lib_AncSR2_trans_mat))
EGKA_AncSR2_start[EGKA_AncSR2i] <- 1

GSKV_AncSR2i <- which(rownames(aa_lib_AncSR2_trans_mat) == "GSKV")

# simulate evolution for 5 steps
AncSR2_5step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_trans_mat, 5)))
AncSR2_w_5step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_wtrans_mat, 5)))
AncSR2_sw_5step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_swtrans_mat, 5)))
AncSR2_nsw_5step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_nswtrans_mat, 5)))

plot(AncSR2_5step, AncSR2_w_5step, 
     main = "AncSR2 EGKA 5 steps", xlab = "Matrix 2 (unweighted)", 
     ylab = "Matrix 3 (codon bias)")
abline(a=0, b=1, col="red")
points(AncSR2_5step[EGKA_AncSR2i], AncSR2_w_5step[EGKA_AncSR2i], col = "magenta", pch=19)
mtext(paste("Pearson's R =", 
            round(cor(AncSR2_5step, AncSR2_w_5step), 2)), 
      side = 3, adj = 0)
plot(AncSR2_w_5step, AncSR2_sw_5step, 
     main = "AncSR2 EGKA 5 steps", xlab = "Matrix 3 (codon bias)", 
     ylab = "Matrix 4 (codon bias + synonymous mutations)")
abline(a=0, b=1, col="red")
points(AncSR2_w_5step[EGKA_AncSR2i], AncSR2_sw_5step[EGKA_AncSR2i], col = "magenta", pch=19)
mtext(paste("Pearson's R =", 
            round(cor(AncSR2_w_5step, AncSR2_sw_5step), 2)), 
      side = 3, adj = 0)
plot(AncSR2_sw_5step, AncSR2_nsw_5step, 
     main = "AncSR2 EGKA 5 steps", xlab = "Matrix 4 (codon bias + synonymous mutations)", 
     ylab = "Matrix 5 (codon bias + synonymous + nonfunctional mutations)")
abline(a=0, b=1, col="red")
points(AncSR2_sw_5step[EGKA_AncSR2i], AncSR2_nsw_5step[EGKA_AncSR2i], col = "magenta", pch=19)
mtext(paste("Pearson's R =", 
            round(cor(AncSR2_sw_5step, AncSR2_nsw_5step), 2)), 
      side = 3, adj = 0)






####### OTHER STUFF ###########

# # what is the average number of nucleotide steps to go from EGKA to GSKV on the AncSR2 network?
# 
# # lindep <- detect.lindep(as.matrix(aa_lib_AncSR2_swtrans_mat[-GSKV_AncSR2i,-GSKV_AncSR2i]) - 
# #                           diag(rep(1,nrow(aa_lib_AncSR2_swtrans_mat)-1)))
# # EGKA_GSKV_nt_steps_syn <- solve((aa_lib_AncSR2_swtrans_mat[-GSKV_AncSR2i,-GSKV_AncSR2i] - 
# #                                    diag(rep(1,nrow(aa_lib_AncSR2_swtrans_mat)-1))), 
# #                                 rep(-1, nrow(aa_lib_AncSR2_swtrans_mat)-1))
# 
# 
# # subset of only AAs that are in direct path between EGKA and GSKV
# subnet <- c("EGKA", "EGKV", "ESKA", "ESKV", "GGKA", "GGKV", "GSKA", "GSKV")
# # subset of AAs that have K at position 3
# # subnet <- which(sapply(strsplit(aa_lib_AncSR2, ""), function(x) x[3] == "K"))
# # subnet <- aa_lib_AncSR2[subnet]
# 
# # syn transition matrix
# EGKA_GSKV_subnet <- aa_lib_AncSR2_swtrans_mat[rownames(aa_lib_AncSR2_swtrans_mat) %in% subnet,
#                                               colnames(aa_lib_AncSR2_swtrans_mat) %in% subnet]
# 
# # no back-mutations
# EGKA_GSKV_subnet_no_back <- EGKA_GSKV_subnet
# EGKA_GSKV_subnet_no_back[lower.tri(EGKA_GSKV_subnet_no_back)] <- 0
# 
# # renormalize
# EGKA_GSKV_subnet <- t(apply(EGKA_GSKV_subnet, 1, function(x) x / sum(x)))
# EGKA_GSKV_subnet_no_back <- t(apply(EGKA_GSKV_subnet_no_back, 1, function(x) x / sum(x)))
# 
# # solve for average time (number of nucleotide mutations) from EGKA to GSKV
# # lindep <- detect.lindep(EGKA_GSKV_subnet[-GSKV_AncSR2_subneti,-GSKV_AncSR2_subneti] - 
# #                           diag(rep(1,nrow(EGKA_GSKV_subnet)-1)))
# # EGKA_GSKV_nt_steps_syn <- solve((EGKA_GSKV_subnet[-GSKV_AncSR2_subneti,-GSKV_AncSR2_subneti] - 
# #                                   diag(rep(1,nrow(EGKA_GSKV_subnet)-1)))[-lindep,-lindep], 
# #                                 rep(-1, nrow(EGKA_GSKV_subnet)-1-length(lindep)))
# EGKA_GSKV_nt_steps_syn <- solve((EGKA_GSKV_subnet[1:3,1:3] - diag(rep(1,3))), rep(-1, 3))[1]
# 
# EGKA_GSKV_nt_steps_syn_no_back <- solve(EGKA_GSKV_subnet_no_back[1:3,1:3] - diag(rep(1,3)), rep(-1, 3))[1]
# 
# 
# ### weighted matrix, no synonymous
# EGKA_GSKV_subnet <- aa_lib_AncSR2_wtrans_mat[rownames(aa_lib_AncSR2_wtrans_mat) %in% subnet,
#                                              colnames(aa_lib_AncSR2_wtrans_mat) %in% subnet]
# # no back-mutations
# EGKA_GSKV_subnet_no_back <- EGKA_GSKV_subnet
# EGKA_GSKV_subnet_no_back[lower.tri(EGKA_GSKV_subnet_no_back)] <- 0
# 
# # renormalize
# EGKA_GSKV_subnet <- t(apply(EGKA_GSKV_subnet, 1, function(x) x / sum(x)))
# EGKA_GSKV_subnet_no_back <- t(apply(EGKA_GSKV_subnet_no_back, 1, function(x) x / sum(x)))
# 
# # solve for average time (number of nucleotide mutations) from EGKA to GSKV
# EGKA_GSKV_nt_steps_nosyn <- solve(EGKA_GSKV_subnet[1:3,1:3] - diag(rep(1,3)), rep(-1, 3))[1]
# EGKA_GSKV_nt_steps_nosyn_no_back <- solve(EGKA_GSKV_subnet_no_back[1:3,1:3] - diag(rep(1,3)), rep(-1, 3))[1]
# 
# 
# ### unweighted submatrix
# EGKA_GSKV_subnet <- aa_lib_AncSR2_trans_mat[rownames(aa_lib_AncSR2_trans_mat) %in% subnet,
#                                             colnames(aa_lib_AncSR2_trans_mat) %in% subnet]
# # no back-mutations
# EGKA_GSKV_subnet_no_back <- EGKA_GSKV_subnet
# EGKA_GSKV_subnet_no_back[lower.tri(EGKA_GSKV_subnet_no_back)] <- 0
# 
# # renormalize
# EGKA_GSKV_subnet <- t(apply(EGKA_GSKV_subnet, 1, function(x) x / sum(x)))
# EGKA_GSKV_subnet_no_back <- t(apply(EGKA_GSKV_subnet_no_back, 1, function(x) x / sum(x)))
# 
# # solve for average time (number of nucleotide mutations) from EGKA to GSKV
# EGKA_GSKV_nt_steps_unw <- solve(EGKA_GSKV_subnet[1:3,1:3] - diag(rep(1,3)), rep(-1, 3))[1]
# EGKA_GSKV_nt_steps_unw_no_back <- solve(EGKA_GSKV_subnet_no_back[1:3,1:3] - diag(rep(1,3)), rep(-1, 3))[1]
# 
# 
# ### weighted + synonymous + nonfunctinal matrix
# EGKA_GSKV_subnet <- aa_lib_AncSR2_nswtrans_mat[rownames(aa_lib_AncSR2_nswtrans_mat) %in% subnet,
#                                                colnames(aa_lib_AncSR2_nswtrans_mat) %in% subnet]
# # no back-mutations
# EGKA_GSKV_subnet_no_back <- EGKA_GSKV_subnet
# EGKA_GSKV_subnet_no_back[lower.tri(EGKA_GSKV_subnet_no_back)] <- 0
# 
# # renormalize
# EGKA_GSKV_subnet <- t(apply(EGKA_GSKV_subnet, 1, function(x) x / sum(x)))
# EGKA_GSKV_subnet_no_back <- t(apply(EGKA_GSKV_subnet_no_back, 1, function(x) x / sum(x)))
# 
# # solve for average time (number of nucleotide mutations) from EGKA to GSKV
# EGKA_GSKV_nt_steps_nsw <- solve(EGKA_GSKV_subnet[1:3,1:3] - diag(rep(1,3)), rep(-1, 3))[1]
# EGKA_GSKV_nt_steps_nsw_no_back <- solve(EGKA_GSKV_subnet_no_back[1:3,1:3] - diag(rep(1,3)), rep(-1, 3))[1]
# 
# 
# 
# 
# ### comparisons of normalized trajectories
# 
# # simulate evolution for EGKA to GSKV steps
# AncSR2_3step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_trans_mat, 3)))
# AncSR2_w_3step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_wtrans_mat, 3)))
# AncSR2_sw_3step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_swtrans_mat, 3)))
# AncSR2_nsw_3step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_nswtrans_mat, 3)))
# AncSR2_sw_41step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_swtrans_mat, 41)))
# AncSR2_nsw_112step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_nswtrans_mat, 112)))
# names(AncSR2_nsw_112step) <- rownames(aa_lib_AncSR2_nswtrans_mat)
# 
# plot(AncSR2_3step, AncSR2_w_3step, 
#      main = "AncSR2 EGKA normalized steps", xlab = "Matrix 2 (unweighted)", 
#      ylab = "Matrix 3 (codon bias)")
# abline(a=0, b=1, col="red")
# points(AncSR2_3step[c(EGKA_AncSR2i, GSKV_AncSR2i)], 
#        AncSR2_w_3step[c(EGKA_AncSR2i, GSKV_AncSR2i)], 
#        col = c("magenta", "green"), pch=19)
# mtext(paste("Pearson's R =", 
#             round(cor(AncSR2_3step, AncSR2_w_3step), 2)), 
#       side = 3, adj = 0)
# plot(AncSR2_w_3step, AncSR2_sw_41step, 
#      main = "AncSR2 EGKA normalized steps", xlab = "Matrix 3 (codon bias)", 
#      ylab = "Matrix 4 (codon bias + synonymous mutations)"
#      # , xlim=c(0, 0.04), ylim=c(0,0.04)
#      )
# abline(a=0, b=1, col="red")
# points(AncSR2_w_3step[c(EGKA_AncSR2i, GSKV_AncSR2i)], 
#        AncSR2_sw_41step[c(EGKA_AncSR2i, GSKV_AncSR2i)], 
#        col = c("magenta", "green"), pch=19)
# mtext(paste("Pearson's R =", 
#             round(cor(AncSR2_w_3step, AncSR2_sw_41step), 2)), 
#       side = 3, adj = 0)
# # plot(AncSR2_w_3step, AncSR2_sw_3step, 
# #      main = "AncSR2 EGKA 3 steps", xlab = "Matrix 3 (codon bias)", 
# #      ylab = "Matrix 4 (codon bias + synonymous mutations)"
# #      # , xlim=c(0, 0.04), ylim=c(0,0.04)
# # )
# # abline(a=0, b=1, col="red")
# # points(AncSR2_w_3step[GSKV_AncSR2i], AncSR2_sw_3step[GSKV_AncSR2i], col = "green", pch=19)
# # mtext(paste("Pearson's R =", 
# #             round(cor(AncSR2_w_3step, AncSR2_sw_3step), 2)), 
# #       side = 3, adj = 0)
# # plot(AncSR2_sw_3step, AncSR2_nsw_3step, 
# #      main = "AncSR2 EGKA 3 steps", xlab = "Matrix 4 (codon bias + synonymous mutations)", 
# #      ylab = "Matrix 4 (codon bias + synonymous + nonfunctional mutations)"
# #      # , xlim=c(0, 0.04), ylim=c(0,0.04)
# # )
# # abline(a=0, b=1, col="red")
# # points(AncSR2_sw_3step[GSKV_AncSR2i], AncSR2_nsw_3step[GSKV_AncSR2i], col = "green", pch=19)
# # mtext(paste("Pearson's R =", 
# #             round(cor(AncSR2_sw_3step, AncSR2_nsw_3step), 2)), 
# #       side = 3, adj = 0)
# plot(AncSR2_sw_41step, AncSR2_nsw_112step, 
#      main = "AncSR2 EGKA normalized steps", xlab = "Matrix 4 (codon bias + synonymous mutations)", 
#      ylab = "Matrix 5 (codon bias + synonymous + nonfunctional mutations)")
# abline(a=0, b=1, col="red")
# points(AncSR2_sw_41step[c(EGKA_AncSR2i, GSKV_AncSR2i)], 
#        AncSR2_nsw_112step[c(EGKA_AncSR2i, GSKV_AncSR2i)], 
#        col = c("magenta", "green"), pch=19)
# mtext(paste("Pearson's R =", 
#             round(cor(AncSR2_sw_41step, AncSR2_nsw_112step), 2)), 
#       side = 3, adj = 0)
# 
# 
# 
# # with back-mutation steps
# AncSR2_9step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_trans_mat, 9)))
# AncSR2_w_10step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_wtrans_mat, 10)))
# AncSR2_sw_80step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_swtrans_mat, 80)))
# AncSR2_nsw_208step <- as.vector(EGKA_AncSR2_start %*% (mat_exp(aa_lib_AncSR2_nswtrans_mat, 208)))
# 
# plot(AncSR2_9step, AncSR2_w_10step, 
#      main = "AncSR2 EGKA normalized steps no back", xlab = "Matrix 2 (unweighted)", 
#      ylab = "Matrix 3 (codon bias)")
# abline(a=0, b=1, col="red")
# points(AncSR2_9step[c(EGKA_AncSR2i, GSKV_AncSR2i)], 
#        AncSR2_w_10step[c(EGKA_AncSR2i, GSKV_AncSR2i)], 
#        col = c("magenta", "green"), pch=19)
# mtext(paste("Pearson's R =", 
#             round(cor(AncSR2_9step, AncSR2_w_10step), 2)), 
#       side = 3, adj = 0)
# 
# plot(AncSR2_w_10step, AncSR2_sw_80step, 
#      main = "AncSR2 EGKA normalized steps no back", xlab = "Matrix 3 (codon bias)", 
#      ylab = "Matrix 4 (codon bias + synonymous mutations)")
# abline(a=0, b=1, col="red")
# points(AncSR2_w_10step[c(EGKA_AncSR2i, GSKV_AncSR2i)], 
#        AncSR2_sw_80step[c(EGKA_AncSR2i, GSKV_AncSR2i)], 
#        col = c("magenta", "green"), pch=19)
# mtext(paste("Pearson's R =", 
#             round(cor(AncSR2_w_10step, AncSR2_sw_80step), 2)), 
#       side = 3, adj = 0)
# 
# 
# 
# 
# AncSR1_nsw_112step <- as.vector(EGKA_AncSR1_start %*% (mat_exp(aa_lib_AncSR1_nswtrans_mat, 112)))
# names(AncSR1_nsw_112step) <- rownames(aa_lib_AncSR1_nswtrans_mat)
# 
# plot(AncSR1_nsw_112step, AncSR2_nsw_112step[names(AncSR1_nsw_112step)])
# abline(a=0, b=1, col="red")
# 
# plot(AncSR1_nsw_112step, aa_lib_AncSR1_nswstat_dist)
# 
