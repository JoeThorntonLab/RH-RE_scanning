#!/usr/bin/env Rscript
#####

setwd("/project/joet1/santiagoherrera/DMS_RH_RE/221202_A00639_1438_AHLYYJDRX2-JT-JP-SHA-1s/results/MutSel_v2")

# LOAD PACKAGES #

packages <- c("tidyr", "MASS", "ggplot2", "Matrix", "stringr", "tibble", 
              "dplyr", "patchwork", "foreach", "doParallel", "matrixStats",
              "Biostrings", "readr","furrr","purrr","igraph","recommenderlab","Matrix")
installed_packages <- packages %in% rownames(installed.packages())

if(any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}

# load packages
invisible(lapply(packages, library, character.only=TRUE))

#####

# LOAD FUNCTIONS #
source("./MC_MutSel_functions.R")

# LOAD DATA #

meanF_data <- readr::read_csv("../meanF_data_corrected_NovaSeq.csv.gz") %>% 
  dplyr::select(AA_var,REBC,bg,RE,type,avg_meanF,padj,sig)


# For each variant:
phenotypes_tbl <- meanF_data %>% filter(type != "control" & sig == "significant") %>%
  group_by(AA_var, bg) %>%
  summarise(n_bound_REs = n(), # how many DNA elements can bind
            meanF_bREs = mean(avg_meanF), # the average meanF across all DNA elements bound
            max_meanF_bREs = max(avg_meanF), # the max meanF across all DNA elements bound
            min_meanF_bREs = min(avg_meanF), # the min meanF across all DNA elements bound 
            specific = ifelse(n_bound_REs == 1, "YES","NO"), # functional binding to only one DNA element?
            specificity = ifelse(n_bound_REs == 1, RE, "Promiscuous"), # Determine the type of specificity
            promiscuity =  ifelse(specific=="NO",n_bound_REs * max_meanF_bREs/sum(avg_meanF),0), # Level of promiscuity: promiscuous (1), specific (0)
            bound_REs = list(RE)) # assign DNA elements bound


#####

# GLOBAL PARAMETERS #

# Reference wild-type ancestral genotypes:
AncSR1_ERE_ref <- meanF_data %>% filter(AA_var == "EGKA" & REBC == "AncSR1_REBC3") %>% pull(avg_meanF)
AncSR2_SRE_ref <- meanF_data %>% filter(AA_var == "GSKV" & REBC == "AncSR2_REBC1") %>% pull(avg_meanF)
Ne = 10e8 # population size (reasonable average for chordates e.g.,https://doi.org/10.7554/eLife.67509)

# STEP FUNCTION PARAMETERS:
# Initial parameters correspond to the same as the logistic function. 'mF_ref' = Minimal meanF for active variants
# Fitness corresponds to the exponential growth rate = N*r
MIN_ACTIVE <- phenotypes_tbl %>% with(min(min_meanF_bREs))
STEP.PARAM = c(MIN_ACTIVE)

# number of cores
N_CORES=detectCores()-1

#####

# BUILD GENOTYPE NETWORKS #

sr1_variants <- phenotypes_tbl %>% filter(bg == "AncSR1") %>% pull(AA_var)
sr2_variants <- phenotypes_tbl %>% filter(bg == "AncSR2") %>% pull(AA_var)

# Weighted networks for matrix A
net_sr1_weights <- build_genotype_network(nodes=sr1_variants,weighted="fraction",cores=N_CORES)
net_sr2_weights <- build_genotype_network(nodes=sr2_variants,weighted="fraction",cores=N_CORES)

# Regular networks
adj_mat_sr1 <- build_mutation_matrix(sr1_variants,type=1,N_CORES)
net_sr1 <- build_genotype_network(build_mat=FALSE,adj_mat=adj_mat_sr1)

adj_mat_sr2 <- build_mutation_matrix(sr2_variants,type=1,N_CORES)
net_sr2 <- build_genotype_network(build_mat=FALSE,adj_mat=adj_mat_sr2)

# BUILD TRANSITION PROBABILITY MATRICES #

# Generate P matrices for drift scenario under three approaches:
# A - Mutation as fraction of codons. Using original functions (see 'MutSel_matricesA.R' script. This is like a positive control)
# B - Mutation as fraction of codons. Using new functions (should be exactly the same as A)
# C - Mutation as number of mut. paths. Using new functions (compare how similar is it to A/B)

# Specify models for each genetic background
MODEL.PARAM_SR1_drift <- list("AncSR1","drift",Ne,STEP.PARAM,TRUE,"max") # AncSR1
MODEL.PARAM_SR2_drift <- list("AncSR2","drift",Ne,STEP.PARAM,TRUE,"max") # AncSR2

# Build matrices
# AncSR1
M_drift_sr1_A <- build_transition_matrix(sr1_variants,net_sr1_weights,phenotypes_tbl,MODEL.PARAM_SR1_drift,N_CORES)

adj_mat_fraction_sr1 <- build_mutation_matrix(sr1_variants,type=2,N_CORES)
M_drift_sr1_B <- build_transition_matrix_v2(sr1_variants,adj_mat_fraction_sr1,phenotypes_tbl,MODEL.PARAM_SR1_drift,N_CORES)

adj_mat_count_sr1 <- build_mutation_matrix(sr1_variants,type=3,N_CORES)
M_drift_sr1_C <- build_transition_matrix_v2(sr1_variants,adj_mat_count_sr1,phenotypes_tbl,MODEL.PARAM_SR1_drift,N_CORES)

# AncSR2
M_drift_sr2_A <- build_transition_matrix(sr2_variants,net_sr2_weights,phenotypes_tbl,MODEL.PARAM_SR2_drift,N_CORES)

adj_mat_fraction_sr2 <- build_mutation_matrix(sr2_variants,type=2,N_CORES)
M_drift_sr2_B <- build_transition_matrix_v2(sr2_variants,adj_mat_fraction_sr2,phenotypes_tbl,MODEL.PARAM_SR2_drift,N_CORES)

adj_mat_count_sr2 <- build_mutation_matrix(sr2_variants,type=3,N_CORES)
M_drift_sr2_C <- build_transition_matrix_v2(sr2_variants,adj_mat_count_sr2,phenotypes_tbl,MODEL.PARAM_SR2_drift,N_CORES)


save.image(file = "./matrix_comparison.RData")

