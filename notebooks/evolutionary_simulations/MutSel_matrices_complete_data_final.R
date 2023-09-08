#!/usr/bin/env Rscript
#####

setwd("/project/joet1/santiagoherrera/DMS_RH_RE/221202_A00639_1438_AHLYYJDRX2-JT-JP-SHA-1s/results/MutSel_final")

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

# Complete data: measured + predicted meanF; only functional variants (see 'results/genotype_phentoype_distributions' folder)
meanF_data <- readr::read_csv(file.path("..","mutation_effects_model","meanF_data_fxnal.csv.gz"))

#####

# GLOBAL PARAMETERS #
AncSR1_ERE_ref <- meanF_data %>% filter(AA_var == "EGKA" & bg == "AncSR1" & RE == "ERE (GT)") %>% pull(avg_meanF)
AncSR2_SRE_ref <- meanF_data %>% filter(AA_var == "GSKV" & bg == "AncSR2" & RE == "SRE1 (AA)") %>% pull(avg_meanF)
Ne = 10e8 # population size (reasonable average for chordates e.g.,https://doi.org/10.7554/eLife.67509)
N_CORES=detectCores()-1 # number of cores

# BUILD PHENOTYPE TABLE #

# For each *functional* variant (those with meanF >= AncSR2_SRE_ref):
phenotypes_tbl <- meanF_data %>% group_by(AA_var, bg) %>%
  summarise(n_bound_REs = n(), # how many DNA elements can bind
            meanF_bREs = mean(avg_meanF), # the average meanF across all DNA elements bound
            max_meanF_bREs = max(avg_meanF), # the max meanF across all DNA elements bound
            min_meanF_bREs = min(avg_meanF), # the min meanF across all DNA elements bound 
            specific = ifelse(n_bound_REs == 1, "YES","NO"), # functional binding to only one DNA element?
            specificity = ifelse(n_bound_REs == 1, RE, "Promiscuous"), # Determine the type of specificity
            promiscuity =  ifelse(specific=="NO",n_bound_REs * max_meanF_bREs/sum(avg_meanF),0), # Level of promiscuity: promiscuous (1), specific (0)
            bound_REs = list(RE)) # assign DNA elements bound

# GLOBAL PARAMETERS FITNESS FUNCTIONS #

# LOGISTIC FUNCTION PARAMETERS:
# Initial parameters of the logistic curve based on maximum likelihood estimates from DMS+ASR on SR phylogeny (fiting deltaF values with respect to AncSR1 reference)
# Fitness corresponds to the exponential growth rate = N*r
L = 2.7915391 # maximum fitness
k = 4.33760376 # logistic grouth rate (steepness)
x_o = -0.54577352 # midpoint
LOG.PARAM = c(L,k,x_o)

# NORMAL FUNCTION PARAMETERS:
# Initial parameters correspond to the mean (and sd) fluorescence of ERE specific variants.
# Fitness corresponds to the exponential growth rate = N*r
u = phenotypes_tbl %>% filter(specific == "YES" & bound_REs == "ERE (GT)") %>% with(mean(max_meanF_bREs)) # mean
sd = phenotypes_tbl %>% filter(specific =="YES" & bound_REs == "ERE (GT)") %>% with(sd(max_meanF_bREs)) # stdev
scale = L # make the maximum fitness match the maximum fitness of logistic function
NORM.PARAM = c(u,sd,scale)

# STEP FUNCTION PARAMETERS:
# Initial parameters correspond to the same as the logistic function. 'mF_ref' = Minimal meanF for active variants
# Fitness corresponds to the exponential growth rate = N*r
MIN_ACTIVE <- phenotypes_tbl %>% with(min(min_meanF_bREs))
STEP.PARAM = c(MIN_ACTIVE)

#####

# BUILD GENOTYPE NETWORKS #

sr1_variants <- phenotypes_tbl %>% filter(bg == "AncSR1") %>% pull(AA_var)
sr2_variants <- phenotypes_tbl %>% filter(bg == "AncSR2") %>% pull(AA_var)

# Regular networks
adj_mat_sr1 <- build_mutation_matrix(sr1_variants,type=1,N_CORES)
net_sr1 <- build_genotype_network(build_mat=FALSE,adj_mat=adj_mat_sr1)

adj_mat_sr2 <- build_mutation_matrix(sr2_variants,type=1,N_CORES)
net_sr2 <- build_genotype_network(build_mat=FALSE,adj_mat=adj_mat_sr2)

#####

# BUILD TRANSITION PROBABILITY MATRICES #

# Generate transition probability matrix under three evolutionary scenarios for AncSR1 and AncSR2 backgrounds (From = rows, to = cols)

# Build model scenarios for each DBD background:
# AncSR1
MODEL.PARAM_SR1_drift <- list("AncSR1","drift",Ne,STEP.PARAM,TRUE,"max")
MODEL.PARAM_SR1_directional <- list("AncSR1","directional",Ne,LOG.PARAM,TRUE,"max")
#MODEL.PARAM_SR1_stabilizing <- list("AncSR1","stabilizing",Ne,NORM.PARAM,TRUE,"max")

# AncSR2
MODEL.PARAM_SR2_drift <- list("AncSR2","drift",Ne,STEP.PARAM,TRUE,"max")
MODEL.PARAM_SR2_directional <- list("AncSR2","directional",Ne,LOG.PARAM,TRUE,"max")
#MODEL.PARAM_SR2_stabilizing <- list("AncSR2","stabilizing",Ne,NORM.PARAM,TRUE,"max")

# Build matrices
# AncSR1
adj_mat_count_sr1 <- build_mutation_matrix(sr1_variants,type=3,N_CORES)
M_drift_sr1 <- build_transition_matrix_v2(sr1_variants,adj_mat_count_sr1,phenotypes_tbl,MODEL.PARAM_SR1_drift,N_CORES,complex = F)
M_dir_sr1 <- build_transition_matrix_v2(sr1_variants,adj_mat_count_sr1,phenotypes_tbl,MODEL.PARAM_SR1_directional,N_CORES,complex = F)

# AncSR2
adj_mat_count_sr2 <- build_mutation_matrix(sr2_variants,type=3,N_CORES)
M_drift_sr2 <- build_transition_matrix_v2(sr2_variants,adj_mat_count_sr2,phenotypes_tbl,MODEL.PARAM_SR2_drift,N_CORES,complex = F)
M_dir_sr2 <- build_transition_matrix_v2(sr2_variants,adj_mat_count_sr2,phenotypes_tbl,MODEL.PARAM_SR2_directional,N_CORES,complex = F)

# convert adjacency matrices of mutation rates to transition probability matrices
adj_mat_count_sr1 <- t(apply(adj_mat_count_sr1, 1, function(x) x / sum(x)))
adj_mat_count_sr1 <- replace(adj_mat_count_sr1,is.nan(adj_mat_count_sr1),0)
adj_mat_count_sr1 <- as(adj_mat_count_sr1, "sparseMatrix")

adj_mat_count_sr2 <- t(apply(adj_mat_count_sr2, 1, function(x) x / sum(x)))
adj_mat_count_sr2 <- replace(adj_mat_count_sr2,is.nan(adj_mat_count_sr2),0)
adj_mat_count_sr2 <- as(adj_mat_count_sr2, "sparseMatrix")

save.image(file = "./MutSel_matrices_final.RData")

