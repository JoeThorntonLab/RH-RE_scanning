



This report was automatically generated with the R package **knitr**
(version 1.43).


```r
---
title: "Evolutionary simulations on empirical GP maps"
author: "Santiago Herrera"
date: "2023-08-23"
output: github_document
editor_options:
  chunk_output_type: console
---
```

```
## Error: <text>:10:0: unexpected end of input
## 8: ---
## 9: 
##   ^
```


This notebook contains the analysis for simulating evolutionary trajectories on empirical GP maps using markov chains. It has three mains sections:

- A theoretical explanation for modeling evolution on empirical GP maps [here](#strong-selection-weak-mutation-regime)
- The simulations on the protein genotype networks [here](#evolutionary-simulations-using-discrete-markov-chains)
- The simulations on the protein-DNA complex genotype networks [here](#evolution-of-protein-dna-complexes)

---

### Reading in the data

*Note:* The transition matrices and genotype networks contained in the `MutSel_matrices_final.RData` file were generated in the Midway3 cluster using the R script `MutSel_matrices_complete_data_final.R`.


```r
# Check whether the matrices have been created already
if(!file.exists(file.path(".", "MutSel_matrices_final.RData"))) {
  # Load complete data from mutation effects model
  meanF_data <- readr::read_csv(file.path("..","..","results","genotype_phenotype_distributions","meanF_data_fxnal.csv.gz"))
  
  # Reference wild-type ancestral genotypes:
  AncSR1_ERE_ref <- meanF_data %>% filter(AA_var == "EGKA" & bg == "AncSR1" & RE == "ERE (GT)") %>% pull(avg_meanF)
  AncSR2_SRE_ref <- meanF_data %>% filter(AA_var == "GSKV" & bg == "AncSR2" & RE == "SRE1 (AA)") %>% pull(avg_meanF)
  
  # BUILD PHENOTYPE TABLE #
  # For each *functional* variant (those with meanF >= AncSR2_SRE_ref):
  phenotypes_tbl <- meanF_data %>% group_by(AA_var, bg) %>%
    summarise(n_bound_REs = n(), # how many DNA elements can bind
              meanF_bREs = mean(avg_meanF), # the average meanF across all DNA elements bound
              max_meanF_bREs = max(avg_meanF), # the max meanF across all DNA elements bound
              min_meanF_bREs = min(avg_meanF), # the min meanF across all DNA elements bound 
              specific = ifelse(n_bound_REs == 1, "YES","NO"), # functional binding to only one DNA element?
              specificity = ifelse(n_bound_REs == 1, RE, "Promiscuous"), # Determine the type of specificity
              bound_REs = list(RE)) # assign DNA elements bound
  
  #phenotypes_tbl_prot <- phenotypes_tbl # save it for later

} else {
  # load matrices if already created
  load("./MutSel_matrices_final.RData")
  phenotypes_tbl_prot <- phenotypes_tbl # save it for later
}

# load functions
source("../MC_MutSel_functions.R")
```


## Strong Selection-Weak Mutation regime

Under the SSWM regime, the mutation rate is low enough such that the time to fixation of a mutation is much lower than the time between subsequent mutations (Gillespie, 1984). Thus, trajectories on a genotype landscape can be modeled as an stepwise, origin-fixation process. Specifically, the rate of fixation from allele $i$ to $j$ is a product of the rate of introduction of allele $j$ in the population and the probability that it goes to fixation, like so 

$$
Q(i,j)=2N_e\mu_{ij} \times P_{\text{fix}}(j;s_{ij},N_e),
$$

where $\mu_ {ij}$ is the mutation rate from allele $i$ to $j$, $s_ {ij}$ is the selection coefficient, and $N_e$ is the effective population size. When $s_ {ij} \neq 0$, the fixation probability, $P_{\text{fix}}(j;s_ {ij},N_e)$ is given by the Kimura equation (Kimura, 1962):

$$
P_{\text{fix}}(j) = \begin{cases}
\frac{1-e^{-2s_{ij}}}{1-e^{-4Ns_{ij}}} & \text{when } s_{ij} \neq 0 \\
\frac{1}{2N_e} & \text{when } s_{ij} = 0
\end{cases}
$$

To precisely model an origin-fixation walk on an empirical genotype-phenotype (GP) map, we need to compute $P(i,j)$, the probability that the *next* mutation will be $j$. That is, we need to account for the local structure of the network around the focal node $i$. Thus, 

$$
P(i,j) = \frac{Q(i,j)}{ \sum_{k \neq i} Q(i,k)} ,
$$

where $k$ are all nodes connected to node $i$ (McCandlish and Stoltzfus, 2014).

**Some assumptions:** Assuming a constant population size ($N_e$), and no mutation bias and reversibility ($\mu_ {ij} = \mu_ {ji} = \mu_ {ik}$), the probability that the next mutation will be $j$ becomes a function of the rescaled fixation probabilities alone, such that 

$$
P(i,j) = \frac{2N_e\mu_{ij} \times P_{\text{fix}}(j)}{\sum_{k \neq i}2N_e\mu_{ik} \times P_{\text{fix}}(k)} \
= (\frac{2N_e\mu_{ij}}{2N_e\mu_{ik}}) \times \frac{P_{\text{fix}}(j)}{\sum_{k \neq i} P_{\text{fix}}(k)} \\
= \frac{P_{\text{fix}}(j)}{\sum_{k \neq i}P_{\text{fix}}(k)}
$$

Although we are assuming no mutation bias and reversibility, we can weight the fixation probablities by a *mutational accessibility rate* (MAR). The MAR aims to capture an important aspect of the mutation process: the mapping from codon-to-amino acid, thus capturing the accessibility of genotype variants as mediated through the structure of the genetic code.

MARs ($\rho_ {ij}$) can be defined in multiple ways, but we use two equivalent definitions:

- All single-step amino acid mutations to functional genotypes given the genetic code are accessible, with probability equal to the *fraction of codons* from amino acid $i$ that can access amino acid $j$ via single nucleotide mutations (i.e., mutational propensity).
- All single-step amino acid mutations to functional genotypes given the genetic code are accessible, with probability proportional to the *number of nucleotide changes* that can encode each amino acid change - accounting for all the possible synonymous backgrounds (i.e., codon bias).

A reasonable assumption in both cases is that a population fixed for a given amino acid genotype, may explore all synonymous codons; i.e., the population is *delocalized* at the nucleotide level but fixed at the amino acid level.

Re-writing the previous equation to include mutational propensities, we have:

$$
P(i,j) = \frac{\rho_{ij} \times P_{\text{fix}}(j)}{\sum_{k \neq i}(\rho_{ik} \times P_{\text{fix}}(k))}
$$

This specification of an origin-fixation model allows to model evolution on a genotype network as a discrete Markov process, where each time step consists of a single *amino acid* substitution.  

### Molecular evolutionary scenarios

We are interested in characterizing the probability distribution of functional variation (PDFV) accessible around a particular genotype. Two types of PDFV can be distinguished: a *variational PDFV*, which refers to the probability distribution of intrinsic tendency of variation, summarizing how the system can change solely by mutation given the structure of the GP map and the genetic code; and an *evolvable PDFV*, determined by the mutation-only PDFV, which summarizes the probability distribution of variation in response to (purifying, directional or stabilizing) selection and drift.

We will simulate walks on the GP map under three selective regimes: 1) Random walk, 2) Directional selection, and 3) Stabilizing selection. These regimes act on the mean fluorescence estimated for each variant. 

- *Random walk on sequence space* (Maynard-Smith, 1970; King and Jukes, 1969): Under this model, proteins can only traverse the sequence space if they form a continuous network without crossing through nonfunctional intermediate nodes. There exists a minimum level of functionality that protein variants must have, below which variants are removed by purifying selection. Variants with a higher level of function are equivalent and thus their evolution is affected solely by genetic drift. (The other two models are built from the basic assumptions of this model - i.e., continuous network of functional variants, and drift+purifying selection are always present).

- *Directional selection* (Gillespie, 1984): $P(i,j) \propto s_{ij}$, the strength of the selection coefficient. There is a positive relationship between phenotype and fitness, thus, selection will act to increase the phenotypic value along the trajectory. However, we won't use a linear fitness function, because it assumes that fitness increases (or decreases) infinitely and linearly with $s_{ij}$. In contrast, we will use a logistic fitness function which asumes that after a certain increase (or decrease) in meanF, the fitness is unlikely to change.

- *Stabilizing selection*: $P(i,j) \propto s_{ij}$. This model assumes that there is an (static) fitness optimum, thus, selection will act to maintain a given phenotypic value along the trajectory. In this case, a bell-shaped function centered around a reference phenotypic value allows to assign low selection coefficients to mutations that decrease or increase too much meanF. 


### Calculating selection coefficients

There are effectively infinite ways of relating a phenotype to fitness. Two key aspects of the phenotype-to-fitness map are 1) the shape of the function and 2) the values of the parameters of the funciton. To find the parameters of the fitness functions, we relied on a previous study that performed ancestral sequence reconstruction on the steroid receptor (SR) phylogeny using experimentally informed models of sequence evolution (Muñiz et al. in prep). Briefly, the study used a DMS library of AncSR1, the ancestor of the SR family, which measured the functional effect of every single amino acid substitution at every site in the protein, and used this information to build a mutation-selection model to optimize branch lengths and reconstruct ancestral sequences on the tree. The method infers the maximum likelihood estimates of the relevant population-level parameters of the model given the experimental data and the tree topology.

We used Muñiz et al.'s dataset to estimate the parameters of two relevant fitness functions: 1) a logistic function ($F(g)_{\text{Log}}$) and 2) a bell-shaped function ($F(g)_{\text{Norm}}$) (see the section above on *Molecular evolutionary scenarios* for details on the selection of these functions). Each function has three parameters that can be estimated via maximum likelihood, and are defined in terms of the effect of mutations on the phenotype relative to the AncSR1 reference, i.e. mean fluorescence ($\Delta F$):

$$
F(g)_{\text{Log}} = \frac{L}{1 + e^{-k(\Delta F - \Delta F_0)}}
$$

where, $L$ is the maximum fitness, $k$ is the steepness, and $\Delta F_0$ is the midpoint.

$$
F(g)_{\text{Norm}} = S .e^{\frac{1}{2}(\frac{(\Delta F - \bar{\Delta F}}{\sigma})^2}
$$

where $S$ is a scaling factor that determines the maximum fitness, $\bar{\Delta F}$ is the average effect on phenotype, and $\sigma$ is the standard deviation of the effect on phenotype.

The phylogenetic mutation-selection model assumes exponential growth, continuous time, overlapping generations, and no density- or frequency-dependent selection. Under these assumptions the fitness corresponds to the 'normalized' growth rate, $N_e \times r_i$, where $N_e$ is population size and $r_i$ is the variant-specific instantaneous growth rate.

We also specified a step function that models the scenario of random walk. In this case, fitness of functional variants are the same:

$$
F(g)_{\text{step}} = \begin{cases}
0 & \text{when } \hat F < F_{\text{ref. genotype}} \\
1 & \text{when } \hat F \geq F_{\text{ref. genotype}}
\end{cases}
$$

For a population with two competing genotypes, $i$ and $j$, with respective fitnesses (i.e., growth rates) $N_er_i$ and $N_er_j$, we can define the selection coefficient $s_ {ij}$ as:

$$
N_e \times s_{ij} = (N_e \times r_j) - (N_e \times r_i) \\ 
= N_e \times (r_j - r_i)
$$

Let's have a look at the shape of the fitness functions. Note that AncSR2/SRE1 (green dashed line) is the reference genotype for clasifying a genotype as functional (green region):


```r
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

mFs <- seq(min(meanF_data$avg_meanF),max(meanF_data$avg_meanF),0.01)

# Logistic fitness function
log.ft.fn <- data.frame(mFs = mFs, fitness = fitness_logistic_ref(mFs,L,k,x_o,AncSR2_SRE_ref)) %>%
  ggplot(aes(x=mFs,y=fitness)) + ggtitle("Purifying + Drift + Directional") +
  annotate(geom = "rect", xmin =  min(meanF_data$avg_meanF), xmax = MIN_ACTIVE, ymin = 0, ymax = 3,
           fill = "#e37f74", alpha = 0.2) +
  annotate(geom = "rect", xmin =  MIN_ACTIVE, xmax = max(meanF_data$avg_meanF), ymin = 0, ymax = 3,
           fill = "#35a831", alpha = 0.2) +
  geom_line(linewidth=1.3) + theme_classic() + 
  xlab("meanF (phenotype)") + ylab("Fitness (Ne*r)") + 
  scale_fill_manual(values=c("#91c9a8","#c99591")) +
  geom_vline(xintercept = AncSR1_ERE_ref,col="#7710b3",linetype="dashed",linewidth=1.2) + # AncSR1/ERE reference phenotype and fitness
  geom_vline(xintercept = AncSR2_SRE_ref,col="#14b009",linetype="dashed",linewidth=1.2) # AncSR2/SRE reference phenotype and fitness


# Normal fitness function
norm.ft.fn <- data.frame(mFs = mFs, fitness = fitness_normal_ref(mFs,u,sd,scale,AncSR2_SRE_ref)) %>%
  ggplot(aes(x=mFs,y=fitness)) + ggtitle("Purifying + Stabilizing + Drift") +
  annotate(geom = "rect", xmin =  min(meanF_data$avg_meanF), xmax = MIN_ACTIVE, ymin = 0, ymax = 3,
           fill = "#e37f74", alpha = 0.2) +
  annotate(geom = "rect", xmin =  MIN_ACTIVE, xmax = max(meanF_data$avg_meanF), ymin = 0, ymax = 3,
           fill = "#35a831", alpha = 0.2) +
  geom_line(linewidth=1.3) + theme_classic() + 
  xlab("meanF (phenotype)") + ylab("Fitness (Ne*r)") + 
  geom_vline(xintercept = AncSR1_ERE_ref,col="#7710b3",linetype="dashed",linewidth=1.2) + # AncSR1/ERE reference phenotype and fitness
  geom_vline(xintercept = AncSR2_SRE_ref,col="#14b009",linetype="dashed",linewidth=1.2) # AncSR2/SRE reference phenotype and fitness

# Step fitness function
step.ft.fn <- data.frame(mFs = mFs, fitness = fitness_purifying(mFs,MIN_ACTIVE)) %>%
  ggplot(aes(x=mFs,y=fitness)) + ggtitle("Purifying + Drift") +
  annotate(geom = "rect", xmin =  min(meanF_data$avg_meanF), xmax = MIN_ACTIVE, ymin = 0, ymax = 3,
           fill = "#e37f74", alpha = 0.2) +
  annotate(geom = "rect", xmin =  MIN_ACTIVE, xmax = max(meanF_data$avg_meanF), ymin = 0, ymax = 3,
           fill = "#35a831", alpha = 0.2) +
  geom_line(linewidth=1.3) + theme_classic() + 
  xlab("meanF (phenotype)") + ylab("Fitness (Ne*r)") + 
  geom_vline(xintercept = AncSR1_ERE_ref,col="#7710b3",linetype="dashed",linewidth=1.2) + # AncSR1/ERE reference phenotype and fitness
  geom_vline(xintercept = AncSR2_SRE_ref,col="#14b009",linetype="dashed",linewidth=1.2) # AncSR2/SRE reference phenotype and fitness

# Visualize fitness functions
step.ft.fn + log.ft.fn #+ norm.ft.fn
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdfitness_functions-1.png" alt="plot of chunk fitness_functions"  />
<p class="caption">plot of chunk fitness_functions</p>
</div>

### Discrete Markov process on a GP map

We are interested in understanding how likely was each of the 16 DNA binding phenotypes to evolve under different biologically relevant evolutionary scenarios. Conceptually, the problem can be stated as follows: Under a given scenario, proteins traverse the sequence space by mutations, the GP map determines the effect of a mutation on phenotype, and the phenotype-fitness map defines its probability of fixation. Thus, the evolution of a DNA binding phenotype can be formulated in terms of its probability of occurring: What is the probability of evolving any given DNA binding phenotype from a specific starting genotype(s) after a number of substitutions, without losing function? 

The evolution of a DNA binding phenotype depends on two features of the genotype-phenotype map: 1) the number of accessible neighbors, at each mutation step, with the new function, and 2) the total number of genotypes that encode the function. Let $\{G\}$ be the complete set of genotypes in the network. Thus, the probability of evolving ${RE}_i$ ($P({RE}_i)$) is conditional on the starting set of genotypes $\{G_0\}$, the length of the mutational trajectory $S$, and the structure of the genotype-phenotye map.

We need to find the probabilty distribution over the states of $\{G\}$ after N steps of a Markov chain. That is, the probability that the substitution process will result in each genotype given the origin-fixation dynamics given by the transition matrix $P$ (where the probability of moving from genotype $i$ to genotype $j$ is $P(i,j)$). The row vector $\pi_{(S)}$ containing the distribution over the states of $\{G\}$ after $S$ substitutions is:

$$
\pi_{(S)} = \pi_{(0)} \times P^S
$$

for $S>0$, where $\pi_{(0)}$ is the vector of state frequencies at time step = 0 for every state in $\{G\}$. Since the transition matrices $P$ can be specified under different selection scenarios, this effectively captures the probabilities of the possible realizations of the process for each biological scenario.

We can use the DMS data to assign a function to each genotype of the vector $\pi_{(S)}$. Thus, the conditional probability of evolving any DNA binding funciton expressed in terms of the Markov chain becomes:

$$
P(RE_i|\pi_{(0)},S,P) = \frac{\sum_{j \in RE_i} \pi_{(S)j}}{\sum_{i=1}^{k} (\sum_{j \in RE_i} \pi_{(S)j})_i}
$$

where the numerator is the sum of the probailities of all the genotypes encoding the DNA binding function ${RE}_i$ after $S$ substitutions, and the denominator is a normalization constant with $k = 16$ such that $\sum_i P({RE}_i|\pi_{(0)},S,P) = 1$.

For $k$ different DNA binding phenotypes, we obtain a probability distribution of functional variation (PDFV), a multinomial probability distribution, around any set of starting genotypes $\{G_0\}$ that quantifies the likelihood that evolution will produce a particular phenotypic outcome given a set of conditions.  


## Evolutionary simulations using discrete markov chains


```r
## GLOBAL PARAMETERS FOR MARKOV CHAINS ##
REF_GENOTYPE = "EGKA" # EGKA (historical genotype) 
REF_GENOTYPE_COMPLEX = "EGKAGT" #EGKA/GT (historical protein-DNA complex)
N_CORES=detectCores()-1 # number of cores for parallel processing
cols <- RColorBrewer::brewer.pal(8, name="Dark2") # color palette
```


### Building genotype networks and probability transition matrices

For every evolutionary scenario, we will create a `P` matrix for each background, where each entry $P[i,j]$ corresponds to the probability $P(i,j)$ under each model. Initially, we will calculate selection coefficients based on the *maximum meanF* across the DNA elements bound by each variant.


```r
# Check whether the matrices have been loaded already (above), if not proceed to build the networks and matrices
if(!file.exists(file.path(".", "MutSel_matrices_final.RData"))) {
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

}
```

### Global properties of ancestral DBD genotype networks

Before delving into the neighborhoods around a focal node, let's take a look at the global properties of the GP map for each ancestral protein. Specifically, we're interested in three questions: 1) what is the phenotypic variation attainable through random mutation?, 2) what is the stationary distribution of the 16 DNA binding funtions?, and 3) are mutation and selection aligned in the GP map?

Answering the first question allows us to understand what is the variational propensity of each DBD background, that is, how likely is each phenotype to arise by random mutation. Answering the second question allows us to understand how the structure of the GP map determines the probability distribution of phenotypic variation at equilibrium (after suficiently long time has passed). Answering the third question can tell us whether the the mutation process is more likely to produce high fitness genotypes due the structure of the GP map. 

To answer these questions, first, we can use the transition probablity (`P`) matrices to determine the probabilities of every mutational trajectory in the GP map and compute the relative probability of evolving each DNA binding funtion. Second, we can use the `P` matrices to compute the stationary distribution of all genotypes (an their phenotypes) in the network. Third, we will use the mutation accessibility matrix (`M`) to ask whether the stationary distributoin of states generated from this matrix correlates with the fitness of those RH variants.

The transition matrices of Markov chains describe the probabilities of change between every state of the system at each step - the probability distribution of the states evolves over time ($X_{t}$). The stationary distribution of a Markov chain describes the distribution of states after suficiently long time has passed ($t\rightarrow\infty$) such that $X_{t}$ does noth change any longer (i.e., the system has reached an equilibrium).

Let's extract the square matrices from the main components of the genotype networks (the biggest connected sub-network) and compute the stationary distributions of the genotypes.


```r
# TRANSITION MATRICES (mut + selection) #
# Drift
P_drift_sr1_ntwrk <- extract_main_ntwrk(net_sr1,M_drift_sr1) # create square matrix of main component of genotype network
```

```
## This graph was created by an old(er) igraph version.
##   Call upgrade_graph() on it to use with the current igraph version
##   For now we convert it on the fly...
```

```r
P_drift_sr1_ntwrk_statdist <- stationary_dist(P_drift_sr1_ntwrk) # compute stationary distribution of genotypes in the main network component

P_drift_sr2_ntwrk <- extract_main_ntwrk(net_sr2,M_drift_sr2) 
```

```
## This graph was created by an old(er) igraph version.
##   Call upgrade_graph() on it to use with the current igraph version
##   For now we convert it on the fly...
```

```r
P_drift_sr2_ntwrk_statdist <- stationary_dist(P_drift_sr2_ntwrk)

# Directional selection
P_dir_sr1_ntwrk <- extract_main_ntwrk(net_sr1,M_dir_sr1) 
P_dir_sr1_ntwrk_statdist <- stationary_dist(P_dir_sr1_ntwrk) 

P_dir_sr2_ntwrk <- extract_main_ntwrk(net_sr2,M_dir_sr2) 
P_dir_sr2_ntwrk_statdist <- stationary_dist(P_dir_sr2_ntwrk)

# MUTATION MATRICES #
M_mat_sr1_ntwrk <- extract_main_ntwrk(net_sr1,adj_mat_count_sr1)
M_mat_sr1_ntwrk_statdist <- stationary_dist(M_mat_sr1_ntwrk)

M_mat_sr2_ntwrk <- extract_main_ntwrk(net_sr2,adj_mat_count_sr2)
M_mat_sr2_ntwrk_statdist <- stationary_dist(M_mat_sr2_ntwrk)
```

We can explore how each molecular evolutionary scenario interacts with the structure of the GP map. Specifically, we can compare the `P` matrices to the `M` matrices, the transition matrices estimated only from the mutation rates. The `M` matrices capture how the mutation process alone can explore the network, accounting for the structure of the genetic code and the structure of the genotype network. The `P` matrices include the effect of selection + drift, so they can give us an idea of the interaction of the structure of the GP map and the evolutionary processes. Overall, we can begin to understand the effect of several GP mappings: 1) codon-to-amino acid, 2) amino acid-to-function, and 3) function-to-fitness.

Because the scenario of purifying selection + drift treats all functional variants as equally fit, this scenario is a purely mutation-driven process and we should expect the stationary distributions of genotypes from the `M` and `P` matrices to be the same. In contrast, directional selection can introduce deviations because of fitness differences between variants.


```r
print("Statonary distribution of genotypes:")
```

```
## [1] "Statonary distribution of genotypes:"
```

```r
# AncSR1
print(paste("Correlation between M and P matrices for AncSR1 bg under drift:",cor(M_mat_sr1_ntwrk_statdist,P_drift_sr1_ntwrk_statdist)))
```

```
## [1] "Correlation between M and P matrices for AncSR1 bg under drift: 1"
```

```r
print(paste("Correlation between M and P matrices for AncSR1 bg under directional sln:",
            cor(M_mat_sr1_ntwrk_statdist,P_dir_sr1_ntwrk_statdist)))
```

```
## [1] "Correlation between M and P matrices for AncSR1 bg under directional sln: 0.717559145086064"
```

```r
a <- data.frame(X=M_mat_sr1_ntwrk_statdist,Y=P_drift_sr1_ntwrk_statdist) %>% ggplot(aes(x=X,y=Y)) + 
  geom_point(fill="black") + geom_abline(slope = 1,intercept = 0,col="red") + theme_classic() + xlab("Genotype prob. (M matrix)") +
  ylab("Genotype prob. (P matrix)") + ggtitle("Drift")

b <- data.frame(X=M_mat_sr1_ntwrk_statdist,Y=P_dir_sr1_ntwrk_statdist) %>% ggplot(aes(x=X,y=Y)) + 
  geom_point(fill="black") + geom_abline(slope = 1,intercept = 0,col="red") + theme_classic() + xlab("Genotype prob. (M matrix)") +
  ylab("Genotype prob. (P matrix)") + ggtitle("Directional sln.")

a+b
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdm_p_matrices-1.png" alt="plot of chunk m_p_matrices"  />
<p class="caption">plot of chunk m_p_matrices</p>
</div>

```r
# AncSR2
print(paste("Correlation between M and P matrices for AncSR2 bg under drift:",cor(M_mat_sr2_ntwrk_statdist,P_drift_sr2_ntwrk_statdist)))
```

```
## [1] "Correlation between M and P matrices for AncSR2 bg under drift: 1"
```

```r
print(paste("Correlation between M and P matrices for AncSR2 bg under directional sln:",
            cor(M_mat_sr2_ntwrk_statdist,P_dir_sr2_ntwrk_statdist)))
```

```
## [1] "Correlation between M and P matrices for AncSR2 bg under directional sln: 0.846831977209567"
```

```r
c <- data.frame(X=M_mat_sr2_ntwrk_statdist,Y=P_drift_sr2_ntwrk_statdist) %>% ggplot(aes(x=X,y=Y)) + 
  geom_point(fill="black") + geom_abline(slope = 1,intercept = 0,col="red") + theme_classic() + xlab("Genotype prob. (M matrix)") +
  ylab("Genotype prob. (P matrix)") + ggtitle("Drift")

d <- data.frame(X=M_mat_sr2_ntwrk_statdist,Y=P_dir_sr2_ntwrk_statdist) %>% ggplot(aes(x=X,y=Y)) + 
  geom_point(fill="black") + geom_abline(slope = 1,intercept = 0,col="red") + theme_classic() + xlab("Genotype prob. (M matrix)") +
  ylab("Genotype prob. (P matrix)") + ggtitle("Directional sln.")

c+d
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdm_p_matrices-2.png" alt="plot of chunk m_p_matrices"  />
<p class="caption">plot of chunk m_p_matrices</p>
</div>

As expected, we see that the stationary distribution of genotypes from `M` and `P` matrices for the random walk scenario are the same. In contrast, directional selection reduces the correlation by ~25-30%.

---

Let's now explore the questions in more detail. For this, we are interested in characterizing the **Probability Distribution of Functional Variation (PDFV)**, that is, the relative probability of accessing any of the 16 DNA binding elements by single mutations. 

For the first question, we are interested in the PDFV that arises through random mutation, the *variational PDFV*. That is, the likelihood of a DNA binding phentoype is equal to the proportion of genotypes encoding each function. This distribution is only computed from the functional genotypes - it incorporates the effects of strong truncation selection - and therefore represents the probability that a random mutation to a functional genotype will produce a given phenotype. By including the amino acid-to-function mapping, this can give us an idea of the variational properties of the system.


```r
# Variational propensities
var.prop_AncSR1_df <- get_PDFV_v2(type="network",Bg="AncSR1",model="Binding")
var.prop_AncSR1_df2 <- get_PDFV_v2(type="network",Bg="AncSR1",model="Specific Binding",specific = TRUE) %>% 
  filter(RE != "Promiscuous") %>% 
  mutate(Norm_F_prob = Norm_F_prob /sum(Norm_F_prob))
var.prop_AncSR2_df <- get_PDFV_v2(type="network",Bg="AncSR2",model="Binding")
var.prop_AncSR2_df2 <- get_PDFV_v2(type="network",Bg="AncSR2",model="Specific Binding",specific = TRUE) %>% 
  filter(RE != "Promiscuous") %>% 
  mutate(Norm_F_prob = Norm_F_prob /sum(Norm_F_prob))

# Plot variational propensities for each DBD
a <- circular_PDFV_v2(list(var.prop_AncSR1_df,var.prop_AncSR1_df2),cols = c("#1D1D35","#00ACBF"),fill=F,title = "AncSR1\nVariational propensity")
b <- circular_PDFV_v2(list(var.prop_AncSR2_df,var.prop_AncSR2_df2),cols = c("#1D1D35","#00ACBF"),fill=F,legend = F,title = "AncSR2\nVariational propensity")

a + b
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdvariational_prop-1.png" alt="plot of chunk variational_prop"  />
<p class="caption">plot of chunk variational_prop</p>
</div>

```r
# Compare the distribution of AncSR2 with the “expected” distribution of AncSR1
exp_sr1 <- phenotypes_tbl %>% filter(specific == "YES" & bg == "AncSR1") %>% ungroup() %>%
  reframe(RE =  REs[[1]],
          count = table(factor(specificity,levels=REs[[1]])),
          prop = count/sum(count))

obs_sr2 <- phenotypes_tbl %>% filter(specific == "YES" & bg == "AncSR2") %>% ungroup() %>%
  reframe(RE =  REs[[1]],
          count = table(factor(specificity,levels=REs[[1]])))

d <- inner_join(exp_sr1,obs_sr2,by="RE") %>% 
  mutate(exp_count = prop * sum(count.y),
         mod_prop = ifelse(prop==0,.Machine$double.xmin,prop), # replace zeroes with .Machine$double.xmin for XNomial test
         mod_exp_count = mod_prop * sum(count.y))

print("P-values from Multinomial Exact test and Chi2 test:")
```

```
## [1] "P-values from Multinomial Exact test and Chi2 test:"
```

```r
xmonte(obs=as.vector(d$count.y),expr = as.vector(d$mod_exp_count))
```

```
## 
## P value (LLR) = 0 ± 0
```

```r
chisq.test(x=as.vector(d$count.y),p=as.vector(d$mod_prop))
```

```
## 
## 	Chi-squared test for given probabilities
## 
## data:  as.vector(d$count.y)
## X-squared = 1.4835e+308, df = 15, p-value < 2.2e-16
```

These figures show the PDFV that arises from random mutation in each DBD background. There are three key results here:

- The AncSR1 background only produces 12/16 phenotypes upon random mutation, while AncSR2 can produce all 16 phenotypes.
- From the 12 phenotypes produced on AncSR1, nine of them are bound specifically, while the remaining three are only bound by promiscuous genotypes. From the 16 phentoypes produced on AncSR2, 14 of them are bound specifically and the remaining two are bound by promiscuous genotypes.
- The PDFV for each DBD background are aligned with the wild-type phenotypes. The most frequently produced phenotype is the wild type-phenotype of each ancestral protein, ERE and SRE1, respectively.
- The shape of the PDFV changed along the phylogenetic trajectory. The variational properties of AncSR1 and AncSR2 are significantly different (Goodness of fit *P-values* << 0.01).

---

Now let's explore the second question. Since we have already computed the stationary distributions of genotypes for each `P` matrix, we can use these distributions to compute the PDFV at equilibrium.


```r
##############################
# STATIONARY PDFVs
# AncSR1
stat_pdfv_drift_sr1 <- get_PDFV_v2(P_drift_sr1_ntwrk_statdist,type="simulated mc",Bg="AncSR1",model="Drift")
stat_pdfv_dir_sr1 <- get_PDFV_v2(P_dir_sr1_ntwrk_statdist,type="simulated mc",Bg="AncSR1",model="Dir. sln.")

p1 <- circular_PDFV_v2(list(var.prop_AncSR1_df,stat_pdfv_drift_sr1,stat_pdfv_dir_sr1),cols = c(cols[2:3],"gray60"),title = "AncSR1: Stationary PDFVs",fill = F)

# AncSR2
stat_pdfv_drift_sr2 <- get_PDFV_v2(P_drift_sr2_ntwrk_statdist,type="simulated mc",Bg="AncSR2",model="Drift")
stat_pdfv_dir_sr2 <- get_PDFV_v2(P_dir_sr2_ntwrk_statdist,type="simulated mc",Bg="AncSR2",model="Dir. sln.")

p2 <- circular_PDFV_v2(list(var.prop_AncSR2_df,stat_pdfv_drift_sr2,stat_pdfv_dir_sr2),cols = c(cols[2:3],"gray60"),title = "AncSR2: Stationary PDFVs",legend = F,fill = F)

p1+p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-RmdGPmapGlobalProperties_b-1.png" alt="plot of chunk GPmapGlobalProperties_b"  />
<p class="caption">plot of chunk GPmapGlobalProperties_b</p>
</div>

From these figures we can see that the stationary PDFVs in AncSR2 approximates the proportional expectation (variational propensity). In contrast, the stationary PDFVs in AncSR1 deviate from the expected variational propensity. This difference is likely due to the GP map of AncSR2 being more connected, and therefore, evolutionary processes can navigate more efficiently the network and explore the phenotypic space.

In fact, we can explore how the network architecture affects the PDFV at equilibrium. The genotype networks, and the PDFVs inferred from them, incorporate the structure of the genetic code. There are two ways in which we can remove the genetic code to understand its effects, First, we can build a network based on hamming distances, where two variants are connected if they differ by a single amino acid change regardless if those amino acids are accessible through nucleotide mutations. Second, we can build a maximally connected network, where all nodes are connected to each other. We expect the maximally connected network to have a PDFV exactly equal to the variational PDFV because evolution can traverse this network with maximal efficiency.


```r
##############################
# STATIONARY PDFVs FOR SIMULATED GP MAPS
# AncSR1
# full graph
adj_mat_fullgraph_sr1 <- simulate_GPmap(net_sr1,type=3,which="mat")
adj_mat_fullgraph_sr1 <- t(apply(adj_mat_fullgraph_sr1, 1, function(x) x / sum(x)))
adj_mat_fullgraph_sr1 <- replace(adj_mat_fullgraph_sr1,is.nan(adj_mat_fullgraph_sr1),0)
adj_mat_fullgraph_sr1 <- as(adj_mat_fullgraph_sr1, "sparseMatrix")

# hamming graph
hamming_sr1 <- simulate_GPmap(net_sr1,type=1,cores=N_CORES,which="both")
net_hamming_sr1 <- hamming_sr1[[1]]
adj_mat_hamming_sr1 <- hamming_sr1[[2]]
adj_mat_hamming_sr1 <- t(apply(adj_mat_hamming_sr1, 1, function(x) x / sum(x)))
adj_mat_hamming_sr1 <- replace(adj_mat_hamming_sr1,is.nan(adj_mat_hamming_sr1),0)
adj_mat_hamming_sr1 <- as(adj_mat_hamming_sr1, "sparseMatrix")
adj_mat_hamming_sr1 <- extract_main_ntwrk(graph = net_hamming_sr1,tr_mat = adj_mat_hamming_sr1)

adj_mat_fullgraph_sr1_statdist <- stationary_dist(adj_mat_fullgraph_sr1)
adj_mat_hamming_sr1_statdist <- stationary_dist(adj_mat_hamming_sr1)

fullgraph_stat_pdfv_sr1 <- get_PDFV_v2(adj_mat_fullgraph_sr1_statdist,type="simulated mc",Bg="AncSR1",model="Full graph")
hamming_stat_pdfv_sr1 <- get_PDFV_v2(adj_mat_hamming_sr1_statdist,type="simulated mc",Bg="AncSR1",model="Hamming")
genetic_code_pdfv_sr1 <- get_PDFV_v2(M_mat_sr1_ntwrk_statdist,type="simulated mc",Bg="AncSR1",model="Genetic code")

# compare distributions
n_functional_prot_genotypes_sr1 <- 164
n_functional_complexes_sr1 <- 189

d <- rbind(var.prop_AncSR1_df,genetic_code_pdfv_sr1,fullgraph_stat_pdfv_sr1,hamming_stat_pdfv_sr1) %>% 
  mutate(Norm_F_prob = ifelse(Norm_F_prob==0,.Machine$double.xmin,Norm_F_prob), # replace zeroes with .Machine$double.xmin for XNomial test)
         Exp_n = Norm_F_prob * n_functional_complexes_sr1) %>% select(RE,model,Exp_n) %>%
  pivot_wider(names_from="model",values_from="Exp_n") %>% inner_join(.,var.prop_AncSR1_df,by="RE") %>%
  mutate(Norm_F_prob = ifelse(Norm_F_prob==0,.Machine$double.xmin,Norm_F_prob))

print("AncSR1: Chi2 test with variational PDFV as reference:")
```

```
## [1] "AncSR1: Chi2 test with variational PDFV as reference:"
```

```r
chisq.test(x=d$`Genetic code`,p=d$Norm_F_prob)
```

```
## 
## 	Chi-squared test for given probabilities
## 
## data:  d$`Genetic code`
## X-squared = 35.872, df = 15, p-value = 0.001844
```

```r
chisq.test(x=d$`Full graph`,p=d$Norm_F_prob)
```

```
## 
## 	Chi-squared test for given probabilities
## 
## data:  d$`Full graph`
## X-squared = 1.7254e-26, df = 15, p-value = 1
```

```r
chisq.test(x=d$Hamming,p=d$Norm_F_prob)
```

```
## 
## 	Chi-squared test for given probabilities
## 
## data:  d$Hamming
## X-squared = 12.215, df = 15, p-value = 0.6627
```

```r
###################

# AncSR2
# full graph
adj_mat_fullgraph_sr2 <- simulate_GPmap(net_sr2,type=3,which="mat")
adj_mat_fullgraph_sr2 <- t(apply(adj_mat_fullgraph_sr2, 1, function(x) x / sum(x)))
adj_mat_fullgraph_sr2 <- replace(adj_mat_fullgraph_sr2,is.nan(adj_mat_fullgraph_sr2),0)
adj_mat_fullgraph_sr2 <- as(adj_mat_fullgraph_sr2, "sparseMatrix")

# hamming graph
hamming_sr2 <- simulate_GPmap(net_sr2,type=1,cores=N_CORES,which="both")
net_hamming_sr2 <- hamming_sr2[[1]]
adj_mat_hamming_sr2 <- hamming_sr2[[2]]
adj_mat_hamming_sr2 <- t(apply(adj_mat_hamming_sr2, 1, function(x) x / sum(x)))
adj_mat_hamming_sr2 <- replace(adj_mat_hamming_sr2,is.nan(adj_mat_hamming_sr2),0)
adj_mat_hamming_sr2 <- as(adj_mat_hamming_sr2, "sparseMatrix")
adj_mat_hamming_sr2 <- extract_main_ntwrk(graph = net_hamming_sr2,tr_mat = adj_mat_hamming_sr2)

adj_mat_fullgraph_sr2_statdist <- stationary_dist(adj_mat_fullgraph_sr2)
adj_mat_hamming_sr2_statdist <- stationary_dist(adj_mat_hamming_sr2)

fullgraph_stat_pdfv_sr2 <- get_PDFV_v2(adj_mat_fullgraph_sr2_statdist,type="simulated mc",Bg="AncSR2",model="Full graph")
hamming_stat_pdfv_sr2 <- get_PDFV_v2(adj_mat_hamming_sr2_statdist,type="simulated mc",Bg="AncSR2",model="Hamming")
genetic_code_pdfv_sr2 <- get_PDFV_v2(M_mat_sr2_ntwrk_statdist,type="simulated mc",Bg="AncSR2",model="Genetic code")

# compare distributions
n_functional_prot_genotypes_sr2 <- 2036
n_functional_complexes_sr2 <- 3983

d2 <- rbind(var.prop_AncSR2_df,genetic_code_pdfv_sr2,fullgraph_stat_pdfv_sr2,hamming_stat_pdfv_sr2) %>% 
  mutate(Norm_F_prob = ifelse(Norm_F_prob==0,.Machine$double.xmin,Norm_F_prob), # replace zeroes with .Machine$double.xmin for XNomial test)
         Exp_n = Norm_F_prob * n_functional_complexes_sr2) %>% select(RE,model,Exp_n) %>%
  pivot_wider(names_from="model",values_from="Exp_n") %>% inner_join(.,var.prop_AncSR2_df,by="RE") %>%
  mutate(Norm_F_prob = ifelse(Norm_F_prob==0,.Machine$double.xmin,Norm_F_prob))

print("AncSR2: Chi2 test with variational PDFV as reference:")
```

```
## [1] "AncSR2: Chi2 test with variational PDFV as reference:"
```

```r
chisq.test(x=d2$`Genetic code`,p=d2$Norm_F_prob)
```

```
## 
## 	Chi-squared test for given probabilities
## 
## data:  d2$`Genetic code`
## X-squared = 152.64, df = 15, p-value < 2.2e-16
```

```r
chisq.test(x=d2$`Full graph`,p=d2$Norm_F_prob)
```

```
## 
## 	Chi-squared test for given probabilities
## 
## data:  d2$`Full graph`
## X-squared = 1.3735e-24, df = 15, p-value = 1
```

```r
chisq.test(x=d2$Hamming,p=d2$Norm_F_prob)
```

```
## 
## 	Chi-squared test for given probabilities
## 
## data:  d2$Hamming
## X-squared = 19.38, df = 15, p-value = 0.197
```

```r
################

# plots
p1 <- circular_PDFV_v2(list(var.prop_AncSR1_df,genetic_code_pdfv_sr1,hamming_stat_pdfv_sr1),cols = cols[4:6],title = "AncSR1: Stationary PDFVs",fill = F)

p2 <- rbind(data.frame(AA_var=names(degree(net_sr1)),degree=degree(net_sr1)/2,network="Genetic code"),
            data.frame(AA_var=names(degree(net_hamming_sr1)),degree=degree(net_hamming_sr1)/2,network="Hamming")) %>% 
  ggplot(aes(x=degree,fill=network)) + geom_histogram(bins=25,color="black") + theme_classic() + ggtitle("Degree distribution")

p3 <- circular_PDFV_v2(list(var.prop_AncSR2_df,genetic_code_pdfv_sr2,hamming_stat_pdfv_sr2),cols = cols[4:6],title = "AncSR2: Stationary PDFVs",fill = F)

p4 <- rbind(data.frame(AA_var=names(degree(net_sr2)),degree=degree(net_sr2)/2,network="Genetic code"),
            data.frame(AA_var=names(degree(net_hamming_sr2)),degree=degree(net_hamming_sr2)/2,network="Hamming")) %>% 
  ggplot(aes(x=degree,fill=network)) + geom_histogram(bins=25,color="black") + theme_classic() + ggtitle("Degree distribution")

p1 + p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-RmdGPmapGlobalProperties_c-1.png" alt="plot of chunk GPmapGlobalProperties_c"  />
<p class="caption">plot of chunk GPmapGlobalProperties_c</p>
</div>

```r
p3 + p4
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-RmdGPmapGlobalProperties_c-2.png" alt="plot of chunk GPmapGlobalProperties_c"  />
<p class="caption">plot of chunk GPmapGlobalProperties_c</p>
</div>

Wee see that the PDFV at equilibrium for the full graph and the hamming network are not significantly different from the variational PDFV. In contrast, the PDFV from the network built accounting for the genetic code differs from the expected PDFV. Since all of these networks have the same amino acid-to-function map, this indicates that the mapping from codon to amino acid introduces an additional layer of complexity to the GP map. For instance, the average number of edges per node (degree) for the networks in AncSR1 are: 81 for the full network; 3.92 for the hamming network; and 1.91 for the genetic code network. In AncSR2 the averages are: 1017 for the full network; 10.86 for the hamming network; and 5.09 for the genetic code network. In both cases, introducing the genetic code reduces by half the average node connectivity relative to the hamming network, and this is enough to devaite the equilibrium dynamics from the expectectation.

---

Finally, we can use the `M` matrices and their stationary distributions to ask whether the structure of the network is 'aligned' with the selection surface. That is, whether the functional genotypes that are more likely to evolve through random mutation have also higher fitness.


```r
##############################
# GLOBAL ALIGNMENT OF MUTATION AND SELECTION

# Compute fitness of genotypes variants of main sub-networks and merge with stationary distributions from M-matrices
Aln_matrix_sr1 <- phenotypes_tbl %>% filter(bg == "AncSR1") %>% filter(AA_var %in% names(M_mat_sr1_ntwrk_statdist)) %>% 
  select(AA_var,max_meanF_bREs) %>% mutate(FL = fitness_logistic(max_meanF_bREs,L,k,x_o)) %>%
  inner_join(.,data.frame(AA_var=names(M_mat_sr1_ntwrk_statdist),m_freq=M_mat_sr1_ntwrk_statdist),by="AA_var") 

cor.test(Aln_matrix_sr1$m_freq,Aln_matrix_sr1$FL)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  Aln_matrix_sr1$m_freq and Aln_matrix_sr1$FL
## t = 5.423, df = 241, p-value = 1.423e-07
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.2127690 0.4374719
## sample estimates:
##       cor 
## 0.3297832
```

```r
p1 <- Aln_matrix_sr1 %>% pivot_longer(cols=3,names_to="fitness_function",values_to="fitness") %>%
  ggplot(aes(x=m_freq,y=fitness)) + geom_point(shape=21,color="black",fill="gray") +
  geom_smooth(method = "loess",color="black") + theme_classic() + scale_fill_manual(values=c("#b00c59","#219e92")) +
  labs(title = "AncSR1",x="M-matrix stationary distribution of RH genotype",y="Fitness")

Aln_matrix_sr2 <- phenotypes_tbl %>% filter(bg == "AncSR2") %>% filter(AA_var %in% names(M_mat_sr2_ntwrk_statdist)) %>% 
  select(AA_var,max_meanF_bREs) %>% mutate(FL = fitness_logistic(max_meanF_bREs,L,k,x_o)) %>%
  inner_join(.,data.frame(AA_var=names(M_mat_sr2_ntwrk_statdist),m_freq=M_mat_sr2_ntwrk_statdist),by="AA_var")

cor.test(Aln_matrix_sr2$m_freq,Aln_matrix_sr2$FL)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  Aln_matrix_sr2$m_freq and Aln_matrix_sr2$FL
## t = 14.089, df = 2367, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.2405906 0.3149121
## sample estimates:
##       cor 
## 0.2781676
```

```r
p2 <- Aln_matrix_sr2 %>% pivot_longer(cols=3,names_to="fitness_function",values_to="fitness") %>%
  ggplot(aes(x=m_freq,y=fitness)) + geom_point(shape=21,color="black",fill="gray") +
  geom_smooth(method = "loess",color="black") + theme_classic() + scale_fill_manual(values=c("#b00c59","#219e92")) +
  labs(title = "AncSR2",x="M-matrix stationary distribution of RH genotype",y="Fitness")

p1 + p2 
```

```
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-RmdGPmapGlobalProperties_d-1.png" alt="plot of chunk GPmapGlobalProperties_d"  />
<p class="caption">plot of chunk GPmapGlobalProperties_d</p>
</div>

The plots show that mutation and selection are not globally aligned in the GPmaps, that is, the structure of the GPmap do not tend to produce higher fitness genotypes at equilibrium.

### Discrete Markov chains

### Evolution around EGKA on AncSR1 and AncSR2 networks

Let's now focus on the historical trajectory. Based on ASR, we know that the RH genotype on AncSR1 was 'EGKA' - an ERE-specific genotype - and the genotype on AncSR2 was 'GSKV' - an SRE-specific genotype. Although the hamming distance in the RH between AncSR1:EGKA and AncSR2:GSKV is 3 amino acid substitutions, this doesn't necessarily mean that the historical trajectory involved exactly 3 mutation steps. The branch length between AncSR1 and AncSR2 is ~1.5 substitutions/site which suggests that multiple substitutions per site may have happened during the phylogenetic trajectory (the branch length after the cephalochordate kSRs split is 0.97 subs/site which still is relatively long). 

We can try to infer the **length** of the historical trajectory by using the starting and ending points of the trajectory, the genotype networks and the `P` matrices. The genotype network can help us to identify the shortest functional mutational path between the historical RH genotypes. With the `P` matrices we can simulate markov chains of increasing lengths and infer the probability of the derived genotype. 

Lastly, we can ask two things about the historical mutational trajectory: 1) Was the derived RH genotype - GSKV - accessible on AncSR1 and AncSR2? and 2) Why did GSKV evolve? Answering these questions can give some insight into the effects of chance and contingency.
 

```r
# Shortest mutational trajectory between EGKA and GSKV in each genotype network
mut_tr_sr1 <- all_shortest_paths(net_sr1,from = which(names(V(net_sr1))=="EGKA"),
                               to = which(names(V(net_sr1))=="GSKV"))$res

mut_tr_sr2 <- all_shortest_paths(net_sr2,from = which(names(V(net_sr2))=="EGKA"),
                               to = which(names(V(net_sr2))=="GSKV"))$res

print("Mutational trajectories between EGKA and GSKV on the AncSR1 network:")
```

```
## [1] "Mutational trajectories between EGKA and GSKV on the AncSR1 network:"
```

```r
mut_tr_sr1
```

```
## list()
```

```r
print("Mutational trajectories between EGKA and GSKV on the AncSR2 network:")
```

```
## [1] "Mutational trajectories between EGKA and GSKV on the AncSR2 network:"
```

```r
mut_tr_sr2
```

```
## [[1]]
## + 4/2390 vertices, named, from cee954f:
## [1] EGKA GGKA GGKV GSKV
```

We can see that:

- GSKV is not accessible from EGKA on the AncSR1 background.
- GSKV becomes accessible from EGKA on the AncSR2 background. The single trajectory consists of 3 mutation steps, and this trajectory recapitulates the one inferred on Anderson, McKeown & Thornton (2015). 

Based on the AncSR2 genotype network, a length of 3 mutation steps is the shortest path length between RH genotypes. We can now use the `P` matrices to track the probability of evolving GSKV from EGKA under different markov chain lengths. We will use a relative likelihood approach to determine the most likely estimate of $S$, the length of the trajectory. The relative likelihood compares the relative plausibilities of different candidate models or of different values of a parameter of a single model. In our case, the model is the same (i.e., $\pi_{(S)} = \pi_{(0)} \times P^S$) but we want to infer the best estimate of $S$ for $P(GSKV)$. 

If the maximum likelihood estimate for $S$ is $\hat{S}$, the relative plausibilities of other $S$ values may be found by comparing the likelihoods of those other values with the likelihood of $\hat{S}$. The relative likelihood of $S$ is defined to be

$$
{\displaystyle {\frac {~{\mathcal {L}}(S \mid x)~}{~{\mathcal {L}}({\hat {S }}\mid x)~}}}
$$

Further, by defining a significance level $\alpha$, we can establish a likelihood interval that captures the set of $S$ estimates that are not significantly different from $\hat{S}$. The likelihood interval ($LI$) is defined:

$$
LI = {\displaystyle \left\{S :{\frac {{\mathcal {L}}(S \mid x)}{{\mathcal {L}}({\hat {S \,}}\mid x)}}\geq {1 - \alpha}\right\}.}
$$

To find the value of $\hat{S}$, we will run a markov chain for multiple values of $S$ to get an estimate of the of the shape of the likelihood function. 


```r
# Track the change in probability of evolving GSKV on AncSR2 --> Identify the mutation step at which P(GSKV) is maximized
# Random walk
p_gskv <- c()
mut_step <- 15
for(i in 1:mut_step){
  tmp_mc <- simulate_markov_chain(REF_GENOTYPE,P_drift_sr2_ntwrk,n_steps = i)
  p_gskv <- c(p_gskv,tmp_mc[names(tmp_mc)=="GSKV"])
}
stat_p_gskv <- data.frame(mut_step = mut_step+1, p_gskv = P_drift_sr2_ntwrk_statdist[names(P_drift_sr2_ntwrk_statdist)=="GSKV"])

# Relative likelihood
d <- data.frame(mut_step = seq(1,mut_step,1), p_gskv = p_gskv) %>% rbind(.,stat_p_gskv) %>%
  mutate(LR = p_gskv/max(p_gskv), # relative likelihood
         LI = ifelse(LR >= 0.95, TRUE,FALSE)) # Likelihood interval (1-a, for a = 0.05)

p1 <- d %>%
  ggplot(aes(x=mut_step,y=p_gskv)) + geom_point(col="black") + geom_line() + 
  xlab("Mutation step") + ylab("P(GSKV)") + geom_vline(data = d %>% filter(LI==T),aes(xintercept=mut_step),col="red",linetype="dashed") +
  ggtitle("AncSR2: Random walk") + 
  scale_x_continuous(breaks = seq(1,mut_step+1,1),labels = c(seq(1,mut_step,1),Inf)) +
  theme_classic()

# Directional selection
p_gskv <- c()
mut_step <- 15
for(i in 1:mut_step){
  tmp_mc <- simulate_markov_chain(REF_GENOTYPE,P_dir_sr2_ntwrk,n_steps = i)
  p_gskv <- c(p_gskv,tmp_mc[names(tmp_mc)=="GSKV"])
}
stat_p_gskv <- data.frame(mut_step = mut_step+1, p_gskv = P_dir_sr2_ntwrk_statdist[names(P_dir_sr2_ntwrk_statdist)=="GSKV"])

# Relative likelihood
d2 <- data.frame(mut_step = seq(1,mut_step,1), p_gskv = p_gskv) %>% rbind(.,stat_p_gskv) %>%
  mutate(LR = p_gskv/max(p_gskv), # relative likelihood
         LI = ifelse(LR >= 0.95, TRUE,FALSE)) # Likelihood interval (1-a, for a = 0.05)

p2 <- d2 %>%
  ggplot(aes(x=mut_step,y=p_gskv)) + geom_point(col="black") + geom_line() + 
  xlab("Mutation step") + ylab("P(GSKV)") + geom_vline(data = d2 %>% filter(LI==T),aes(xintercept=mut_step),col="red",linetype="dashed") +
  ggtitle("AncSR2: Directional sln") + 
  scale_x_continuous(breaks = seq(1,mut_step+1,1),labels = c(seq(1,mut_step,1),Inf)) +
  theme_classic()

p1 + p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdmut_trajectory2-1.png" alt="plot of chunk mut_trajectory2"  />
<p class="caption">plot of chunk mut_trajectory2</p>
</div>

The figures shows that the probability of evolving GSKV is in fact maximized when the markov chain is ran for 3 mutation steps (i.e., $\hat{S} = 3$), and it steadly declines over time. Further, the $LI_{\alpha = 0.05}$ (red lines) is $3 \leq S \leq 5$. Based on this, we can conclude that the most likely mutational trajectory between the ancestral and derived RH genotypes involved 3 amino acid substitutions.


```r
# Set 'PATH_LENGTH' as a global variable for downstream analyses
PATH_LENGTH = 3 # path length (neighborhood size) to find mutational trajectoties in the protein network
```

We now know that `P(GSKV)` is maximized for a mutational path length of 3, however, we still don't know whether GSKV was the *most* likely genotype to evolve amogst all the accessible genotypes. For this, we will compute the probabilities of every accessible genotype in the mutational neighborhood around AncSR2:EGKA (for a path length of 3 mutation steps).


```r
# Plot probability of every genotype accessible after 3 amino acid substitutions (from markov chains) --> only for AncSR2
# to compare probability of GSKV
# Random walk
mc_Drift_ref_genotype_sr2 <- simulate_markov_chain(REF_GENOTYPE,tr_mat = P_drift_sr2_ntwrk,n_steps = PATH_LENGTH)

df_mut_traj_probs_sr2_drift <- data.frame(AA_var=names(mc_Drift_ref_genotype_sr2),prob=mc_Drift_ref_genotype_sr2,
                                          scenario="Random walk") %>% mutate(rank = dense_rank(prob))

p1 <- ggplot(df_mut_traj_probs_sr2_drift,aes(x=rank,y=prob)) + geom_point(col="black") + 
  geom_point(data=df_mut_traj_probs_sr2_drift %>% filter(AA_var=="GSKV"),color="red",size=3) +
  geom_text(data = df_mut_traj_probs_sr2_drift %>% filter(AA_var=="GSKV"), 
            aes(x = rank, y = prob + 0.002, label = "GSKV"),angle=45,hjust=0) +
  xlab("Genotype rank") + ylab("Probability") + ggtitle("AncSR2: Random walk") +
  theme_classic()

# Directional selection
mc_DirSln_ref_genotype_sr2 <- simulate_markov_chain(REF_GENOTYPE,tr_mat = P_dir_sr2_ntwrk,n_steps = PATH_LENGTH)

df_mut_traj_probs_sr2_dir <- data.frame(AA_var=names(mc_DirSln_ref_genotype_sr2),prob=mc_DirSln_ref_genotype_sr2,
                                          scenario="Dir Sln") %>% mutate(rank = dense_rank(prob))

p2 <- ggplot(df_mut_traj_probs_sr2_dir,aes(x=rank,y=prob)) + geom_point(col="black") + 
  geom_point(data=df_mut_traj_probs_sr2_dir %>% filter(AA_var=="GSKV"),color="red",size=3) +
  geom_text(data = df_mut_traj_probs_sr2_dir %>% filter(AA_var=="GSKV"), 
            aes(x = rank, y = prob + 0.002, label = "GSKV"),angle=45,hjust=0) +
  xlab("Genotype rank") + ylab("Probability") + ggtitle("Directional sln") +
  theme_classic()

p1 + p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdmut_trajectory3-1.png" alt="plot of chunk mut_trajectory3"  />
<p class="caption">plot of chunk mut_trajectory3</p>
</div>

While GSKV was accessible on AncSR2, it was not amongst the most likely genotypes to evolve after 3 mutation steps. Under random walk, DGKA, an ATRE-specific binder, was the most likely to evolve; under directional selection was EGRA, a promiscuous genotype that binds ERE, GGRE and TTRE.

Thus, the evolution of GSKV was highly contingent on change events.



---

Now that we have determined the relevant size of the mutational neighborhood around EGKA, we can ask: What is the phenotypic variation *accessible* to the ancestral EGKA genotype after 3 mutation steps? We can do this on both ancestral DBD backgrounds to further evaluate the effect of the non-RH substitutions that fixed along the phylogenetic interval.


```r
# run a discrete Markov chain from 'REF_GENOTYPE' of length 'PATH_LENGTH' for different population genetic scenarios
# REF_GENOTYPE' and 'PATH_LENGTH' are global variables
mc_Drift_ref_genotype_sr1 <- simulate_markov_chain(REF_GENOTYPE,P_drift_sr1_ntwrk,n_steps = PATH_LENGTH)
mc_DirSln_ref_genotype_sr1 <- simulate_markov_chain(REF_GENOTYPE,P_dir_sr1_ntwrk,n_steps = PATH_LENGTH)

mc_Drift_ref_genotype_sr2 <- simulate_markov_chain(REF_GENOTYPE,P_drift_sr2_ntwrk,n_steps = PATH_LENGTH)
mc_DirSln_ref_genotype_sr2 <- simulate_markov_chain(REF_GENOTYPE,P_dir_sr2_ntwrk,n_steps = PATH_LENGTH)

# Compute the PDFV from the Markov chains
pdfv_Drift_ref_genotype_sr1 <- get_PDFV_v2(mc_Drift_ref_genotype_sr1,Bg = "AncSR1",
                                           model = "Random walk",specific = TRUE,type="simulated mc")
pdfv_DirSln_ref_genotype_sr1 <- get_PDFV_v2(mc_DirSln_ref_genotype_sr1,Bg = "AncSR1",
                                            model = "Dir. Sln.",specific = TRUE,type="simulated mc")

pdfv_Drift_ref_genotype_sr2 <- get_PDFV_v2(mc_Drift_ref_genotype_sr2,Bg = "AncSR2",
                                           model = "Random walk",specific = TRUE,type="simulated mc")
pdfv_DirSln_ref_genotype_sr2 <- get_PDFV_v2(mc_DirSln_ref_genotype_sr2,Bg = "AncSR2",
                                            model = "Dir. Sln.",specific = TRUE,type="simulated mc")

# plots including promiscuous as extra phenotype
p1 <- circular_PDFV_v2(list(pdfv_Drift_ref_genotype_sr1,pdfv_DirSln_ref_genotype_sr1),
                       cols = cols[2:3],title = "AncSR1:EGKA",fill = F)
p2 <- circular_PDFV_v2(list(pdfv_Drift_ref_genotype_sr2,pdfv_DirSln_ref_genotype_sr2),
                       cols = cols[2:3],title = "AncSR2:EGKA",legend = F,fill = F)
#p1 + p2

###################

# plot without promiscuous (including all binders)
pdfv_Drift_ref_genotype_sr1_b <- get_PDFV_v2(mc_Drift_ref_genotype_sr1,Bg = "AncSR1",
                                           model = "Random walk",specific = F,type="simulated mc")
pdfv_DirSln_ref_genotype_sr1_b <- get_PDFV_v2(mc_DirSln_ref_genotype_sr1,Bg = "AncSR1",
                                            model = "Dir. Sln.",specific = F,type="simulated mc")

pdfv_Drift_ref_genotype_sr2_b <- get_PDFV_v2(mc_Drift_ref_genotype_sr2,Bg = "AncSR2",
                                           model = "Random walk",specific = F,type="simulated mc")
pdfv_DirSln_ref_genotype_sr2_b <- get_PDFV_v2(mc_DirSln_ref_genotype_sr2,Bg = "AncSR2",
                                            model = "Dir. Sln.",specific = F,type="simulated mc")

p3 <- circular_PDFV_v2(list(pdfv_Drift_ref_genotype_sr1_b,pdfv_DirSln_ref_genotype_sr1_b),
                       cols = cols[2:3],title = "AncSR1:EGKA",fill = F)
p4 <- circular_PDFV_v2(list(pdfv_Drift_ref_genotype_sr2_b,pdfv_DirSln_ref_genotype_sr2_b),
                       cols = cols[2:3],title = "AncSR2:EGKA",legend = F,fill = F)
p3 + p4
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdmc_chains_ref_genotype_a-1.png" alt="plot of chunk mc_chains_ref_genotype_a"  />
<p class="caption">plot of chunk mc_chains_ref_genotype_a</p>
</div>

The figures above showed that after 3 mutation steps, the ancestral RH genotype, EGKA, had very different PDFVs on each GP map. From AncSR1 to AncSR2 the probability of retaining the ancestral function decreased, and the probability of evolving SRE1 binding increased. Interestingly, the probability of evolving ATRE binding, a non-historical phenotype, was higher than that of SRE1 on the AncSR2 background. 

---

Let's now explore the temporal dynamics of pehenotypic variation: how the probability of each DNA binding phenotype changes at each step along the Markov chain. We're particularly interested in tracking the change in probability of the derived historical function SRE1 on both DBD backgrounds. We will focus for now on the `P` matrix from the Drift scenario.


```r
# Number of iterations to run the Markov chain
mc_iter <- 10
SPECIFICITY=TRUE
PROMISCUOUS=TRUE

# Simulations under random walk, directional selection, and stabilizing selection
pdfv_mc_multistep_sr1 <- simulate_markov_chain_multistep(REF_GENOTYPE,P_drift_sr1_ntwrk,mc_iter,"AncSR1",specific = SPECIFICITY)
pdfv_mc_multistep_sr2 <- simulate_markov_chain_multistep(REF_GENOTYPE,P_drift_sr2_ntwrk,mc_iter,"AncSR2",specific = SPECIFICITY)

# re-normalize functions excluding promiscuous
if(!PROMISCUOUS){
  pdfv_mc_multistep_sr1 <- lapply(pdfv_mc_multistep_sr1, function(x) x %>% filter(RE != "Promiscuous") %>% 
                                    mutate(Norm_F_prob = Norm_F_prob /sum(Norm_F_prob)))

  pdfv_mc_multistep_sr2 <- lapply(pdfv_mc_multistep_sr2, function(x) x %>% filter(RE != "Promiscuous") %>% 
                                    mutate(Norm_F_prob = Norm_F_prob /sum(Norm_F_prob)))
}

# Plot the trajectory: Include stationary distributoins as "extra mutation steps"
stationary_PDFV_sr1 <- get_PDFV_v2(P_drift_sr1_ntwrk_statdist,Bg = "AncSR1",model = mc_iter+1,specific = T,type="simulated mc")
stationary_PDFV_sr2 <- get_PDFV_v2(P_drift_sr2_ntwrk_statdist,Bg = "AncSR2",model = mc_iter+1,specific = T,type="simulated mc")
#s0 <- data.frame(RE=REs[[1]],Norm_F_prob=c(rep(0,2),1,rep(0,13)),model=0)
s0 <- data.frame(RE=REs[[2]],Norm_F_prob=c(rep(0,2),1,rep(0,14)),model=0)
#RE_COLS <- hex_RE_colors
RE_COLS <- c(hex_RE_colors,"Promiscuous"="black")

col_df <- do.call(rbind,pdfv_mc_multistep_sr1) %>% rbind(.,s0,stationary_PDFV_sr1) %>% 
 filter(RE=="SRE1 (AA)")
p1 <- do.call(rbind,pdfv_mc_multistep_sr1) %>% rbind(.,s0,stationary_PDFV_sr1) %>%
  ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
  geom_line(data=col_df,aes(x=model,y=Norm_F_prob),col="black",linewidth=3) + 
  geom_point() + geom_line() + scale_color_manual(values = RE_COLS) + 
  theme_classic() +
  labs(x="Mutation step",y="Probability",title = paste("AncSR1:",REF_GENOTYPE)) +
  scale_x_continuous(breaks=seq(0,mc_iter+1,1),labels=c(seq(0,mc_iter,1),Inf)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.position = "none")

col_df <- do.call(rbind,pdfv_mc_multistep_sr2) %>% rbind(.,s0,stationary_PDFV_sr2) %>% 
  filter(RE=="SRE1 (AA)")
p2 <- do.call(rbind,pdfv_mc_multistep_sr2) %>% rbind(.,s0,stationary_PDFV_sr2) %>%
  ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
  geom_line(data=col_df,aes(x=model,y=Norm_F_prob),col="black",linewidth=3) + 
  geom_point() + geom_line() + scale_color_manual(values = RE_COLS) + 
  theme_classic() +
  labs(x="Mutation step",y="Probability",title = paste("AncSR2:",REF_GENOTYPE)) +
  scale_x_continuous(breaks=seq(0,mc_iter+1,1),labels=c(seq(0,mc_iter,1),Inf)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))

p1 + p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdmc_chains_ref_genotype_b-1.png" alt="plot of chunk mc_chains_ref_genotype_b"  />
<p class="caption">plot of chunk mc_chains_ref_genotype_b</p>
</div>

```r
## heatmap showing the probability of each phenotype per mutation step
#pdfv_mc_multistep_sr1_mat <- do.call(rbind,pdfv_mc_multistep_sr1) %>% rbind(.,s0) %>%
#  pivot_wider(.,names_from="RE",values_from="Norm_F_prob") %>% dplyr::arrange(., model,model) %>%
#  column_to_rownames(var="model") %>%  select(order(names(.))) %>% relocate(`SRE1 (AA)`) %>% relocate(`SRE2 (GA)`, .after=`SRE1 (AA)`) %>% 
#  relocate(`ERE (GT)`, .after=`SRE2 (GA)`) %>% 
#  with(as.matrix(.))
#
#pdfv_mc_multistep_sr2_mat <- do.call(rbind,pdfv_mc_multistep_sr2) %>% rbind(.,s0) %>%
#  pivot_wider(.,names_from="RE",values_from="Norm_F_prob") %>% dplyr::arrange(., model,model) %>%
#  column_to_rownames(var="model") %>% select(order(names(.))) %>% relocate(`SRE1 (AA)`) %>% relocate(`SRE2 (GA)`, .after=`SRE1 (AA)`) %>% 
#  relocate(`ERE (GT)`, .after=`SRE2 (GA)`) %>%
#  with(as.matrix(.))
#
#breaksList <- seq(-1,4,0.1)
#pheatmap(t(pdfv_mc_multistep_sr1_mat[4,]),cluster_rows = F,cluster_cols = F,na_col = "black",
#         border_color = "black",
#         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
#         breaks = breaksList, legend_labels = seq(-1,4,0.1),scale="row")
#
#pheatmap(t(pdfv_mc_multistep_sr2_mat[4,]),cluster_rows = F,cluster_cols = F,na_col = "black",
#         border_color = "black",
#         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
#         breaks = breaksList, legend_labels = seq(-1,4,0.1),scale="row")
```

The figure shows the change in the probability of evolving each phenotype over time when evolution begins from the ancestral EKGA genotype. We see that SRE1 (green+black background) is inaccesible during the first 10 mutation steps on the AncSR1 background, whereas it becomes quickly accessible on the AncSR2 background. However, it becomes the most accessible specificity only after 6 mutation steps. By the third mutation step, the most likely phenotypic outcome is ATRE, whereas ERE and SRE1 are equally likely. This implies that, under the historical conditions, the evolution of SRE1 was not the most likey outcome and, therefore, there was an important role for chance (*results are robust if we use the directional sln. `P` matrix). 

---

We also know that the AncSR2 network has more promiscuous genotypes, so we can ask what is the effect of accessing promiscuous genotypes for phenotypic evolution.


```r
## Tracking the number of specific vs promiscuous genotypes at each step of the MC
genotype_types_sr1 <- list()
for(i in 1:mc_iter){
  mc_tmp <- simulate_markov_chain(REF_GENOTYPE,P_drift_sr1_ntwrk,n_steps = i)
  genotype_types_sr1[[i]] <- genotype_type_mc(mc_tmp,Bg="AncSR1",model=i)
}

genotype_types_sr2 <- list()
for(i in 1:mc_iter){
  mc_tmp <- simulate_markov_chain(REF_GENOTYPE,P_drift_sr2_ntwrk,n_steps = i)
  genotype_types_sr2[[i]] <- genotype_type_mc(mc_tmp,Bg="AncSR2",model=i)
}

s0_type <- data.frame(type=c("Specific","Promiscuous"),count=c(1,0),fraction=c(1,0),model=0)

p1 <- do.call(rbind,genotype_types_sr1) %>% rbind(.,s0_type) %>% 
  ggplot(aes(x=model,y=fraction,color=type)) + geom_vline(xintercept = 3,linetype="dashed",col="gray",linewidth=1.2) +
  geom_point() + geom_line() + theme_classic() + scale_color_manual(values = c("black","#c21a0e")) + 
  labs(x="Mutation step",y="Fraction of genotypes",title = paste("AncSR1:",REF_GENOTYPE)) +
  scale_x_continuous(breaks=seq(0,10,1),labels=c(seq(0,10,1))) + theme(legend.position = "none")

p2 <- do.call(rbind,genotype_types_sr2) %>% rbind(.,s0_type) %>% 
  ggplot(aes(x=model,y=fraction,color=type)) + geom_vline(xintercept = 3,linetype="dashed",col="gray",linewidth=1.2) +
  geom_point() + geom_line() + theme_classic() + scale_color_manual(values = c("black","#c21a0e")) + 
  labs(x="Mutation step",y="Fraction of genotypes",title = paste("AncSR2:",REF_GENOTYPE)) +
  scale_x_continuous(breaks=seq(0,10,1),labels=c(seq(0,10,1)))

p1 + p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdmc_chains_ref_genotype_c-1.png" alt="plot of chunk mc_chains_ref_genotype_c"  />
<p class="caption">plot of chunk mc_chains_ref_genotype_c</p>
</div>

```r
## Exploration of the phenotypic space 
# Data from 'pdfv_mc_multistep_sr1' and 'pdfv_mc_multistep_sr2', with option 'specific=FALSE'
step_phenotypes <- data.frame(bg=rep(c("AncSR1","AncSR2"),each=3),
                              step=rep(c(1,2,3),times=2),
                              n_phenotypes = c(1,1,2,7,11,15),
                              n_phenotypes_prom = c(0,0,1,4,8,8),
                              n_phenotypes_spec = c(1,1,1,3,3,7)) %>% 
  pivot_longer(cols = 4:5, names_to="type",values_to="n_pheno_split") %>%
  mutate(legend = ifelse(type=="n_phenotypes_prom","Promiscuous","Specific"))

p3 <- step_phenotypes %>% filter(bg=="AncSR1") %>% 
  ggplot(aes(x=step,y=n_pheno_split,fill=legend)) + 
  geom_bar(stat="identity",color="black") +
  scale_fill_manual(name="Genotype type",values = c("#f0e92e","gray")) +
  theme_classic() + scale_y_continuous(limits = c(0,16),breaks = c(4,8,12,16)) +
  labs(x="Mutation step",y="No. of phenotypes") +
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title=element_text(size=12,face="bold"),
        legend.position = "none") +
  geom_text(aes(x=step,label=n_phenotypes,y=n_phenotypes),vjust=-0.2)

p4 <- step_phenotypes %>% filter(bg=="AncSR2") %>% 
  ggplot(aes(x=step,y=n_pheno_split,fill=legend)) + 
  geom_bar(stat="identity",color="black") +
  scale_fill_manual(name="Genotype type",values = c("#f0e92e","gray")) +
  theme_classic() + scale_y_continuous(limits = c(0,16),breaks = c(4,8,12,16)) +
  labs(x="Mutation step",y="No. of phenotypes") +
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title=element_text(size=12,face="bold")) +
  geom_text(aes(x=step,label=n_phenotypes,y=n_phenotypes),vjust=-0.2)

p3 + p4
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdmc_chains_ref_genotype_c-2.png" alt="plot of chunk mc_chains_ref_genotype_c"  />
<p class="caption">plot of chunk mc_chains_ref_genotype_c</p>
</div>

We see that EGKA can access more promiscuous genotypes on AncSR2, and also faster - within few mutations it gains access to more promiscuous genotypes. The second figure shows the effects of accessing promiscuous genotypes for how the ancestral genotype can explore the phenotypic space. First, we see that EGKA quickly explores the phenotypic space on the AncSR2 netwrok: by the third mutation step it has explord 15/16 phenotypes, compared to only 2/16 on the AncSR1 network. Second, we can see that at each step on the AncSR2 network, at least 50% of the new phenotypes are due to access to promiscuous genotypes.


### Phenotypic transitions

Besides asking what phenotypic variation was accessible to the ancestral genotype - EGKA -, we can also ask which phenotypic transitions where more likely to happen at each point in the evolutionary trajectory. That is, we know that the transition ERE-->SRE is what happened during evolution, but we still can ask: 1) what was the likelihood of this transition? and 2) What was the likelihood of all possible phenotypic transitions?

To answer these questions, we can specify Markov chains to simulate evolution on each ancestral genotype network. The initial states for each chain will be the genotypes that bind *specifically* to each DNA element **and** are part of the main network. This specification allows us to estimate the probability of each phenotypic transition averaged over all possible initial states.


```r
# Phenotypic transitions only between specific genotypes --> change in specificity
pheno_transition_sr1_spec <- phenotypic_transitions(from=REs[[1]],to=REs[[1]],tr_mat=P_drift_sr1_ntwrk,bg = "AncSR1",n_steps = PATH_LENGTH,specific = T)
pheno_transition_sr2_spec <- phenotypic_transitions(from=REs[[1]],to=REs[[1]],tr_mat=P_drift_sr2_ntwrk,bg = "AncSR2",n_steps = PATH_LENGTH,specific = T)

# optional: compute global transition probs --> from all specific genotype in main network to every specific genotype
# AncSR1 network
spec_nodes_sr1 <- phenotypes_tbl %>% filter(specific=="YES" & bg=="AncSR1") %>% pull(AA_var)
ntwrk_transition_sr1 <- phenotypic_transitions(from=NULL,from_nodes = spec_nodes_sr1,tr_mat=P_drift_sr1_ntwrk,
                                               bg = "AncSR1",n_steps = PATH_LENGTH,specific = T)
ntwrk_transition_sr1 <- apply(ntwrk_transition_sr1,2,sum,na.rm=T)
ntwrk_transition_sr1 <- as.matrix(t(ntwrk_transition_sr1/sum(ntwrk_transition_sr1))); rownames(ntwrk_transition_sr1) <- "Network"
pheno_transition_sr1_spec <- rbind(ntwrk_transition_sr1,pheno_transition_sr1_spec)

# AncSR2 network
spec_nodes_sr2 <- phenotypes_tbl %>% filter(specific=="YES" & bg=="AncSR2") %>% pull(AA_var)
ntwrk_transition_sr2 <- phenotypic_transitions(from=NULL,from_nodes = spec_nodes_sr2,tr_mat=P_drift_sr2_ntwrk,
                                               bg = "AncSR2",n_steps = PATH_LENGTH,specific = T)
ntwrk_transition_sr2 <- apply(ntwrk_transition_sr2,2,sum,na.rm=T)
ntwrk_transition_sr2 <- as.matrix(t(ntwrk_transition_sr2/sum(ntwrk_transition_sr2))); rownames(ntwrk_transition_sr2) <- "Network"
pheno_transition_sr2_spec <- rbind(ntwrk_transition_sr2,pheno_transition_sr2_spec)

# scaled probabilities
pheno_transition_sr1_spec_scaled <- t(apply(pheno_transition_sr1_spec,1,scale)); colnames(pheno_transition_sr1_spec_scaled) <- REs[[1]]
pheno_transition_sr2_spec_scaled <- t(apply(pheno_transition_sr2_spec,1,scale)); colnames(pheno_transition_sr2_spec_scaled) <- REs[[1]]

breaksList <- seq(-1,4,0.1)
pheatmap(pheno_transition_sr1_spec_scaled,cluster_rows = F,cluster_cols = F,na_col = "black",
         border_color = "black",
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, legend_labels = seq(-1,4,0.2),main="AncSR1 transitions")
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdpheno_transitions-1.png" alt="plot of chunk pheno_transitions"  />
<p class="caption">plot of chunk pheno_transitions</p>
</div>

```r
pheatmap(pheno_transition_sr2_spec_scaled,cluster_rows = F,cluster_cols = F,na_col = "black",
         border_color = "black",
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, legend_labels = seq(-1,4,0.2),main="AncSR2 transitions")
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdpheno_transitions-2.png" alt="plot of chunk pheno_transitions"  />
<p class="caption">plot of chunk pheno_transitions</p>
</div>

These plots show two very contrasting patterns: On the one hand, on AncSR1 genetic background, SRE1 was only more likely to evolve if evolution had started from GA-specific genotyes. In most cases, the most likely event was to remain on the same phenotype. On the other hand, on AncSR2 genetic background, SRE1 is amogst the most likely phenotypes to evolve from almost all phenotypes. An interesting similarity between both backgrounds, is that when integrating the evolutionary trajectory over the entire network, SRE1 is the most likely phenotype to evolve. However, the historical phenotypic transition (ERE --> SRE1), was very unlikely on AncSR1 but becomes a high-probability transition on AncSR2, suggesting that the GP map and its features evolved along the phylogenetic interval (*Note*: black rows represent inexistent phenotypes or phenotypes that are not bound specifically).

We can also see that in almost every case (on both backgrounds), the most likely event is to remain with the same DNA binding specificity. This suggests that robustness is an intrinsic feature of these GP maps, and need not be a "selected" feature.

Overall, the AncSR2 genetic background has two interesting patterns: First, we observed before pattern of *phenotype bias* where SRE1 is the most likely phenotype to arise by random mutation; second, we can now observe a pattern of *phenotypic transition bias* where the evolutionary transition X --> SRE1 is the most likely to occur. This highlights the importance of the structure of the GP map on phenotypic evolution (* all the results are robust if we use the directional sln. `P` matrix).


## Evolution of protein-DNA complexes

Thus far we have only considered the GP map from the perspective of the protein. However, in reality, the SR regulatory module evolves as a protein-DNA complex, where both interacting partners can mutate. In this scenario, each node of the GP map coresponds to a protein-DNA complex genotype (instead of a protein genotype), and two complexes are neighbors if they are one mutation appart at the protein OR DNA side (note that there are no "promiscuous" genotypes anymore, because each complex is a unique genotype).

This co-evolutionary process can have important implications for the evolution of novel regulatory phenotypes, because it might drastically change the dymaics of the mutation process on the genotype network.

### Reading in the data

*Note:* The transition matrices and genotype networks contained in the `MutSel_matrices_complexes_final.RData` file were generated in the Midway3 cluster using the R script `MutSel_matrices_complete_data_complexes_final.R`.


```r
# Check whether the matrices have been created already
if(!file.exists(file.path(".", "MutSel_matrices_complexes_final.RData"))) {
  # Load complete data from mutation effects model
  meanF_data <- readr::read_csv(file.path("..","..","results","genotype_phenotype_distributions","meanF_data_fxnal.csv.gz"))

  # Add column of protein-DNA complexes
  meanF_data <- meanF_data %>% mutate(RE_mod = case_when(RE == "ERE (GT)" ~ "GT",
                                                         RE == "SRE1 (AA)" ~ "AA",
                                                         RE == "SRE2 (GA)" ~ "GA",
                                                         TRUE ~ RE),
                                      complex = paste(AA_var,RE_mod,sep=""))

  
  # Reference wild-type ancestral genotypes:
  AncSR1_ERE_ref <- meanF_data %>% filter(AA_var == "EGKA" & bg == "AncSR1" & RE == "ERE (GT)") %>% pull(avg_meanF)
  AncSR2_SRE_ref <- meanF_data %>% filter(AA_var == "GSKV" & bg == "AncSR2" & RE == "SRE1 (AA)") %>% pull(avg_meanF)
  
  # BUILD PHENOTYPE TABLE #
  # For each *functional* complex (those with meanF >= AncSR2_SRE_ref):
  phenotypes_tbl <- meanF_data %>% group_by(complex, bg) %>%
    summarise(n_bound_REs = n(), # how many DNA elements can bind
              meanF_bREs = avg_meanF, # the average meanF for the prot-DNA complex
              specific = ifelse(n_bound_REs == 1, "YES","NO"), # functional binding to only one DNA element?
              specificity = ifelse(n_bound_REs == 1, RE, "Promiscuous"), # Determine the type of specificity
              bound_REs = list(RE)) # assign DNA elements bound
  
  #phenotypes_tbl_complex <- phenotypes_tbl # save it for later

} else {
  # load matrices if already created
  load("./MutSel_matrices_complexes_final.RData")
  phenotypes_tbl_complex <- phenotypes_tbl # save it for later
}
```

### Building protein-DNA genotype networks and probability transition matrices

For every evolutionary scenario, we will create a `P` matrix for each background, using the meanF for each protein-DNA complex.


```r
# Check whether the matrices have been loaded already (above), if not proceed to build the networks and matrices
if(!file.exists(file.path(".", "MutSel_matrices_complexes_final.RData"))) {
  
  # BUILD GENOTYPE NETWORKS #
  sr1_complexes <- phenotypes_tbl %>% filter(bg == "AncSR1") %>% pull(complex)
  sr2_complexes <- phenotypes_tbl %>% filter(bg == "AncSR2") %>% pull(complex)
  
  # Regular networks
  net_sr1_complex <- build_genotype_network(nodes=sr1_complexes, type=4, cores=N_CORES)
  net_sr2_complex <- build_genotype_network(nodes=sr2_complexes, type=4, cores=N_CORES)
  
  #####
  
  # BUILD TRANSITION PROBABILITY MATRICES #
  
  # Generate transition probability matrix under three evolutionary scenarios for AncSR1 and AncSR2 backgrounds (From =   rows, to = cols)
  
  # Build model scenarios for each DBD background:
  # AncSR1
  MODEL.PARAM_SR1_drift <- list("AncSR1","drift",Ne,STEP.PARAM,TRUE,"mean")
  MODEL.PARAM_SR1_directional <- list("AncSR1","directional",Ne,LOG.PARAM,TRUE,"mean")
  #MODEL.PARAM_SR1_stabilizing <- list("AncSR1","stabilizing",Ne,NORM.PARAM,TRUE,"max")
  
  # AncSR2
  MODEL.PARAM_SR2_drift <- list("AncSR2","drift",Ne,STEP.PARAM,TRUE,"mean")
  MODEL.PARAM_SR2_directional <- list("AncSR2","directional",Ne,LOG.PARAM,TRUE,"mean")
  #MODEL.PARAM_SR2_stabilizing <- list("AncSR2","stabilizing",Ne,NORM.PARAM,TRUE,"max")
  
  # Build matrices
  # AncSR1
  adj_mat_complex_sr1 <- build_mutation_matrix(sr1_complexes,type=4,N_CORES)
  M_drift_sr1 <- build_transition_matrix_v2(sr1_complexes,adj_mat_complex_sr1,phenotypes_tbl,MODEL.PARAM_SR1_drift,N_CORES, complex=TRUE)
  M_dir_sr1 <- build_transition_matrix_v2(sr1_complexes,adj_mat_complex_sr1,phenotypes_tbl,MODEL.PARAM_SR1_directional,N_CORES,complex=TRUE)
  
  # AncSR2
  adj_mat_complex_sr2 <- build_mutation_matrix(sr2_complexes,type=4,N_CORES)
  M_drift_sr2 <- build_transition_matrix_v2(sr2_complexes,adj_mat_complex_sr2,phenotypes_tbl,MODEL.PARAM_SR2_drift,N_CORES,complex=TRUE)
  M_dir_sr2 <- build_transition_matrix_v2(sr2_complexes,adj_mat_complex_sr2,phenotypes_tbl,MODEL.PARAM_SR2_directional,N_CORES,complex=TRUE)
  
  # convert adjacency matrices of mutation rates to transition probability matrices
  adj_mat_complex_sr1 <- t(apply(adj_mat_complex_sr1, 1, function(x) x / sum(x)))
  adj_mat_complex_sr1 <- replace(adj_mat_complex_sr1,is.nan(adj_mat_complex_sr1),0)
  adj_mat_complex_sr1 <- as(adj_mat_complex_sr1, "sparseMatrix")
  
  adj_mat_complex_sr2 <- t(apply(adj_mat_complex_sr2, 1, function(x) x / sum(x)))
  adj_mat_complex_sr2 <- replace(adj_mat_complex_sr2,is.nan(adj_mat_complex_sr2),0)
  adj_mat_complex_sr2 <- as(adj_mat_complex_sr2, "sparseMatrix")
}
```


### Global properties of ancestral protein-DNA genotype networks

Like before, let's explore some global properties of the networks of protein-DNA genotypes. First, let's evaluate the PDFVs. The variational PDFV should be the same as the *specific binding PDFV* from the protein GP maps, because the functional complexes are still the same. However, the stationary PDFV under each evolutionary scenario might be different because structure of the genotype networks is different. 


```r
## Extract main components of each network and build the corresponding square P matrices
# Drift
P_drift_sr1_ntwrk_complex <- extract_main_ntwrk(as.undirected(net_sr1_complex),M_drift_sr1)
```

```
## This graph was created by an old(er) igraph version.
##   Call upgrade_graph() on it to use with the current igraph version
##   For now we convert it on the fly...
```

```r
P_drift_sr1_ntwrk_statdist_complex <- stationary_dist(P_drift_sr1_ntwrk_complex) 

P_drift_sr2_ntwrk_complex <- extract_main_ntwrk(as.undirected(net_sr2_complex),M_drift_sr2) 
```

```
## This graph was created by an old(er) igraph version.
##   Call upgrade_graph() on it to use with the current igraph version
##   For now we convert it on the fly...
```

```r
P_drift_sr2_ntwrk_statdist_complex <- stationary_dist(P_drift_sr2_ntwrk_complex)

# Directional selection
P_dir_sr1_ntwrk_complex <- extract_main_ntwrk(as.undirected(net_sr1_complex),M_dir_sr1) 
P_dir_sr1_ntwrk_statdist_complex <- stationary_dist(P_dir_sr1_ntwrk_complex) 

P_dir_sr2_ntwrk_complex <- extract_main_ntwrk(as.undirected(net_sr2_complex),M_dir_sr2) 
P_dir_sr2_ntwrk_statdist_complex <- stationary_dist(P_dir_sr2_ntwrk_complex)

##########

## Variational properties of the complexes should be the same as the "binding PDFV" of the protein-only GP map
var.prop_AncSR1_complex_df <- get_PDFV_v2(type="network",Bg="AncSR1",model="Binding (complex)")
var.prop_AncSR2_complex_df <- get_PDFV_v2(type="network",Bg="AncSR2",model="Binding (complex)")

# Plot variational propensities for each DBD
a <- circular_PDFV_v2(list(var.prop_AncSR1_complex_df,var.prop_AncSR1_df),
                      cols = c("#1D1D35",cols[1]),fill=F,title = "AncSR1\nVariational propensity")
b <- circular_PDFV_v2(list(var.prop_AncSR2_complex_df,var.prop_AncSR2_df),
                      cols = c("#1D1D35",cols[1]),fill=F,legend = F,title = "AncSR2\nVariational propensity")

a + b
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdglobal_prop_complexes-1.png" alt="plot of chunk global_prop_complexes"  />
<p class="caption">plot of chunk global_prop_complexes</p>
</div>

```r
# Multinomial Exact test
n_functional_prot_genotypes_sr1 <- 164
n_functional_complexes_sr1 <- 189
n_functional_prot_genotypes_sr2 <- 2036
n_functional_complexes_sr2 <- 3983

var.prop_AncSR1 <- inner_join(var.prop_AncSR1_complex_df,var.prop_AncSR1_df,by="RE") %>%
  mutate(Norm_F_prob.x = ifelse(Norm_F_prob.x==0,.Machine$double.xmin,Norm_F_prob.x),
         Norm_F_prob.y = ifelse(Norm_F_prob.y==0,.Machine$double.xmin,Norm_F_prob.y),
         obs_complexes = Norm_F_prob.x * n_functional_complexes_sr1, exp_prot = Norm_F_prob.y * n_functional_prot_genotypes_sr1)

var.prop_AncSR2 <- inner_join(var.prop_AncSR2_complex_df,var.prop_AncSR2_df,by="RE") %>%
  mutate(Norm_F_prob.x = ifelse(Norm_F_prob.x==0,.Machine$double.xmin,Norm_F_prob.x),
         Norm_F_prob.y = ifelse(Norm_F_prob.y==0,.Machine$double.xmin,Norm_F_prob.y),
         obs_complexes = Norm_F_prob.x * n_functional_complexes_sr1, exp_prot = Norm_F_prob.y * n_functional_prot_genotypes_sr1)

print("MET - Variational properties AncSR1:")
```

```
## [1] "MET - Variational properties AncSR1:"
```

```r
xmonte(obs = var.prop_AncSR1$obs_complexes, expr = var.prop_AncSR1$exp_prot)
```

```
## 
## P value (LLR) = 1 ± 0
```

```r
print("MET - Variational properties AncSR2:")
```

```
## [1] "MET - Variational properties AncSR2:"
```

```r
xmonte(obs = var.prop_AncSR2$obs_complexes, expr = var.prop_AncSR2$exp_prot)
```

```
## 
## P value (LLR) = 1 ± 0
```

```r
##########

# STATIONARY PDFVs
# AncSR1
stat_pdfv_complex_drift_sr1 <- get_PDFV_v2(P_drift_sr1_ntwrk_statdist_complex,type="simulated mc",
                                           Bg="AncSR1",model="Drift (complex)",complex = TRUE)
stat_pdfv_complex_dir_sr1 <- get_PDFV_v2(P_dir_sr1_ntwrk_statdist_complex,type="simulated mc",
                                         Bg="AncSR1",model="Dir. sln. (complex)",complex = TRUE)

p1 <- circular_PDFV_v2(list(var.prop_AncSR1_complex_df,stat_pdfv_complex_drift_sr1,stat_pdfv_complex_dir_sr1),
                       cols = c(cols[2:3],"gray60"),title = "AncSR1: Stationary PDFVs",fill = F)

# AncSR2
stat_pdfv_complex_drift_sr2 <- get_PDFV_v2(P_drift_sr2_ntwrk_statdist_complex,type="simulated mc",
                                           Bg="AncSR2",model="Drift (complex)",complex = TRUE)
stat_pdfv_complex_dir_sr2 <- get_PDFV_v2(P_dir_sr2_ntwrk_statdist_complex,type="simulated mc",
                                         Bg="AncSR2",model="Dir. sln. (complex)",complex = TRUE)

p2 <- circular_PDFV_v2(list(var.prop_AncSR2_complex_df,stat_pdfv_complex_drift_sr2,stat_pdfv_complex_dir_sr2),
                       cols = c(cols[2:3],"gray60"),title = "AncSR2: Stationary PDFVs",legend = F,fill = F)

p1+p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdglobal_prop_complexes-2.png" alt="plot of chunk global_prop_complexes"  />
<p class="caption">plot of chunk global_prop_complexes</p>
</div>

```r
##########

# Comparison of stationary PDFVs: protein vs. prot-DNA networks.
p3 <- circular_PDFV_v2(list(stat_pdfv_drift_sr1,stat_pdfv_complex_drift_sr1), cols = cols[1:2],
                       title = "AncSR1: Random walk PDFVs",fill = F)

p4 <- circular_PDFV_v2(list(stat_pdfv_dir_sr1,stat_pdfv_complex_dir_sr1), cols = cols[3:4],
                       title = "AncSR1: Directional sln PDFVs",fill = F)

p5 <- circular_PDFV_v2(list(stat_pdfv_drift_sr2,stat_pdfv_complex_drift_sr2), cols = cols[1:2],
                       title = "AncSR2: Random walk PDFVs",fill = F,legend = F)

p6 <- circular_PDFV_v2(list(stat_pdfv_dir_sr2,stat_pdfv_complex_dir_sr2), cols = cols[3:4],
                       title = "AncSR2: Directional sln PDFVs",fill = F,legend = F)


p3 + p5
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdglobal_prop_complexes-3.png" alt="plot of chunk global_prop_complexes"  />
<p class="caption">plot of chunk global_prop_complexes</p>
</div>

```r
p4 + p6
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdglobal_prop_complexes-4.png" alt="plot of chunk global_prop_complexes"  />
<p class="caption">plot of chunk global_prop_complexes</p>
</div>

```r
# Multinomial Exact test (MET)
rdm_walk_sr1_df <- inner_join(stat_pdfv_complex_drift_sr1,stat_pdfv_drift_sr1,by="RE") %>% 
  mutate(Norm_F_prob.x = ifelse(Norm_F_prob.x==0,.Machine$double.xmin,Norm_F_prob.x),
         Norm_F_prob.y = ifelse(Norm_F_prob.y==0,.Machine$double.xmin,Norm_F_prob.y),
         obs_complexes = Norm_F_prob.x * n_functional_complexes_sr1, exp_prot = Norm_F_prob.y * n_functional_prot_genotypes_sr1)

dirsln_sr1_df <- inner_join(stat_pdfv_complex_dir_sr1,stat_pdfv_dir_sr1,by="RE") %>% 
  mutate(Norm_F_prob.x = ifelse(Norm_F_prob.x==0,.Machine$double.xmin,Norm_F_prob.x),
         Norm_F_prob.y = ifelse(Norm_F_prob.y==0,.Machine$double.xmin,Norm_F_prob.y),
         obs_complexes = Norm_F_prob.x * n_functional_complexes_sr1, exp_prot = Norm_F_prob.y * n_functional_prot_genotypes_sr1)

rdm_walk_sr2_df <- inner_join(stat_pdfv_complex_drift_sr2,stat_pdfv_drift_sr2,by="RE") %>% 
  mutate(Norm_F_prob.x = ifelse(Norm_F_prob.x==0,.Machine$double.xmin,Norm_F_prob.x),
         Norm_F_prob.y = ifelse(Norm_F_prob.y==0,.Machine$double.xmin,Norm_F_prob.y),
         obs_complexes = Norm_F_prob.x * n_functional_complexes_sr2, exp_prot = Norm_F_prob.y * n_functional_prot_genotypes_sr2)

dirsln_sr2_df <- inner_join(stat_pdfv_complex_dir_sr2,stat_pdfv_dir_sr2,by="RE") %>% 
  mutate(Norm_F_prob.x = ifelse(Norm_F_prob.x==0,.Machine$double.xmin,Norm_F_prob.x),
         Norm_F_prob.y = ifelse(Norm_F_prob.y==0,.Machine$double.xmin,Norm_F_prob.y),
         obs_complexes = Norm_F_prob.x * n_functional_complexes_sr2, exp_prot = Norm_F_prob.y * n_functional_prot_genotypes_sr2)

print("MET - Drift AncSR1:")
```

```
## [1] "MET - Drift AncSR1:"
```

```r
xmonte(obs = rdm_walk_sr1_df$obs_complexes, expr = rdm_walk_sr1_df$exp_prot)
```

```
## 
## P value (LLR) = 0 ± 0
```

```r
print("MET - Directional sln AncSR1:")
```

```
## [1] "MET - Directional sln AncSR1:"
```

```r
xmonte(obs = rdm_walk_sr2_df$obs_complexes, expr = rdm_walk_sr2_df$exp_prot)
```

```
## 
## P value (LLR) = 0 ± 0
```

```r
print("MET - Drift AncSR2:")
```

```
## [1] "MET - Drift AncSR2:"
```

```r
xmonte(obs = dirsln_sr1_df$obs_complexes, expr = dirsln_sr1_df$exp_prot)
```

```
## 
## P value (LLR) = 0 ± 0
```

```r
print("MET - Directional sln AncSR2:")
```

```
## [1] "MET - Directional sln AncSR2:"
```

```r
xmonte(obs = dirsln_sr2_df$obs_complexes, expr = dirsln_sr2_df$exp_prot)
```

```
## 
## P value (LLR) = 0 ± 0
```

These results show that while the variational PDFVs are exactly the same between the protein and protein-DNA GP maps, the stationary distributions under drift and directional selection are significantly different. This suggests that protein-DNA co-evolution results in different equilibrium dynamics for the mutation-selection process, and affects the relative accessibility of new DNA binding phenotypes.

### Discrete Markov chains

### Evolution around EGKA on AncSR1 and AncSR2 protein-DNA networks

let's consider how protein-DNA co-evolution affects the mutational trajectory between the two historical complexes: EGKA/GT and GSKV/AA. We know that the derived genotype was not accessible on the AncSR1 background, so we will only consider the mutational trajectories on the AncSR2 background. 

Like before, we want to infer the length of the mutational trajectory between ancestral and derived complexes. Now the hamming distance between complexes is 5, but it doesn't mean that the trajectory was exactly 5 substitutions. 


```r
# Mutational trajectory on the protein genotype network
print("Mutational trajectories on the protein genotype network:")
```

```
## [1] "Mutational trajectories on the protein genotype network:"
```

```r
all_shortest_paths(net_sr2,from = which(names(V(net_sr2))=="EGKA"),
                               to = which(names(V(net_sr2))=="GSKV"))$res
```

```
## [[1]]
## + 4/2390 vertices, named, from cee954f:
## [1] EGKA GGKA GGKV GSKV
```

```r
# Mutational trajectory on the protein-DNA genotype network
print("Mutational trajectories on the protein-DNA genotype network:")
```

```
## [1] "Mutational trajectories on the protein-DNA genotype network:"
```

```r
all_shortest_paths(as.undirected(net_sr2_complex),
                                       from = which(names(V(net_sr2_complex))=="EGKAGT"),
                                       to = which(names(V(net_sr2_complex))=="GSKVAA"))$res
```

```
## [[1]]
## + 15/5032 vertices, named, from e6f47fd:
##  [1] EGKAGT EGKSGT EGKIGT EGKKGT EAKKGT EAKKGA QAKKGA QAKMGA QAKMAA LAKMAA
## [11] WAKMAA WGKMAA GGKMAA GSKMAA GSKVAA
```

We see four important things:

- There is only one mutational trajectory on the protein network and on the complex network.
- The length of the mutational trajectories under co-evolution is longer. It requires 14 mutation steps to reach the derived complex, vs. 3 mutation steps on the protein network.
- The identity of the intermediate protein genotypes is different. On the protein network, the trajectory involves only intermediate genotypes with the ancestral and derived states (i.e., GGKA and GGKV); on the complex network, it involves protein genotypes with alternative (non-historical) amino acid states states, i.e., E, Q, L or W at the 1st position; A at the 2nd position and S,M,I or K at the 4th positiion.
- The evolution of the derived specificity can be achieved within 1 mutation step on the protein network (GGKA is SRE1-specific), whereas in the trajectory of the complex network, the evolutionof the derived specificity requires passing through intermediate GARE-binding genotypes.

Unlike before, the length of the mutational trjectory on the protein-DNA network does not have the same length as the number of differences between genotypes. Let's now track the probability of the genotyoe complex GSKV:AA over time to get an idea of the likelihood function.


```r
# Track the change in probability of evolving GSKVAA on AncSR2 --> Identify the mutation step at which P(GSKV:AA) is maximized
p_gskvaa <- c()
mut_step <- 20
for(i in 1:mut_step){
  tmp_mc <- simulate_markov_chain(REF_GENOTYPE_COMPLEX,P_drift_sr2_ntwrk_complex,n_steps = i)
  p_gskvaa <- c(p_gskvaa,tmp_mc[names(tmp_mc)=="GSKVAA"])
}
stat_p_gskvaa <- data.frame(mut_step = mut_step+1, p_gskvaa = P_drift_sr2_ntwrk_statdist_complex[names(P_drift_sr2_ntwrk_statdist_complex)=="GSKVAA"])

# Relative likelihood
d <- data.frame(mut_step = seq(1,mut_step,1), p_gskvaa = p_gskvaa) %>% rbind(.,stat_p_gskvaa) %>%
  mutate(LR = p_gskvaa/max(p_gskvaa), # relative likelihood
         LI = ifelse(LR >= 0.95, TRUE,FALSE)) # Likelihood interval (1-a, for a = 0.05)

d %>%
  ggplot(aes(x=mut_step,y=p_gskvaa)) + geom_point(col="black") + geom_line() + 
  xlab("Mutation step") + ylab("Pr(GSKV:AA)") + 
  geom_vline(data = d %>% filter(LI==T),aes(xintercept=mut_step),col="red",linetype="dashed") +
  ggtitle("AncSR2: Random walk (complex)") + 
  scale_x_continuous(breaks = seq(1,mut_step+1,1),labels = c(seq(1,mut_step,1),Inf)) +
  theme_classic()
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdmut_trajectory_complex2-1.png" alt="plot of chunk mut_trajectory_complex2"  />
<p class="caption">plot of chunk mut_trajectory_complex2</p>
</div>

Also unlike before, we see that the likelihood function for $P(GSKV:AA)$ continuously (monotonically) increases, and is maximized when the path length $\approx \infty$ (i.e., $\hat{S} = \infty$). Note also that the $P(GSKV:AA)$ by 5 mutation steps is zero. Based on this we also know that the $LI_{\alpha = 0.05}$ must be $S > 20$. (By 50 steps the relative likelihood is 0.92). We will use the shortest path length indicated by the network, with the caveat that the relative likelihood for $S = 14$ is 0.42.


```r
PATH_LENGTH_COMPLEX = 14 # path length (neighborhood size) to find mutational trajectoties in the prot-DNA network
```

### Local neighborhoods in protein-DNA genotype-phenotype maps

Let's now ask what is the effect of co-evolution for the structure of the local neighborhood around the ancestral genotype (EGKA). To visualize the effects, let's focus on AncSR1:EGKA. We know that the historical trajectory began with the complex EGKA/GT (ERE) and evolved to GSKV/AA (SRE1). This trajectory involves 3 differences on the protein side and 2 on the DNA side (for a hamming distance of 5). 

We will plot the 3-mutational neighborhood around EGKA on the protein and protein-DNA networks, and also the 5-mutational neighborhood on the protein-DNA network to account for the changes in the complex as a whole (Note that 5 mutation steps is just for visualization purposes, we established before that will use a path length of 8).


```r
# Protein GP map
subntwrk_sr1 <- make_ego_graph(net_sr1,order=PATH_LENGTH,nodes = "EGKA")[[1]]
#phenotypes_tbl_prot %>% filter(AA_var %in% names(V(subntwrk_sr1)) & bg == "AncSR1") # identify phenotypes of the genotypes

# assign vertex attributes: DNA binding phenotypes
subntwrk_sr1 <- set.vertex.attribute(subntwrk_sr1, 'phenotype', V(subntwrk_sr1)[-7], 'ERE')
subntwrk_sr1 <- set.vertex.attribute(subntwrk_sr1, 'phenotype', V(subntwrk_sr1)[7], 'Promiscuous')

plot(subntwrk_sr1,layout=layout_with_fr,edge.curved=.1,vertex.color = c(rep(hex_RE_colors[3],6),"white",rep(hex_RE_colors[3],10)),
     vertex.label.color="black", vertex.label.dist=2,main="local neighborhood EGKA: 3 mutations")
legend(x=-1.5, y=0, c("ERE","ERE/GG"), pch=21, pt.bg=c(hex_RE_colors[3],"white"), pt.cex=2, cex=.8, bty="n", ncol=1)
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdlocal_neighborhood_complex-1.png" alt="plot of chunk local_neighborhood_complex"  />
<p class="caption">plot of chunk local_neighborhood_complex</p>
</div>

```r
##########

# Protein-DNA GP map: 3 steps
subntwrk_sr1_complex <- make_ego_graph(as.undirected(net_sr1_complex),order=PATH_LENGTH,nodes = "EGKAGT")[[1]]
#phenotypes_tbl_complex %>% filter(complex %in% names(V(subntwrk_sr1_complex)) & bg == "AncSR1")

# assign vertex attributes: DNA binding phenotypes
subntwrk_sr1_complex <- set.vertex.attribute(subntwrk_sr1_complex, 'phenotype', V(subntwrk_sr1_complex), 'ERE')

plot(subntwrk_sr1_complex,layout=layout_with_fr,edge.curved=.1,vertex.color = c(rep(hex_RE_colors[3],11)),
     vertex.label.color="black", vertex.label.dist=2,main="local neighborhood EGKA: 3 mutations (complex)")
legend(x=-1.5, y=0, c("ERE"), pch=21, pt.bg=c(hex_RE_colors[3]), pt.cex=2, cex=.8, bty="n", ncol=1)
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdlocal_neighborhood_complex-2.png" alt="plot of chunk local_neighborhood_complex"  />
<p class="caption">plot of chunk local_neighborhood_complex</p>
</div>

```r
##########

# Protein-DNA GP map: 5 steps
subntwrk_sr1_complex2 <- make_ego_graph(as.undirected(net_sr1_complex),order=5,nodes = "EGKAGT")[[1]]
#phenotypes_tbl_complex %>% filter(complex %in% names(V(subntwrk_sr1_complex2)) & bg == "AncSR1")

# assign vertex attributes: DNA binding phenotypes
subntwrk_sr1_complex2 <- set.vertex.attribute(subntwrk_sr1_complex2, 'phenotype', V(subntwrk_sr1_complex2)[-c(4,13,16)], 'ERE')
subntwrk_sr1_complex2 <- set.vertex.attribute(subntwrk_sr1_complex2, 'phenotype', V(subntwrk_sr1_complex2)[c(13,16)], 'GG')
subntwrk_sr1_complex2 <- set.vertex.attribute(subntwrk_sr1_complex2, 'phenotype', V(subntwrk_sr1_complex2)[-4], 'AT')

plot(subntwrk_sr1_complex2,layout=layout_with_fr,edge.curved=.1,
     vertex.color = c(rep(hex_RE_colors[3],3),hex_RE_colors[6],rep(hex_RE_colors[3],8),hex_RE_colors[12],rep(hex_RE_colors[3],2),
                      hex_RE_colors[12],rep(hex_RE_colors[3],21)),
     vertex.label.color="black", vertex.label.dist=2,vertex.label.cex=.7,main="local neighborhood EGKA: 5 mutations (complex)")
legend(x=-1.5, y=0, c("ERE","GG","AT"), pch=21, pt.bg=c(hex_RE_colors[3],hex_RE_colors[12],hex_RE_colors[6]), 
       pt.cex=2, cex=.8, bty="n", ncol=1)
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdlocal_neighborhood_complex-3.png" alt="plot of chunk local_neighborhood_complex"  />
<p class="caption">plot of chunk local_neighborhood_complex</p>
</div>

These plots show some interesting consequences of introducing protein-DNA co-evolution to the GP map. We saw before that co-evolution changes the probabilities of evolving each phenotype at equilibrium, but we can also zoom-in to inspect what is happenning around a focal genotype. We can observe three things:

1. The 3-mutational neighborhood on the protein network contains a promiscuous genotype (EARS) that can bind to ERE and GGRE. However, the 3-mutational neighborhood on the protein-DNA network only contains ERE-specific genotypes.
2. EARS/GG is four mutations away from EGKA/GT, three in the protein and one in the DNA, and this makes it innaccesible within 3 mutations in the protein-DNA network. EARS/GG becomes accessible in the 5-mutational neighborhood.
3. The 5-mutational neighborhood on the protein-DNA network contain the phenotypes (ERE, GG and AT), suggesting that co-evolution increases the length of mutational trajectories to reach new phenotypes.

---

We saw one specific example of the effects of introducing co-evolution to the genotype network. However, we want to ask more systematically what are the effects of co-evolution. Here we will compute the PDFV accesible to the ancestral complex. For now, we will focus for now on the `P` matrix from the Drift scenario.


```r
# run a discrete Markov chain from 'REF_GENOTYPE' of length 'PATH_LENGTH' for different population genetic scenarios
mc_Drift_ref_geno_14m_sr1_complex <- simulate_markov_chain(REF_GENOTYPE_COMPLEX,P_drift_sr1_ntwrk_complex,n_steps = PATH_LENGTH_COMPLEX) # 5-mutational neighborhood
mc_Drift_ref_geno_14m_sr2_complex <- simulate_markov_chain(REF_GENOTYPE_COMPLEX,P_drift_sr2_ntwrk_complex,n_steps = PATH_LENGTH_COMPLEX) 

# Compute the PDFV from the Markov chains
pdfv_Drift_ref_geno_14m_sr1_complex <- get_PDFV_v2(mc_Drift_ref_geno_14m_sr1_complex,Bg = "AncSR1",
                                           model = "Random walk (14m, complex)",specific = F,type="simulated mc",complex = T)

pdfv_Drift_ref_geno_14m_sr2_complex <- get_PDFV_v2(mc_Drift_ref_geno_14m_sr2_complex,Bg = "AncSR2",
                                           model = "Random walk (14m, complex)",specific = F,type="simulated mc",complex = T)


# plots
p1 <- circular_PDFV_v2(list(pdfv_Drift_ref_genotype_sr1_b,pdfv_Drift_ref_geno_14m_sr1_complex),
                       cols = cols[2:4],title = "AncSR1:EGKA or EGKA:GT",fill = F)
p2 <- circular_PDFV_v2(list(pdfv_Drift_ref_genotype_sr2_b,pdfv_Drift_ref_geno_14m_sr2_complex),
                       cols = cols[2:4],title = "AncSR2:EGKA or EGKA:GT",legend = F,fill = F)
p1 + p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdmc_chain_ref_genotype_complex-1.png" alt="plot of chunk mc_chain_ref_genotype_complex"  />
<p class="caption">plot of chunk mc_chain_ref_genotype_complex</p>
</div>

We see that the PDFVs around the local neighborhood of the EGKA/GT complex have important differences with that of the protein network. For instance, on the AncSR2 background, the historical derived phenotype (SRE1) is less accessible and ATRE is more accessible when we incorporate co-evolution.

---

Let's now explore the temporal dynamics of pehenotypic variation: how the probability of each DNA binding phenotype changes at each step along the Markov chain. We're particularly interested in tracking the change in probability of the derived historical function SRE1 on both DBD backgrounds. 


```r
# Number of iterations to run the Markov chain
mc_iter <- 20
SPECIFICITY=FALSE

# Simulations under random walk, directional selection, and stabilizing selection
pdfv_mc_multistep_sr1_complex <- simulate_markov_chain_multistep(REF_GENOTYPE_COMPLEX,P_drift_sr1_ntwrk_complex,mc_iter,
                                                                 "AncSR1",specific = SPECIFICITY,complex = TRUE)
pdfv_mc_multistep_sr2_complex <- simulate_markov_chain_multistep(REF_GENOTYPE_COMPLEX,P_drift_sr2_ntwrk_complex,mc_iter,
                                                                 "AncSR2",specific = SPECIFICITY,complex = TRUE)

# Plot the trajectory: Include stationary distributoins as "extra mutation steps"
stationary_PDFV_sr1_complex <- get_PDFV_v2(P_drift_sr1_ntwrk_statdist_complex,Bg = "AncSR1",
                                           model = mc_iter+1,specific = T,type="simulated mc",complex = TRUE)
stationary_PDFV_sr2_complex <- get_PDFV_v2(P_drift_sr2_ntwrk_statdist_complex,Bg = "AncSR2",
                                   model = mc_iter+1,specific = T,type="simulated mc",complex = TRUE)

s0 <- data.frame(RE=REs[[2]],Norm_F_prob=c(rep(0,2),1,rep(0,14)),model=0)
RE_COLS <- hex_RE_colors

col_df <- do.call(rbind,pdfv_mc_multistep_sr1_complex) %>% rbind(.,s0,stationary_PDFV_sr1_complex) %>% 
 filter(RE=="SRE1 (AA)")
p1 <- do.call(rbind,pdfv_mc_multistep_sr1_complex) %>% rbind(.,s0,stationary_PDFV_sr1_complex) %>%
  ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
  geom_line(data=col_df,aes(x=model,y=Norm_F_prob),col="black",linewidth=3) + 
  geom_point() + geom_line() + scale_color_manual(values = RE_COLS) + 
  theme_classic() +
  labs(x="Mutation step",y="Probability",title = paste("AncSR1:",REF_GENOTYPE_COMPLEX)) +
  scale_x_continuous(breaks=seq(0,mc_iter+1,1),labels=c(seq(0,mc_iter,1),Inf)) +
  theme(axis.text.x = element_text(size=10,angle = 45),
        axis.text.y = element_text(size=15),
        legend.position = "none")

col_df <- do.call(rbind,pdfv_mc_multistep_sr2_complex) %>% rbind(.,s0,stationary_PDFV_sr2_complex) %>% 
  filter(RE=="SRE1 (AA)")
p2 <- do.call(rbind,pdfv_mc_multistep_sr2_complex) %>% rbind(.,s0,stationary_PDFV_sr2_complex) %>%
  ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
  geom_line(data=col_df,aes(x=model,y=Norm_F_prob),col="black",linewidth=3) + 
  geom_point() + geom_line() + scale_color_manual(values = RE_COLS) + 
  theme_classic() +
  labs(x="Mutation step",y="Probability",title = paste("AncSR2:",REF_GENOTYPE_COMPLEX)) +
  scale_x_continuous(breaks=seq(0,mc_iter+1,1),labels=c(seq(0,mc_iter,1),Inf)) +
  theme(axis.text.x = element_text(size=10,angle = 45),
        axis.text.y = element_text(size=15))

p1 + p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdmc_chain_ref_genotype_complex_b-1.png" alt="plot of chunk mc_chain_ref_genotype_complex_b"  />
<p class="caption">plot of chunk mc_chain_ref_genotype_complex_b</p>
</div>

The figure shows the change in the probability of evolving each phenotype over time when evolution begins from the ancestral EKGA/GT complex genotype. We see that SRE1 (green+black background) is inaccesible during the first 20 mutation steps on the AncSR1 background. On the AncSR2 background, SRE1 becomes more likely over time, but it takes 16 mutations to be as likely as ERE (as opposed to 3 mutations in the protein network), and more than 20 mutations to become the most likely phenotypic outcome. In contrast, ATRE becomes accessible within fewer mutations, and becomes the most likely outcome after 12 mutations. 

By the 14th mutation step, accounting for the changes that must have occurred in the complex to swith from the ancestral to the derived protein-DNA genotype, we see that SRE1's probability is zero on AncSR1 and is the third most likely outcome on AncSR2, implying again that, under the historical conditions, the evolution of SRE1 involved a substantial role for chance.

Another important difference regarding the dynamics of the evolutionary process between the protein and protein-DNA networks, is the persistence of the ancestral phenotype (ERE) as the most likely outcome. We can quantify this by the number of mutations it takes for the ancestral phenotype to reduce it's probabiltity by 50% (i.e., $L_{50}$). On the AncSR1 background, the $L_{50} \approx 9$ on the protein network, whereas for the complex network the $L_{50} \approx \infty$ (the probability of ERE at equilibrium is 0.524). On the AncSR2 background, the $L_{50} \approx 1$ on the protein network, whereas for the complex network the $L_{50} \approx 6$. These results imply that the protein-DNA network is more robust than the protein network, and co-evolution can substantially delay the evolution of novel phentoypes.

### Phenotypic transitions

Finally, let's explore the effects of co-evolution for the transitions between phenotypes. Since we're considering the evolution of the SR complex, we will determine the probability of each phenotypic transition after 5 mutation steps. 


```r
# Phenotypic transitions only between specific genotypes --> change in specificity
pheno_transition_sr1_complex <- phenotypic_transitions(from=REs[[1]],to=REs[[1]],tr_mat=P_drift_sr1_ntwrk_complex,
                                                    bg = "AncSR1",n_steps = PATH_LENGTH_COMPLEX,specific = T,complex = T)
pheno_transition_sr2_complex <- phenotypic_transitions(from=REs[[1]],to=REs[[1]],tr_mat=P_drift_sr2_ntwrk_complex,
                                                    bg = "AncSR2",n_steps = PATH_LENGTH_COMPLEX,specific = T,complex = T)

# optional: compute global transition probs --> from all specific genotype in main network to every specific genotype
# AncSR1 network
spec_nodes_sr1 <- phenotypes_tbl %>% filter(specific=="YES" & bg=="AncSR1") %>% pull(complex)
ntwrk_transition_sr1 <- phenotypic_transitions(from=NULL,from_nodes = spec_nodes_sr1,tr_mat=P_drift_sr1_ntwrk_complex,
                                               bg = "AncSR1",n_steps = PATH_LENGTH_COMPLEX,specific = T,complex=T)
ntwrk_transition_sr1 <- apply(ntwrk_transition_sr1,2,sum,na.rm=T)
ntwrk_transition_sr1 <- as.matrix(t(ntwrk_transition_sr1/sum(ntwrk_transition_sr1))); rownames(ntwrk_transition_sr1) <- "Network"
pheno_transition_sr1_complex <- rbind(ntwrk_transition_sr1,pheno_transition_sr1_complex)

# AncSR2 network
spec_nodes_sr2 <- phenotypes_tbl %>% filter(specific=="YES" & bg=="AncSR2") %>% pull(complex)
ntwrk_transition_sr2 <- phenotypic_transitions(from=NULL,from_nodes = spec_nodes_sr2,tr_mat=P_drift_sr2_ntwrk_complex,
                                               bg = "AncSR2",n_steps = PATH_LENGTH_COMPLEX,specific = T,complex=T)
ntwrk_transition_sr2 <- apply(ntwrk_transition_sr2,2,sum,na.rm=T)
ntwrk_transition_sr2 <- as.matrix(t(ntwrk_transition_sr2/sum(ntwrk_transition_sr2))); rownames(ntwrk_transition_sr2) <- "Network"
pheno_transition_sr2_complex <- rbind(ntwrk_transition_sr2,pheno_transition_sr2_complex)

# plot
# scaled probabilities
pheno_transition_sr1_complex_scaled <- t(apply(pheno_transition_sr1_complex,1,scale)); colnames(pheno_transition_sr1_complex_scaled) <- REs[[1]]
pheno_transition_sr2_complex_scaled <- t(apply(pheno_transition_sr2_complex,1,scale)); colnames(pheno_transition_sr2_complex_scaled) <- REs[[1]]

breaksList <- seq(-1,4,0.1)
pheatmap(pheno_transition_sr1_complex_scaled,cluster_rows = F,cluster_cols = F,na_col = "black",
         border_color = "black",
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, legend_labels = seq(-1,4,0.2),main="AncSR1 transitions (complex)")
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdpheno_transitions_complex-1.png" alt="plot of chunk pheno_transitions_complex"  />
<p class="caption">plot of chunk pheno_transitions_complex</p>
</div>

```r
pheatmap(pheno_transition_sr2_complex_scaled,cluster_rows = F,cluster_cols = F,na_col = "black",
         border_color = "black",
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, legend_labels = seq(-1,4,0.2),main="AncSR2 transitions (complex)")
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdpheno_transitions_complex-2.png" alt="plot of chunk pheno_transitions_complex"  />
<p class="caption">plot of chunk pheno_transitions_complex</p>
</div>

These figures show contrasting results with the phenotypic transitions computed from the protein network. Specifically, we see that on AncSR1, SRE1 is the most likely outcome from CARE-binding genotypes, instead of GARE-binding phenotypes, which are not even part of the main complex network. The phenotypic transition bias towards SRE1 on AncSR2 that we observed before is now less clear, and in fact, SRE1 and ATRE have similar probabilities. As for the historical transition, ERE --> SRE1, we see that it is highly unlikely in both backgrounds. In fact, on AncSR1, ERE is unlikely to change to other phenotypes, and on AncSR2, ERE is more likely to transition to ATRE. 

When considering the entire networks, we see that SRE1, ERE and ATRE are the most likely outcomes. ATRE-binding is the second most common phenotype on AncSR1 and the third most common phenotype produced by mutation on AncSR2 (see the variational PDFVs). Thus, we now see a new pattern phenotypic transition bias when accounting for co-evolution.

We also see more clearly the robustness of the protein-DNA networks. The probabilities along the diagonals are substantially higher than the off-diaginals.

---

To better understand the new patterns, we can use the transition profiles of each phenotype to detect whether some phenotypes share more similar profiles. This can reveal some underlying structure to the landscape of phentoypic transitions. Specifically, we can cluster the rows in the heatmaps according to their patterns of transition; we will use hierarchical clustering and PCA to find the best number of clusters in each ancestral background.


```r
# Find clusters using gap statistics:
# Gap statistics compares the total within intra-cluster variation for different values of K with their expected values 
# under null reference distribution of the data.
pheno_transition_sr1_complex_scaled_filt <- data.frame(na.omit(pheno_transition_sr1_complex_scaled)[-1,]) # remove NAs and Network
pheno_transition_sr2_complex_scaled_filt <- data.frame(pheno_transition_sr2_complex_scaled[-1,])

gap.stat_sr1 <- clusGap(pheno_transition_sr1_complex_scaled_filt, FUN = kmeans, K.max = 6)
gap.stat_sr2 <- clusGap(pheno_transition_sr2_complex_scaled_filt, FUN = kmeans, K.max = 6)

# visualize gap statistics
p1 <- fviz_nbclust(pheno_transition_sr1_complex_scaled_filt, kmeans, method = "wss", k.max = 6) + theme_minimal() + ggtitle("AncSR1 transitions (SSs per cluster)")
p2 <- fviz_nbclust(pheno_transition_sr2_complex_scaled_filt, kmeans, method = "wss", k.max = 15) + theme_minimal() + ggtitle("AncSR2 transitions (SSs per cluster)")

p3 <- fviz_gap_stat(gap.stat_sr1) + ggtitle("AncSR1 transitions (gap stat)") + theme_minimal()
p4 <- fviz_gap_stat(gap.stat_sr2) + ggtitle("AncSR2 transitions (gap stat)") + theme_minimal()

p1 + p3
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdclustered_heatmaps-1.png" alt="plot of chunk clustered_heatmaps"  />
<p class="caption">plot of chunk clustered_heatmaps</p>
</div>

```r
p2 + p4
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdclustered_heatmaps-2.png" alt="plot of chunk clustered_heatmaps"  />
<p class="caption">plot of chunk clustered_heatmaps</p>
</div>

```r
##########################

# Find clusters using PCA
pca_sr1 <- prcomp(pheno_transition_sr1_complex_scaled_filt)
pca_sr2 <- prcomp(pheno_transition_sr2_complex_scaled_filt)

# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
#fviz_eig(pca_sr1)
#fviz_eig(pca_sr2)

# Graph of individuals. Individuals with a similar profile are grouped together.
# Graph of variables. Contributions and directions of the loadings of each variable
p1 <- fviz_pca_ind(pca_sr1,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) + ggtitle("AncSR1 transitions (complex) - PCA")
p2 <- fviz_pca_var(pca_sr1, col.var = "contrib",gradient.cols = c("white", "blue", "red"),repel = T) + 
  scale_x_continuous(limits = c(-2.5,1.5)) + scale_y_continuous(limits = c(-1,1.5))
p1 + p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdclustered_heatmaps-3.png" alt="plot of chunk clustered_heatmaps"  />
<p class="caption">plot of chunk clustered_heatmaps</p>
</div>

```r
p3 <- fviz_pca_ind(pca_sr2,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) + ggtitle("AncSR2 transitions (complex) - PCA")
p4 <- fviz_pca_var(pca_sr2, col.var = "contrib",gradient.cols = c("white", "blue", "red"),repel = T) + 
  scale_x_continuous(limits = c(-1,1.75)) + scale_y_continuous(limits = c(-1,1))
p3 + p4
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdclustered_heatmaps-4.png" alt="plot of chunk clustered_heatmaps"  />
<p class="caption">plot of chunk clustered_heatmaps</p>
</div>

The gap statistics and PCA for AncSR1 recover between 3-4 clusters. The first is `[SRE1+CA]` that transition more likely to SRE1; the second is `[ERE+GG+SRE2]` that transition more likely to ERE; the third is `[AT+TT+AG]` that transition more likely to ATRE; the fourth cluster would be only `TA` that transitions equally likely to ERE and ATRE.

The gap statistics and PCA for AncSR2 recover 3 clusters. The first is `[SRE2+GC]` that transition more likely to SRE2; the second is `[CG+CA+SRE1+TA+CC]` that transition more likely to SRE1 and to CARE; the third is `[ERE+GG+AC+AG+TG+TT+CT+AT+CC]` that transition more likely to ATRE.

We can see the clustering on the heatmaps:


```r
# Clustered heatmaps
pheatmap(pheno_transition_sr1_complex_scaled_filt,cluster_rows = T,cluster_cols = F,na_col = "black",
         border_color = "black",
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, legend_labels = seq(-1,4,0.2),main="AncSR1 (complex)",
         cutree_rows = 4,clustering_distance_rows="correlation")
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdclustered_heatmaps2-1.png" alt="plot of chunk clustered_heatmaps2"  />
<p class="caption">plot of chunk clustered_heatmaps2</p>
</div>

```r
pheatmap(pheno_transition_sr2_complex_scaled_filt,cluster_rows = T,cluster_cols = F,na_col = "black",
         border_color = "black",
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList, legend_labels = seq(-1,4,0.2),main="AncSR2 (complex)",
         cutree_rows = 3,clustering_distance_rows="correlation")
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdclustered_heatmaps2-2.png" alt="plot of chunk clustered_heatmaps2"  />
<p class="caption">plot of chunk clustered_heatmaps2</p>
</div>

This analysis reveals some features of the underlying structure of the protein-DNA GP map and has important implications for phenotypic evolution. First, these clusters broadly correspond to transitions to the most frequent DNA binding phenotypes produced in each background (ERE, SRE1 and ATRE, and CARE), and the members within clusters are *generally* one mutation away in the DNA from the most likely phenotypic transition(s). Second, while in the protein genotype network SRE1 was the most likely phentypic outcome on the AncSR2 network, we now see that in the protein-DNA network the evolution of SRE1 was highly contingent on the starting specificity. Specifically, SRE1 was a likely phenotypic outcome only from one of the three clusters (cluster #2: `[CG+CA+SRE1+TA+CC]`). 

Regarding the historical transition, ERE --> ATRE is more likely than ERE --> SRE1 because ATRE is mutationally closer to ERE (GT) than SRE1 (AA): access to ATRE requires only one mutation on the DNA, whereas access to SRE1 requires two. In fact, we see that ERE belongs to the "ATRE" clutser on the AncSR2 network.

These results show how the frequency by which a phenotype arises through random mutation and the co-evolutionary accessibility interact to shape the landscape of phenotypic transitions. 

---

We can also compare directly the relative total probabilities of evolving each phenotype under the different genotype networks. 


```r
total_prob_phenotype_sr1 <- apply(abs(pheno_transition_sr1_spec_scaled[-1,]),2,sum,na.rm=T) # remove whole network prob.
total_prob_phenotype_sr1_complex <- apply(abs(pheno_transition_sr1_complex_scaled[-1,]),2,sum,na.rm=T)

df_total_pheno_prob_sr1 <- rbind(data.frame(RE=names(total_prob_phenotype_sr1),prob=total_prob_phenotype_sr1,type="protein"),
                                 data.frame(RE=names(total_prob_phenotype_sr1_complex),prob=total_prob_phenotype_sr1_complex,type="protein-DNA"))

total_prob_phenotype_sr2 <- apply(abs(pheno_transition_sr2_spec_scaled[-1,]),2,sum,na.rm=T)
total_prob_phenotype_sr2_complex <- apply(abs(pheno_transition_sr2_complex_scaled[-1,]),2,sum,na.rm=T)

df_total_pheno_prob_sr2 <- rbind(data.frame(RE=names(total_prob_phenotype_sr2),prob=total_prob_phenotype_sr2,type="protein"),
                                 data.frame(RE=names(total_prob_phenotype_sr2_complex),prob=total_prob_phenotype_sr2_complex,type="protein-DNA"))

# plot
p1 <- df_total_pheno_prob_sr1 %>% mutate(RE = factor(RE,levels =REs[[1]])) %>%
  ggplot(aes(x=RE,y=prob,group=type,color=type)) + 
  geom_point() + geom_line() +
  scale_color_manual(values = c("#b01796","#1796b0")) +
  theme_classic() + ggtitle("AncSR1") +
  labs(x="DNA element",y="Relative total transition probability") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 10),
        axis.text.y = element_text(size = 13),
        axis.title=element_text(size=13,face="bold"),
        legend.position = "none")

p2 <- df_total_pheno_prob_sr2 %>% mutate(RE = factor(RE,levels =REs[[1]])) %>%
  ggplot(aes(x=RE,y=prob,group=type,color=type)) + 
  geom_point() + geom_line() + 
  scale_color_manual(values = c("#b01796","#1796b0")) +
  theme_classic() + ggtitle("AncSR2") +
  labs(x="DNA element",y="Relative total transition probability") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 10),
        axis.text.y = element_text(size = 13),
        axis.title=element_text(size=13,face="bold"))

p1 + p2
```

<div class="figure" style="text-align: center">
<img src="figure/MutSel-complete-data-final-Rmdcomparison_pheno_tr-1.png" alt="plot of chunk comparison_pheno_tr"  />
<p class="caption">plot of chunk comparison_pheno_tr</p>
</div>

The inclusion of co-evolution drastically changes the landscape of phenotypic transitions, and creates new phenotypic transition biases. Specifically, ATRE is now a high probability phenotype on both ancestral backgrounds, and also is amongst the most frequently produced phenotypes. 
---

<!-- We saw that co-evolution leads to a different pattern of phenotypic transition bias, where ATRE binding is now amongst the most likely phenotypes to evolve. This new pattern can arise because ATRE is mutationally closer to ERE (GT) than SRE1 (AA): access to ATRE requires only one mutation on the DNA, whereas access to SRE1 (AA) requires two. Thus, introducing protein-DNA co-evolution can substatially change the landscape of phenotypic transitions.

We can make sense of these patterns by exploring the proximity between neutral networks, that is, how close in sequence space are the phenotypes. For this, we can compute the number of direct links between genotypes in each neutral network; this represents the the number of direct mutational paths that exist between phenotypes, and can provide an approximation of how easy it is for evolution to access new phenotypes. -->




### References

1. Gillespie J. Molecular Evolution Over the Mutational Landscape. Evolution (N Y). 1984;38: 1116–1129.
2. Mccandlish DM, Stoltzfus A. Modeling Evolution Using the Probability of Fixation: History and Implications. Q Rev Biol. 2014;89: 225–252. 
3. Maynard-Smith J. Natural Selection and the Concept of a Protein Space. Nature. 1970;225: 726–734. 
4. King JL, Jukes TH. Non-Darwinian evolution. Science. 1969. pp. 788–798. doi:10.1126/science.164.3881.788

```

The R session information (including the OS info, R version and all
packages used):


```r
sessionInfo()
```

```
## Warning: 'timedatectl' indicates the non-existent timezone name 'n/a'
```

```
## unable to deduce timezone name from 'America/Chicago'
```

```
## R version 4.3.1 (2023-06-16)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 8
## 
## Matrix products: default
## BLAS/LAPACK: /software/openblas-0.3.13-el8-x86_64/lib/libopenblas_skylakexp-r0.3.13.so;  LAPACK version 3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: NA
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] factoextra_1.0.7     cluster_2.1.4        reshape2_1.4.4      
##  [4] XNomial_1.0.4        pheatmap_1.0.12      ggradar_0.2         
##  [7] recommenderlab_1.0.4 proxy_0.4-27         arules_1.7-6        
## [10] igraph_1.5.0.1       purrr_1.0.1          furrr_0.3.1         
## [13] future_1.33.0        readr_2.1.4          Biostrings_2.68.1   
## [16] GenomeInfoDb_1.36.3  XVector_0.40.0       IRanges_2.34.1      
## [19] S4Vectors_0.38.1     BiocGenerics_0.46.0  matrixStats_1.0.0   
## [22] doParallel_1.0.17    iterators_1.0.14     foreach_1.5.2       
## [25] patchwork_1.1.2      dplyr_1.1.2          tibble_3.2.1        
## [28] stringr_1.5.0        Matrix_1.6-0         ggplot2_3.4.2       
## [31] MASS_7.3-60          tidyr_1.3.0          knitr_1.43          
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.2.0        farver_2.1.1            bitops_1.0-7           
##  [4] RCurl_1.98-1.12         recosystem_0.5.1        digest_0.6.33          
##  [7] lifecycle_1.0.3         magrittr_2.0.3          compiler_4.3.1         
## [10] rlang_1.1.1             tools_4.3.1             utf8_1.2.3             
## [13] ggsignif_0.6.4          labeling_0.4.2          plyr_1.8.8             
## [16] RColorBrewer_1.1-3      abind_1.4-5             registry_0.5-1         
## [19] withr_2.5.0             grid_4.3.1              fansi_1.0.4            
## [22] ggpubr_0.6.0            colorspace_2.1-0        globals_0.16.2         
## [25] scales_1.2.1            cli_3.6.1               crayon_1.5.2           
## [28] generics_0.1.3          rstudioapi_0.15.0       tzdb_0.4.0             
## [31] splines_4.3.1           zlibbioc_1.46.0         vctrs_0.6.3            
## [34] carData_3.0-5           car_3.1-2               hms_1.1.3              
## [37] rstatix_0.7.2           ggrepel_0.9.3           irlba_2.3.5.1          
## [40] listenv_0.9.0           glue_1.6.2              parallelly_1.36.0      
## [43] codetools_0.2-19        stringi_1.7.12          gtable_0.3.3           
## [46] munsell_0.5.0           pillar_1.9.0            float_0.3-1            
## [49] GenomeInfoDbData_1.2.10 R6_2.5.1                evaluate_0.21          
## [52] lattice_0.21-8          highr_0.10              backports_1.4.1        
## [55] broom_1.0.5             Rcpp_1.0.11             nlme_3.1-162           
## [58] mgcv_1.8-42             xfun_0.39               forcats_1.0.0          
## [61] pkgconfig_2.0.3
```

```r
Sys.time()
```

```
## [1] "2023-09-08 15:25:26 CDT"
```
