Evolutionary modeling on empirical GP maps
================
Santiago Herrera
2024-08-13

This notebook contains the analysis for modeling evolutionary trajectories on empirical GP maps using markov chains. It has three mains sections:

-   A theoretical explanation for modeling evolution on empirical GP maps [here](#specification-of-evolutionary-model)
-   Evolutionary modeling on the protein genotype networks [here](#evolutionary-modeling-using-discrete-markov-chains)
-   Evolutionary modeling on the protein-DNA genotype networks [here](#evolution-of-protein-dna-complexes)

------------------------------------------------------------------------

## Specification of evolutionary model

Under the strong selection-Weak mutation (SSWM) regime, the mutation rate is low enough such that the time to fixation of a mutation is much lower than the time between subsequent mutations (Gillespie, 1984). Thus, trajectories on a genotype landscape can be modeled as an stepwise, origin-fixation process. Specifically, the rate of fixation from allele *i* to *j* is a product of the rate of introduction of allele *j* in the population and the probability that it goes to fixation, like so

*Q*(*i*, *j*)=2*N*<sub>*e*</sub>*μ*<sub>*i**j*</sub> × *P*<sub>fix</sub>(*j*; *s*<sub>*i**j*</sub>, *N*<sub>*e*</sub>),

where *μ*<sub>*i**j*</sub> is the mutation rate from allele *i* to *j*, *s*<sub>*i**j*</sub> is the selection coefficient, and *N*<sub>*e*</sub> is the effective population size. When *s*<sub>*i**j*</sub> ≠ 0, the fixation probability, *P*<sub>fix</sub>(*j*; *s*<sub>*i**j*</sub>, *N*<sub>*e*</sub>) is given by the Kimura equation (Kimura, 1962):

$$
P\_{\\text{fix}}(j) = \\begin{cases}
\\frac{1-e^{-2s\_{ij}}}{1-e^{-4Ns\_{ij}}} & \\text{when } s\_{ij} \\neq 0 \\\\
\\frac{1}{2N\_e} & \\text{when } s\_{ij} = 0
\\end{cases}
$$

Note, however, that because we were interested in isolating the effect of the GP map’s structure on evolution, we considered a scenario in which the fixation process is unbiased, i.e., the fixation probability is only affected by drift. Under this scenario the selection coefficient *s*<sub>*i**j*</sub> = 0 and therefore *P*<sub>*f**i**x*</sub>(*j*)=1/2*N*<sub>*e*</sub> in every case.

The fitness landscape thus looks like this:

``` r
# STEP FUNCTION PARAMETERS:
# Initial parameters correspond to the same as the logistic function. 'mF_ref' = Minimal meanF for active variants
MIN_ACTIVE <- -4
STEP.PARAM = c(MIN_ACTIVE)

# create a vector of fluorescences
mFs <- seq(-4.5,-2.5,0.01)

# Step fitness function
data.frame(mFs = mFs, fitness = fitness_purifying(mFs,MIN_ACTIVE)) %>%
  ggplot(aes(x=mFs,y=fitness)) + ggtitle("Purifying sln + Drift") +
  annotate(geom = "rect", xmin =  -4.5, xmax = MIN_ACTIVE, ymin = 0, ymax = 3,
           fill = "#e37f74", alpha = 0.2) +
  annotate(geom = "rect", xmin =  MIN_ACTIVE, xmax = -2.5, ymin = 0, ymax = 3,
           fill = "#35a831", alpha = 0.2) +
  geom_line(linewidth=1.3) + theme_classic() + 
  xlab("Fluorescence (phenotype)") + ylab("Fitness") + 
  geom_vline(xintercept = -3.98,col="#7710b3",linetype="dashed",linewidth=1.2) + # reference phenotype and fitness
  annotate("text", x = -4.25, y = 2.5, label = "Non-Functional\ncomplexes") + 
  annotate("text", x = -3.5, y = 2.5, label = "Functional\ncomplexes") 
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/fitness_function-1.png)

To precisely model an origin-fixation walk on an empirical genotype-phenotype (GP) map, we need to compute *P*(*i*, *j*), the probability that the *next* mutation will be *j*. That is, we need to account for the local structure of the network around the focal node *i*. Thus,

$$
P(i,j) = \\frac{Q(i,j)}{ \\sum\_{k \\neq i} Q(i,k)} ,
$$

where *k* are all nodes connected to node *i* (McCandlish and Stoltzfus, 2014).

**Some assumptions:** Assuming a constant population size (*N*<sub>*e*</sub>), and no mutation bias and reversibility (*μ*<sub>*i**j*</sub> = *μ*<sub>*j**i*</sub> = *μ*<sub>*i**k*</sub>), the probability that the next mutation will be *j* becomes a function of the rescaled fixation probabilities alone:

$$
P(i,j) = \\frac{2N\_e\\mu\_{ij} \\times P\_{\\text{fix}}(j)}{\\sum\_{k \\neq i}2N\_e\\mu\_{ik} \\times P\_{\\text{fix}}(k)} \\
= (\\frac{2N\_e\\mu\_{ij}}{2N\_e\\mu\_{ik}}) \\times \\frac{P\_{\\text{fix}}(j)}{\\sum\_{k \\neq i} P\_{\\text{fix}}(k)} \\\\
= \\frac{P\_{\\text{fix}}(j)}{\\sum\_{k \\neq i}P\_{\\text{fix}}(k)}
$$

such that *P*(*i*, *j*) is a function of a single population-level parameter: the effective population size (*N*<sub>*e*</sub>). Since *P*(*i*, *j*) is a ratio of probabilities, the exact value of *N*<sub>*e*</sub> is unimportant; however, we chose a reasonable estimate of the average population size for chordates of 10<sup>6</sup> (Buffalo, 2021).

While we assumed no mutation bias at the nucleotide level, we did take into account the effect of the genetic code on mutation rates at the amino acid level. We weighted the fixation probability by a mutational propensity rate (*ρ*<sub>*i**j*</sub>) between alleles *i* and *j*. The mutational propensity rate represents the total number of mutational paths that exist between two protein sequences that differ at a single amino acid; it captures the mapping from codon-to-amino acid and its effect on the accessibility of genotype variants as mediated through the structure of the genetic code. We assume that a population fixed for a given amino acid genotype may explore all synonymous codons, i.e., the population is delocalized at the nucleotide level but fixed at the amino acid level. Thus, for each amino acid mutation in the genotype network, we computed ij as the number of nucleotide changes that can encode the amino acid change, while accounting for all the possible synonymous codon “backgrounds” at the other 3 invariant sites in which the mutation can occur:

*ρ*<sub>*i**j*</sub> = *η*<sub>*i*, *j*</sub><sup>*k*<sup>\*</sup></sup> × ∏<sub>*k* ≠ *k*<sup>\*</sup></sub>*c*<sub>*k*</sub>,

where *k* indexes the site in the amino acid sequence, *k*<sup>\*</sup> is the site at which the amino acid change occurs, *η*<sub>*i**j*</sub> is the number of single nucleotide changes between codons of amino acids *i* and *j*, and *c*<sub>*k*</sub> is the number of codons that encode each of the invariant amino acid states at the other sites (we tested a couple different ways of computing the mutational propensity, and the results are consistent between different specifications - see `scripts/mutation_matrix_comparison/matrix_comparison.md`).

Finally, we can rewrite the probability that the next mutation will be j including mutational propensities as:

$$
P(i,j) = \\frac{\\rho\_{ij} \\times P\_{\\text{fix}}(j)}{\\sum\_{k \\neq i}(\\rho\_{ik} \\times P\_{\\text{fix}}(k))}
$$

This specification of an origin-fixation model allows to model evolution on a genotype network as a discrete Markov process, where each time step consists of a single *amino acid* substitution.

### Discrete Markov process on a GP map

We were interested in understanding how likely each of the 16 DNA specificity phenotypes is to evolve. Conceptually, the process can be stated as follows: proteins traverse the sequence space by single amino acid substitutions, the genotype network determines the relative accessibility of an amino acid change, the fitness landscape defines its probability of fixation, and the GP map defines which specificity phenotype is encoded by each genotype. Thus, we can reformulate the problem in terms of the probability that a specific phenotypic outcome will happen: What is the probability of evolving any given specificity phenotype from a given starting genotype after a given number of substitutions, without losing function?

This process can be modeled by a discrete Markov chain, which can be fully specified with a transition probability matrix (*P*). The transition matrix *P* contains all possible *P*(*i*, *j*)’s between every pair of functional genotypes in the network. We specified the entries of the *P* matrix as follows:

$$
P(i,j) = \\begin{cases}
\\frac{\\rho\_{ij} \\times P\_{\\text{fix}}(j)}{\\sum\_{k \\neq i}(\\rho\_{ik} \\times P\_{\\text{fix}}(k))} & \\text{when } d = 1 \\\\
0 & \\text{when } d &gt;  0 \\\\
0 & \\text{when } d = 0
\\end{cases}
$$

where *d* is the number of edges between two genotypes in the network. Thus, *P*(*i*, *j*) is non-zero only when the genotypes are single-step neighbors. To ensure that every time step in the Markov chain represents an amino acid substitution we further constrained the probability of staying in the same genotype *P*(*i*, *i*) to zero (i.e., *d* = 0).

To compute the probability of evolving *R**E*<sub>*i*</sub>, we must first obtain the probability distribution over the complete set of functional genotypes in the network after S steps of the Markov chain. The vector containing the probabilities of the possible realizations of the process across every genotype is:

*π*<sub>(*S*)</sub> = *π*<sub>(0)</sub> × *P*<sup>*S*</sup>

where *P* is the transition matrix, *S* &gt; 0, and *π*<sub>(0)</sub> is the vector of state frequencies at time step *S* = 0. If the Markov chain starts from a single genotype *g*<sub>*i*</sub>, we set the entry *f*(*g*<sub>*i*</sub>) in the vector *π*<sub>(0)</sub> to 1 and the rest of the states have an initial frequency of zero. If the Markov chain starts from a set of genotypes \[*g*\], the initial frequencies across g are drawn from a uniform distribution, and the rest of the states have an initial frequency of zero.

Finally, we can use the experimental data to assign a phenotype to each probability state of *π*<sub>(*S*)</sub>. The conditional probability of evolving REi expressed in terms of the Markov chain becomes:

$$
P(RE\_i|\\pi\_{(0)},S,P) = \\frac{\\sum\_{j \\in RE\_i} \\pi\_{(S)j}}{\\sum\_{i=1}^{k} (\\sum\_{j \\in RE\_i} \\pi\_{(S)j})\_i}
$$

where the numerator is the sum of the probailities of all the genotypes encoding the DNA binding function *R**E*<sub>*i*</sub> after *S* substitutions, and the denominator is a normalization constant with *k* = 16 such that ∑<sub>*i*</sub>*P*(*R**E*<sub>*i*</sub>|*π*<sub>(0)</sub>, *S*, *P*)=1.

For *k* different specificity phenotypes, we obtain a multinomial probability distribution around any set of starting genotypes that quantifies the likelihood that evolution on the genotype network will produce a particular phenotypic outcome given a set of conditions – i.e., the spectrum of evolutionary outcomes.

### References

1.  Gillespie J. Molecular Evolution Over the Mutational Landscape. Evolution (N Y). 1984;38: 1116–1129.
2.  Mccandlish DM, Stoltzfus A. Modeling Evolution Using the Probability of Fixation: History and Implications. Q Rev Biol. 2014;89: 225–252.
3.  Maynard-Smith J. Natural Selection and the Concept of a Protein Space. Nature. 1970;225: 726–734.
4.  King JL, Jukes TH. Non-Darwinian evolution. Science. 1969. pp. 788–798. <doi:10.1126/science.164.3881.788>

------------------------------------------------------------------------

## Evolutionary modeling using discrete markov chains

First, we will build the transition probability matrices for the protein genotypes networks in AncSR1 and AncSR2 backgrounds. As a phenotype, we will use the maximum mean fluorescence across the DNA elements bound by each variant protein variant.

``` r
# Check whether the matrices have been created already
if(!file.exists(file.path("..","..","results","evolutionary_modeling","MutSel_matrices_final.RData"))) {
  # Load complete data from mutation effects model
  meanF_data <- readr::read_csv(file.path("..","..","results","classifying_functional_variants","meanF_data_fxnal.csv.gz"))
  
  # BUILD PHENOTYPE TABLE #
  # For each *functional* variant (those with meanF >= AncSR1_ERE_ref or meanF >= AncSR2_SRE_ref):
  phenotypes_tbl <- meanF_data %>% filter(functional==T) %>% group_by(AA_var, bg) %>%
    reframe(n_bound_REs = n(), # how many DNA elements can bind
              meanF_bREs = mean(avg_meanF), # the average meanF across all DNA elements bound
              max_meanF_bREs = max(avg_meanF), # the max meanF across all DNA elements bound
              min_meanF_bREs = min(avg_meanF), # the min meanF across all DNA elements bound 
              specific = ifelse(n_bound_REs == 1, "YES","NO"), # functional binding to only one DNA element?
              specificity = ifelse(n_bound_REs == 1, RE, "Promiscuous"), # Determine the RE specificity
              promiscuity =  ifelse(specific=="NO",n_bound_REs * max_meanF_bREs/sum(avg_meanF),0), # Level of promiscuity: promiscuous (1), specific (0)
              bound_REs = list(RE)) # assign DNA elements bound
  
  # BUILD GENOTYPE NETWORKS #
  # select the set of funcitonal variants for each DBD background
  sr1_variants <- phenotypes_tbl %>% filter(bg == "AncSR1") %>% pull(AA_var)
  sr2_variants <- phenotypes_tbl %>% filter(bg == "AncSR2") %>% pull(AA_var)
  
  # buld genotype netwotk under Maynard-Smith's model of sequence space
  net_sr1 <- build_genotype_network(nodes=sr1_variants, type=1, cores=N_CORES)
  net_sr2 <- build_genotype_network(nodes=sr2_variants, type=1, cores=N_CORES)
  
  # BUILD TRANSITION MATRICES #
  Ne = 10e6 # population size (reasonable average for chordates e.g.,https://doi.org/10.7554/eLife.67509)
  N_CORES=detectCores()-1 # number of cores
  
  # Set-up parameters for the fitness function:
  # Fitness corresponds to the exponential growth rate = N*r
  MIN_ACTIVE <- phenotypes_tbl %>% with(min(min_meanF_bREs))
  STEP.PARAM = c(MIN_ACTIVE)

  # Build model scenarios for each DBD background:
  MODEL.PARAM_SR1_drift <- list("AncSR1","drift",Ne,STEP.PARAM,TRUE,"max")
  MODEL.PARAM_SR2_drift <- list("AncSR2","drift",Ne,STEP.PARAM,TRUE,"max")
  
  # Build transition matrices
  adj_mat_count_sr1 <- build_mutation_matrix(sr1_variants,type=3,N_CORES) # adjacency matrix with mutational propensties
  M_drift_sr1 <- build_transition_matrix_v2(sr1_variants,adj_mat_count_sr1,phenotypes_tbl,MODEL.PARAM_SR1_drift,N_CORES,complex = F) # compute P(i,j) for each connected pair
  
  adj_mat_count_sr2 <- build_mutation_matrix(sr2_variants,type=3,N_CORES) # adjacency matrix with mutational propensties
  M_drift_sr2 <- build_transition_matrix_v2(sr2_variants,adj_mat_count_sr2,phenotypes_tbl,MODEL.PARAM_SR2_drift,N_CORES,complex = F) # compute P(i,j) for each connected pair
  
} else {
  # load matrices if already created
  load(file.path("..","..","results","evolutionary_modeling","MutSel_matrices_final.RData"))
  phenotypes_tbl_prot <- phenotypes_tbl # save it for later
  rm(phenotypes_tbl) # clear up space
  N_CORES=detectCores()-1 # re-set cores if running locally
  load(file.path("..","..","results","evolutionary_modeling","GPmap_evol_model_sr1.RData")) # load modeling results for AncSR1
  load(file.path("..","..","results","evolutionary_modeling","GPmap_evol_model_sr2.RData")) # load modeling results for AncSR2
}
```

Let's check how many genotypes are part of the main network component in each background

``` r
knitr::kable(rbind(data.frame(Bg="AncSR1",total_fxnal_vars = length(phenotypes_tbl_prot %>% filter(bg == "AncSR1") %>% pull(AA_var)), fxnal_vars_in_net = length(extract_main_ntwrk(net_sr1,nodes = T,tr_mat = NULL))),
                   data.frame(Bg="AncSR2",total_fxnal_vars = length(phenotypes_tbl_prot %>% filter(bg == "AncSR2") %>% pull(AA_var)), fxnal_vars_in_net = length(extract_main_ntwrk(net_sr2,nodes = T,tr_mat = NULL)))))
```

| Bg     |  total\_fxnal\_vars|  fxnal\_vars\_in\_net|
|:-------|-------------------:|---------------------:|
| AncSR1 |                 107|                    97|
| AncSR2 |                2407|                  2402|

Since some genotypes in the network are disconnected, and would not be accessible via single steps, we will only consider the genotypes in the biggest component of the genotype networks (and the corresponding *P* matrices) to model the evolutionary process. This ensures that the Markov chain is irreducible – that is, every state is accessible from every other state – and therefore has a unique stationary distribution.

``` r
## Extract main components of each network and build the corresponding square P matrices
P_drift_sr1_ntwrk <- extract_main_ntwrk(net_sr1,M_drift_sr1) # square P matrix
P_drift_sr1_ntwrk_statdist <- stationary_dist(P_drift_sr1_ntwrk) # compute stationary distribution

P_drift_sr2_ntwrk <- extract_main_ntwrk(net_sr2,M_drift_sr2) # square P matrix
P_drift_sr2_ntwrk_statdist <- stationary_dist(P_drift_sr2_ntwrk) # compute stationary distribution
```

We will now compute the production spectrum, equilibrium outcome spectrum, global bias and equilibrium bias for each GP map.

``` r
## GLOBAL OBJECTS ##

REF_GENOTYPE = "EGKA" # EGKA (historical genotype) 

# Global production spectra (PDFV = Probability Distribution of Functional Variation)
var.prop_AncSR1_df <- get_PDFV_v2(type="production",Bg="AncSR1",model="production_sr1",pheno_tbl = phenotypes_tbl_prot, specific = F) # Binding
var.prop_spec_AncSR1_df <- get_PDFV_v2(type="production",Bg="AncSR1",model="production_sr1",pheno_tbl = phenotypes_tbl_prot, specific = T) # Spec. Binding

var.prop_AncSR2_df <- get_PDFV_v2(type="production",Bg="AncSR2",model="production_sr2",pheno_tbl = phenotypes_tbl_prot, specific = F) # Binding
var.prop_spec_AncSR2_df <- get_PDFV_v2(type="production",Bg="AncSR2",model="production_sr2",pheno_tbl = phenotypes_tbl_prot, specific = T) # Spec. Binding

# Equilibrium outcome spectra
stationary_PDFV_sr1 <- get_PDFV_v2(P_drift_sr1_ntwrk_statdist,type="simulated mc",Bg="AncSR1",model="Stationary",pheno_tbl = phenotypes_tbl_prot, specific = F) # Binding
stationary_PDFV_spec_sr1 <- get_PDFV_v2(P_drift_sr1_ntwrk_statdist,type="simulated mc",Bg="AncSR1",model="Stationary", pheno_tbl = phenotypes_tbl_prot, specific = T) # Specific Binding

stationary_PDFV_sr2 <- get_PDFV_v2(P_drift_sr2_ntwrk_statdist,type="simulated mc",Bg="AncSR2",model="Stationary",pheno_tbl = phenotypes_tbl_prot, specific = F) # Binding
stationary_PDFV_spec_sr2 <- get_PDFV_v2(P_drift_sr2_ntwrk_statdist,type="simulated mc",Bg="AncSR2",model="Stationary",pheno_tbl = phenotypes_tbl_prot, specific = T) # Specific Binding

# Phenotype frequencies in main network
var.prop_AncSR1_net_df <- get_PDFV_v2(type="main network",Bg="AncSR1",model="Binding Network",graph = net_sr1,pheno_tbl = phenotypes_tbl_prot, specific = F) # Binding 
var.prop_AncSR1_net_spec_df <- get_PDFV_v2(type="main network",Bg="AncSR1",model="Binding Network",graph = net_sr1, pheno_tbl = phenotypes_tbl_prot, specific = T) # Specific Binding

var.prop_AncSR2_net_df <- get_PDFV_v2(type="main network",Bg="AncSR2",model="Binding Network",graph = net_sr2,pheno_tbl = phenotypes_tbl_prot, specific = F) # Binding
var.prop_AncSR2_net_spec_df <- get_PDFV_v2(type="main network",Bg="AncSR2",model="Binding Network",graph = net_sr2, pheno_tbl = phenotypes_tbl_prot, specific = T) # Specific Binding

# Global bias
GlobalB_sr1 <- get_bias(var.prop_spec_AncSR1_df %>% remove_promiscuous(norm = T),Bg="AncSR1",input_freqs = TRUE)
GlobalB_sr2 <- get_bias(var.prop_spec_AncSR2_df %>% remove_promiscuous(norm = T),Bg="AncSR2",input_freqs = TRUE)

# Equilibrium bias
StatB_sr1 <- get_bias(stationary_PDFV_spec_sr1 %>% remove_promiscuous(norm = T),Bg="AncSR1",input_freqs = TRUE)
StatB_sr2 <- get_bias(stationary_PDFV_spec_sr2 %>% remove_promiscuous(norm = T),Bg="AncSR2",input_freqs = TRUE)

# Summary tables: Whether each phenotype is encoded in the GP map, whether it is encoded by specific variants, whether it is in the main network component, 
# and whether it is encoded by specific variants in the main network component.
# AncSR1 summary
summary_sr1 <- foreach(i = 1:length(REs[[1]]), .combine = 'rbind') %do% {
  encoded <- is.encoded(REs[[1]][i],graph = net_sr1,pheno_df = phenotypes_tbl_prot,Bg = "AncSR1")
  specific <- is.bound.specific(REs[[1]][i],graph = net_sr1,pheno_df = phenotypes_tbl_prot,Bg = "AncSR1")
  in_net <- is.in.ntwrk(REs[[1]][i],graph = net_sr1,pheno_df = phenotypes_tbl_prot,Bg = "AncSR1")
  in_net_spec <- is.in.ntwrk(REs[[1]][i],graph = net_sr1,pheno_df = phenotypes_tbl_prot,Bg = "AncSR1",specific = T)
  v <- c(encoded,specific,in_net,in_net_spec)
} %>% `rownames<-`(REs[[1]]) %>% `colnames<-`(c("Encoded","Specific","In_Net","In_Net_Spec"))

# AncSR2 summary
summary_sr2 <- foreach(i = 1:length(REs[[1]]), .combine = 'rbind') %do% {
  encoded <- is.encoded(REs[[1]][i],graph = net_sr2,pheno_df = phenotypes_tbl_prot,Bg = "AncSR2")
  specific <- is.bound.specific(REs[[1]][i],graph = net_sr2,pheno_df = phenotypes_tbl_prot,Bg = "AncSR2")
  in_net <- is.in.ntwrk(REs[[1]][i], graph = net_sr2,pheno_df = phenotypes_tbl_prot,Bg = "AncSR2")
  in_net_spec <- is.in.ntwrk(REs[[1]][i],graph = net_sr2,pheno_df = phenotypes_tbl_prot,Bg = "AncSR2",specific = T)
  v <- c(encoded,specific,in_net,in_net_spec)
} %>% `rownames<-`(REs[[1]]) %>% `colnames<-`(c("Encoded","Specific","In_Net","In_Net_Spec")) 


############
# Other useful global objects

# Functional variants per background
func_vars_sr1 <- phenotypes_tbl_prot %>% filter(bg == "AncSR1") %>% pull(AA_var)
func_vars_sr2 <- phenotypes_tbl_prot %>% filter(bg == "AncSR2") %>% pull(AA_var)

# Nodes in main network component
nodes_in_ntwrk_sr1 <- extract_main_ntwrk(net_sr1,nodes=T,tr_mat = NULL)
nodes_in_ntwrk_sr2 <- extract_main_ntwrk(net_sr2,nodes=T,tr_mat = NULL)

# Nodes in main network component that are specific
nodes_in_ntwrk_sr1_spec <- data.frame(AA_var=nodes_in_ntwrk_sr1,bg="AncSR1") %>% inner_join(.,phenotypes_tbl_prot,by=c("bg","AA_var")) %>% filter(specific=="YES") %>% pull(AA_var)
nodes_in_ntwrk_sr1_spec <- nodes_in_ntwrk_sr1_spec[nodes_in_ntwrk_sr1_spec %in% nodes_in_ntwrk_sr1]
nodes_in_ntwrk_sr2_spec <- data.frame(AA_var=nodes_in_ntwrk_sr2,bg="AncSR2") %>% inner_join(.,phenotypes_tbl_prot,by=c("bg","AA_var")) %>% filter(specific=="YES") %>% pull(AA_var)
nodes_in_ntwrk_sr2_spec <- nodes_in_ntwrk_sr2_spec[nodes_in_ntwrk_sr2_spec %in% nodes_in_ntwrk_sr2]

# arrange phenotypes in descending order according to the production spectrum of each background
REs[[3]] <- factor(rownames(summary_sr1),levels = as.character(remove_promiscuous(var.prop_spec_AncSR1_df) %>% arrange(desc(Norm_F_prob)) %>% pull(RE)))
REs[[4]] <- factor(rownames(summary_sr2),levels = as.character(remove_promiscuous(var.prop_spec_AncSR2_df) %>% arrange(desc(Norm_F_prob)) %>% pull(RE)))

REs[[5]] <- factor(c(rownames(summary_sr1),"Promiscuous"),
                     levels = c(as.character(remove_promiscuous(var.prop_spec_AncSR1_df) %>% arrange(desc(Norm_F_prob)) %>% pull(RE)),"Promiscuous"))
REs[[6]] <- factor(c(rownames(summary_sr2),"Promiscuous"),
                     levels = c(as.character(remove_promiscuous(var.prop_spec_AncSR2_df) %>% arrange(desc(Norm_F_prob)) %>% pull(RE)),"Promiscuous"))
```

------------------------------------------------------------------------

## AncSR1 GP map

### Structure of the AncSR1 GP map

We will first investigate the structure of the AncSR1 GP map. The GP map would have no effect in evolutionary outcomes only if it were isotropic (producing all possible phenotypes with uniform frequency) and homogeneous (producing the same distribution of phenotypes from each possible starting genotype). We will now assess each of these two aspects of the map.

Isotropic production of phenotypic variation means that all 16 phenotypes should be equally frequent (6.25%), and the bias would be *B* = 0. Deviations from this expectation would imply that production is anisotropic. A maximally anisotropic GPmap - where only a single phenotype is encoded - would have a bias of *B* = 1. Let's see how does the AncSR1's global production spectrum look like.

``` r
####################
# AncSR1 global production spectrum
###################
SPECIFIC = T

if(SPECIFIC){
  # Only specific binding
  global_spec_sr1_plot <- phenotypes_tbl_prot %>% ungroup() %>% filter(bg=="AncSR1") %>% filter(specific == "YES") %>% 
    reframe(RE = REs[[3]],
            count = table(factor(specificity,levels=REs[[3]]))) %>% 
    ggplot(aes(x=RE,y=count,fill=RE)) + 
    geom_bar(stat="identity",position = "dodge",color="black") +
    scale_fill_manual(values = hex_RE_colors(1)) + 
    theme_classic() + 
    labs(x="DNA element",y="Number of\nprotein variants") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size=11)) +
    guides(fill="none") +
    geom_hline(yintercept = round(91*(1/16)),linetype="dashed",col="black")
} else{
  # With promiscuous
    global_spec_sr1_plot <- phenotypes_tbl_prot %>% filter(bg=="AncSR1") %>% unnest(bound_REs) %>% group_by(specific) %>% 
      reframe(RE = REs[[3]],
              count = table(factor(bound_REs,levels=REs[[3]]))) %>% 
    ggplot(aes(x=RE,y=count,fill=RE,alpha=specific)) + 
    geom_bar(stat="identity",color="black") +
    scale_fill_manual(values = hex_RE_colors(1)) + 
    scale_alpha_manual(values =c(0,1)) +
    theme_classic() + 
    labs(x="DNA element",y="Number of\nprotein variants") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size=11)) +
    guides(fill="none") +
    geom_hline(yintercept = round(107*(1/16)),linetype="dashed",col="black")
}

global_spec_sr1_plot + annotate("text",x=3,y=30,label = paste("B =",round(GlobalB_sr1,2)))
```

    ## Don't know how to automatically pick scale for object of type <table>.
    ## Defaulting to continuous.

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/gpmap_sr1_structure_1-1.png)

We can see that the AncSR1 GP map is strongly anisotropic (with a bias of *B* = 0.42). ERE specificitiy is the most common phenotype, followed by SRE. And 7 specificity phenotypes are never encoded by any protein variant. We call this nonuniform production of variation global bias.

Global bias implies that the AncSR1 protein background has an intrinsic propensity to produce certain phenotypes more often by mutations in the RH. But, since evolution is a genetic process, we also need to account for the connectivity of functional variants for that would affect the mutational access to different phenotypes. Let's first investigate how the functional variants are distributed in the space of 160,000 possible RH variants.

``` r
####################
# Logo plots functional variants
###################

func_vars_sr1_plot <- list(func_vars_sr1)
names(func_vars_sr1_plot) <- "Functional variants"
ggseqlogo::ggseqlogo(func_vars_sr1,ncol=3,method='p',seq_type='aa') + theme(legend.position = "none")
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/gpmap_sr1_structure_2-1.png)

``` r
####################
# Connectivity of AncSR1 network
####################

# Compare the connectivity of the AcSR1 network to that of random networks with the same number of variants 
cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

# number of functional genotypes in AncSR1
N_VAR_SR1 <- phenotypes_tbl_prot %>% filter(bg=="AncSR1") %>% pull(AA_var) %>% length()

set.seed(1246)
# Random networks: Summary statistics for 200 random networks of the same size as AncSR1 network
connect_sr1_null <- foreach(i = 1:200, .combine = 'rbind', .errorhandling="remove") %dopar% {
  net_tmp <- simulate_GPmap(type=4,cores=N_CORES,graph = NULL,n_sample = N_VAR_SR1) %>% as.undirected()
  n_edges <- length(E(net_tmp))
  size_main_ntwrk <- ifelse(n_edges==0,1,length(extract_main_ntwrk(net_tmp,nodes=T,tr_mat = NULL)))
  data.frame(bg="AncSR1",n_nodes=N_VAR_SR1,size_main_ntwrk=size_main_ntwrk,n_edges=n_edges,connect=mean(degree(net_tmp)))
}

stopCluster(cl)
```

``` r
# Non-parametric bootstrap test: p-value is the fraction of random networks with a connectivity as high as the one observed in the AncSR1 network.
obs_connectivity <- mean(degree(net_sr1))

#descdist(connect_sr1_null$connect, discrete = F) # possible candidate distributions for permuted dataset
lnorm_fit_permut <- fitdistrplus::fitdist(connect_sr1_null$connect+1, distr = "lnorm",discrete = F) # fit log-normal distribution to permuted dataset
#hist(connect_sr1_null$connect+1,xlim = c(1,1.15),freq = F)
#lines(x=seq(1,1.15,length.out = 100),y=dlnorm(x=seq(1,1.15,length.out = 100), meanlog = lnorm_fit_permut$estimate[1], sdlog = lnorm_fit_permut$estimate[2]),col="red")

pval_connect <- plnorm(q = obs_connectivity, meanlog =  lnorm_fit_permut$estimate[1],sdlog = lnorm_fit_permut$estimate[2],lower.tail = F)
print(paste("Observed connectivity of AncSR1 network: ", round(obs_connectivity,2),". P-value = ",pval_connect))
```

    ## [1] "Observed connectivity of AncSR1 network:  3.5 . P-value =  0"

These results show that the 107 functional nodes are more similar to each other than expected; the set of functional gentoypes is clumped in the possible sequence space, with an average of 3.5 edges per node, much higher than expected if they were scattered randomly in sequence space.

Now, we will look at the distribution of phenotypes in the network. For this, we will identify genotype clusters in the main network component and then compute the frequency of each phenotype within each cluster.

``` r
####################
# Detect genotype clusters in network
####################

# Genotype clusters
cl_sr1 <- cluster_edge_betweenness(net_sr1)

# Associate individual genotypes to each cluster
modules_sr1 <- data.frame(AA_var=names(membership(cl_sr1)),cluster=as.integer(membership(cl_sr1)),bg="AncSR1") %>%
  inner_join(.,phenotypes_tbl_prot,by=c("bg","AA_var")) %>% select(bg,AA_var,cluster,specificity,bound_REs) %>%
  mutate(module = case_when(cluster == 1 ~ "A",
                            cluster == 2 ~ "B",
                            cluster == 4 ~ "C",
                            cluster == 5 ~ "D",
                            cluster == 6 ~ "E",
                            cluster == 9 ~ "F")) %>%
  inner_join(.,data.frame(AA_var=nodes_in_ntwrk_sr1,type="net"),by="AA_var") %>% arrange(cluster) %>% as_tibble()

# Frequency of phenotypes per cluster
phenotypes_sr1 <- data.frame(summary_sr1) %>% rownames_to_column(var = "RE") %>% filter(Encoded == TRUE & Specific == TRUE) %>% pull(RE)
pheno_freqs_modules_sr1 <- modules_sr1 %>% group_by(module) %>% 
  reframe(RE = phenotypes_sr1, 
          RE = factor(RE, levels(REs[[3]])),
          count_t = n(),
          count_p = table(factor(specificity,levels=phenotypes_sr1)),
          freq_p = count_p / count_t) %>% acast(module~RE, value.var="freq_p")

# Compare to the global production spectrum
global_freqs_sr1 <- remove_promiscuous(var.prop_spec_AncSR1_df,norm = T) %>% mutate(model="Global") %>% 
  filter(Norm_F_prob > 0) %>% mutate(RE = factor(RE,levels(REs[[3]]))) %>% acast(model~RE, value.var="Norm_F_prob")

pheno_freqs_modules_sr1 <- rbind(pheno_freqs_modules_sr1,global_freqs_sr1) %>% as.matrix(.)
pheno_freqs_modules_sr1 <- t(apply(pheno_freqs_modules_sr1, 1, function(x) x/sum(x))) # renormalize frequencies within each cluster

# Detect enrichment of phenotypes in each cluster using Fisher's exact test
module_sr1_ids <- unique(modules_sr1$module)
REs_in_net_sr1 <- unique(modules_sr1 %>% filter(specificity != "Promiscuous") %>% pull(specificity)) 

m_enrichment <- matrix(NA,ncol=length(REs_in_net_sr1),nrow=length(module_sr1_ids))
rownames(m_enrichment) <- module_sr1_ids
colnames(m_enrichment) <- REs_in_net_sr1

for(i in 1:length(module_sr1_ids)){
  #modules_sr1 <- modules_sr1 %>% unnest(bound_REs)
  for(j in 1:length(REs_in_net_sr1)){
    a <- modules_sr1 %>% filter(module == module_sr1_ids[i] & specificity == REs_in_net_sr1[j]) %>% reframe(c=n()) %>% pull(c) # RE j in module i
    b <- modules_sr1 %>% filter(module == module_sr1_ids[i] & specificity != REs_in_net_sr1[j]) %>% reframe(c=n()) %>% pull(c) # Other REs in module i
    c <- modules_sr1 %>% filter(module != module_sr1_ids[i] & specificity == REs_in_net_sr1[j]) %>% reframe(c=n()) %>% pull(c) # RE j in other modules 
    d <- modules_sr1 %>% filter(module != module_sr1_ids[i] & specificity != REs_in_net_sr1[j]) %>% reframe(c=n()) %>% pull(c) # Other REs in other modules
    
    # perform fisher's exact test on counts
    x <- fisher.test( matrix(c(a,b,c,d),ncol=2,nrow=2,byrow = T),alternative = "greater")
    m_enrichment[i,j] <- x$p.value
    
    # Perform fsher's exact test using hypergeometric distribution (one-sided p): see 'https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/'
    #m <- a + b      # Nodes in cluster i
    #n <- c + d      # Nodes NOT IN cluster i
    #k <- a + c      # Node hits, that is, from RE j

    #m_enrichment[i,j] <- phyper(a, m, n, k, lower.tail = FALSE,log.p = FALSE) # probability of observing 'a' nodes or more with REj in module i
  }
}

# correct for multiple testting (Bonferroni)
m_enrichment <- matrix(p.adjust(m_enrichment,method="bonferroni"),nrow = dim(m_enrichment)[1],
                       ncol = dim(m_enrichment)[2],byrow = F,dimnames = list(rownames(m_enrichment),colnames(m_enrichment)))
m_enrichment <- m_enrichment[,order(match(colnames(m_enrichment),REs[[1]]))] # order columns

knitr::kable(as.tibble(m_enrichment,rownames=NA) %>% rownames_to_column(var="Cluster") %>%
               rowwise() %>% mutate(n_enriched_phenotypes = sum(c_across(`SRE (AA)`:TA)<0.05),
                                    enriched_phenotypes = toString(names(pick(`SRE (AA)` : TA))[c_across(`SRE (AA)` : TA) < 0.05])))
```

    ## Warning: `as.tibble()` was deprecated in tibble 2.0.0.
    ## ℹ Please use `as_tibble()` instead.
    ## ℹ The signature and semantics have changed, see `?as_tibble`.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

| Cluster |  SRE (AA)|         GA|   ERE (GT)|   AC|       AT|   TA|  n\_enriched\_phenotypes| enriched\_phenotypes |
|:--------|---------:|----------:|----------:|----:|--------:|----:|------------------------:|:---------------------|
| A       |         1|  1.0000000|  1.0000000|    1|  1.0e+00|    0|                        1| TA                   |
| B       |         0|  1.0000000|  1.0000000|    1|  1.0e+00|    1|                        1| SRE (AA)             |
| C       |         1|  0.0007081|  1.0000000|    1|  1.0e+00|    1|                        1| GA                   |
| D       |         1|  1.0000000|  0.5124773|    1|  1.0e+00|    1|                        0|                      |
| E       |         1|  1.0000000|  1.0000000|    1|  9.5e-06|    1|                        1| AT                   |
| F       |         1|  1.0000000|  0.0000128|    1|  1.0e+00|    1|                        1| ERE (GT)             |

The table above shows the corrected p-values for the enrichment test. 5 out of 6 clusters are significantly enriched for a single specificity phenotype (*p* &lt; 0.05). Now let's visualize the frequencies of each phentoype within each cluster.

``` r
####################
# Phenotype distribution per cluster and bias
###################

# compute cluster bias
c_bias <- foreach(i = 1:dim(pheno_freqs_modules_sr1)[1], .combine = 'c') %do% {
  round(get_bias(pheno_freqs_modules_sr1[i,],Bg = "AncSR1",input_freqs = TRUE),2)
}

# phenotype frequencies per cluster
bp_cl <- as_tibble(pheno_freqs_modules_sr1,rownames = NA) %>% rownames_to_column(var="cluster") %>%
  reshape2::melt() %>% 
  arrange(value) %>% 
  ggplot(aes(x=cluster,y=value,fill=variable)) + geom_bar(stat="identity",width=0.95) + 
  scale_fill_manual(values = hex_RE_colors(1)) +
  theme_classic() +
  labs(x="Genotype cluster",y="Phenotype frequency",fill="") +
  theme(axis.title = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13,angle=45,hjust = 1,vjust = 1),
        legend.position="bottom")
```

    ## Using cluster as id variables

``` r
# strength of bias per cluster
b_cl <- data.frame(cl = rownames(pheno_freqs_modules_sr1),bias = c_bias) %>%
  ggplot(aes(x = cl,y=bias)) + geom_bar(stat="identity",width=0.95,fill="gray70") +
  theme_classic() + labs(y="Bias (B)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank()) +
  geom_text(aes(label=bias), vjust=1.5, size=2.5)

(b_cl / bp_cl) + plot_layout(ncol=1,nrow=2,heights = c(1,4))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
####################
# Logo plots per cluster
###################

vars_clusterA <- modules_sr1 %>% filter(module == "A") %>% pull(AA_var)
vars_clusterB <- modules_sr1 %>% filter(module == "B") %>% pull(AA_var)
vars_clusterC <- modules_sr1 %>% filter(module == "C") %>% pull(AA_var)
vars_clusterD <- modules_sr1 %>% filter(module == "D") %>% pull(AA_var)
vars_clusterE <- modules_sr1 %>% filter(module == "E") %>% pull(AA_var)
vars_clusterF <- modules_sr1 %>% filter(module == "F") %>% pull(AA_var)

logo_vars <- list(vars_clusterA,vars_clusterB,vars_clusterC,vars_clusterD,vars_clusterE,vars_clusterF)
names(logo_vars) <- c("A","B","C","D","E","F")

cl_logo <- ggseqlogo::ggseqlogo(logo_vars,ncol=3,method='p',seq_type='aa') + theme(legend.position = "none")

####################
# Logo plots per phenotype
###################

vars_SRE <- modules_sr1 %>% filter(specificity == "SRE (AA)") %>% pull(AA_var)
vars_GA <- modules_sr1 %>% filter(specificity == "GA") %>% pull(AA_var)
vars_ERE <- modules_sr1 %>% filter(specificity == "ERE (GT)") %>% pull(AA_var)
vars_AC <- modules_sr1 %>% filter(specificity == "AC") %>% pull(AA_var)
vars_AT <- modules_sr1 %>% filter(specificity == "AT") %>% pull(AA_var)
vars_TA <- modules_sr1 %>% filter(specificity == "TA") %>% pull(AA_var)

logo_vars <- list(vars_ERE,vars_SRE,vars_GA,vars_AC,vars_AT,vars_TA)
names(logo_vars) <- c("ERE","SRE","GA","AC","AT","TA")

pheno_logo <- ggseqlogo::ggseqlogo(logo_vars,ncol=3,method='p',seq_type='aa') + theme(legend.position = "none")

cl_logo + pheno_logo
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-6-1.png)

We can see that clusters are strongly biased, more than the global production spectrum, and that the enriched phenotypes differ between clusters. As a result, amino acid variants enriched within each cluster also correspond with the genetic determinants of specificity. These results show that the map is strongly heterogeneous: different regions of the genotype network have different phenotype distributions, a phenomenon we call local bias.

Now, let's get a finer view of the local bias, by assessing the distribution of phenotypes accessible to every variant within its one-mutational neighborhood.

``` r
###################################
# local bias: One-mutation neighborhood
###################################

cl <- makeCluster(N_CORES,type = "FORK",outfile="")
registerDoParallel(cl)

# Number of phenotypes accessible in the one-mutational neighborhood
one_mut_nei_sr1 <- foreach(i = 1:length(nodes_in_ntwrk_sr1_spec), .combine = 'rbind') %dopar% {
  # focal variant and phenotype
  f_var <- phenotypes_tbl_prot %>% filter(bg=="AncSR1") %>% filter(AA_var %in% nodes_in_ntwrk_sr1_spec[i]) %>% select(AA_var,specificity)
  
  # neighbor variants and phenotype
  n_i <- names(neighbors(graph = net_sr1,v = nodes_in_ntwrk_sr1_spec[i]))
  p_n_i <- phenotypes_tbl_prot %>% ungroup() %>% filter(bg=="AncSR1") %>% filter(AA_var %in% n_i) %>% select(AA_var,specificity) %>%
    mutate(focal_var = f_var$AA_var, focal_pheno = f_var$specificity)
}

# local bias computation of the one-mutational neighborhood
local_bias_sr1 <- foreach(i = 1:length(nodes_in_ntwrk_sr1_spec), .combine = 'rbind') %dopar% {
  # extract neighbors of focal genotype
  n_i <- names(neighbors(graph = net_sr1,v = nodes_in_ntwrk_sr1_spec[i]))
  # frequencies of phenotypes among one-mutant neighbors
  p_n_i <- phenotypes_tbl_prot %>% ungroup() %>% filter(bg=="AncSR1") %>% filter(AA_var %in% n_i) %>%
    reframe(RE = REs[[1]],
            count = table(factor(specificity,levels=REs[[1]])),
            fr = count/length(n_i))
  # local bias
  local_b_i <- get_bias(as.vector(p_n_i$fr),input_freqs = T)
  data.frame(g=nodes_in_ntwrk_sr1_spec[i],b=local_b_i)}

stopCluster(cl)
```

``` r
# Number of accessible and *new* phenotypes per genotype
accessible_phenotypes_singleStep_sr1 <- one_mut_nei_sr1 %>% filter(specificity != "Promiscuous") %>% 
  mutate(new_pheno = ifelse(focal_pheno == specificity,FALSE,TRUE)) %>%
  group_by(focal_var) %>% 
  reframe(n_acc_pheno = n_distinct(specificity),
          n_new_pheno = n_distinct(specificity) - any(specificity == focal_pheno)) %>% 
  ungroup() 

# plot number of new phenotypes accessible per specific RH variant   
p1 <- accessible_phenotypes_singleStep_sr1 %>%
  ggplot(aes(x=n_new_pheno)) + 
  geom_bar(aes(y = after_stat(count))) +
  labs(x="Number of new\nphenotypes accessible",y="Number of specific RH genotypes") +
  geom_text(aes(label = scales::percent(after_stat(count)/sum(after_stat(count))),y= after_stat(count)), stat= "count", vjust = -.5,size=3) +
  theme_classic() + 
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=12))
  
# Local bias per genotype vs. global bias
p2 <- local_bias_sr1 %>% ggplot(aes(x="AncSR1",y=b)) +
  geom_violin() + 
  #geom_point(alpha=0.2,position = position_jitter(width = 0.1)) + 
  geom_pointrange(stat = "summary",size=0.4,col="red") +
  geom_hline(yintercept = GlobalB_sr1,col="red",linetype="dashed",linewidth=1.5) +
  labs(x="",y="Phenotype Bias (B)") +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25)) +
  theme_classic() +
  theme(axis.title = element_text(size=11),
        axis.text.y = element_text(size=13),
        axis.text.x = element_blank())

p1 + p2
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
print(paste("Ratio of mean local bias in one-mutant neighborhood to global bias: ", round(local_bias_sr1 %>% mutate(r = b/GlobalB_sr1) %>% reframe(m_r = mean(r)) %>% pull(m_r),2)))
```

    ## [1] "Ratio of mean local bias in one-mutant neighborhood to global bias:  2.18"

We can see that ~70% of specific genotypes cannot acces any new genotype within one mutation, and when they do, they mostly can access a single new phenotype. As a result, the local bias of the one-mutant neighborhood is on average 2X higher than the global bias. This also suggests that genotypes are surrounded mostly by neutral neighbors. We will check this now.

``` r
###################################
# fraction of neighbor nodes with same or different phenotype 
###################################

cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

# compute fraction of neighbors with same or different phenotype as focal variant (only for specific genotypes)
neighbors_sr1 <- foreach(i = 1:length(nodes_in_ntwrk_sr1_spec), .combine = 'rbind') %dopar% {
  fraction_of_neighbors_per_genotype(from_g = nodes_in_ntwrk_sr1_spec[i],Bg = "AncSR1",graph = net_sr1,pheno_tbl = phenotypes_tbl_prot)}

# To assess whether the avearge fraction of neutral neighbors is higher than expected, perform a non-parametric bootstrap test by 
# randomizing phenotypic annotations for RH variants, but keeping same network structure.
# Compute mean fraction of neutral neighbors for each scrambled GPmap and compute p-value as the fraction of bootstrap samples
# that are higher than the observed value.
N_BOOT <- 50

mean_prop_same_neighbors_permut_sr1 <- foreach(i = 1:N_BOOT) %:%
  foreach(j = 1:length(nodes_in_ntwrk_sr1_spec), .combine = 'rbind') %dopar% {
    pheno_tbl_permut <- phenotypes_tbl_prot %>% filter(bg=="AncSR1") %>% transform(.,specificity = sample(specificity,replace = F))
    fraction_of_neighbors_per_genotype(from_g = nodes_in_ntwrk_sr1_spec[j],Bg = "AncSR1",
                                       graph = net_sr1,pheno_tbl = pheno_tbl_permut)} %>%
  map_dbl(.f = function(x) x %>% filter(pheno_match == T) %>% with(mean(prop)))

stopCluster(cl)
```

``` r
# plot proportion of neutral neighbors per specific genotype
neighbors_sr1 %>% filter(pheno_match == T) %>%
  ggplot(aes(x=factor("AncSR1"),y=prop)) +
  geom_violin() +
  geom_pointrange(stat = "summary",size=0.4,col="red") + theme_classic() +
  labs(x="",y="Proportion of neutral neighbors per genotype") +
  scale_y_continuous(limits=c(0,1)) +
  theme(axis.title = element_text(size=10),
        axis.text.y = element_text(size=13),
        axis.text.x = element_blank())
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
# compute p-value
obs_mean <- neighbors_sr1 %>% filter(pheno_match == T) %>% with(mean(prop))
norm_fit_permut <- fitdistrplus::fitdist(mean_prop_same_neighbors_permut_sr1, "norm") # fit normal distribution to permuted dataset
#hist(mean_prop_same_neighbors_permut_sr1,xlim = c(0.2,0.8),freq = F,ylim = c(0,12))
#lines(x=seq(0.2,0.8,0.01),y=dnorm(x=seq(0.2,0.8,0.01),mean = norm_fit_permut$estimate[1],sd = norm_fit_permut$estimate[2]),col="red")

pval_neutral <- pnorm(obs_mean, mean = norm_fit_permut$estimate[1], sd = norm_fit_permut$estimate[2], lower.tail = FALSE, log.p = FALSE) # p-value
#sum(mean_prop_same_neighbors_permut_sr1>=obs_mean)/length(mean_prop_same_neighbors_permut_sr1) # empirical p-val

print(paste("Mean fraction of neutral neighbors in AncSR1 network: ", round(obs_mean,2),". P-value = ",pval_neutral))
```

    ## [1] "Mean fraction of neutral neighbors in AncSR1 network:  0.79 . P-value =  3.01738223898129e-20"

As expected, we see that 79% of single mutations result in genotypes with the same phenotype. That is, conservation is the most likely outcome of single mutations in the AncSr1 GP map.

Finally, let's investgate the distance between specificity phenotypes in the network. This gives us an idea of how likely are different phenotypic transitions to occur. We will assess to two types of distance: direct contact between phenotypes (whether phenotypes can be transitioned directly via a single amino acid change), and mean mutational distance (average number of amino acid changes between nodes encoding two phenotypes).

``` r
###################################
# Distance between neutral networks
###################################

RE_combos <- expand.grid(REs[[1]],REs[[1]])

# Number of direct links between phenotypes
direct_links_sr1 <- RE_combos %>% mutate(links = map2_dbl(.x=as.character(Var1),.y=as.character(Var2),
                                                    pairwise_neutral_network_proximity,Bg="AncSR1",graph=net_sr1,type=2,
                                                    pheno_df=phenotypes_tbl_prot,from_specific=T,to_specific=T)) %>% acast(Var1~Var2, value.var="links")

direct_links_sr1 <- direct_links_sr1[rowSums(is.na(direct_links_sr1)) != ncol(direct_links_sr1), ] # remove all NAs-rows
direct_links_sr1 <- direct_links_sr1[, colSums(is.na(direct_links_sr1)) != nrow(direct_links_sr1)] # remove all NAs-cols
direct_links_sr1 <- direct_links_sr1[order(match(rownames(direct_links_sr1),REs[[1]])),order(match(colnames(direct_links_sr1),REs[[1]]))] # order rows and columns
```

``` r
# heatmap of direct links
diag(direct_links_sr1) <- NA
direct_links_sr1 <- apply(direct_links_sr1,1,function(x) ifelse(x > 0,1,0))
pheatmap(direct_links_sr1,cluster_rows = F,cluster_cols = F,angle_col = 45, na_col = "white",
         legend = T,border_color = "black",color = c("white","gray30"),cellwidth = 20,cellheight = 20,main="Direct transitions")
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
# Average path length between phenotypes
pairwise_dists_sr1 <- RE_combos %>% mutate(links = map2_dbl(.x=as.character(Var1),.y=as.character(Var2),
                                                    pairwise_neutral_network_proximity,Bg="AncSR1",graph=net_sr1,type=3,
                                                    pheno_df=phenotypes_tbl_prot,from_specific=T,to_specific=T)) %>% acast(Var1~Var2, value.var="links")

pairwise_dists_sr1 <- pairwise_dists_sr1[rowSums(is.na(pairwise_dists_sr1)) != ncol(pairwise_dists_sr1), ] # remove all NAs-rows
pairwise_dists_sr1 <- pairwise_dists_sr1[, colSums(is.na(pairwise_dists_sr1)) != nrow(pairwise_dists_sr1)] # remove all NAs-cols
pairwise_dists_sr1 <- pairwise_dists_sr1[order(match(rownames(pairwise_dists_sr1),REs[[1]])),order(match(colnames(pairwise_dists_sr1),REs[[1]]))] # order rows and columns


# To assess whether distances are shorter or longer than expected, perform a non-parametric bootstrap test by 
# randomizing phenotypic annotations for RH variants, but keeping same network structure.
# Compute pairwise phenotype distances for each scrambled GPmap and compute p-value as the fraction of bootstrap samples
# that are shorter or longer than the observed value.
cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

N_BOOT <- 100
set.seed(3546)

# Compute phenotype distances for 100 bootstrap samples
mean_dist_pheno_matrix_permut <- foreach(i = 1:N_BOOT) %dopar% {
  pheno_tbl_permut <- phenotypes_tbl_prot %>% filter(bg=="AncSR1") %>% transform(.,specificity = sample(specificity))
  tmp_mat <- RE_combos %>% mutate(links = map2_dbl(.x=as.character(Var1),.y=as.character(Var2),
                                                    pairwise_neutral_network_proximity,Bg="AncSR1",graph=net_sr1,type=3,
                                                    pheno_df=pheno_tbl_permut,from_specific=T,to_specific=T)) %>%
    acast(Var1~Var2,value.var="links")
  tmp_mat <- tmp_mat[rowSums(is.na(tmp_mat)) != ncol(tmp_mat), ] # remove all NAs-rows
  tmp_mat <- tmp_mat[, colSums(is.na(tmp_mat)) != nrow(tmp_mat)] # remove all NAs-cols
  tmp_mat <- tmp_mat[order(match(rownames(tmp_mat),REs[[1]])),order(match(colnames(tmp_mat),REs[[1]]))] # order rows and columns
  tmp_mat <- tmp_mat[rownames(tmp_mat) != "CA",colnames(tmp_mat) != "CA"] # remove genotypes with CA specificity as this phenotype is not part of the main network
}

stopCluster(cl)

# remove bootstrap samples that do no have the same phenotypes as the emprical GPmap
mean_dist_pheno_matrix_permut <- mean_dist_pheno_matrix_permut[-which(lapply(mean_dist_pheno_matrix_permut,dim) %>% map(.,function(x) all(x[1] != c(6,6))) %>% unlist() == T)]
```

``` r
# empirical p-values for each phenotype pair
pval_matrix <- matrix(NA,nrow = 6,ncol = 6,dimnames = list(rownames(pairwise_dists_sr1),colnames(pairwise_dists_sr1)))

for(i in 1:nrow(pval_matrix)){
  for(j in 1:ncol(pval_matrix)){
    if(i != j){
      perm_vals <- map_dbl(mean_dist_pheno_matrix_permut, .f = function(x) x[i,j])
      if(pairwise_dists_sr1[i,j] > mean(perm_vals)) p <- sum(perm_vals>=pairwise_dists_sr1[i,j])/length(perm_vals)
      else p <- sum(perm_vals<=pairwise_dists_sr1[i,j])/length(perm_vals)
      pval_matrix[i,j] <- p
    }
  }
}

fdr_pvals <- p.adjust(pval_matrix,method = "bonferroni")
fdr_pval_matrix <- matrix(fdr_pvals,byrow = T,nrow = 6,ncol = 6,
                      dimnames = list(rownames(pairwise_dists_sr1),colnames(pairwise_dists_sr1)))

# heatmap
diag(pairwise_dists_sr1) <- NA
pheatmap(pairwise_dists_sr1,cluster_rows = F,cluster_cols = F,angle_col = 45, na_col = "white",
         legend = T,border_color = "black",color = rev(hcl.colors(10,palette = "viridis")),cellwidth = 20,cellheight = 20,
         breaks = seq(3.5,8,0.5),legend_breaks = seq(3.5,8,0.5),legend_labels = seq(3.5,8,0.5),main = "Average distances")
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
# phenotype pairs that are significantly closer/farther than expected
knitr::kable(inner_join(as.data.frame.table(fdr_pval_matrix), as.data.frame.table(pairwise_dists_sr1),by=c("Var1","Var2")) %>% filter(Freq.x < 0.05) %>%
               mutate_at(vars(Var1,Var2),list(as.character)) %>% rowwise() %>% mutate(pair=paste(min(Var1,Var2),max(Var1,Var2),sep="_")) %>% ungroup() %>%
               distinct(.,pair,.keep_all = TRUE) %>% select(Var1,Var2,Freq.y) %>% dplyr::rename(Phenotype_A=Var1,Phenotype_B=Var2,avg_mut_distance=Freq.y))
```

| Phenotype\_A | Phenotype\_B |  avg\_mut\_distance|
|:-------------|:-------------|-------------------:|
| GA           | SRE (AA)     |            3.523810|
| ERE (GT)     | SRE (AA)     |            6.552288|
| TA           | SRE (AA)     |            5.694444|
| ERE (GT)     | GA           |            7.176471|
| TA           | ERE (GT)     |            7.474265|

The figures above show that not all direct phenotypic transitions are possible, and that phenotypes are nonrandomly distributed in the network. In many cases, phenotypes are closer together or farther apart than expected due to the genetic determinants of specificity (*p* &lt; 0.05).

### Phenotypic evolution on the AncSR1 GP map

Since we have established that the AncSR1 GPmap is anisotrpic and heterogeneous, we will now focus on investigating the effects of this structure on evolutionary outcomes. For this, we will model evolution on the genotype network as a discrete Markov chain.

First, we will assess the timescales at which global and local bias affect evolutionary outcomes. The equilibrium outcome spectrum provides an approximation of the most likely evolutionary outcomes when there has been enough sequence evolution, such that the outcomes no longer depend on the specific starting genotype. By comparing the equilibrium spectrum to the global production spectrum, we can determine whether evolutionary outcomes are still biased at long timescales.

``` r
# Production spectrum and equilibrium outcome spectrum
df_long_term_sr1 <- inner_join(remove_promiscuous(var.prop_spec_AncSR1_df,norm = T),remove_promiscuous(stationary_PDFV_spec_sr1,norm = T),by="RE")

# correlation between spectra
cor <- df_long_term_sr1 %>% with(cor(Norm_F_prob.x,Norm_F_prob.y)^2)

# plot correlation
p1 <- df_long_term_sr1 %>%
  ggplot(aes(x=Norm_F_prob.x,y=Norm_F_prob.y,fill=RE)) + 
  geom_point(shape=21,col="black",size=4) + scale_fill_manual(values = hex_RE_colors(1)) +
  scale_y_continuous(limits = c(0,0.45)) + scale_x_continuous(limits = c(0,0.45)) + 
  geom_abline(intercept = 0,slope = 1,linetype="dashed",col="black") +
  labs(x="Global production spectrum",y="Equilibrium outcome spectrum") + 
  theme_classic() +
  theme(axis.text = element_text(size = 10),axis.title = element_text(size=11)) +
  #annotate("text",x=0.05,y=0.42,label=expression(paste("r"^2,"= ",.(round(cor,2)))))
  annotate("text", x = 0.05, y = 0.42, label = bquote(r^2 == .(round(cor,2))))

# plot barplot per spectrum
bias_stat_global <- data.frame(model=c("Production bias","Equilibrium"),b=c(GlobalB_sr1,StatB_sr1),RE="TA",y=0.3)
p2 <- bars_PDFV(list(remove_promiscuous(var.prop_spec_AncSR1_df,norm = T) %>% mutate(model="Production bias"),
               remove_promiscuous(stationary_PDFV_spec_sr1,norm = T) %>% mutate(model="Equilibrium")),
          arrange_by = "Production bias",fill_by_RE = T,hex_id = 1) +
  geom_text(data=bias_stat_global,aes(x=RE,y=y,label=paste("B =",round(b,2))))


p1 + p2
```

    ## Warning in is.na(x): is.na() applied to non-(list or vector) of type 'language'

    ## Don't know how to automatically pick scale for object of type <table>.
    ## Defaulting to continuous.

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/timescales_sr1_1-1.png)

We see that the equilibrium outcome spectrum is highly correlated with the global production spectrum, suggesting that the GPmap biases evolutionary outcomes in the long term. However, we also see that the two spectra are not identical, which captures the effect of the network connectivity in the long-term equilibrium occupancy of genotypes and phenotypes. Thus, both global bias and heterogeneity shape evolutionary outcomes even over infinitely timescales.

Now we will evaluate the effect of the GPmap on shorter timescales. For this, we will compute the bias in evolutionary outcomes from every genotype in the network and will track its dynamics a function of the number of amino acid substitutions.

``` r
###################################
# Change in outcome bias over time
###################################

cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

MAX_STEP = 150
change_bias_df_sr1 <- data.frame(step=NULL,bias=NULL,AA_var=NULL)

for(s in c(1,seq(5,MAX_STEP,5))){
  # Compute outcome spectra for each genotype per step
  tmp_mat_i <- foreach(i = 1:length(nodes_in_ntwrk_sr1), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr1[i],tr_mat = P_drift_sr1_ntwrk,n_steps = s) # run markov chain from genotype i
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR1",graph=net_sr1,model = nodes_in_ntwrk_sr1[i],specific = T,pheno_tbl = phenotypes_tbl_prot) %>% # compute outcome spectrum
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[2]]))]
  }
  rownames(tmp_mat_i) <- nodes_in_ntwrk_sr1
  tmp_mat_i <- tmp_mat_i[, which(colSums(tmp_mat_i) != 0)]
  
  # compute bias per genotype per step
  tmp_mat_i <- t(apply(tmp_mat_i[,1:ncol(tmp_mat_i)-1],1,function(x) x/sum(x))) # remove promiscuous phenotype and renormalize probs
  B_g <- apply(tmp_mat_i, 1, get_bias,Bg="AncSR1",input_freqs=T) # compute bias
  df_tmp <- data.frame(step=s,bias=B_g,AA_var=nodes_in_ntwrk_sr1)
  change_bias_df_sr1 <- rbind(change_bias_df_sr1,df_tmp)
}

stopCluster(cl)
```

``` r
# Number of substitutions for mean bias in outcomes to be within 0.05 units of equilibrium value
forget_lbias_sr1 <- change_bias_df_sr1 %>% group_by(AA_var) %>% filter(bias <= StatB_sr1+0.05) %>% 
  distinct(.,AA_var,.keep_all = TRUE) %>% with(mean(step))

print(paste("Number of substitutions for mean bias in outcomes to be within 0.05 units of equilibrium value = ",round(forget_lbias_sr1,2)))
```

    ## [1] "Number of substitutions for mean bias in outcomes to be within 0.05 units of equilibrium value =  41.91"

``` r
# Number of substitutions for bias in outcomes from EGKA to be within 0.05 units of equilibrium value
forget_lbias_EGKA_sr1 <- change_bias_df_sr1 %>% filter(AA_var==REF_GENOTYPE) %>% filter(bias <= StatB_sr1+0.05) %>% head(., 1) %>% pull(step)

print(paste("Number of substitutions for bias in outcomes from EGKA to be within 0.05 units of equilibrium value = ",forget_lbias_EGKA_sr1))
```

    ## [1] "Number of substitutions for bias in outcomes from EGKA to be within 0.05 units of equilibrium value =  80"

``` r
# change in outcome bias over time from ancestral genotype EGKA
egka_bias_sr1 <- change_bias_df_sr1 %>% filter(step<=80) %>% filter(AA_var=="EGKA")

# average change in outcome bias over time across genotypes
mean_bias_sr1 <- change_bias_df_sr1 %>% filter(step<=80) %>% group_by(step) %>% reframe(mean_B=mean(bias,na.rm=T))

mut_to_brlen <- function(x) x / 4 # branch length units
MAX_STEP = 150
change_bias_df_sr1 %>% filter(step<=80) %>%
  ggplot(aes(x=step,y=bias)) + geom_point(color="gray60",size=0.5) + 
  xlab("Substitution step") + ylab("Outcome Bias (B)") +
  geom_line(data=mean_bias_sr1,aes(x=step,y=mean_B),color="orange",size=1.5) +
  geom_hline(yintercept = GlobalB_sr1,col="red",linetype="dashed",linewidth=0.7) +
  geom_hline(yintercept = StatB_sr1,col="orange",linetype="dashed",linewidth=0.7) +
  geom_vline(xintercept = forget_lbias_sr1,col="gray80",linewidth=1.5) +
  geom_line(data=egka_bias_sr1,aes(x=step,y=bias),color="black",size=0.5) +
  theme_classic() +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0,MAX_STEP,10),
                     sec.axis = dup_axis(trans = ~mut_to_brlen(.), name = "Branch length (subs / RH site)",breaks = seq(0, 20, by = 2))) +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        legend.position = "none")
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-15-1.png)

This shows that local bias affects evolutionary outcomes at short-to-moderate timescales. Over time, the outcome bias decays to its equilibrium expectation, but it takes ~42 substitutions on average to do so, much longer than relevant substitution timescales for natural proteins.

Since we know that there are 3 amino acid differences between the ancestral and derived RH genotypes (EGKA to GSKV), we can further assess how local biases affect the number of phenotypes that can evolve from each genotype after 3 substitutions.

``` r
cl <- makeCluster(N_CORES,type = "FORK", outfile="")
registerDoParallel(cl)

N_MUTATIONS <- 3 # ~ 0.75 subs per site
SPEC_BINDING <- TRUE # whether to consider only specific binders (promiscuous as a separate "phenotype")
NORM_PROB <- TRUE # whether to remove promiscuous end-points and renormalize probabilities (only if SPEC_BINDING = TRUE)

if(SPEC_BINDING){
  # if considering phenotype as specific binding
  markov_chain_sr1_3s <- foreach(i = 1:length(nodes_in_ntwrk_sr1), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr1[i],tr_mat = P_drift_sr1_ntwrk,n_steps = N_MUTATIONS) # run markov chain
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR1",graph=net_sr1,model = nodes_in_ntwrk_sr1[i],specific = T,pheno_tbl = phenotypes_tbl_prot) %>% # compute outcome spectrum
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[2]]))]
  }
  rownames(markov_chain_sr1_3s) <- nodes_in_ntwrk_sr1
  markov_chain_sr1_3s <- markov_chain_sr1_3s[, which(colSums(markov_chain_sr1_3s) != 0)]
  
  if(NORM_PROB){
    # remove promiscuous and renormalize probabilities
    markov_chain_sr1_3s <- markov_chain_sr1_3s[,1:ncol(markov_chain_sr1_3s)-1]
    markov_chain_sr1_3s <- t(apply(markov_chain_sr1_3s, 1, function(x) x/sum(x)))
  }
} else{
  # if considering phenotype as binding
  markov_chain_sr1_3s <- foreach(i = 1:length(nodes_in_ntwrk_sr1), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr1[i],tr_mat = P_drift_sr1_ntwrk,n_steps = N_MUTATIONS)
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR1",graph=net_sr1,model = nodes_in_ntwrk_sr1[i],specific = F, pheno_tbl = phenotypes_tbl_prot) %>%
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
  }
  rownames(markov_chain_sr1_3s) <- nodes_in_ntwrk_sr1
  markov_chain_sr1_3s <- markov_chain_sr1_3s[, which(colSums(markov_chain_sr1_3s) != 0)]
}

stopCluster(cl)

# Distribution of bias in outcomes after 3 steps
lb_3s_df <- data.frame(AA_var = rownames(markov_chain_sr1_3s),b = apply(markov_chain_sr1_3s,1,get_bias,input_freqs=T))

# compute mean outcome bias
mean_lb_3s <- lb_3s_df %>% with(mean(b))

# plot distribution
lb_3s_df %>%
  ggplot(aes(x=b)) + geom_histogram(binwidth = 0.01,col="black",fill="gray65") + 
  geom_vline(xintercept = mean_lb_3s,linetype="dashed") +
  theme_classic() +
  labs(x = "Bias in outcomes after 3 steps", y= "Frequency") +
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=9))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
print(paste("Mean bias in outcomes across genotypes after 3 steps = ",round(mean_lb_3s,2)))
```

    ## [1] "Mean bias in outcomes across genotypes after 3 steps =  0.79"

``` r
# Distribution of novel phenotypes that can evolve from each starting genotype after 3 substitutions
new_pheno_3s_df <- phenotypes_tbl_prot %>% filter(bg=="AncSR1") %>% select(AA_var,specificity) %>% 
  inner_join(.,as_tibble(markov_chain_sr1_3s,rownames = NA) %>% rownames_to_column(var="AA_var"),by="AA_var") %>%
  reshape2::melt() %>% 
  # keep phenotypes that have non-zero probability to evolve
  filter(value > 0) %>%
  # check whether phenotype is different from initial phenotype
  mutate(new_pheno = ifelse((specificity != variable & value > 0),1,0)) %>%
  group_by(AA_var) %>% reframe(n_new_pheno = sum(new_pheno))
```

    ## Using AA_var, specificity as id variables

``` r
# compute mean number of novel phenotypes
mean_new_pheno_3s <- new_pheno_3s_df %>% with(mean(n_new_pheno))

# plot distribution
new_pheno_3s_df %>%
  ggplot(aes(x=n_new_pheno)) + geom_histogram(binwidth = 1,col="black",fill="gray65") + 
  geom_vline(xintercept = mean_new_pheno_3s,linetype="dashed") +
  theme_classic() +
  labs(x = "Number of new phenotypes per genotype", y= "Frequency") +
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=9))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-16-2.png)

``` r
print(paste("Mean number of new phenotypes per genotype = ",round(mean_new_pheno_3s,2)))
```

    ## [1] "Mean number of new phenotypes per genotype =  2.43"

We can see that most genotypes can only evolve around 2 or 3 novel phenotypes after 3 substitutions. Thus, local bias exerts a strong influence in the outcomes of phenotypic diversification at short-to-moderate timescales.

Now that we have established the timescales at which local and global bias matter, let's investigate what is the effect of local bias in the direction of phenotypic diversification. For this, we will analyze the spectrum of outcomes after 8 steps, for which the influence of local bias is still high, but several phenotypes become more likely to evolve.

``` r
###################################
# Compute outcome spectra from every genotype after 8 substitutions
###################################
cl <- makeCluster(N_CORES,type = "FORK",outfile="")
registerDoParallel(cl)

N_MUTATIONS <- 8 # ~ 2 subs per site
SPEC_BINDING <- TRUE # whether to consider only specific binders (promiscuous as a separate "phenotype")
NORM_PROB <- TRUE # whether to remove promiscuous end-points and renormalize probanilities (only if SPEC_BINDING = TRUE)

# Run Markov chain from every genotype in the network for 8 amino acid substitutions and compute the outcome spectra
if(SPEC_BINDING){
  # if considering phenotype as specific binding
  markov_chain_sr1 <- foreach(i = 1:length(nodes_in_ntwrk_sr1), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr1[i],tr_mat = P_drift_sr1_ntwrk,n_steps = N_MUTATIONS) # run markov chain for genotype i
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR1",graph=net_sr1,model = nodes_in_ntwrk_sr1[i],specific = T,pheno_tbl = phenotypes_tbl_prot) %>% # compute outcome spectrum
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[2]]))]
  }
  rownames(markov_chain_sr1) <- nodes_in_ntwrk_sr1
  markov_chain_sr1 <- markov_chain_sr1[, which(colSums(markov_chain_sr1) != 0)] # remove phenotypes not encoded
  
  if(NORM_PROB){
    # remove promiscuous and renormalize probabilities
    markov_chain_sr1 <- markov_chain_sr1[,1:ncol(markov_chain_sr1)-1]
    markov_chain_sr1 <- t(apply(markov_chain_sr1, 1, function(x) x/sum(x)))
  }
} else{
  # if considering phenotype as binding
  markov_chain_sr1 <- foreach(i = 1:length(nodes_in_ntwrk_sr1), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr1[i],tr_mat = P_drift_sr1_ntwrk,n_steps = N_MUTATIONS)
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR1",graph=net_sr1,model = nodes_in_ntwrk_sr1[i],specific = F,pheno_tbl = phenotypes_tbl_prot) %>%
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
  }
  rownames(markov_chain_sr1) <- nodes_in_ntwrk_sr1
  markov_chain_sr1 <- markov_chain_sr1[, which(colSums(markov_chain_sr1) != 0)]
}

stopCluster(cl)
```

``` r
# Run Markov chain for 8 amino acid substitutions averaged across all genotypes
N_MUTATIONS <- 8 # ~ 2 subs per site
if(SPEC_BINDING){
  mc_avg <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr1,tr_mat = P_drift_sr1_ntwrk,n_steps = N_MUTATIONS)
  outcome_spec_avg_all_sr1 <- get_PDFV_v2(mc_avg,type="simulated mc", Bg="AncSR1",graph=net_sr1,model = "All nodes",specific = T,pheno_tbl = phenotypes_tbl_prot) %>%
        acast(model~RE,value.var="Norm_F_prob")
  outcome_spec_avg_all_sr1 <- outcome_spec_avg_all_sr1[,order(match(colnames(outcome_spec_avg_all_sr1),REs[[3]]))]
  outcome_spec_avg_all_sr1 <- outcome_spec_avg_all_sr1[which(outcome_spec_avg_all_sr1 != 0)]
  if(NORM_PROB){
    outcome_spec_avg_all_sr1 <- outcome_spec_avg_all_sr1[1:length(outcome_spec_avg_all_sr1)-1]
    outcome_spec_avg_all_sr1 <- outcome_spec_avg_all_sr1/sum(outcome_spec_avg_all_sr1)
  }
} else {
  mc_avg <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr1,tr_mat = P_drift_sr1_ntwrk,n_steps = N_MUTATIONS)
  outcome_spec_avg_all_sr1 <- get_PDFV_v2(mc_avg,type="simulated mc", Bg="AncSR1",graph=net_sr1,model = "All nodes",specific = F,pheno_tbl = phenotypes_tbl_prot) %>%
        acast(model~RE,value.var="Norm_F_prob")
  outcome_spec_avg_all_sr1 <- outcome_spec_avg_all_sr1[,order(match(colnames(outcome_spec_avg_all_sr1),REs[[3]]))]
  outcome_spec_avg_all_sr1 <- outcome_spec_avg_all_sr1[which(outcome_spec_avg_all_sr1 != 0)]
}

# Run Markov chain for 8 amino acid substitutions averaged across non-ERE genotypes
nonERE_genotypes_ntwrk_sr1 <- phenotypes_tbl_prot %>% filter(bg =="AncSR1" & specificity != "ERE (GT)") %>% pull(AA_var)
nonERE_genotypes_ntwrk_sr1 <- nonERE_genotypes_ntwrk_sr1[nonERE_genotypes_ntwrk_sr1 %in% nodes_in_ntwrk_sr1]

if(SPEC_BINDING){
  mc_avg_2 <- simulate_markov_chain(states_0 = nonERE_genotypes_ntwrk_sr1,tr_mat = P_drift_sr1_ntwrk,n_steps = N_MUTATIONS)
  outcome_spec_avg_all_sr1_2 <- get_PDFV_v2(mc_avg_2,type="simulated mc", Bg="AncSR1",graph=net_sr1,model = "non-ERE nodes",specific = T,pheno_tbl = phenotypes_tbl_prot) %>%
        acast(model~RE,value.var="Norm_F_prob")
  outcome_spec_avg_all_sr1_2 <- outcome_spec_avg_all_sr1_2[,order(match(colnames(outcome_spec_avg_all_sr1_2),REs[[3]]))]
  outcome_spec_avg_all_sr1_2 <- outcome_spec_avg_all_sr1_2[which(outcome_spec_avg_all_sr1_2 != 0)]
  if(NORM_PROB){
    outcome_spec_avg_all_sr1_2 <- outcome_spec_avg_all_sr1_2[1:length(outcome_spec_avg_all_sr1_2)-1]
    outcome_spec_avg_all_sr1_2 <- outcome_spec_avg_all_sr1_2/sum(outcome_spec_avg_all_sr1_2)
  }
} else {
  mc_avg_2 <- simulate_markov_chain(states_0 = nonERE_genotypes_ntwrk_sr1,tr_mat = P_drift_sr1_ntwrk,n_steps = N_MUTATIONS)
  outcome_spec_avg_all_sr1_2 <- get_PDFV_v2(mc_avg,type="simulated mc", Bg="AncSR1",graph=net_sr1,model = "non-ERE nodes",specific = F,pheno_tbl = phenotypes_tbl_prot) %>%
        acast(model~RE,value.var="Norm_F_prob")
  outcome_spec_avg_all_sr1_2 <- outcome_spec_avg_all_sr1_2[,order(match(colnames(outcome_spec_avg_all_sr1_2),REs[[3]]))]
  outcome_spec_avg_all_sr1_2 <- outcome_spec_avg_all_sr1_2[which(outcome_spec_avg_all_sr1_2 != 0)]
}

# "Admixture plot"-like figure
# Associate genotypes to clusters
clusters_sr1 <- modules_sr1 %>% select(AA_var,specificity,module) %>% filter(AA_var %in% nodes_in_ntwrk_sr1) %>%
  dplyr::rename("RE" = "specificity")

# reorder RH variants to match clusters_sr1 data frame
markov_chain_sr1_plot <- as_tibble(markov_chain_sr1) %>% mutate(AA_var=nodes_in_ntwrk_sr1) %>% inner_join(.,clusters_sr1,by="AA_var") %>%
  column_to_rownames(var="AA_var") %>% select(-c(RE,module)) %>% as.matrix(.)

# Phenotype of starting genotype
starting_DNAspec_df <- clusters_sr1 %>% mutate(variable=RE,value = -0.05)

# build data frame to plot
markov_chain_sr1_plot <- as_tibble(markov_chain_sr1_plot,rownames=NA) %>% rownames_to_column(var="AA_var") %>% inner_join(.,clusters_sr1,by="AA_var") %>%
  reshape2::melt() %>%
  rbind(.,starting_DNAspec_df) %>%
  arrange(desc(factor(variable,levels = levels(REs[[3]])))) %>% 
  arrange(factor(RE,levels = levels(REs[[5]]))) %>% arrange(module) %>%
  mutate(variable = forcats::as_factor(variable),
         module = forcats::as_factor(module),
         AA_var = forcats::as_factor(AA_var),
         RE = forcats::as_factor(RE))
```

    ## Using AA_var, RE, module as id variables

``` r
# plot from every genotype
admx_plot_sr1 <- markov_chain_sr1_plot %>%
  ggplot(aes(x=AA_var,y=value,fill=variable,group=module)) + geom_bar(stat="identity",width=1) + 
  scale_fill_manual(values = hex_RE_colors(2)) +
  theme_classic() + facet_grid(cols=vars(module),scales="free_x",space="free_x") + 
  geom_hline(yintercept = 0,linewidth=2,color="white") +
  labs(x="Starting genotype",y="Pr(Ending phenotype)",fill="") +
  theme(axis.title = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.text.x = element_blank(),
        legend.position="bottom")

# plot from average outcomes across genotypes
admx_plot_avgs_sr1 <- rbind(data.frame(model = "All",pheno = names(outcome_spec_avg_all_sr1), prob = outcome_spec_avg_all_sr1),
                       data.frame(model = "non ERE",pheno = names(outcome_spec_avg_all_sr1_2), prob = outcome_spec_avg_all_sr1_2)) %>%
  arrange(factor(pheno,levels = levels(REs[[3]]))) %>% 
  mutate(pheno = forcats::as_factor(pheno)) %>%
  ggplot(aes(x=model,y=prob,fill=pheno)) + geom_bar(stat="identity",width=0.9) + 
  scale_fill_manual(values = hex_RE_colors(2)) +
  theme_classic() +
  geom_hline(yintercept = 0,linewidth=2,color="white") +
  labs(x="",y="Pr(Ending phenotype)",fill="") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        legend.position="none")

admx_plot_sr1 + admx_plot_avgs_sr1 + plot_layout(widths = c(6,0.5))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-18-1.png)

This plot shows the probability of evolving each phenotypic outcome from every starting genotype after 8 substitutions. We can see that genotypes in the same cluster tend to evolve more similar phenotypic outcomes, which captures the effect of the heterogeneity of the GPmap.

Since phenotypes are strongly associated with specific clusters, we can also assess how does the probability of particular phenotypic transitions vary with the starting phenotype.

``` r
###################################
# Probabilities of transition aggregated by phenotype
###################################
cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

# phenotypes in network
REs_spec_sr1 <- data.frame(summary_sr1) %>% rownames_to_column(var="RE") %>% filter(Encoded == T & In_Net_Spec == T) %>% pull(RE)

# compute phenotypic transitions: outcome spectra aggregated for RH variants encoding the same phenotype
pt_sr1 <- foreach(i = 1:length(REs_spec_sr1), .combine = 'rbind') %dopar% {
  # select RH variants with same phenotype
  phenotype_vars <- phenotypes_tbl_prot %>% filter(bg == "AncSR1" & specificity == REs_spec_sr1[i]) %>% pull(AA_var)
  # run markov chain from all neutral variants and compute outcome spectra
  tmp_mc <- simulate_markov_chain(states_0 = phenotype_vars,tr_mat = P_drift_sr1_ntwrk,n_steps = N_MUTATIONS)
  tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR1",graph=net_sr1,model = REs_spec_sr1[i],specific = T,pheno_tbl = phenotypes_tbl_prot) %>%
    acast(model~RE,value.var="Norm_F_prob")
  tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
}
rownames(pt_sr1) <- REs_spec_sr1
pt_sr1 <- pt_sr1[, which(colSums(pt_sr1) != 0)]
pt_sr1 <- pt_sr1[,1: ncol(pt_sr1)-1]
pt_sr1 <- t(apply(pt_sr1, 1, function(x) x/sum(x)))

parallel::stopCluster(cl)
```

``` r
# Admixture plot-like figure
as_tibble(pt_sr1,rownames=NA) %>% rownames_to_column(var="Start_Pheno") %>%
  reshape2::melt() %>%
  rbind(.,data.frame(Start_Pheno = rownames(pt_sr1),variable=rownames(pt_sr1),value = -0.05)) %>%
  arrange(desc(factor(variable, levels = levels(REs[[3]])))) %>%
  arrange(factor(Start_Pheno, levels = levels(REs[[3]]))) %>%
  mutate(variable = factor(variable, levels = levels(REs[[3]])),
         Start_Pheno = factor(Start_Pheno, levels = levels(REs[[3]]))) %>%
  ggplot(aes(x=Start_Pheno,y=value,fill=variable)) + 
  geom_bar(stat="identity",width=0.95) + 
  scale_fill_manual(values = hex_RE_colors(1)) +
  theme_classic() +
  geom_hline(yintercept = 0,linewidth=2,color="white") +
  labs(x="Starting phenotype",y="Pr(Ending phenotype)",fill="") +
  theme(axis.title = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13,angle=45,hjust = 1,vjust = 1),
        legend.position="bottom")
```

    ## Using Start_Pheno as id variables

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-20-1.png)

The likelihood of particular phenotypic transitions also vary across starting phenotypes. These results establish that the structure of the GP map - local bias, in particular - introduces lineage-specificity in phenotypic evolution: the likely outcomes of evolution depend on the starting genotype and phenotype.

Heterogeneity and mutational connectivity also imply that as genotypes evolve in the RH network, the spectrum of outcomes will change. Thus, we can quantity the effect of heterogeneity in the map by assessing how much of the variation in the similariry of evolutionary outcomes is explained by cluster membership, and also assess the effect of mutational distance on the variation in similarity of phenotypic outcomes.

``` r
###################################
# Compute similarity of phenotypic outcomes - Pairwise correlations of outcome spectra between genotypes
###################################

# pairwise r2 of outcome spectra between genotypes
pairwise_corr_pdfvs_sr1 <- foreach(i = 1:dim(markov_chain_sr1)[1],.combine='rbind',.errorhandling="pass") %:%
  foreach(j = (i+1):dim(markov_chain_sr1)[1],.combine='rbind') %do% {
    g_i <- rownames(markov_chain_sr1)[i]
    g_j <- rownames(markov_chain_sr1)[j]
    c <- cor(markov_chain_sr1[i,],markov_chain_sr1[j,])^2
    data.frame(g_i = g_i, g_j = g_j, r2 = c)
  }
pairwise_corr_pdfvs_sr1 <- pairwise_corr_pdfvs_sr1 %>% filter(g_i != g_j) # remove same-genotype comparisons (if any)

# assign genotype clusters of starting gentoypes
cluster_g_i <- data.frame(AA_var=unique(pairwise_corr_pdfvs_sr1$g_i),bg="AncSR1") %>% inner_join(.,modules_sr1,by=c("bg","AA_var")) %>% 
  select(AA_var,module) %>% dplyr::rename("g_i"="AA_var")
cluster_g_j <- data.frame(AA_var=unique(pairwise_corr_pdfvs_sr1$g_j),bg="AncSR1") %>% inner_join(.,modules_sr1,by=c("bg","AA_var")) %>% 
  select(AA_var,module) %>% dplyr::rename("g_j"="AA_var")

# assign pairwise mutational distance between starting genotypes
dist_tbl_sr1 <- distances(net_sr1)
assign_dist <- function(var1,var2,d_tbl){
  return(d_tbl[which(rownames(d_tbl)==var1),which(colnames(d_tbl)==var2)])}

plan("multisession", workers = N_CORES)
df_pairwise_corr_pdfvs_sr1 <- pairwise_corr_pdfvs_sr1 %>%
  inner_join(.,cluster_g_i,by="g_i") %>% inner_join(.,cluster_g_j,by="g_j") %>%
  mutate(cluster = ifelse(module.x == module.y,"SAME","DIFF"),
         m_dist = furrr::future_map2_int(.x=g_i,.y=g_j,.f=assign_dist,d_tbl=dist_tbl_sr1))
plan("sequential")

# Average pairwise similarity in outcomes for genotypes in the same or different clusters
knitr::kable(df_pairwise_corr_pdfvs_sr1 %>% group_by(cluster) %>% reframe(mean_similarity = mean(r2)))
```

| cluster |  mean\_similarity|
|:--------|-----------------:|
| DIFF    |         0.2021386|
| SAME    |         0.8017876|

``` r
# Linear model to assess the effect of mutational distance and cluster membership on outsome similarity
full_lm_mod <- lm(df_pairwise_corr_pdfvs_sr1$r2~df_pairwise_corr_pdfvs_sr1$m_dist + df_pairwise_corr_pdfvs_sr1$cluster)
anova(full_lm_mod)
```

    ## Analysis of Variance Table
    ## 
    ## Response: df_pairwise_corr_pdfvs_sr1$r2
    ##                                      Df Sum Sq Mean Sq F value    Pr(>F)    
    ## df_pairwise_corr_pdfvs_sr1$m_dist     1 209.25 209.254  2924.7 < 2.2e-16 ***
    ## df_pairwise_corr_pdfvs_sr1$cluster    1  88.98  88.983  1243.7 < 2.2e-16 ***
    ## Residuals                          4653 332.91   0.072                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# plot similarity in outcomes as a function of mutational distance between genotypes
# average similarity per mutational distance
meanC_per_step <- df_pairwise_corr_pdfvs_sr1 %>% group_by(m_dist) %>% reframe(mean_cor = mean(r2))

mut_dist_plot <- df_pairwise_corr_pdfvs_sr1 %>% ggplot(aes(x=m_dist,y=r2)) +
  geom_point(col="gray60",size=0.8,alpha=0.5) + 
  geom_line(data=meanC_per_step,aes(x=m_dist,y=mean_cor),color="orange",linewidth=1.5) +
  labs(x="Mutational distance",y=expression(paste("Similarity in evolutionary outcomes (r"^2,")"))) +
  scale_x_continuous(limits = c(1,11),breaks = seq(1,11,1)) +
  theme_classic() + 
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=11))

# plot distribution of similarity in outcomes by cluster membership
hist_cor <- df_pairwise_corr_pdfvs_sr1 %>%
  ggplot(aes(x=r2,fill=cluster)) + geom_histogram(color="black", position = 'stack',binwidth = 0.02) + 
  scale_fill_manual(values=c("#ff6db6", "#00ACFF")) +
  labs(x=expression(paste("Similarity in evolutionary outcomes (r"^2,")")), y="Frequency",fill="Genotype cluster") + 
  theme_classic() + 
  theme(axis.title.x = element_text(size=13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=9),
        axis.text.y = element_blank(),
        legend.position=c(0.6,0.5),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7)) +
  coord_flip()

mut_dist_plot + hist_cor
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/outcome_similarity_sr1-1.png)

These analyses show that as genotypes diverge from each other, the similarity in outcomes declines rapidly. And as expected, genotypes in the same cluster have on average more similar outcome spectra (mean *r*<sup>2</sup> = 0.8) compared to genotypes in different clusters (mean *r*<sup>2</sup> = 0.2). Thus, heterogeneity in the GPmap affect the direction of phenotypic diversification, and evolving lineages diverge rapidly in their most likely outcomes of evolution.

Since we saw before that local bias tends to favor phenotypic conservation, because most single mutations are neutral, we can also ask whether phenotypic conservation is the most likely outcome of evolution.

``` r
###################################
# Probability of conservation and transition per genotype
###################################

#associate genotypes to phenotypes
geno_pheno_sr1 <- phenotypes_tbl_prot %>% filter(bg=="AncSR1") %>% select(AA_var,specificity) %>% dplyr::rename("RE" = "specificity")
 
# extract probabilities of the most likely phenotypic outcome
max_prob_pheno <- apply(markov_chain_sr1,1,function(x) colnames(markov_chain_sr1)[which(x == max(x))]) # phenotype with the maximum transition probability
max_prob <- apply(markov_chain_sr1,1,max) # max prob.

# Most likely phenotypic outcome per *specific* genotype (conservation or transition)
highest_Lik_outcomes_sr1_df <- inner_join(geno_pheno_sr1,data.frame(AA_var=names(max_prob),max_prob=max_prob)) %>% 
  inner_join(.,data.frame(AA_var=names(max_prob_pheno),max_prob_pheno=max_prob_pheno)) %>%
  mutate(hl_outcome = ifelse(RE == max_prob_pheno,"Conservation","Transition")) %>% filter(RE != "Promiscuous")
```

    ## Joining with `by = join_by(AA_var)`
    ## Joining with `by = join_by(AA_var)`

``` r
# Extract probability of conservation and probability of phenotypic transition per genotype
if(SPEC_BINDING && NORM_PROB){
  # only including specificity phenotypes
  pheno_outcomes_sr1_df <- as_tibble(markov_chain_sr1,rownames = NA) %>% 
    rownames_to_column(var="AA_var") %>% 
    inner_join(.,geno_pheno_sr1,by="AA_var") %>%
    pivot_longer(cols=2:7,names_to="Pheno",values_to="Prob") %>% 
    filter(RE != "Promiscuous") %>%
    group_by(AA_var) %>% filter(RE == Pheno) %>%
    reframe(p_cons = Prob,
            p_tr = 1-p_cons) %>%
    pivot_longer(cols=2:3,names_to="outcome",values_to="prob")
} else {
  # including specificity phenotypes and promiscuity as an extra phenotype
  pheno_outcomes_sr1_df <- as_tibble(markov_chain_sr1,rownames = NA) %>% 
    rownames_to_column(var="AA_var") %>% 
    inner_join(.,geno_pheno_sr1,by="AA_var") %>%
    pivot_longer(cols=2:8,names_to="Pheno",values_to="Prob") %>% 
    filter(RE != "Promiscuous") %>%
    group_by(AA_var) %>% filter(RE == Pheno) %>%
    reframe(p_cons = Prob,
            p_tr = 1-p_cons) %>%
    pivot_longer(cols=2:3,names_to="outcome",values_to="prob")
}

# average probability of conservation across specific genotypes
avg_pr_cons <- pheno_outcomes_sr1_df %>% filter(outcome == "p_cons") %>% with(mean(prob))
print(paste("Average probability of conservation = ",round(avg_pr_cons,2)))
```

    ## [1] "Average probability of conservation =  0.61"

``` r
# percentaje of genotypes for which conservation or transition is the most likely outcome
knitr::kable(highest_Lik_outcomes_sr1_df %>% ungroup() %>% group_by(hl_outcome) %>%
               reframe(count = n()) %>% ungroup() %>% mutate(total = sum(count), frac = count/total) %>% select(hl_outcome,frac))
```

| hl\_outcome  |       frac|
|:-------------|----------:|
| Conservation |  0.8313253|
| Transition   |  0.1686747|

``` r
inner_join(pheno_outcomes_sr1_df,highest_Lik_outcomes_sr1_df,by="AA_var") %>%
  filter(outcome == "p_cons") %>%
  ggplot(aes(x=prob,fill=hl_outcome)) + geom_histogram(position = "stack",col="black") + 
  scale_fill_manual(labels = c("Conservation (83%)","Transition (17%)"), values = c("white","gray40")) +
  geom_vline(xintercept = avg_pr_cons,linetype="dashed") +
  theme_classic() + scale_x_continuous(limits = c(0,1)) +
  labs(x="P(conservation) per genotype",y="Frequency",fill="Highest likelihood\noutcome") +
  theme(axis.title = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust = 1,vjust = 1),
        legend.position=c(0.3,0.8))
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 4 rows containing missing values (`geom_bar()`).

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-21-1.png)

Conservation is indeed the most likely outcome for 83% of genotypes, and the mean probability of conservation is also high (P = 0.61). However, individual genotypes vary substantially in their propensities to evolve new phenotypes as observed from the variation along the x-axis. All these analyses show that local bias introduces lineage specificity both in the propensity that a genotype will evolve a novel phenotype, and in the phenotype that will most likely evolve from it.

Finally, we can assess the likely outcomes of evolution from the ancestral RH genotype EGKA. For this, we will track the outcome spectra over time.

``` r
###################################
# Dynamics of evolutionary outcomes over time from EGKA 
###################################

mc_iter <- 25 # number of substitutions to track
SPEC_BINDING <- TRUE

# Run Markov chain for 25 substitutions and track the outcome spectra at every time point
pdfv_mc_multistep_sr1 <- simulate_markov_chain_multistep(REF_GENOTYPE,P_drift_sr1_ntwrk,mc_iter,"AncSR1",specific = SPEC_BINDING)
```

``` r
mc_iter <- 25 # number of substitutions to track
SPEC_BINDING <- TRUE
REMOVE_PROMISCUOUS <- TRUE
NORMALIZE <- TRUE

if(!SPEC_BINDING){
  # outcome spectra at time step 0
  s0 <- data.frame(RE=REs[[1]],Norm_F_prob=c(rep(0,2),1,rep(0,13)),model=0)
} else{
  s0 <- data.frame(RE=REs[[2]],Norm_F_prob=c(rep(0,2),1,rep(0,14)),model=0)
  if(REMOVE_PROMISCUOUS){
    if(NORMALIZE){
      pdfv_mc_multistep_sr1 <- lapply(pdfv_mc_multistep_sr1, remove_promiscuous,norm=T)
    }
    else{
      pdfv_mc_multistep_sr1 <- lapply(pdfv_mc_multistep_sr1, remove_promiscuous,norm=F)
    }
    s0 <- data.frame(RE=REs[[1]],Norm_F_prob=c(rep(0,2),1,rep(0,13)),model=0)
  }
}

## plot the dynamics of outcome spectra over time
df_sims_sr1 <- do.call(rbind,pdfv_mc_multistep_sr1) %>% rbind(.,s0)
ERE_SRE_ONLY <- F # whether to color lines only for ERE and SRE

if(ERE_SRE_ONLY){
  ere_sre <- df_sims_sr1 %>% filter(RE %in% c("ERE (GT)","SRE (AA)"))
  df_sims_sr1 %>%
    ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
    geom_line(linewidth=1.3) + 
    geom_line(data=ere_sre,linewidth=1.8,aes(color=RE)) + 
    scale_color_manual(values = hex_RE_colors(3)) + 
    theme_classic() +
    labs(x="Substitution step",y="Probability",title = paste("AncSR1:",REF_GENOTYPE),color="DNA specificity",fill="DNA specificity") +
    scale_x_continuous(breaks=seq(0,mc_iter,5),labels=seq(0,mc_iter,5)) +
    theme(axis.title = element_text(size=15),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13)) +
    guides(color = guide_legend(override.aes = list(size = 5)))
} else{
  df_sims_sr1 %>%
    ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
    geom_line(linewidth=1.3,aes(color=RE)) + 
    scale_color_manual(values = hex_RE_colors(1)) +
    theme_classic() +
    labs(x="Substitution step",y="Probability",title = paste("AncSR1:",REF_GENOTYPE),color="DNA specificity",fill="DNA specificity") +
    scale_x_continuous(breaks=seq(0,mc_iter,5),labels=seq(0,mc_iter,5)) +
    theme(axis.title = element_text(size=15),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13)) +
    guides(color = guide_legend(override.aes = list(size = 5)))
}
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
knitr::kable(df_sims_sr1 %>% filter(model %in% c(3,8,25) & RE %in% c("ERE (GT)","SRE (AA)")) %>% dplyr::rename(Probability=Norm_F_prob,step=model))
```

| RE       |  Probability|  step|
|:---------|------------:|-----:|
| ERE (GT) |    1.0000000|     3|
| SRE (AA) |    0.0000000|     3|
| ERE (GT) |    0.9000039|     8|
| SRE (AA) |    0.0008318|     8|
| ERE (GT) |    0.5699111|    25|
| SRE (AA) |    0.0376695|    25|

Even after 25 amino acid substitutions (6.26 subs / RH site), ERE specificity is the most likely phenotypic outcome. This shows that the probability of evolving a new phenotype from the historical ancestor was extremely low. After 3 steps no new phenotype is likely to evolve, and after 8 steps the probability of evolving a new phenotype is just 0.1 - the probability of SRE in particular is only 0.0008.

The structure of the AncSR1 GP map (global and local biases) is congruent with the patterns of ERE conservation observed in the steroid receptor phylogeny. But it also implies that the evolution of SRE specificity was an extremely unlikely event. As we will see below, the structure of the GP map in AncSR2 can help explain SRE's evolution in the kSR lineage.

------------------------------------------------------------------------

## AncSR2 GP map

### Structure of the AncSR2 GP map (and evolution of the GP map)

We will now turn to analyze the AncSR2 GP map. Following a similar structure as before, we fill first assess whether the AncSR2 map is isotropic and homogeneous, and compare its structure to that of the AncSR1 map.

First, we can see that the AncSR2 GP map has &gt;20X more functional variants than the AncSR1 GP map. The amino acid profiles of these functional variants are very similar to those in AncSR1, but the profile is less stringent in AncSR2.

``` r
# Number of functional variants per GP map
knitr::kable(rbind(data.frame(Bg="AncSR1",total_fxnal_vars = length(phenotypes_tbl_prot %>% filter(bg == "AncSR1") %>% pull(AA_var)), fxnal_vars_in_net = length(extract_main_ntwrk(net_sr1,nodes = T,tr_mat = NULL))),
                   data.frame(Bg="AncSR2",total_fxnal_vars = length(phenotypes_tbl_prot %>% filter(bg == "AncSR2") %>% pull(AA_var)), fxnal_vars_in_net = length(extract_main_ntwrk(net_sr2,nodes = T,tr_mat = NULL)))))
```

| Bg     |  total\_fxnal\_vars|  fxnal\_vars\_in\_net|
|:-------|-------------------:|---------------------:|
| AncSR1 |                 107|                    97|
| AncSR2 |                2407|                  2402|

``` r
# Logo plots functional variants
# AncSR1
names(func_vars_sr1_plot) <- "Functional variants AncSR1"

# AncSR2
func_vars_sr2_plot <- list(func_vars_sr2)
names(func_vars_sr2_plot) <- "Functional variants AncSR2"

p1 <- ggseqlogo::ggseqlogo(func_vars_sr1_plot,ncol=3,method='p',seq_type='aa') + theme(legend.position = "none")
p2 <- ggseqlogo::ggseqlogo(func_vars_sr2_plot,ncol=3,method='p',seq_type='aa') + theme(legend.position = "none")

p1 + p2
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/gpmap_sr2_structure_1-1.png)

``` r
# Correlation between amino acid profiles of funcitonal variants between backgrounds
# Amino acid frequencies per site for functional variants
freq_matrix_sr1_func <- shannon_entropy(func_vars_sr1,which="freq")
freq_matrix_sr2_func <- shannon_entropy(func_vars_sr2,which="freq")

# compute inner product correlation for matrices
inner_product_corr <- function(matrix1, matrix2){
  # Compute inner product correlation between two matrices (https://angeloyeo.github.io/2019/08/20/correlation_and_inner_product_en.html)
  mat_crossprod <- function(matrix1, matrix2){matrix1 %*% t(matrix2)}
  sum(mat_crossprod(matrix1,matrix2)) / (sqrt(sum(mat_crossprod(matrix1,matrix1))) * sqrt(sum(mat_crossprod(matrix2,matrix2))))
}

obs_genetic_similarity_mat <- inner_product_corr(freq_matrix_sr1_func,freq_matrix_sr2_func)
print(paste("Similarity of amino acid profiles for functional variants (r) = ", round(obs_genetic_similarity_mat,2)))
```

    ## [1] "Similarity of amino acid profiles for functional variants (r) =  0.89"

Second, the connectivity between functional genotypes increased by 3X. Since the functional variants in AncSR2 are similar to those in AncSR1, the genotypes that remain functional along the AncSR1-AncSr2 branch gain an average of 11 neighbors.

``` r
# connectivity: number of neighbors per genoytpe
knitr::kable(data.frame(bg=c("AncSR1","AncSR2"),avg_connectivity=c(mean(degree(net_sr1)),mean(degree(net_sr2)))))
```

| bg     |  avg\_connectivity|
|:-------|------------------:|
| AncSR1 |           3.495327|
| AncSR2 |          10.737017|

``` r
rbind(data.frame(bg="AncSR1",n_neighbors=degree(net_sr1)),
            data.frame(bg="AncSR2",n_neighbors=degree(net_sr2))) %>% 
  ggplot(aes(x=bg,y=n_neighbors)) +
  geom_violin() +
  geom_pointrange(stat = "summary",size=0.4,color="red") + 
  theme_classic() + labs(x="Network",y="Number of neighbors\nper genotype") +
  theme(axis.title = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=13,angle=45,hjust = 1,vjust = 1))
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
wilcox.test(degree(net_sr2),degree(net_sr1))
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  degree(net_sr2) and degree(net_sr1)
    ## W = 240663, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0

``` r
# Increase in connectivity
# For genotypes that remain functional from AncSR1 to AncSR2, compute the number of *new* neighbors gained 
# functional RH variants in both backgrounds
func_vars_both <- func_vars_sr2[func_vars_sr2 %in% func_vars_sr1]

neighbors_gained <- foreach(i = 1:length(func_vars_both), .combine = "rbind") %do% {
  n_neigh_sr1 <- names(neighbors(graph = net_sr1,v = func_vars_both[i])) # neighbors in AncSR1 map
  n_neigh_sr2 <- names(neighbors(graph = net_sr2,v = func_vars_both[i])) # neighbors in AncSR2 map
  data.frame(AA_var = func_vars_both[i], new_neighbors = length(n_neigh_sr2[!(n_neigh_sr2 %in% n_neigh_sr1)])) # number of new neighbors
}

# plot
neighbors_gained %>% ggplot(aes(x=new_neighbors)) +
  geom_histogram(color="black",fill="gray60",binwidth = 2) + 
  geom_vline(xintercept = mean(neighbors_gained$new_neighbors),linetype="dashed") +
  labs(x="Number of neighbors\ngained in AncSR2", y="Number of RH genotypes in AncSR1") +
  theme_classic() +
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=13))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-23-2.png)

``` r
print(paste("Average new neighbors gained: ", round(mean(neighbors_gained$new_neighbors)),2))
```

    ## [1] "Average new neighbors gained:  11 2"

Third, we see that the AncSR2 GP map is also strongly anisotropic (with a bias of *B* = 0.37). However, 14 specificity phenotypes are now encoded on the map (twice as many as in AncSR1), and the direction of the bias of the production spectrum also changed from favoring ERE to SRE.

``` r
####################
# AncSR2 global production spectrum
###################
SPECIFIC = T

if(SPECIFIC){
  # Only specific binding
  global_spec_sr2_plot <- phenotypes_tbl_prot %>% ungroup() %>% filter(bg=="AncSR2") %>% filter(specific == "YES") %>% 
    reframe(RE = REs[[4]],
            count = table(factor(specificity,levels=REs[[4]]))) %>%
    ggplot(aes(x=RE,y=count,fill=RE)) + 
    geom_bar(stat="identity",position = "dodge",color="black") +
    scale_fill_manual(values = hex_RE_colors(1)) + 
    theme_classic() + 
    labs(x="DNA element",y="Number of\nprotein variants") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size=11)) +
    guides(fill="none") +
    geom_hline(yintercept = round(1321*(1/16)),linetype="dashed",col="black")
} else{
  # With promiscuous
    global_spec_sr2_plot <- phenotypes_tbl_prot %>% filter(bg=="AncSR2") %>% unnest(bound_REs) %>% group_by(specific) %>% 
    reframe(RE = REs[[4]],
            count = table(factor(bound_REs,levels=REs[[4]]))) %>% 
    ggplot(aes(x=RE,y=count,fill=RE,alpha=specific)) + 
    geom_bar(stat="identity",color="black") +
    scale_fill_manual(values = hex_RE_colors(1)) + 
    scale_alpha_manual(values =c(0,1)) +
    theme_classic() + 
    labs(x="DNA element",y="Number of\nprotein variants") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size=11)) +
    guides(fill="none") +
    geom_hline(yintercept = round(4840*(1/16)),linetype="dashed",col="black")
}

global_spec_sr2_plot + annotate("text",x=3,y=450,label = paste("B =",round(GlobalB_sr2,2)))
```

    ## Don't know how to automatically pick scale for object of type <table>.
    ## Defaulting to continuous.

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-24-1.png)

The increase in the number of specificity phenotypes encoded in the AncSR2 map is due entirely to the gain of function of genotypes that were non-functional in the AncSR1 background - there are no instances of where an RH genotype swithces its specificity or where a former promiscuous genotype gains specificity that was not present in the AncSR1 map.

``` r
# Functional effect of 31 background substitutions
possible_AAvar <- do.call(paste,expand.grid(AA_STANDARD,AA_STANDARD,AA_STANDARD,AA_STANDARD)) %>% gsub(" ","",.)
not_in_sr1 <- possible_AAvar[!(possible_AAvar %in% func_vars_sr1)] # non-functional genotypes in AncSR1
not_in_sr2 <- possible_AAvar[!(possible_AAvar %in% func_vars_sr2)] # non-functional genotypes in AncSR2

# Genotypes that remain functional from AncSR1 to AncSR2, and whether they change specificity
functional_effect_df <- full_join(rbind(phenotypes_tbl_prot %>% filter(bg=="AncSR1") %>%
                                          select(AA_var,bg,specificity) %>% mutate(category = ifelse(specificity == "Promiscuous","Promiscuous","Specific")),
                                        data.frame(AA_var=not_in_sr1,bg="AncSR1",specificity="Non_functional",category = "Non_functional")),
                                  rbind(phenotypes_tbl_prot %>% filter(bg=="AncSR2") %>%
                                          select(AA_var,bg,specificity) %>% mutate(category = ifelse(specificity == "Promiscuous","Promiscuous","Specific")),
                                        data.frame(AA_var=not_in_sr2,bg="AncSR2",specificity="Non_functional",category = "Non_functional")), by="AA_var") %>%
  mutate(bg.x = ifelse(is.na(bg.x),"AncSR1",bg.x),
         bg.y = ifelse(is.na(bg.y),"AncSR2",bg.y),
         # assign funcitonal effect along the branch
         functional_effect = case_when((specificity.x == specificity.y) & (specificity.x != "Non_functional" & specificity.y != "Non_functional") ~ "Same specificity",
                                       (specificity.x != specificity.y) & (specificity.x != "Promiscuous" & specificity.y != "Promiscuous") &
                                         (specificity.x != "Non_functional" & specificity.y != "Non_functional") ~ "Switch specificity",
                                       specificity.x == "Promiscuous" & specificity.y != "Promiscuous" & specificity.y != "Non_functional" ~ "Gain specificity",
                                       specificity.x != "Promiscuous" & specificity.x != "Non_functional" & specificity.y == "Promiscuous" ~ "Loss specificity",
                                       specificity.x == "Non_functional" & specificity.y != "Non_functional" ~ "Gain function",
                                       specificity.x != "Non_functional" & specificity.y == "Non_functional" ~ "Loss function")) %>% ungroup() %>%
  filter(!is.na(functional_effect)) # Only keep variants that are functional in either background

# Summary of functional effects along the branch
knitr::kable(functional_effect_df %>% 
  group_by(category.x,category.y,functional_effect) %>%
  reframe(count = n()) %>%
  mutate(total = sum(count),
         prop = count/total) %>%
  arrange(desc(prop)) %>% select(-c(count,total)) %>% dplyr::rename(category_AncSR1 = category.x,category_AncSR2 = category.y))
```

| category\_AncSR1 | category\_AncSR2 | functional\_effect |       prop|
|:-----------------|:-----------------|:-------------------|----------:|
| Non\_functional  | Specific         | Gain function      |  0.5374122|
| Non\_functional  | Promiscuous      | Gain function      |  0.4183547|
| Specific         | Promiscuous      | Loss specificity   |  0.0248036|
| Specific         | Specific         | Same specificity   |  0.0078545|
| Promiscuous      | Promiscuous      | Same specificity   |  0.0057875|
| Specific         | Non\_functional  | Loss function      |  0.0049607|
| Promiscuous      | Specific         | Gain specificity   |  0.0008268|

``` r
## barplot for promiscuous + specific genotypes on AncSR2 separated by new/old phenotypes
REs_spec_sr1 <- data.frame(summary_sr1) %>% rownames_to_column(var="RE") %>% filter(Encoded == T & Specific == T) %>% pull(RE) # phenotypes encoded in AncSR1 map
REs_spec_sr2 <- data.frame(summary_sr2) %>% rownames_to_column(var="RE") %>% filter(Encoded == T & Specific == T) %>% pull(RE) # phenotypes encoded in AncSR2 map
new_REs_sr2 <- REs_spec_sr2[!(REs_spec_sr2 %in% REs_spec_sr1)] # Phenotypes in AncSR2 and not in AncSR1 (new phenotypes)


functional_effect_df %>% filter(category.y != "Non_functional") %>%
  mutate(phenotype = case_when(specificity.y %in% new_REs_sr2 ~ "New phenotype",
                                   !(specificity.y %in% new_REs_sr2) & specificity.y != "Promiscuous" ~ "Existing phenotype",
                                   specificity.y == "Promiscuous" ~ "Promiscuous")) %>%
  group_by(category.x,category.y,functional_effect,phenotype) %>%
  reframe(count = n()) %>%
  mutate(total = sum(count),
         prop = count/total,
         functional_effect = factor(functional_effect, levels = c("Same specificity","Switch specificity","Gain specificity","Loss specificity","Gain function","Loss function"))) %>%
  arrange(factor(phenotype, levels = c("Existing phenotype","New phenotype","Promiscuous"))) %>%
  mutate(category.x = forcats::as_factor(category.x)) %>% ungroup() %>%
  group_by(phenotype) %>%
  mutate(total_prop = sum(prop),
         norm_prop = prop/total_prop) %>%
  ggplot(aes(x=phenotype,y=norm_prop,fill=category.x)) +
  scale_fill_manual(values=c("gray45","blue","red")) +
  labs(x="Phenotype class in AncSR2",y="Fraction of RH\ngenotypes in AncSR1",fill="") +
  geom_rect(xmin = 0.5, xmax = 2.5,   ymin = -0.02, ymax = 1.1, fill = NA,color="black",linetype=2,size=0.2) +
  annotate("text", x = 1.5, y = 1.05, label = "Specific") +
  geom_bar(stat="identity",color="black",size=0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(size=10,angle = 45,vjust = 1,hjust = 1))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-25-1.png)

``` r
# How are the encoded specificities in the AncSR2 map being generated?
knitr::kable(functional_effect_df %>% filter(specificity.y %in% union(REs_spec_sr1,REs_spec_sr2)) %>%
  group_by(specificity.y) %>% 
  reframe(functional_category = c("Same specificity","Switch specificity","Gain specificity","Loss specificity","Gain function","Loss function"),
          count = table(factor(functional_effect,levels = c("Same specificity","Switch specificity","Gain specificity","Loss specificity","Gain function","Loss function")))) %>%
  ungroup() %>% mutate(new_phenotype = ifelse(specificity.y %in% new_REs_sr2,T,F)) %>%
  filter(count > 0) %>% dplyr::rename(phenotype_AncSR2 = specificity.y,n_RH_vars = count) %>% arrange(new_phenotype))
```

| phenotype\_AncSR2 | functional\_category |  n\_RH\_vars| new\_phenotype |
|:------------------|:---------------------|------------:|:---------------|
| AC                | Same specificity     |            1| FALSE          |
| AC                | Gain function        |           68| FALSE          |
| AT                | Gain function        |           57| FALSE          |
| CA                | Same specificity     |            3| FALSE          |
| CA                | Gain function        |          163| FALSE          |
| ERE (GT)          | Same specificity     |           12| FALSE          |
| ERE (GT)          | Gain function        |           26| FALSE          |
| GA                | Same specificity     |            2| FALSE          |
| GA                | Gain function        |          353| FALSE          |
| SRE (AA)          | Gain specificity     |            2| FALSE          |
| SRE (AA)          | Gain function        |          513| FALSE          |
| TA                | Same specificity     |            1| FALSE          |
| TA                | Gain function        |           26| FALSE          |
| AG                | Gain function        |           54| TRUE           |
| CG                | Gain function        |            9| TRUE           |
| CT                | Gain function        |           11| TRUE           |
| GC                | Gain function        |            9| TRUE           |
| GG                | Gain function        |            1| TRUE           |
| TG                | Gain function        |            2| TRUE           |
| TT                | Gain function        |            8| TRUE           |

Fourth, while genotype clusters are still present in the AncSR2 map, a lower proportion of clusters are enriched for single phenotypes. As a consequence, bias within clusters is on average weaker than in the AncSR1 map, &gt;50% of genotypes can access between 1-4 novel phenotypes in a single mutation, the local bias of the one-mutational neighborhoods are lower, and a higher fraction of direct phenotypic transitions are possible than in the AncSR1 map.

``` r
####################
# Detect genotype clusters in network
####################

cl_sr2 <- cluster_edge_betweenness(net_sr2)
```

``` r
# assocate genotypes with clusters
modules_sr2 <- data.frame(AA_var=names(membership(cl_sr2)),cluster=as.integer(membership(cl_sr2)),bg="AncSR2") %>%
  inner_join(.,phenotypes_tbl_prot,by=c("bg","AA_var")) %>% select(bg,AA_var,cluster,specificity,bound_REs) %>%
  inner_join(.,data.frame(AA_var=nodes_in_ntwrk_sr2,type="net"),by="AA_var") %>% arrange(cluster) %>% as_tibble()

# select largest clusters that contain 90% of all variants in the network
modules_sr2 %>% group_by(cluster) %>% reframe(count=n()) %>% arrange(desc(count)) %>% ungroup() %>% mutate(total = sum(count)) %>% 
  group_by(cluster) %>% mutate(frac = count/total) %>% ungroup() %>% mutate(cum_frac = cumsum(frac),cluster=factor(cluster,levels = cluster)) %>% 
  ggplot(aes(x=cluster,y=cum_frac)) + geom_point() + geom_hline(yintercept = 0.9) + labs(y="Cumulative fraction of variants in network") + theme_classic()
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
# # Associate individual genotypes to each cluster for clusters that capture 90% of variants in the network
clusters_90_sr2 <- modules_sr2 %>% group_by(cluster) %>% reframe(count=n()) %>% arrange(desc(count)) %>% ungroup() %>% mutate(total = sum(count)) %>% 
  group_by(cluster) %>% mutate(frac = count/total) %>% ungroup() %>% mutate(cum_frac = cumsum(frac)) %>%
  mutate(cluster_new = rank(factor(cluster,levels=cluster)),
         cluster_new = ifelse(cum_frac <= 0.9, cluster_new,26),
         cluster_new = map_chr(cluster_new,.f=function(i) LETTERS[i])) %>% 
  inner_join(.,modules_sr2,by="cluster") %>% select(-c(count,total,cum_frac,frac))

# Frequency of phenotypes per cluster
phenotypes_sr2 <- data.frame(summary_sr2) %>% rownames_to_column(var = "RE") %>% filter(Encoded == TRUE & Specific == TRUE) %>% pull(RE)
pheno_freqs_modules_sr2 <- clusters_90_sr2 %>% group_by(cluster_new) %>% 
  reframe(RE = phenotypes_sr2, 
          RE = factor(RE, levels(REs[[4]])),
          count_t = n(),
          count_p = table(factor(specificity,levels=phenotypes_sr2)),
          freq_p = count_p / count_t) %>% acast(cluster_new~RE, value.var="freq_p")

# Compare to global production spectrum
freq_in_net_sr2 <- remove_promiscuous(var.prop_spec_AncSR2_df,norm = T) %>% mutate(RE, levels(REs[[1]]),model="Global") %>% 
  filter(Norm_F_prob > 0) %>% mutate(RE = factor(RE,levels(REs[[4]]))) %>% acast(model~RE, value.var="Norm_F_prob")

pheno_freqs_modules_sr2 <- rbind(pheno_freqs_modules_sr2,freq_in_net_sr2) %>% as.matrix(.)
pheno_freqs_modules_sr2 <- t(apply(pheno_freqs_modules_sr2, 1, function(x) x/sum(x))) # renormalize frequencies

# Detect enrichment of phenotypes per cluster using Fisher's exact test
module_sr2_ids <- unique(clusters_90_sr2$cluster_new)
REs_in_net_sr2 <- unique(clusters_90_sr2 %>% filter(specificity != "Promiscuous") %>% pull(specificity))

m_enrichment_sr2 <- matrix(NA,ncol=length(REs_in_net_sr2),nrow=length(module_sr2_ids))
rownames(m_enrichment_sr2) <- module_sr2_ids
colnames(m_enrichment_sr2) <- REs_in_net_sr2

for(i in 1:length(module_sr2_ids)){
  #modules_sr2 <- modules_sr2 %>% unnest(bound_REs)
  for(j in 1:length(REs_in_net_sr2)){
    a <- clusters_90_sr2 %>% filter(cluster_new == module_sr2_ids[i] & specificity == REs_in_net_sr2[j]) %>% reframe(c=n()) %>% pull(c) # RE j in module i
    b <- clusters_90_sr2 %>% filter(cluster_new == module_sr2_ids[i] & specificity != REs_in_net_sr2[j]) %>% reframe(c=n()) %>% pull(c) # Other REs in module i
    c <- clusters_90_sr2 %>% filter(cluster_new != module_sr2_ids[i] & specificity == REs_in_net_sr2[j]) %>% reframe(c=n()) %>% pull(c) # RE j in other modules modules 
    d <- clusters_90_sr2 %>% filter(cluster_new != module_sr2_ids[i] & specificity != REs_in_net_sr2[j]) %>% reframe(c=n()) %>% pull(c) # Other REs in other modules modules
    
    # perform fisher's exact test on counts
    x <- fisher.test( matrix(c(a,b,c,d),ncol=2,nrow=2,byrow = T),alternative = "greater")
    m_enrichment_sr2[i,j] <- x$p.value
    
    # Perform fsher's exact test using hypergeometric distribution (one-sided p): see 'https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/'
    #m <- a + b      # Nodes in cluster
    #n <- c + d      # Nodes NOT IN cluster
    #k <- a + c      # Node hits, that is, from RE j

    #m_enrichment_sr2[i,j] <- phyper(a, m, n, k, lower.tail = FALSE,log.p = FALSE) # probability of observing 'a' nodes or more with REj in module i
  }
}

# correct for multiple testting (Bonferroni)
m_enrichment_sr2 <- matrix(p.adjust(m_enrichment_sr2,method="bonferroni"),nrow = dim(m_enrichment_sr2)[1],
                       ncol = dim(m_enrichment_sr2)[2],byrow = F, dimnames = list(rownames(m_enrichment_sr2),colnames(m_enrichment_sr2)))

m_enrichment_sr2 <- m_enrichment_sr2[,order(match(colnames(m_enrichment_sr2),REs[[1]]))] # order columns

# compare the number of enriched clusters (for a single phenotype; p<0.05) between AncSR1 and AncSR2
knitr::kable(rbind(data.frame(cluster = rownames(m_enrichment), n_enriched_phenos = apply(m_enrichment,1,function(x) length(x[x<0.05])),bg="AncSR1",n_clusters=nrow(m_enrichment)),
      data.frame(cluster = rownames(m_enrichment_sr2), n_enriched_phenos = apply(m_enrichment_sr2,1,function(x) length(x[x<0.05])),bg="AncSR2",n_clusters=nrow(m_enrichment_sr2)-1)) %>% 
  filter(n_enriched_phenos == 1) %>% group_by(bg) %>% reframe(frac_enriched_clusters = n()/n_clusters) %>% distinct(.,bg,.keep_all = TRUE))
```

| bg     |  frac\_enriched\_clusters|
|:-------|-------------------------:|
| AncSR1 |                 0.8333333|
| AncSR2 |                 0.5000000|

``` r
####################
# Phenotype distribution per cluster and bias
###################

# Compute bias for all clusters
c_bias_sr2 <- foreach(i = 1:dim(pheno_freqs_modules_sr2)[1], .combine = 'c') %do% {
  round(get_bias(pheno_freqs_modules_sr2[i,],Bg = "AncSR2",input_freqs = TRUE),2)
} 

bp_cl <- as_tibble(pheno_freqs_modules_sr2,rownames = NA) %>% rownames_to_column(var="cluster") %>%
  reshape2::melt() %>% 
  arrange(value) %>% 
  mutate(cluster = factor(cluster, levels = c(LETTERS[1:14],LETTERS[26],"Global"))) %>%
  ggplot(aes(x=cluster,y=value,fill=variable)) + geom_bar(stat="identity",width=0.95) + 
  scale_fill_manual(values = hex_RE_colors(1)) +
  theme_classic() +
  labs(x="Genotype cluster",y="Phenotype frequency",fill="") +
  theme(axis.title = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13,angle=45,hjust = 1,vjust = 1),
        legend.position="bottom")
```

    ## Using cluster as id variables

``` r
b_cl <- data.frame(cl = rownames(pheno_freqs_modules_sr2),bias = c_bias_sr2) %>%
  mutate(cl = factor(cl, levels = c(LETTERS[1:14],LETTERS[26],"Global"))) %>%
  ggplot(aes(x = cl,y=bias)) + geom_bar(stat="identity",width=0.95,fill="gray70") +
  theme_classic() + labs(y="Bias (B)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank()) +
  geom_text(aes(label=bias), vjust=1.5, size=2.5) +
  scale_y_continuous(limits = c(0,1))

(b_cl / bp_cl) + plot_layout(ncol=1,nrow=2,heights = c(1,4))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-27-2.png)

``` r
# mean cluster Bias for each GP map
knitr::kable(data.frame(bg=c("AncSR1","AncSR2"),mean_cluster_B = c(mean(c_bias),mean(c_bias_sr2))))
```

| bg     |  mean\_cluster\_B|
|:-------|-----------------:|
| AncSR1 |         0.7842857|
| AncSR2 |         0.5812500|

``` r
###################################
# local bias: One-mutation neighborhood
###################################

cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

# Number of phenotypes accessible per genoytpe
one_mut_nei_sr2 <- foreach(i = 1:length(nodes_in_ntwrk_sr2_spec), .combine = 'rbind') %dopar% {
  # focal variant and phenotype
  f_var <- phenotypes_tbl_prot %>% filter(bg=="AncSR2") %>% filter(AA_var %in% nodes_in_ntwrk_sr2_spec[i]) %>% select(AA_var,specificity)
  
  # neighbor variants and their phenotypes
  n_i <- names(neighbors(graph = net_sr2,v = nodes_in_ntwrk_sr2_spec[i]))
  p_n_i <- phenotypes_tbl_prot %>% ungroup() %>% filter(bg=="AncSR2") %>% filter(AA_var %in% n_i) %>% select(AA_var,specificity) %>%
    mutate(focal_var = f_var$AA_var, focal_pheno = f_var$specificity)
}

# local bias computation
local_bias_sr2 <- foreach(i = 1:length(nodes_in_ntwrk_sr2_spec), .combine = 'rbind') %dopar% {
  # extract neighbors
  n_i <- names(neighbors(graph = net_sr2,v = nodes_in_ntwrk_sr2_spec[i]))
  # frequency of phenotypes in one-mutational neighborg¡hood
  p_n_i <- phenotypes_tbl_prot %>% ungroup() %>% filter(bg=="AncSR2") %>% filter(AA_var %in% n_i) %>%
    reframe(RE = REs[[1]],
            count = table(factor(specificity,levels=REs[[1]])),
            fr = count/length(n_i))
  # local bias
  local_b_i <- get_bias(as.vector(p_n_i$fr),input_freqs = T)
  data.frame(g=nodes_in_ntwrk_sr2_spec[i],b=local_b_i)}

stopCluster(cl)
```

``` r
# Numberof accessible and *new* phenotypes per genotype
accessible_phenotypes_singleStep_sr2 <- one_mut_nei_sr2 %>% filter(specificity != "Promiscuous") %>% 
  mutate(new_pheno = ifelse(focal_pheno == specificity,FALSE,TRUE)) %>%
  group_by(focal_var) %>% 
  reframe(n_acc_pheno = n_distinct(specificity),
          n_new_pheno = n_distinct(specificity) - any(specificity == focal_pheno)) %>% 
  ungroup() 

# plot number of new phenotypes accessible per genotype
p1 <- accessible_phenotypes_singleStep_sr2 %>%
  ggplot(aes(x=n_new_pheno)) + 
  geom_bar(aes(y = after_stat(count))) +
  labs(x="Number of new\nphenotypes accessible",y="Number of specific RH genotypes") +
  geom_text(aes(label = scales::percent(after_stat(count)/sum(after_stat(count))),y= after_stat(count)), stat= "count", vjust = -.5,size=3) +
  theme_classic() + 
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=12))
  
# Local bias per genotype comparison
local_b_one_mut_neighborhood <- rbind(local_bias_sr1 %>% mutate(bg="AncSR1"),local_bias_sr2 %>% mutate(bg="AncSR2"))

p2 <- local_b_one_mut_neighborhood %>%
  ggplot(aes(x=bg,y=b)) +
  geom_violin() +
  geom_pointrange(stat = "summary",size=0.4,col="red") +
  labs(x="Network",y="Phenotype Bias (B)") +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25)) +
  theme_classic() +
  theme(axis.title = element_text(size=10),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13,hjust = 1,vjust = 1,angle = 45))

p1 + p2
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-29-1.png)

``` r
# average local bias of one-mutant neighborhood is lower in AncSR2
wilcox.test(local_b_one_mut_neighborhood$b~local_b_one_mut_neighborhood$bg)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  local_b_one_mut_neighborhood$b by local_b_one_mut_neighborhood$bg
    ## W = 78194, p-value = 3.865e-11
    ## alternative hypothesis: true location shift is not equal to 0

``` r
###################################
# Distance between neutral networks
###################################

RE_combos <- expand.grid(REs[[1]],REs[[1]])

# Number of direct links between phenotypes
direct_links_sr2 <- RE_combos %>% mutate(links = map2_dbl(.x=as.character(Var1),.y=as.character(Var2),
                                                    pairwise_neutral_network_proximity,Bg="AncSR2",graph=net_sr1,type=2,
                                                    pheno_df=phenotypes_tbl_prot,from_specific=T,to_specific=T)) %>% acast(Var1~Var2, value.var="links")

direct_links_sr2 <- direct_links_sr2[rowSums(is.na(direct_links_sr2)) != ncol(direct_links_sr2), ] # remove all NAs-rows
direct_links_sr2 <- direct_links_sr2[, colSums(is.na(direct_links_sr2)) != nrow(direct_links_sr2)] # remove all NAs-cols
direct_links_sr2 <- direct_links_sr2[order(match(rownames(direct_links_sr2),REs[[1]])),order(match(colnames(direct_links_sr2),REs[[1]]))] # order rows and columns
diag(direct_links_sr2) <- NA

# Average distance between phenotypes
pairwise_dists_sr2 <- RE_combos %>% mutate(links = map2_dbl(.x=as.character(Var1),.y=as.character(Var2),
                                                    pairwise_neutral_network_proximity,Bg="AncSR2",graph=net_sr2,type=3,
                                                    pheno_df=phenotypes_tbl_prot,from_specific=T,to_specific=T)) %>% acast(Var1~Var2, value.var="links")

pairwise_dists_sr2 <- pairwise_dists_sr2[rowSums(is.na(pairwise_dists_sr2)) != ncol(pairwise_dists_sr2), ] # remove all NAs-rows
pairwise_dists_sr2 <- pairwise_dists_sr2[, colSums(is.na(pairwise_dists_sr2)) != nrow(pairwise_dists_sr2)] # remove all NAs-cols
pairwise_dists_sr2 <- pairwise_dists_sr2[order(match(rownames(pairwise_dists_sr2),REs[[1]])),order(match(colnames(pairwise_dists_sr2),REs[[1]]))] # order rows and columns
diag(pairwise_dists_sr2) <- NA
```

``` r
# Fraction of direct phenotypic transitions that are possible (comparison between GP maps)
knitr::kable(rbind(as.data.frame.table(direct_links_sr1) %>% filter(Var1 != Var2) %>% mutate_at(vars(Var1,Var2), list(as.character)) %>%
        rowwise() %>% mutate(combo = paste(min(Var1,Var2),max(Var1,Var2),sep="_"),bg="AncSR1") %>% distinct(.,combo,.keep_all = TRUE) %>%
        ungroup() %>% mutate(total_pairs = n(), possible_transitions = sum(Freq)),
      as.data.frame.table(direct_links_sr2) %>% filter(Var1 != Var2) %>% mutate_at(vars(Var1,Var2), list(as.character)) %>%
        rowwise() %>% mutate(combo = paste(min(Var1,Var2),max(Var1,Var2),sep="_"),bg="AncSR2") %>% distinct(.,combo,.keep_all = TRUE) %>%
        ungroup() %>% mutate(total_pairs = n(), possible_transitions = sum(Freq))) %>% 
  distinct(.,bg,.keep_all = TRUE) %>% select(bg:possible_transitions) %>% mutate(frac_possible_transitions = possible_transitions/total_pairs))
```

| bg     |  total\_pairs|  possible\_transitions|  frac\_possible\_transitions|
|:-------|-------------:|----------------------:|----------------------------:|
| AncSR1 |            15|                      6|                    0.4000000|
| AncSR2 |            91|                     55|                    0.6043956|

Finally, all these changes in the accessibility of phenotypic variation can be attributed to a higher fraction of non-neutral neighbors in the AncSR2 map. While the most likely outcome of single mutations is still conservation (genotypes still have more neutral neighbors than expected), the average fraction of neutral neighbors is lower than in AncSR1: decreases from 79% to 54%.

``` r
###################################
# fraction of neighbor nodes with same or different phenotype 
###################################

cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

# compute fraction of neighbors with same or different phenotype as focal variant
neighbors_sr2 <- foreach(i = 1:length(nodes_in_ntwrk_sr2_spec), .combine = 'rbind') %dopar% {
  fraction_of_neighbors_per_genotype(from_g = nodes_in_ntwrk_sr2_spec[i],Bg = "AncSR2",graph = net_sr2,pheno_tbl = phenotypes_tbl_prot)}

# To assess whether the avearge fraction of neutral neighbors is higher than expected, perform a non-parametric bootstrap test by 
# randomizing phenotypic annotations for RH variants, but keeping same network structure.
# Compute mean fraction of neutral neighbors for each scrambled GPmap and compute p-value as the fraction of bootstrap samples
# that are higher than the observed value.
N_BOOT <- 50

mean_prop_same_neighbors_permut_sr2 <- foreach(i = 1:N_BOOT) %:%
  foreach(j = 1:length(nodes_in_ntwrk_sr2_spec), .combine = 'rbind') %dopar% {
    pheno_tbl_permut <- phenotypes_tbl_prot %>% filter(bg=="AncSR2") %>% transform(.,specificity = sample(specificity,replace = F))
    fraction_of_neighbors_per_genotype(from_g = nodes_in_ntwrk_sr2_spec[j],Bg = "AncSR2",
                                       graph = net_sr2,pheno_tbl = pheno_tbl_permut)} %>%
  map_dbl(.f = function(x) x %>% filter(pheno_match == T) %>% with(mean(prop)))

stopCluster(cl)
```

``` r
# compute p-value
obs_mean <- neighbors_sr2 %>% filter(pheno_match == T) %>% with(mean(prop))
norm_fit_permut <- fitdistrplus::fitdist(mean_prop_same_neighbors_permut_sr2, "norm") # fit normal distribution to permuted dataset
#hist(mean_prop_same_neighbors_permut_sr2,xlim = c(0.3,0.6),freq = F,ylim = c(0,100))
#lines(x=seq(0.3,0.6,0.01),y=dnorm(x=seq(0.3,0.6,0.01),mean = norm_fit_permut$estimate[1],sd = norm_fit_permut$estimate[2]),col="red")

pval_neutral <- pnorm(obs_mean, mean = norm_fit_permut$estimate[1], sd = norm_fit_permut$estimate[2], lower.tail = FALSE, log.p = FALSE) # p-value
#sum(mean_prop_same_neighbors_permut>=obs_mean)/length(mean_prop_same_neighbors_permut) # empirical p-val

print(paste("Mean fraction of neutral neighbors in AncSR2 network: ", round(obs_mean,2),". P-value = ",pval_neutral))
```

    ## [1] "Mean fraction of neutral neighbors in AncSR2 network:  0.54 . P-value =  2.89084929865249e-212"

``` r
# Compare fraction of neutral neighbors between GP maps
neutral_neighbors_df <- rbind(neighbors_sr1 %>% filter(pheno_match == T) %>% mutate(bg="AncSR1"),
                             neighbors_sr2 %>% filter(pheno_match == T) %>% mutate(bg="AncSR2"))

knitr::kable(neutral_neighbors_df %>% group_by(bg) %>% reframe(mean_prop_nonneutral_neighbors = mean(1-prop)))
```

| bg     |  mean\_prop\_nonneutral\_neighbors|
|:-------|----------------------------------:|
| AncSR1 |                          0.2106295|
| AncSR2 |                          0.4615453|

``` r
neutral_neighbors_df %>%
  ggplot(aes(x=bg,y=prop)) + 
  geom_violin() +
  geom_pointrange(stat = "summary",size=0.4,col="red") + theme_classic() + 
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,0.25)) +
  labs(x="Network",y="Proportion of neutral neighbors per genotype") +
  #geom_hline(yintercept = 0.5, linetype="dashed") +
  theme(axis.title = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust = 1,vjust = 1))
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-33-1.png)

``` r
wilcox.test(neutral_neighbors_df$prop~neutral_neighbors_df$bg)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  neutral_neighbors_df$prop by neutral_neighbors_df$bg
    ## W = 70018, p-value = 5.492e-15
    ## alternative hypothesis: true location shift is not equal to 0

Overall, these results show that the background substitutions that occurred along the AncSR1-AncSR2 lineage, dramatically changed the structure of the RH GP map: they increased the size of the network, its navigability, the number of encoded phenotypes, and changed the direction of the global bias.

### Phenotypic evolution on the AncSR2 GP map

Now, let's investigate how these changes in the GP map affected phenotypic evolution. The main consequences we observed are threefold: 1) SRE is the most likely phenotypic outcome at equilibrium; 2) SRE is the most liekly phenotypic outcome at moderate timescales across the majority of genotypes and phenotypes; and 3) the probability of phenotypic diversication is higher.

``` r
##################
# Production spectrum and equilibrium outcome spectrum
##################
df_long_term_sr2 <- inner_join(remove_promiscuous(var.prop_spec_AncSR2_df,norm = T),remove_promiscuous(stationary_PDFV_spec_sr2,norm = T),by="RE")

# correlation between spectra
cor <- df_long_term_sr2 %>% with(cor(Norm_F_prob.x,Norm_F_prob.y)^2)

# plot correlation of spectra
p1 <- df_long_term_sr2 %>%
  ggplot(aes(x=Norm_F_prob.x,y=Norm_F_prob.y,fill=RE)) + 
  geom_point(shape=21,col="black",size=4) + scale_fill_manual(values = hex_RE_colors(1)) +
  scale_y_continuous(limits = c(0,0.45)) + scale_x_continuous(limits = c(0,0.45)) + 
  geom_abline(intercept = 0,slope = 1,linetype="dashed",col="black") +
  labs(x="Global production spectrum",y="Equilibrium outcome spectrum") + 
  theme_classic() +
  theme(axis.text = element_text(size = 10),axis.title = element_text(size=11)) +
  annotate("text", x = 0.05, y = 0.42, label = bquote(r^2 == .(round(cor,2))))

# plot barplot per spectra
bias_stat_global <- data.frame(model=c("Production bias","Equilibrium"),b=c(GlobalB_sr2,StatB_sr2),RE="TA",y=0.3)
p2 <- bars_PDFV(list(remove_promiscuous(var.prop_spec_AncSR2_df,norm = T) %>% mutate(model="Production bias"),
               remove_promiscuous(stationary_PDFV_spec_sr2,norm = T) %>% mutate(model="Equilibrium")),
          arrange_by = "Production bias",fill_by_RE = T,hex_id = 1) +
  geom_text(data=bias_stat_global,aes(x=RE,y=y,label=paste("B =",round(b,2))))


p1 + p2
```

    ## Warning in is.na(x): is.na() applied to non-(list or vector) of type 'language'

    ## Don't know how to automatically pick scale for object of type <table>.
    ## Defaulting to continuous.

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-34-1.png)

``` r
###################################
# Compute outcome spectra from every genotype after 8 substitutions
###################################
cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

N_MUTATIONS <- 8 # ~ 2 subs per site
SPEC_BINDING <- TRUE # whether to consider only specific binders (promiscuous as a separate "phenotype")
NORM_PROB <- TRUE # whether to remove promiscuous end-points and renormalize probanilities (only if SPEC_BINDING = TRUE)

# Run Markov chain from every genotype in the network for 8 amino acid substitutions and compute the outcome spectra
if(SPEC_BINDING){
  markov_chain_sr2 <- foreach(i = 1:length(nodes_in_ntwrk_sr2), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr2[i],tr_mat = P_drift_sr2_ntwrk,n_steps = N_MUTATIONS)
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR2",graph=net_sr2,model = nodes_in_ntwrk_sr2[i],specific = T) %>%
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[2]]))]
  }
  rownames(markov_chain_sr2) <- nodes_in_ntwrk_sr2
  markov_chain_sr2 <- markov_chain_sr2[, which(colSums(markov_chain_sr2) != 0)]
  
  if(NORM_PROB){
    markov_chain_sr2 <- markov_chain_sr2[,1:ncol(markov_chain_sr2)-1]
    markov_chain_sr2 <- t(apply(markov_chain_sr2, 1, function(x) x/sum(x)))
  }
} else{
  markov_chain_sr2 <- foreach(i = 1:length(nodes_in_ntwrk_sr2), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr2[i],tr_mat = P_drift_sr2_ntwrk,n_steps = N_MUTATIONS)
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR2",graph=net_sr2,model = nodes_in_ntwrk_sr2[i],specific = F) %>%
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
  }
  rownames(markov_chain_sr2) <- nodes_in_ntwrk_sr2
  markov_chain_sr2 <- markov_chain_sr2[, which(colSums(markov_chain_sr2) != 0)]
}

stopCluster(cl)
```

``` r
# Run Markov chain for 8 amino acid substitutions averaged across all genotypes
N_MUTATIONS <- 8 # ~ 2 subs per site
if(SPEC_BINDING){
  mc_avg <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr2,tr_mat = P_drift_sr2_ntwrk,n_steps = N_MUTATIONS)
  outcome_spec_avg_all_sr2 <- get_PDFV_v2(mc_avg,type="simulated mc", Bg="AncSR2",graph=net_sr2,model = "All nodes",specific = T,pheno_tbl = phenotypes_tbl_prot) %>%
        acast(model~RE,value.var="Norm_F_prob")
  outcome_spec_avg_all_sr2 <- outcome_spec_avg_all_sr2[,order(match(colnames(outcome_spec_avg_all_sr2),REs[[3]]))]
  outcome_spec_avg_all_sr2 <- outcome_spec_avg_all_sr2[which(outcome_spec_avg_all_sr2 != 0)]
  if(NORM_PROB){
    outcome_spec_avg_all_sr2 <- outcome_spec_avg_all_sr2[1:length(outcome_spec_avg_all_sr2)-1]
    outcome_spec_avg_all_sr2 <- outcome_spec_avg_all_sr2/sum(outcome_spec_avg_all_sr2)
  }
} else {
  mc_avg <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr2,tr_mat = P_drift_sr2_ntwrk,n_steps = N_MUTATIONS)
  outcome_spec_avg_all_sr2 <- get_PDFV_v2(mc_avg,type="simulated mc", Bg="AncSR2",graph=net_sr2,model = "All nodes",specific = F,pheno_tbl = phenotypes_tbl_prot) %>%
        acast(model~RE,value.var="Norm_F_prob")
  outcome_spec_avg_all_sr2 <- outcome_spec_avg_all_sr2[,order(match(colnames(outcome_spec_avg_all_sr2),REs[[3]]))]
  outcome_spec_avg_all_sr2 <- outcome_spec_avg_all_sr2[which(outcome_spec_avg_all_sr2 != 0)]
}

# Run Markov chain for 8 amino acid substitutions averaged across non-SRE genotypes
nonSRE_genotypes_ntwrk_sr2 <- phenotypes_tbl_prot %>% filter(bg =="AncSR2" & specificity != "SRE (AA)") %>% pull(AA_var)
nonSRE_genotypes_ntwrk_sr2 <- nonSRE_genotypes_ntwrk_sr2[nonSRE_genotypes_ntwrk_sr2 %in% nodes_in_ntwrk_sr2]

if(SPEC_BINDING){
  mc_avg_2 <- simulate_markov_chain(states_0 = nonSRE_genotypes_ntwrk_sr2,tr_mat = P_drift_sr2_ntwrk,n_steps = N_MUTATIONS)
  outcome_spec_avg_all_sr2_2 <- get_PDFV_v2(mc_avg_2,type="simulated mc", Bg="AncSR2",graph=net_sr2,model = "non-SRE nodes",specific = T,pheno_tbl = phenotypes_tbl_prot) %>%
        acast(model~RE,value.var="Norm_F_prob")
  outcome_spec_avg_all_sr2_2 <- outcome_spec_avg_all_sr2_2[,order(match(colnames(outcome_spec_avg_all_sr2_2),REs[[3]]))]
  outcome_spec_avg_all_sr2_2 <- outcome_spec_avg_all_sr2_2[which(outcome_spec_avg_all_sr2_2 != 0)]
  if(NORM_PROB){
    outcome_spec_avg_all_sr2_2 <- outcome_spec_avg_all_sr2_2[1:length(outcome_spec_avg_all_sr2_2)-1]
    outcome_spec_avg_all_sr2_2 <- outcome_spec_avg_all_sr2_2/sum(outcome_spec_avg_all_sr2_2)
  }
} else {
  mc_avg_2 <- simulate_markov_chain(states_0 = nonSRE_genotypes_ntwrk_sr2,tr_mat = P_drift_sr2_ntwrk,n_steps = N_MUTATIONS)
  outcome_spec_avg_all_sr2_2 <- get_PDFV_v2(mc_avg,type="simulated mc", Bg="AncSR2",graph=net_sr2,model = "non-SRE nodes",specific = F,pheno_tbl = phenotypes_tbl_prot) %>%
        acast(model~RE,value.var="Norm_F_prob")
  outcome_spec_avg_all_sr2_2 <- outcome_spec_avg_all_sr2_2[,order(match(colnames(outcome_spec_avg_all_sr2_2),REs[[3]]))]
  outcome_spec_avg_all_sr2_2 <- outcome_spec_avg_all_sr2_2[which(outcome_spec_avg_all_sr2_2 != 0)]
}

# "Admixture plot"-like figure
# Associate genotypes to clusters
clusters_sr2 <- clusters_90_sr2 %>% select(AA_var,specificity,cluster_new) %>% filter(AA_var %in% nodes_in_ntwrk_sr2) %>%
  dplyr::rename("RE" = "specificity")

# reorder RH variants to match clusters_sr1 data frame
markov_chain_sr2_plot <- as_tibble(markov_chain_sr2) %>% mutate(AA_var=nodes_in_ntwrk_sr2) %>% inner_join(.,clusters_sr2,by="AA_var") %>%
  column_to_rownames(var="AA_var") %>% select(-c(RE,cluster_new)) %>% as.matrix(.)

# phenotype of starting genotype
starting_DNAspec_df <- clusters_sr2 %>% mutate(variable=RE,value = -0.05) 

# data frame for plotting
markov_chain_sr2_plot <- as_tibble(markov_chain_sr2_plot,rownames=NA) %>% rownames_to_column(var="AA_var") %>% inner_join(.,clusters_sr2,by="AA_var") %>%
  reshape2::melt() %>%
  rbind(.,starting_DNAspec_df) %>%
  arrange(desc(factor(variable,levels = levels(REs[[4]])))) %>% 
  arrange(factor(RE,levels = levels(REs[[6]]))) %>% arrange(cluster_new) %>%
  mutate(variable = forcats::as_factor(variable),
         cluster_new = forcats::as_factor(cluster_new),
         AA_var = forcats::as_factor(AA_var),
         RE = forcats::as_factor(RE))
```

    ## Using AA_var, RE, cluster_new as id variables

``` r
# plot from every genotype
admx_plot_sr2 <- markov_chain_sr2_plot %>%
  ggplot(aes(x=AA_var,y=value,fill=variable,group=cluster_new)) + geom_bar(stat="identity",width=1) + 
  scale_fill_manual(values = hex_RE_colors(2)) +
  theme_classic() + facet_grid(cols=vars(cluster_new),scales="free_x",space="free_x") + 
  geom_hline(yintercept = 0,linewidth=2,color="white") +
  labs(x="Starting genotype",y="Pr(Ending phenotype)",fill="") +
  theme(axis.title = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.text.x = element_blank(),
        legend.position="bottom")

# plot from average outcomes across genotypes
admx_avgs_sr2 <- rbind(data.frame(model = "All",pheno = names(outcome_spec_avg_all_sr2), prob = outcome_spec_avg_all_sr2),
                       data.frame(model = "non SRE",pheno = names(outcome_spec_avg_all_sr2_2), prob = outcome_spec_avg_all_sr2_2)) %>%
  arrange(factor(pheno,levels = levels(REs[[4]]))) %>% 
  mutate(pheno = forcats::as_factor(pheno)) %>%
  ggplot(aes(x=model,y=prob,fill=pheno)) + geom_bar(stat="identity",width=0.9) + 
  scale_fill_manual(values = hex_RE_colors(2)) +
  theme_classic() +
  geom_hline(yintercept = 0,linewidth=2,color="white") +
  labs(x="",y="Pr(Ending phenotype)",fill="") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        legend.position="none")

admx_plot_sr2 + admx_avgs_sr2 + plot_layout(widths = c(6,0.5))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-36-1.png)

``` r
# Phenotypic outcome with highest likelood across genotypes in the map
print(paste("Most likely outcome of evolution on average: ", names(outcome_spec_avg_all_sr2)[which(outcome_spec_avg_all_sr2 == max(outcome_spec_avg_all_sr2))]))
```

    ## [1] "Most likely outcome of evolution on average:  SRE (AA)"

``` r
print(paste("Most likely outcome of evolution from non-SRE genotypes on average: ", names(outcome_spec_avg_all_sr2_2)[which(outcome_spec_avg_all_sr2_2 == max(outcome_spec_avg_all_sr2_2))]))
```

    ## [1] "Most likely outcome of evolution from non-SRE genotypes on average:  SRE (AA)"

``` r
###################################
# Probabilities of transition aggregated by phenotype
###################################
cl <- makeCluster(N_CORES,type = "FORK",outfile="")
registerDoParallel(cl)

REs_spec_sr2 <- data.frame(summary_sr2) %>% rownames_to_column(var="RE") %>% filter(Encoded == T & In_Net_Spec == T) %>% pull(RE)

# compute phenotypic transitions: outcome spectra aggregated for RH variants encoding the same phenotype
pt_sr2 <- foreach(i = 1:length(REs_spec_sr2), .combine = 'rbind') %dopar% {
  phenotype_vars <- phenotypes_tbl_prot %>% filter(bg == "AncSR2" & specificity == REs_spec_sr2[i]) %>% pull(AA_var)
  tmp_mc <- simulate_markov_chain(states_0 = phenotype_vars,tr_mat = P_drift_sr2_ntwrk,n_steps = N_MUTATIONS)
  tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR2",graph=net_sr2,model = REs_spec_sr2[i],specific = T,pheno_tbl = phenotypes_tbl_prot) %>%
    acast(model~RE,value.var="Norm_F_prob")
  tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
}
rownames(pt_sr2) <- REs_spec_sr2
pt_sr2 <- pt_sr2[, which(colSums(pt_sr2) != 0)]
pt_sr2 <- pt_sr2[,1: ncol(pt_sr2)-1]
pt_sr2 <- t(apply(pt_sr2, 1, function(x) x/sum(x)))

stopCluster(cl)
```

``` r
# Admixture plot-like figure
as_tibble(pt_sr2,rownames=NA) %>% rownames_to_column(var="Start_Pheno") %>%
  reshape2::melt() %>%
  rbind(.,data.frame(Start_Pheno = rownames(pt_sr2),variable=rownames(pt_sr2),value = -0.05)) %>%
  arrange(desc(factor(variable, levels = levels(REs[[4]])))) %>%
  arrange(factor(Start_Pheno, levels = levels(REs[[4]]))) %>%
  mutate(variable = factor(variable, levels = levels(REs[[4]])),
         Start_Pheno = factor(Start_Pheno, levels = levels(REs[[4]]))) %>%
  ggplot(aes(x=Start_Pheno,y=value,fill=variable)) + 
  geom_bar(stat="identity",width=0.95) + 
  scale_fill_manual(values = hex_RE_colors(1)) +
  theme_classic() +
  geom_hline(yintercept = 0,linewidth=2,color="white") +
  labs(x="Starting phenotype",y="Pr(Ending phenotype)",fill="") +
  theme(axis.title = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13,angle=45,hjust = 1,vjust = 1),
        legend.position="bottom")
```

    ## Using Start_Pheno as id variables

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-38-1.png)

``` r
# Most likely phenotypic transition per phenotype
knitr::kable(as.data.frame.table(pt_sr2) %>% filter(Var1 != Var2) %>% group_by(Var1) %>% filter(Freq == max(Freq)) %>% select(Var1,Var2) %>%
  dplyr::rename("Starting phenotype"=Var1,"Most likely transition"=Var2))
```

| Starting phenotype | Most likely transition |
|:-------------------|:-----------------------|
| GA                 | SRE (AA)               |
| ERE (GT)           | SRE (AA)               |
| AC                 | SRE (AA)               |
| AG                 | SRE (AA)               |
| AT                 | SRE (AA)               |
| CA                 | SRE (AA)               |
| CT                 | SRE (AA)               |
| TA                 | SRE (AA)               |
| TG                 | SRE (AA)               |
| TT                 | SRE (AA)               |
| GC                 | GA                     |
| GG                 | GA                     |
| SRE (AA)           | CA                     |
| CG                 | CA                     |

``` r
###################################
# Probability of conservation and transition per genotype
###################################

#associate genotypes to phenotypes
geno_pheno_sr2 <- phenotypes_tbl_prot %>% filter(bg=="AncSR2") %>% select(AA_var,specificity) %>% dplyr::rename("RE" = "specificity")

# extract probabilities of the most likely phenotypic outcome
max_prob <- apply(markov_chain_sr2,1,max) # max prob.
max_prob_pheno <- apply(markov_chain_sr2,1,function(x) colnames(markov_chain_sr2)[which(x == max(x))]) # phenotype with the maximum transition probability

# Most likely phenotypic outcome per *specific* genotype (conservation or transition)
highest_Lik_outcomes_sr2_df <- inner_join(geno_pheno_sr2,data.frame(AA_var=names(max_prob),max_prob=max_prob)) %>% 
  inner_join(.,data.frame(AA_var=names(max_prob_pheno),max_prob_pheno=max_prob_pheno)) %>%
  mutate(hl_outcome = ifelse(RE == max_prob_pheno,"Conservation","Transition")) %>% filter(RE != "Promiscuous")
```

    ## Joining with `by = join_by(AA_var)`
    ## Joining with `by = join_by(AA_var)`

``` r
# Extract probability of conservation and probability of phenotypic transition per genotype
if(SPEC_BINDING && NORM_PROB){
  pheno_outcomes_sr2_df <- as_tibble(markov_chain_sr2,rownames = NA) %>% 
    rownames_to_column(var="AA_var") %>% 
    inner_join(.,geno_pheno_sr2,by="AA_var") %>%
    pivot_longer(cols=2:15,names_to="Pheno",values_to="Prob") %>% 
    filter(RE != "Promiscuous" & Pheno != "Promiscuous") %>%
    group_by(AA_var) %>% filter(RE == Pheno) %>%
    reframe(p_cons = Prob,
            p_tr = 1-p_cons) %>%
    pivot_longer(cols=2:3,names_to="outcome",values_to="prob")
} else{
  pheno_outcomes_sr2_df <- as_tibble(markov_chain_sr2,rownames = NA) %>% 
    rownames_to_column(var="AA_var") %>% 
    inner_join(.,geno_pheno_sr2,by="AA_var") %>%
    pivot_longer(cols=2:16,names_to="Pheno",values_to="Prob") %>% 
    filter(RE != "Promiscuous" & Pheno != "Promiscuous") %>%
    group_by(AA_var) %>% filter(RE == Pheno) %>%
    reframe(p_cons = Prob,
            p_tr = 1-p_cons) %>%
    pivot_longer(cols=2:3,names_to="outcome",values_to="prob")
}

# average probability of conservation across specific genotypes
avg_pr_cons <- pheno_outcomes_sr2_df %>% filter(outcome == "p_cons") %>% with(mean(prob))
print(paste("Average probability of conservation = ",round(avg_pr_cons,2)))
```

    ## [1] "Average probability of conservation =  0.47"

``` r
# percentaje of genotypes for which conservation or transition is the most likely outcome (comparison to AncSR1)
knitr::kable(rbind(highest_Lik_outcomes_sr1_df %>% mutate(bg="AncSR1") %>% ungroup() %>% group_by(bg,hl_outcome) %>%
               reframe(count = n()) %>% ungroup() %>% mutate(total = sum(count), frac = count/total) %>% select(bg,hl_outcome,frac),
             highest_Lik_outcomes_sr2_df %>% mutate(bg="AncSR2") %>% ungroup() %>% group_by(bg,hl_outcome) %>%
               reframe(count = n()) %>% ungroup() %>% mutate(total = sum(count), frac = count/total) %>% select(bg,hl_outcome,frac)))
```

| bg     | hl\_outcome  |       frac|
|:-------|:-------------|----------:|
| AncSR1 | Conservation |  0.8313253|
| AncSR1 | Transition   |  0.1686747|
| AncSR2 | Conservation |  0.7667173|
| AncSR2 | Transition   |  0.2332827|

``` r
# plot
inner_join(pheno_outcomes_sr2_df,highest_Lik_outcomes_sr2_df,by="AA_var") %>%
  filter(outcome == "p_cons") %>%
  ggplot(aes(x=prob,fill=hl_outcome)) + geom_histogram(position = "stack",col="black") + 
  scale_fill_manual(labels = c("Conservation (77%)","Transition (23%)"), values = c("white","gray40")) +
  geom_vline(xintercept = avg_pr_cons,linetype="dashed") +
  theme_classic() + 
  scale_x_continuous(limits = c(0,1)) +
  labs(x="P(conservation) per genotype",y="Frequency",fill="Highest likelihood\noutcome") +
  theme(axis.title = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust = 1,vjust = 1),
        legend.position=c(0.3,0.8))
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 4 rows containing missing values (`geom_bar()`).

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-39-1.png)

All these consequences occur because of the increased connectivity of the genotype network and the reorientation of the global bias. First, with greater connectivity of the network, genotypes should navigate more efficieintly the GP map and therefore gain access more rapidly to more genotypes and phenotypes. To show this effect, we can track the average change in the oucome bias across genotypes in the AncSR2 network.

``` r
###################################
# Change in outcome bias over time
###################################

cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

set.seed(2468) # Make analysis replicable
n_sample = 100
MAX_STEP = 150

# sample 100 genotypes randomly to start evolutionary trajectories
random_starting_nodes <- c(sample(nodes_in_ntwrk_sr2,size = n_sample,replace = F), REF_GENOTYPE)

change_bias_df_sr2 <- data.frame(step=NULL,bias=NULL,AA_var=NULL)

for(s in c(1,seq(5,MAX_STEP,5))){
  # Compute PDFVs for each genotype per step
  tmp_mat_i <- foreach(i = 1:length(random_starting_nodes), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = random_starting_nodes[i],tr_mat = P_drift_sr2_ntwrk,n_steps = s)
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR2",graph=net_sr2,model = random_starting_nodes[i],specific = T) %>%
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[2]]))]
  }
  rownames(tmp_mat_i) <- random_starting_nodes
  tmp_mat_i <- tmp_mat_i[, which(colSums(tmp_mat_i) != 0)]
  
  # compute entropy per genotype per step (without promiscuous)
  tmp_mat_i <- t(apply(tmp_mat_i[,1:ncol(tmp_mat_i)-1],1,function(x) x/sum(x))) # renormalize probs
  B_g <- apply(tmp_mat_i, 1, get_bias,Bg="AncSR2",input_freqs=T)
  df_tmp <- data.frame(step=s,bias=B_g,AA_var=random_starting_nodes)
  change_bias_df_sr2 <- rbind(change_bias_df_sr2,df_tmp)
}

stopCluster(cl)
```

``` r
# Number of substitutions for mean bias in outcomes to be within 0.05 units of equilibrium value
forget_lbias_sr2 <- change_bias_df_sr2 %>% group_by(AA_var) %>% filter(bias <= StatB_sr2+0.05) %>% 
  distinct(.,AA_var,.keep_all = TRUE) %>% with(mean(step))

print(paste("Number of substitutions for mean bias in outcomes to be within 0.05 units of equilibrium value = ",round(forget_lbias_sr2,2)))
```

    ## [1] "Number of substitutions for mean bias in outcomes to be within 0.05 units of equilibrium value =  13.71"

``` r
# Number of substitutions for bias in outcomes from EGKA to be within 0.05 units of equilibrium value
forget_lbias_EGKA_sr2 <- change_bias_df_sr2 %>% filter(AA_var==REF_GENOTYPE) %>% filter(bias <= StatB_sr2+0.05) %>% head(., 1) %>% pull(step)

print(paste("Number of substitutions for bias in outcomes from EGKA to be within 0.05 units of equilibrium value = ",forget_lbias_EGKA_sr2))
```

    ## [1] "Number of substitutions for bias in outcomes from EGKA to be within 0.05 units of equilibrium value =  5"

``` r
# change in outcome bias over time from ancestral genotype EGKA
egka_bias_sr2 <- change_bias_df_sr2 %>% filter(step<=50) %>% filter(AA_var=="EGKA")

# average change in outcome bias over time across genotypes
mean_bias_sr2 <- change_bias_df_sr2 %>% filter(step<=50) %>% group_by(step) %>% reframe(mean_B=mean(bias,na.rm=T))

mut_to_brlen <- function(x) x / 4 # branch length units
MAX_STEP = 150
change_bias_df_sr2 %>%  filter(step<=50) %>%
  ggplot(aes(x=step,y=bias)) + geom_point(color="gray60",size=0.5) + 
  xlab("Substitution step") + ylab("Phenotype Bias") +
  geom_hline(yintercept = GlobalB_sr2,col="red",linetype="dashed",linewidth=0.7) +
  geom_hline(yintercept = StatB_sr2,col="orange",linetype="dashed",linewidth=0.7) +
  geom_vline(xintercept = forget_lbias_sr2,col="gray80",linewidth=1.5) +
  geom_line(data=egka_bias_sr2,aes(x=step,y=bias),color="black",size=0.5) +
  geom_line(data=mean_bias_sr2,aes(x=step,y=mean_B),color="orange",size=1.5) +
  theme_classic() +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0,MAX_STEP,10),
                     sec.axis = dup_axis(trans = ~mut_to_brlen(.), name = "Branch length (subs / RH site)",breaks = seq(0, 12, by = 2))) +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        legend.position = "none")
```

    ## Warning: Removed 2 rows containing missing values (`geom_point()`).

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-41-1.png)

Consistent with our prediction, we see that the mean outcome bias acorss genotypes converges to the equilibrium bias more rapidly (from ~42 substitutions on AncSR1 to ~14 in AncSR2). This also explains why the probability of diversification at moderate timescales is higher: local bias favors conservations and since its effect is lost more rapidly, the probability of conservation also decreases at moderate timescales.

Second, because genotypes converge more rapidly to the equilibrium, the similarity of phenotypic outcomes at moderate timecsales across genotypes is higher in the AncSR2 map. We can compare the average similarity in outcomes between AncSR1 and AncSR2.

``` r
###################################
# Compute similarity of phenotypic outcomes - Pairwise correlations of outcome spectra between genotypes
###################################

cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

pairwise_corr_pdfvs_sr2 <- foreach(i = 1:dim(markov_chain_sr2)[1],.combine='rbind',.errorhandling="pass") %:%
  foreach(j = (i+1):dim(markov_chain_sr2)[1],.combine='rbind') %dopar% {
    g_i <- rownames(markov_chain_sr2)[i]
    g_j <- rownames(markov_chain_sr2)[j]
    c <- cor(markov_chain_sr2[i,],markov_chain_sr2[j,])^2
    data.frame(g_i = g_i, g_j = g_j, r2 = c)
  }

pairwise_corr_pdfvs_sr2 <- pairwise_corr_pdfvs_sr2 %>% filter(g_i != g_j) # remove same-genotype comparisons (if any)

stopCluster(cl)
```

``` r
# Similarity of phenotypic outcomes - comparison between maps
phen_sim_comparison <- rbind(pairwise_corr_pdfvs_sr1 %>% mutate(bg="AncSR1") %>% select(bg,r2),
            pairwise_corr_pdfvs_sr2 %>% mutate(bg="AncSR2") %>% select(bg,r2))

knitr::kable(phen_sim_comparison %>% group_by(bg) %>% reframe(avg_similarity_outcomes=mean(r2)))
```

| bg     |  avg\_similarity\_outcomes|
|:-------|--------------------------:|
| AncSR1 |                  0.3091636|
| AncSR2 |                  0.8253026|

``` r
phen_sim_comparison %>% ggplot(aes(x=bg,y=r2)) +
  geom_violin() + geom_pointrange(stat = "summary",size=0.4,col="red") + 
  theme_classic() + labs(x="Background",y=expression(paste("Similarity in evolutionary outcomes (r"^2,")"))) +
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,0.25)) +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust = 1,vjust = 1))
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-42-1.png)

``` r
wilcox.test(phen_sim_comparison$r2~phen_sim_comparison$bg)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  phen_sim_comparison$r2 by phen_sim_comparison$bg
    ## W = 2338256548, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0

Finally, because genotypes tend to evolve more similar outcomes across the network and their outcome spectra resemble more the equilibrium spectra, this results in SRE being the most likely phenotype to evolve as a whole and the most frequent novel phenotype from the majority of genotypes.

``` r
###################################
# Most likely phenotypic transition per genotype
###################################

as_tibble(markov_chain_sr2,rownames = NA) %>% 
    rownames_to_column(var="AA_var") %>% 
    inner_join(.,geno_pheno_sr2,by="AA_var") %>%
    pivot_longer(cols=2:15,names_to="Pheno",values_to="Prob") %>%
    group_by(AA_var) %>% 
  # remove conservation cases
  filter(RE != Pheno) %>% 
    filter(Pheno != "Promiscuous") %>%
    group_by(AA_var) %>% 
  # extract transitions with the highest probability 
  filter(Prob == max(Prob)) %>%
  ungroup() %>%
  reframe(RE = REs[[1]],
          count = table(factor(Pheno,levels=REs[[1]])),
          total = sum(count),
          prop = count / total,
          color = case_when(RE == "ERE (GT)"~"A",
                            RE == "SRE (AA)"~"B",
                            TRUE ~ "C")) %>%
  ggplot(aes(x=reorder(RE,count,decreasing = T),y=prop,fill=RE)) + 
  geom_bar(stat = "identity",color="black") +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25),labels = scales::percent) +
  theme_classic() +
  labs(x="Most likely novel DNA specificity\nper genotype",y="% starting genotypes") +
  scale_fill_manual(values = hex_RE_colors(1)) +
  theme(axis.title = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust = 1,vjust = 1)) + guides(fill="none")
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-43-1.png)

SRE also becomes the most likely phenotype to evolve from the ancestral RH genotype EGKA. Let's see the dymamics of the outcome spectrum over time from EGKA.

``` r
###################################
# Dynamics of evolutionary outcomes over time from EGKA 
###################################

mc_iter <- 25 # number of substitutions to track
SPEC_BINDING <- TRUE

# Run Markov chain for 25 substitutions and track the outcome spectra at every time point
pdfv_mc_multistep_sr2 <- simulate_markov_chain_multistep(REF_GENOTYPE,P_drift_sr2_ntwrk,mc_iter,"AncSR1",specific = SPEC_BINDING)
```

``` r
mc_iter <- 25 # number of substitutions to track
SPEC_BINDING <- TRUE
REMOVE_PROMISCUOUS <- TRUE
NORMALIZE <- TRUE

if(!SPEC_BINDING){
  # outcome spectra at time step 0
  s0 <- data.frame(RE=REs[[1]],Norm_F_prob=c(rep(0,2),1,rep(0,13)),model=0)
} else{
  s0 <- data.frame(RE=REs[[2]],Norm_F_prob=c(rep(0,2),1,rep(0,14)),model=0)
  if(REMOVE_PROMISCUOUS){
    if(NORMALIZE){
      pdfv_mc_multistep_sr2 <- lapply(pdfv_mc_multistep_sr2, remove_promiscuous,norm=T)
    }
    else{
      pdfv_mc_multistep_sr2 <- lapply(pdfv_mc_multistep_sr2, remove_promiscuous,norm=F)
    }
    s0 <- data.frame(RE=REs[[1]],Norm_F_prob=c(rep(0,2),1,rep(0,13)),model=0)
  }
}

## plot the dynamics of outcome spectra over time
df_sims_sr2 <- do.call(rbind,pdfv_mc_multistep_sr2) %>% rbind(.,s0)
ERE_SRE_ONLY <- F # whether to color lines only for ERE and SRE

if(ERE_SRE_ONLY){
  ere_sre <- df_sims_sr2 %>% filter(RE %in% c("ERE (GT)","SRE (AA)"))
  df_sims_sr2 %>%
    ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
    geom_line(linewidth=1.3) + 
    geom_line(data=ere_sre,linewidth=1.8,aes(color=RE)) + 
    scale_color_manual(values = hex_RE_colors(3)) + 
    theme_classic() +
    labs(x="Substitution step",y="Probability",title = paste("AncSR2:",REF_GENOTYPE),color="DNA specificity",fill="DNA specificity") +
    scale_x_continuous(breaks=seq(0,mc_iter,5),labels=seq(0,mc_iter,5)) +
    theme(axis.title = element_text(size=15),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13)) +
    guides(color = guide_legend(override.aes = list(size = 5)))
} else{
  df_sims_sr2 %>%
    ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
    geom_line(linewidth=1.3,aes(color=RE)) + 
    scale_color_manual(values = hex_RE_colors(1)) +
    theme_classic() +
    labs(x="Substitution step",y="Probability",title = paste("AncSR2:",REF_GENOTYPE),color="DNA specificity",fill="DNA specificity") +
    scale_x_continuous(breaks=seq(0,mc_iter,5),labels=seq(0,mc_iter,5)) +
    theme(axis.title = element_text(size=15),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13)) +
    guides(color = guide_legend(override.aes = list(size = 5)))
}
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-44-1.png)

``` r
knitr::kable(df_sims_sr2 %>% filter(model %in% c(3,6) & RE %in% c("ERE (GT)","SRE (AA)")) %>% dplyr::rename(Probability=Norm_F_prob,step=model))
```

| RE       |  Probability|  step|
|:---------|------------:|-----:|
| ERE (GT) |    0.2204272|     3|
| SRE (AA) |    0.2078410|     3|
| ERE (GT) |    0.1059626|     6|
| SRE (AA) |    0.3309883|     6|

After just 3 substitutions, SRE becomes as likely to evolve as ERE, and after 6 substitutions SRE is the most likely phenoypic outcome to evolve from the historical genotype.

Overall, we see that the background substitutions changed the RH GP map along the kSR lineage and made the evolution of SRE - the phenotype that historically evolved - the most likely outcome.

------------------------------------------------------------------------

## Evolution of protein-DNA complexes

So far we have only considered the GP map from the perspective of the protein. We also considered a scenario where the SR regulatory module evolves as a protein-DNA complex, where both interacting partners can mutate. In this scenario, each node of the GP map coresponds to a protein-DNA complex genotype (instead of a protein genotype), and two complexes are neighbors if they are one mutation appart at the protein OR DNA side (note that there are no "promiscuous" genotypes anymore, because each complex is a unique genotype).

Like before, let's first build the joint genotype networks and their transition probability matrices.

``` r
# Check whether the matrices have been created already
if(!file.exists(file.path("..","..","results","evolutionary_modeling","MutSel_matrices_complexes_final.RData"))) {
  # Load complete data from mutation effects model
  meanF_data <- readr::read_csv(file.path("..","..","results","classifying_functional_variants","meanF_data_fxnal.csv.gz"))
  
  # BUILD PHENOTYPE TABLE #
  # Add column of protein-DNA complexes
  meanF_data <- meanF_data %>% mutate(RE_mod = case_when(RE == "ERE (GT)" ~ "GT",
                                           RE == "SRE (AA)" ~ "AA",
                                           TRUE ~ RE),
                                      complex = paste(AA_var,RE_mod,sep=""))
  
  # For each *functional* complex (those with meanF >= AncSR1_ERE_ref or meanF >= AncSR2_SRE_ref):
  phenotypes_tbl <- meanF_data %>% filter(functional==T) %>% group_by(complex, bg) %>%
    reframe(AA_var = AA_var,
            meanF_bREs = avg_meanF, # the average meanF for the prot-DNA complex
            bound_REs = list(RE), # assign DNA elements bound
            RE = RE) %>% 
    ungroup() %>% group_by(AA_var,bg) %>%
    reframe(complex = complex,
            meanF_bREs = meanF_bREs,
            bound_REs = bound_REs,
            n_bound_REs = n(), # how many DNA elements can bind the RH
            specific = ifelse(n_bound_REs == 1, "YES","NO"), # functional binding to only one DNA element?
            specificity = ifelse(n_bound_REs == 1, RE, "Promiscuous")) # Determine the type of specificity
  
  # BUILD GENOTYPE NETWORKS #
  sr1_complexes <- phenotypes_tbl %>% filter(bg == "AncSR1") %>% pull(complex)
  sr2_complexes <- phenotypes_tbl %>% filter(bg == "AncSR2") %>% pull(complex)
  
  net_sr1_complex <- build_genotype_network(nodes=sr1_complexes, type=4, cores=N_CORES)
  net_sr2_complex <- build_genotype_network(nodes=sr2_complexes, type=4, cores=N_CORES)

  # BUILD TRANSITION MATRICES #
  Ne = 10e6 # population size (reasonable average for chordates e.g.,https://doi.org/10.7554/eLife.67509)
  N_CORES=detectCores()-1 # number of cores

  # Set-up parameters for the fitness function:
  # Fitness corresponds to the exponential growth rate = N*r
  MIN_ACTIVE <- phenotypes_tbl %>% with(min(meanF_bREs))
  STEP.PARAM = c(MIN_ACTIVE)

  # Build model scenarios for each DBD background:
  MODEL.PARAM_SR1_drift <- list("AncSR1","drift",Ne,STEP.PARAM,TRUE,"mean")
  MODEL.PARAM_SR2_drift <- list("AncSR2","drift",Ne,STEP.PARAM,TRUE,"mean")
  
  # Build transition matrices
  adj_mat_complex_sr1 <- build_mutation_matrix(sr1_complexes,type=4,N_CORES)
  M_drift_sr1 <- build_transition_matrix_v2(sr1_complexes,adj_mat_complex_sr1,phenotypes_tbl,MODEL.PARAM_SR1_drift,N_CORES,complex=TRUE)
  
  adj_mat_complex_sr2 <- build_mutation_matrix(sr2_complexes,type=4,N_CORES)
  M_drift_sr2 <- build_transition_matrix_v2(sr2_complexes,adj_mat_complex_sr2,phenotypes_tbl,MODEL.PARAM_SR2_drift,N_CORES,complex=TRUE)
  
} else {
  # load matrices if already created
  load(file.path("..","..","results","evolutionary_modeling","MutSel_matrices_complexes_final.RData"))
  phenotypes_tbl_complex <- phenotypes_tbl # save it for later
  M_drift_sr1_complex <- M_drift_sr1
  M_drift_sr2_complex <- M_drift_sr2
  rm(list = ls()[ls() %in% c("phenotypes_tbl","M_drift_sr1","M_drift_sr2")]) # clear up space
  N_CORES=detectCores()-1 # re-set cores if running locally
  load(file.path("..","..","results","evolutionary_modeling","GPmap_coevol_model_sr1.RData")) # load modeling results for AncSR1
  load(file.path("..","..","results","evolutionary_modeling","GPmap_coevol_model_sr2.RData")) # load modeling results for AncSR2
}
```

Let's check how many genotypes are part of the main network component in each background

``` r
knitr::kable(rbind(data.frame(Bg="AncSR1",total_fxnal_complexes = length(phenotypes_tbl_complex %>% filter(bg == "AncSR1") %>% pull(complex)), fxnal_complexes_in_net = length(extract_main_ntwrk(net_sr1_complex,nodes = T,tr_mat = NULL))),
                   data.frame(Bg="AncSR2",total_fxnal_complexes = length(phenotypes_tbl_complex %>% filter(bg == "AncSR2") %>% pull(complex)), fxnal_complexes_in_net = length(extract_main_ntwrk(net_sr2_complex,nodes = T,tr_mat = NULL)))))
```

| Bg     |  total\_fxnal\_complexes|  fxnal\_complexes\_in\_net|
|:-------|------------------------:|--------------------------:|
| AncSR1 |                      135|                        112|
| AncSR2 |                     4840|                       4829|

Again, we see that some genotypes are disconnected from the main component so we will now extract the corresponding square *P* matrices from each netwrok and compute their stationary distribution.

``` r
## Extract main components of each network and build the corresponding square P matrices
P_drift_sr1_ntwrk_complex <- extract_main_ntwrk(net_sr1_complex,M_drift_sr1_complex) # square P matrix
P_drift_sr1_ntwrk_complex_statdist <- stationary_dist(P_drift_sr1_ntwrk_complex) # compute stationary distribution

P_drift_sr2_ntwrk_complex <- extract_main_ntwrk(net_sr2_complex,M_drift_sr2_complex) # square P matrix
P_drift_sr2_ntwrk_complex_statdist <- stationary_dist(P_drift_sr2_ntwrk_complex) # compute stationary distribution
```

We will now compute the equilibrium spectrum and its bias for each DBD background. Note that the production spectrum (and the global bias) does not change because the number of protein-DNA complexes remains the same. In contrast, the stationary distribution does change because it directly depends on the structure of the genotype network.

``` r
## GLOBAL OBJECTS ##

REF_GENOTYPE_COMPLEX = "EGKAGT" # EGKA:GT (historical genotype complex) 

# Equilibrium outcome spectra
stationary_PDFV_sr1_complex <- get_PDFV_v2(P_drift_sr1_ntwrk_complex_statdist,type="simulated mc",Bg="AncSR1",model="Stationary",
                                           pheno_tbl = phenotypes_tbl_complex, specific = T,complex = T) # Specific binding
stationary_PDFV_sr2_complex <- get_PDFV_v2(P_drift_sr2_ntwrk_complex_statdist,type="simulated mc",Bg="AncSR2",model="Stationary",
                                           pheno_tbl = phenotypes_tbl_complex, specific = T,complex = T) # Specific binding

# Equilibrium bias
StatB_sr1_complex <- get_bias(stationary_PDFV_sr1_complex %>% remove_promiscuous(norm = T),Bg="AncSR1",input_freqs = TRUE)
StatB_sr2_complex <- get_bias(stationary_PDFV_sr2_complex %>% remove_promiscuous(norm = T),Bg="AncSR2",input_freqs = TRUE)

# Phenotype frequencies in main network
var.prop_AncSR1_net_spec_complex_df <- get_PDFV_v2(type="main network",Bg="AncSR1",model="Binding Network",graph = net_sr1_complex, pheno_tbl = phenotypes_tbl_complex, specific = T, complex = T) # Specific Binding
var.prop_AncSR2_net_spec_complex_df <- get_PDFV_v2(type="main network",Bg="AncSR2",model="Binding Network",graph = net_sr2_complex, pheno_tbl = phenotypes_tbl_complex, specific = T, complex = T) # Specific Binding

# Summary tables: Whether each phenotype is encoded in the GP map, whether it is encoded by specific variants, whether it is in the main network component, 
# and whether it is encoded by specific variants in the main network component.
# AncSR1 summary
summary_sr1_complex <- foreach(i = 1:length(REs[[1]]), .combine = 'rbind') %do% {
  encoded <- is.encoded(REs[[1]][i],graph = net_sr1_complex,pheno_df = phenotypes_tbl_complex,Bg = "AncSR1",complex=T)
  specific <- is.bound.specific(REs[[1]][i],graph = net_sr1_complex,pheno_df = phenotypes_tbl_complex,Bg = "AncSR1",complex = T)
  in_net <- is.in.ntwrk(REs[[1]][i], graph = net_sr1_complex,pheno_df = phenotypes_tbl_complex,Bg = "AncSR1",complex=T)
  in_net_spec <- is.in.ntwrk(REs[[1]][i],graph = net_sr1_complex,pheno_df = phenotypes_tbl_complex,Bg = "AncSR1",specific = T,complex = T)
  v <- c(encoded,specific,in_net,in_net_spec)
} %>% `rownames<-`(REs[[1]]) %>% `colnames<-`(c("Encoded","Specific","In_Net","In_Net_Spec"))

# AncSR2 summary
summary_sr2_complex <- foreach(i = 1:length(REs[[1]]), .combine = 'rbind') %do% {
  encoded <- is.encoded(REs[[1]][i],graph = net_sr2_complex,pheno_df = phenotypes_tbl_complex,Bg = "AncSR2",complex=T)
  specific <- is.bound.specific(REs[[1]][i],graph = net_sr2_complex,pheno_df = phenotypes_tbl_complex,Bg = "AncSR2",complex = T)
  in_net <- is.in.ntwrk(REs[[1]][i], graph = net_sr2_complex,pheno_df = phenotypes_tbl_complex,Bg = "AncSR2",complex=T)
  in_net_spec <- is.in.ntwrk(REs[[1]][i],graph = net_sr2_complex,pheno_df = phenotypes_tbl_complex,Bg = "AncSR2",specific = T,complex = T)
  v <- c(encoded,specific,in_net,in_net_spec)
} %>% `rownames<-`(REs[[1]]) %>% `colnames<-`(c("Encoded","Specific","In_Net","In_Net_Spec"))

############
# Other useful global objects

# Nodes in networks
nodes_in_ntwrk_sr1_complex <- extract_main_ntwrk(net_sr1_complex,nodes=T,tr_mat = NULL)
nodes_in_ntwrk_sr2_complex <- extract_main_ntwrk(net_sr2_complex,nodes=T,tr_mat = NULL)

# Nodes in network that are specific (complexes where the RH binds a single DNA element)
nodes_in_ntwrk_sr1_spec_complex <- data.frame(complex=nodes_in_ntwrk_sr1_complex,bg="AncSR1") %>% inner_join(.,phenotypes_tbl_complex,by=c("bg","complex")) %>% filter(specific=="YES") %>% pull(complex)
nodes_in_ntwrk_sr1_spec_complex <- nodes_in_ntwrk_sr1_spec_complex[nodes_in_ntwrk_sr1_spec_complex %in% nodes_in_ntwrk_sr1_complex]
nodes_in_ntwrk_sr2_spec_complex <- data.frame(complex=nodes_in_ntwrk_sr2_complex,bg="AncSR2") %>% inner_join(.,phenotypes_tbl_complex,by=c("bg","complex")) %>% filter(specific=="YES") %>% pull(complex)
nodes_in_ntwrk_sr2_spec_complex <- nodes_in_ntwrk_sr2_spec_complex[nodes_in_ntwrk_sr2_spec_complex %in% nodes_in_ntwrk_sr2_complex]

# Color palette for plotting
net_colors <- c("orange","#006ddb")
```

### Modeling phenotypic outcomes under coevolution

To understand what are the effects of coevolution, we will compare the outcome spectra in the coevolution scenario to those in the protein-evolution scenario.

Let's first compare the equilibrium probabilities of each phenotype to those in the protein-evolution networks.

``` r
# equilibrium spectra coevolution vs. protein evolution for each DBD background
df_long_term_sr1_2 <- inner_join(remove_promiscuous(stationary_PDFV_spec_sr1,norm = T),remove_promiscuous(stationary_PDFV_sr1_complex,norm = T),by="RE")
df_long_term_sr2_2 <- inner_join(remove_promiscuous(stationary_PDFV_spec_sr2,norm = T),remove_promiscuous(stationary_PDFV_sr2_complex,norm = T),by="RE")
df_long_term_comp <- rbind(df_long_term_sr1_2 %>% mutate(bg="AncSR1"),df_long_term_sr2_2 %>% mutate(bg="AncSR2"))

# correlation between equilibrium spectra
cor_df <- df_long_term_comp %>% group_by(bg) %>% reframe(cor = cor(Norm_F_prob.x,Norm_F_prob.y)^2) %>% 
  mutate(Norm_F_prob.x = 0.1, Norm_F_prob.y = 0.4, RE=NA, label = paste("r2 =",round(cor,3)))

# plot correlation of spectra
df_long_term_comp %>%
  ggplot(aes(x=Norm_F_prob.x,y=Norm_F_prob.y,fill=RE)) + 
  geom_point(shape=21,col="black",size=4) + scale_fill_manual(values = hex_RE_colors(1)) +
  scale_y_continuous(limits = c(0,0.45)) + scale_x_continuous(limits = c(0,0.45)) + 
  geom_abline(intercept = 0,slope = 1,linetype="dashed",col="black") +
  labs(x="Equilibrium outcome spectrum (protein ntwrk)",y="Equilibrium outcome spectrum (coevol ntwrk)") + 
  theme_classic() +
  theme(axis.text = element_text(size = 10),axis.title = element_text(size=11)) +
  facet_grid(~ bg) +
  geom_text(data=cor_df,aes(label=label))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-46-1.png)

We see that at infinitely long timescales, the equilibrium outcome spectra are almost identical. This makes sense because both the protein and protein-DNA networks have the same phenotypes in almost the same frequencies, thus, if evolutionary trajectories are unbounded they should eventually converge to almost the same stationary distribution.

``` r
# phenotype frequencies in protein and coevolution networks
pheno_freqs_in_net_sr1_comp <- inner_join(remove_promiscuous(var.prop_AncSR1_net_spec_df,norm = T),remove_promiscuous(var.prop_AncSR1_net_spec_complex_df,norm = T),by="RE")
pheno_freqs_in_net_sr2_comp <- inner_join(remove_promiscuous(var.prop_AncSR2_net_spec_df,norm = T),remove_promiscuous(var.prop_AncSR2_net_spec_complex_df,norm = T),by="RE")
df_pheno_freqs_in_net_comp <- rbind(pheno_freqs_in_net_sr1_comp %>% mutate(bg="AncSR1"),pheno_freqs_in_net_sr2_comp %>% mutate(bg="AncSR2"))

# compute correlation of phenotype frequencies in the networks 
knitr::kable(df_pheno_freqs_in_net_comp %>% group_by(bg) %>% reframe(r2_phenotype_freqs_in_netwrk = cor(Norm_F_prob.x,Norm_F_prob.y)^2))
```

| bg     |  r2\_phenotype\_freqs\_in\_netwrk|
|:-------|---------------------------------:|
| AncSR1 |                         0.9654795|
| AncSR2 |                         0.9999533|

Now, let's see what is the effect of coevolution at short-to-moderate timescales. We need to account for the fact that a mutation in the coevolution network can occur in the RH or the RE. To compare mutational trajectories between networks, we need to consider only mutations in the RH for the coevolution networks.

To do this, we will first compute the fraction of one-mutant neighbors that result in a change in the RH for every genotype in the coevolution networks. Using this information, we can then estimate the average number of amino acid mutations per protein-DNA step in the coevolution networks.

``` r
########################
# Estimate avg number of amino acid substitutions per amino acid step
########################

cl <- makeCluster(N_CORES,type = "FORK", outfile="")
registerDoParallel(cl)

# AncSR1
# fraction of RH neighbors per protein-DNA complex
fr_aa_subs_sr1 <- frac_aa_subs_per_g(net_sr1_complex_uwud,nodes_in_ntwrk_sr1_complex)

# run Markov chain for 50 prot-DNA steps --> at each step, compute the average number of amino acid changes across genotypes
aa_subs_mc_sr1_complex <- foreach(i = 1:50, .combine = 'rbind') %do% {
  step <- i
  aa_subs <- compute_aa_subs_from_protDNA_steps(nodes_in_ntwrk_sr1_complex,i,P_drift_sr1_ntwrk_complex,fr_aa_subs_sr1)
  data.frame(step=step,avg_aa_subs_at_step=aa_subs)
} %>% mutate(cum_aa_subs = cumsum(avg_aa_subs_at_step))

# AncSR2
# fraction of RH neighbors per protein-DNA complex
fr_aa_subs_sr2 <- frac_aa_subs_per_g(net_sr2_complex_uwud,nodes_in_ntwrk_sr2_complex)

# run Markov chain for 50 prot-DNA steps --> at each step, compute the average number of amino acid changes across genotypes
aa_subs_mc_sr2_complex <- foreach(i = 1:50, .combine = 'rbind') %dopar% {
  step <- i
  aa_subs <- compute_aa_subs_from_protDNA_steps(nodes_in_ntwrk_sr2_complex,i,P_drift_sr2_ntwrk_complex,fr_aa_subs_sr2)
  data.frame(step=step,avg_aa_subs_at_step=aa_subs)
} %>% mutate(cum_aa_subs = cumsum(avg_aa_subs_at_step))

stopCluster(cl)
```

``` r
# plot number of prot-DNA steps vs average number of amino acid steps (dashed line is x=y line)
rbind(aa_subs_mc_sr1_complex %>% mutate(bg="AncSR1"), aa_subs_mc_sr2_complex %>% mutate(bg="AncSR2")) %>%
  ggplot(aes(x=step,y=cum_aa_subs)) +
  geom_line(color="red") + 
  geom_abline(intercept = 0,slope=1,linetype=2,color="gray50") +
  labs(x="protein-DNA step",y="Average number of RH changes across genotypes") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),axis.title = element_text(size=11)) +
  facet_grid(~ bg)
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-48-1.png)

``` r
# Fit linear models to data
# AncSR1
mod1 <- lm(aa_subs_mc_sr1_complex$cum_aa_subs~aa_subs_mc_sr1_complex$step)
AVG_AA_SUBS_PER_STEP_SR1 <- mod1$coefficients[2] # avg number of amino acid changes per protein-DNA step

# AncSR2
mod2 <- lm(aa_subs_mc_sr2_complex$cum_aa_subs~aa_subs_mc_sr2_complex$step)
AVG_AA_SUBS_PER_STEP_SR2 <- mod2$coefficients[2] # avg number of amino acid changes per protein-DNA step

knitr::kable(data.frame(bg=c("AncSR1","AncSR2"),avg_N_AA_subs_per_protDNA_step=c(unname(AVG_AA_SUBS_PER_STEP_SR1),unname(AVG_AA_SUBS_PER_STEP_SR2))))
```

| bg     |  avg\_N\_AA\_subs\_per\_protDNA\_step|
|:-------|-------------------------------------:|
| AncSR1 |                             0.8284577|
| AncSR2 |                             0.8385281|

We see that we require ~1.2 protein-DNA steps to make one amino acid change on average for both backgrounds. We will use this information to correct evolutionary trajectories in order to compare evolutionary outcomes on both types of networks.

We will now compare the dynamics of the outcome bias over short-to-moderate timescales. Like before, we will track the change in the outcome bias across genotypes and compute how long does it take for the bias to decay to its equilibrium expectation.

``` r
###################################
# Change in outcome bias over time - coevolution networks
###################################

cl <- makeCluster(N_CORES,type = "FORK", outfile="")
registerDoParallel(cl)

# AncSR1 - protein-DNA ntwrk
MAX_STEP_B = round(200/AVG_AA_SUBS_PER_STEP_SR1) # number of protein-DNA steps to make 200 RH steps
change_bias_df_sr1_complex <- data.frame(step=NULL,bias=NULL,g=NULL)

for(s in c(1,seq(5,MAX_STEP_B,5))){
  # Compute PDFVs for each genotype per step
  tmp_mat_i <- foreach(i = 1:length(nodes_in_ntwrk_sr1_complex), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr1_complex[i],tr_mat = P_drift_sr1_ntwrk_complex,n_steps = s)
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR1",graph=net_sr1_complex_uwud,model = nodes_in_ntwrk_sr1_complex[i],specific = T,complex = T,pheno_tbl = phenotypes_tbl_complex) %>%
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[2]]))]
  }
  rownames(tmp_mat_i) <- nodes_in_ntwrk_sr1_complex
  tmp_mat_i <- tmp_mat_i[, which(colSums(tmp_mat_i) != 0)]
  
  # compute entropy per genotype per step (without promiscuous)
  tmp_mat_i <- t(apply(tmp_mat_i[,1:ncol(tmp_mat_i)-1],1,function(x) x/sum(x))) # renormalize probs
  B_g <- apply(tmp_mat_i, 1, get_bias,Bg="AncSR1",input_freqs=T)
  df_tmp <- data.frame(step=s,bias=B_g,g=nodes_in_ntwrk_sr1_complex)
  change_bias_df_sr1_complex <- rbind(change_bias_df_sr1_complex,df_tmp)
}


# AncSR2 - protein-DNA ntwrk
MAX_STEP_C = round(100/AVG_AA_SUBS_PER_STEP_SR2) # number of protein-DNA steps to make 100 RH steps
n_sample = 150 # sample 150 genotypes to run the analysis
set.seed(2468) # Make analysis replicable
random_starting_nodes <- c(sample(nodes_in_ntwrk_sr2_complex,size = n_sample,replace = F), REF_GENOTYPE_COMPLEX)
change_bias_df_sr2_complex <- data.frame(step=NULL,bias=NULL,g=NULL)

for(s in c(1,seq(5,MAX_STEP_C,5))){
  # Compute PDFVs for each genotype per step
  tmp_mat_i <- foreach(i = 1:length(random_starting_nodes), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = random_starting_nodes[i],tr_mat = P_drift_sr2_ntwrk_complex,n_steps = s)
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR2",graph=net_sr2_complex_uwud,model = random_starting_nodes[i],specific = T,complex = T,pheno_tbl = phenotypes_tbl_complex) %>%
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[2]]))]
  }
  rownames(tmp_mat_i) <- random_starting_nodes
  tmp_mat_i <- tmp_mat_i[, which(colSums(tmp_mat_i) != 0)]
  
  # compute entropy per genotype per step (without promiscuous)
  tmp_mat_i <- t(apply(tmp_mat_i[,1:ncol(tmp_mat_i)-1],1,function(x) x/sum(x))) # renormalize probs
  B_g <- apply(tmp_mat_i, 1, get_bias,Bg="AncSR2",input_freqs=T)
  df_tmp <- data.frame(step=s,bias=B_g,g=random_starting_nodes)
  change_bias_df_sr2_complex <- rbind(change_bias_df_sr2_complex,df_tmp)
}

stopCluster(cl)
```

``` r
# Number of substitutions for mean bias in outcomes to be within 0.05 units of equilibrium value
forget_lbias_sr1_complex <- change_bias_df_sr1_complex %>% group_by(g) %>% filter(bias <= StatB_sr1_complex+0.05) %>% 
  distinct(.,g,.keep_all = TRUE) %>% with(mean(step))

forget_lbias_sr2_complex <- change_bias_df_sr2_complex %>% group_by(g) %>% filter(bias <= StatB_sr2_complex+0.05) %>% 
  distinct(.,g,.keep_all = TRUE) %>% with(mean(step))

# complare results between networks and backgrounds
# compute avg number of subtitutions for local bias to decay to equilibrium
forget_lbias_mean_df <- data.frame(bg=rep(c("AncSR1","AncSR2"),2),ntwrk_type=rep(c("protein","protein-DNA"),each=2),
                              avg_subs_to_equilibrium = c(forget_lbias_sr1,forget_lbias_sr2,
                                                          forget_lbias_sr1_complex * AVG_AA_SUBS_PER_STEP_SR1,   # correct to count only RH changes
                                                          forget_lbias_sr2_complex * AVG_AA_SUBS_PER_STEP_SR2))  # correct to count only RH changes
knitr::kable(forget_lbias_mean_df)
```

| bg     | ntwrk\_type |  avg\_subs\_to\_equilibrium|
|:-------|:------------|---------------------------:|
| AncSR1 | protein     |                    41.90722|
| AncSR2 | protein     |                    13.71287|
| AncSR1 | protein-DNA |                    71.60242|
| AncSR2 | protein-DNA |                    30.32029|

``` r
# average change in outcome bias over time across genotypes
mean_bias_sr1 <- change_bias_df_sr1 %>% group_by(step) %>% reframe(mean_B=mean(bias,na.rm=T))
mean_bias_sr2 <- change_bias_df_sr2 %>% group_by(step) %>% reframe(mean_B=mean(bias,na.rm=T))
mean_bias_sr1_complex <- change_bias_df_sr1_complex %>% group_by(step) %>% reframe(mean_B=mean(bias,na.rm=T)) %>% mutate(step= step * AVG_AA_SUBS_PER_STEP_SR1)
mean_bias_sr2_complex <- change_bias_df_sr2_complex %>% group_by(step) %>% reframe(mean_B=mean(bias,na.rm=T)) %>% mutate(step= step * AVG_AA_SUBS_PER_STEP_SR2)

# references: Global bias and equilibrium bias
stat_bias_df <- data.frame(b=c(StatB_sr1,StatB_sr2,StatB_sr1_complex,StatB_sr2_complex),ntwrk_type=c("protein","protein","protein-DNA","protein-DNA"),
                           bg=c("AncSR1","AncSR2","AncSR1","AncSR2"))

global_bias_df <- data.frame(b=c(GlobalB_sr1,GlobalB_sr2,GlobalB_sr1,GlobalB_sr2),ntwrk_type=c("protein","protein","protein-DNA","protein-DNA"),
                           bg=c("AncSR1","AncSR2","AncSR1","AncSR2"))

# plot
rbind(mean_bias_sr1 %>% mutate(ntwrk_type="protein",bg="AncSR1"),
      mean_bias_sr2 %>% mutate(ntwrk_type="protein",bg="AncSR2"),
      mean_bias_sr1_complex %>% mutate(ntwrk_type="protein-DNA",bg="AncSR1"),
      mean_bias_sr2_complex %>% mutate(ntwrk_type="protein-DNA",bg="AncSR2")) %>%
  filter(step<=80) %>%
  ggplot(aes(x=step,y=mean_B,group=ntwrk_type)) +
  xlab("Substitution step") + ylab("Outcome Bias") +
  geom_hline(data = stat_bias_df, aes(yintercept = b,col=ntwrk_type),linetype="dashed",linewidth=0.7) +
  geom_hline(data = global_bias_df, aes(yintercept = b),col="red",linetype="dashed",linewidth=0.7) +
  geom_vline(data = forget_lbias_mean_df, aes(xintercept = avg_subs_to_equilibrium,col=ntwrk_type),linewidth=1.2,linetype=3) +
  geom_line(aes(col=ntwrk_type),size=1.5) +
  scale_color_manual(values = net_colors) +
  facet_wrap(nrow = 2,vars(bg)) +
  theme_classic() +
  scale_y_continuous(limits = c(0,1))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-50-1.png)

These results show that coevolution increases the time required for the average local bias to decay to equilibrium. This delay implies that the effect of local bias lasts for an extended period of time and as a consequence the probability of conservation should be higher over moderate timescales. To check this we will first run a markov chain on the coevolution networks.

``` r
###################################
# Compute phenotypic probabilities from every genotype
###################################

N_MUTATIONS_A <- round(8 / AVG_AA_SUBS_PER_STEP_SR1) # Number of protein-DNA stepts to produce an average of 8 AA substitutions
N_MUTATIONS_B <- round(8 / AVG_AA_SUBS_PER_STEP_SR2) # Number of protein-DNA stepts to produce an average of 8 AA substitutions
SPEC_BINDING <- TRUE
NORM_PROB <- TRUE # whether to remove promiscuous end-points and renormalize probanilities (only if SPEC_BINDING = TRUE)

cl <- makeCluster(N_CORES,type = "FORK")
registerDoParallel(cl)

# AncSR1 protein-DNA ntwrk
if(SPEC_BINDING){
 markov_chain_sr1_complex <- foreach(i = 1:length(nodes_in_ntwrk_sr1_complex), .combine = 'rbind') %dopar% {
   tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr1_complex[i],tr_mat = P_drift_sr1_ntwrk_complex,n_steps = N_MUTATIONS_A)
   tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR1",graph=net_sr1_complex_uwud,model = nodes_in_ntwrk_sr1_complex[i],
                          specific = T,complex = T,pheno_tbl = phenotypes_tbl_complex) %>%
     acast(model~RE,value.var="Norm_F_prob")
   tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
  }
  rownames(markov_chain_sr1_complex) <- nodes_in_ntwrk_sr1_complex
  markov_chain_sr1_complex <- markov_chain_sr1_complex[, which(colSums(markov_chain_sr1_complex) != 0)]
  
  if(NORM_PROB){
    markov_chain_sr1_complex <- markov_chain_sr1_complex[,1:ncol(markov_chain_sr1_complex)-1]
    markov_chain_sr1_complex <- t(apply(markov_chain_sr1_complex, 1, function(x) x/sum(x)))
  }
} else{
 markov_chain_sr1_complex <- foreach(i = 1:length(nodes_in_ntwrk_sr1_complex), .combine = 'rbind') %dopar% {
   tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr1_complex[i],tr_mat = P_drift_sr1_ntwrk_complex,n_steps = N_MUTATIONS)
   tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR1",graph=net_sr1_complex_uwud,model = nodes_in_ntwrk_sr1_complex[i],
                          specific = F,complex = T,pheno_tbl = phenotypes_tbl_complex) %>%
     acast(model~RE,value.var="Norm_F_prob")
   tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
  }
  rownames(markov_chain_sr1_complex) <- nodes_in_ntwrk_sr1_complex
  markov_chain_sr1_complex <- markov_chain_sr1_complex[, which(colSums(markov_chain_sr1_complex) != 0)]
  
}

# AncSR2 protein-DNA ntwrk
if(SPEC_BINDING){
  markov_chain_sr2_complex <- foreach(i = 1:length(nodes_in_ntwrk_sr2_complex), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr2_complex[i],tr_mat = P_drift_sr2_ntwrk_complex,n_steps = N_MUTATIONS_B)
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR2",graph=net_sr2_complex_uwud,model = nodes_in_ntwrk_sr2_complex[i],
                            specific = T,complex = T,pheno_tbl = phenotypes_tbl_complex) %>%
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
  }
  rownames(markov_chain_sr2_complex) <- nodes_in_ntwrk_sr2_complex
  markov_chain_sr2_complex <- markov_chain_sr2_complex[, which(colSums(markov_chain_sr2_complex) != 0)]
  
  if(NORM_PROB){
    markov_chain_sr2_complex <- markov_chain_sr2_complex[,1:ncol(markov_chain_sr2_complex)-1]
    markov_chain_sr2_complex <- t(apply(markov_chain_sr2_complex, 1, function(x) x/sum(x)))
  }
} else{
  markov_chain_sr2_complex <- foreach(i = 1:length(nodes_in_ntwrk_sr2_complex), .combine = 'rbind') %dopar% {
    tmp_mc <- simulate_markov_chain(states_0 = nodes_in_ntwrk_sr2_complex[i],tr_mat = P_drift_sr2_ntwrk_complex,n_steps = N_MUTATIONS)
    tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR2",graph=net_sr2_complex_uwud,model = nodes_in_ntwrk_sr2_complex[i],
                            specific = F,complex = T,pheno_tbl = phenotypes_tbl_complex) %>%
      acast(model~RE,value.var="Norm_F_prob")
    tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
  }
  rownames(markov_chain_sr2_complex) <- nodes_in_ntwrk_sr2_complex
  markov_chain_sr2_complex <- markov_chain_sr2_complex[, which(colSums(markov_chain_sr2_complex) != 0)]
  
}

stopCluster(cl)
```

Now we'll compute the average probability of conservation across genotypes:

``` r
###################################
# Probability of consrvation and transition per genotype
###################################

# AncSR1 protein-DNA ntwrk
# associate genotypes to phenotypes
geno_pheno_sr1_complex <- phenotypes_tbl_complex %>% filter(bg=="AncSR1") %>% select(complex,specificity) %>% dplyr::rename("RE" = "specificity")

# Probability of conservation and transition per *specific* genotype 
pheno_outcomes_sr1_complex_df <- as_tibble(markov_chain_sr1_complex,rownames = NA) %>% 
    rownames_to_column(var="complex") %>% 
    inner_join(.,geno_pheno_sr1_complex,by="complex") %>%
    pivot_longer(cols=2:6,names_to="Pheno",values_to="Prob") %>% 
    filter(RE != "Promiscuous") %>%
    group_by(complex) %>% filter(RE == Pheno) %>%
    reframe(p_cons = Prob,
            p_tr = 1-p_cons) %>%
    pivot_longer(cols=2:3,names_to="outcome",values_to="prob")


# AncSR2 protein-DNA ntwrk
# associate genotypes to phenotypes
geno_pheno_sr2_complex <- phenotypes_tbl_complex %>% filter(bg=="AncSR2") %>% select(complex,specificity) %>% dplyr::rename("RE" = "specificity")

# Probability of conservation and transition per *specific* genotype 
pheno_outcomes_sr2_complex_df <- as_tibble(markov_chain_sr2_complex,rownames = NA) %>% 
    rownames_to_column(var="complex") %>% 
    inner_join(.,geno_pheno_sr2_complex,by="complex") %>%
    pivot_longer(cols=2:15,names_to="Pheno",values_to="Prob") %>% 
    filter(RE != "Promiscuous") %>%
    group_by(complex) %>% filter(RE == Pheno) %>%
    reframe(p_cons = Prob,
            p_tr = 1-p_cons) %>%
    pivot_longer(cols=2:3,names_to="outcome",values_to="prob")

# compare average probability of conservation between netwroks for each GP map
knitr::kable(rbind(pheno_outcomes_sr1_df %>% mutate(bg="AncSR1",ntwrk_type="protein") %>% filter(outcome=="p_cons") %>% 
                     group_by(bg,ntwrk_type) %>% reframe(avg_P_conservation = mean(prob)),
                   pheno_outcomes_sr2_df %>% mutate(bg="AncSR2",ntwrk_type="protein") %>% filter(outcome=="p_cons") %>% 
                     group_by(bg,ntwrk_type) %>% reframe(avg_P_conservation = mean(prob)),
                   pheno_outcomes_sr1_complex_df %>% mutate(bg="AncSR1",ntwrk_type="protein-DNA") %>% filter(outcome=="p_cons") %>% 
                     group_by(bg,ntwrk_type) %>% reframe(avg_P_conservation = mean(prob)),
                   pheno_outcomes_sr2_complex_df %>% mutate(bg="AncSR2",ntwrk_type="protein-DNA") %>% filter(outcome=="p_cons") %>% 
                     group_by(bg,ntwrk_type) %>% reframe(avg_P_conservation = mean(prob))))
```

| bg     | ntwrk\_type |  avg\_P\_conservation|
|:-------|:------------|---------------------:|
| AncSR1 | protein     |             0.6101274|
| AncSR2 | protein     |             0.4701205|
| AncSR1 | protein-DNA |             0.8307546|
| AncSR2 | protein-DNA |             0.8030094|

As expected, the average probability of conservation is higher on the coevolution networks. We can also see this effect by investigating the probability of conservation per phenotype.

``` r
########################
# Probability of conservation of phenotypes
########################

cl <- makeCluster(N_CORES,type = "FORK",outfile="")
registerDoParallel(cl)

# AncSR1 - protein-DNA
REs_spec_sr1_complex <- data.frame(summary_sr1_complex) %>% rownames_to_column(var="RE") %>% filter(Encoded == T & In_Net_Spec == T) %>% pull(RE)

pt_sr1_complex <- foreach(i = 1:length(REs_spec_sr1_complex), .combine = 'rbind') %dopar% {
  phenotype_vars <- phenotypes_tbl_complex %>% filter(bg == "AncSR1" & specificity == REs_spec_sr1_complex[i]) %>% pull(complex) # complexes encoding specificity i
  tmp_mc <- simulate_markov_chain(states_0 = phenotype_vars,tr_mat = P_drift_sr1_ntwrk_complex,n_steps = N_MUTATIONS_A) # run markov chain from genotype set i
  tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR1",graph=net_sr1_complex_uwud,model = REs_spec_sr1_complex[i],specific = T,complex = T,pheno_tbl = phenotypes_tbl_complex) %>%
    acast(model~RE,value.var="Norm_F_prob") # compute outcome spectrum from genotype set i
  tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
}
rownames(pt_sr1_complex) <- REs_spec_sr1_complex
pt_sr1_complex <- pt_sr1_complex[, which(colSums(pt_sr1_complex) != 0)]
pt_sr1_complex <- pt_sr1_complex[,1: ncol(pt_sr1_complex)-1]
pt_sr1_complex <- t(apply(pt_sr1_complex, 1, function(x) x/sum(x)))

# AncSR2 - protein-DNA
REs_spec_sr2_complex <- data.frame(summary_sr2_complex) %>% rownames_to_column(var="RE") %>% filter(Encoded == T & In_Net_Spec == T) %>% pull(RE)

pt_sr2_complex <- foreach(i = 1:length(REs_spec_sr2_complex), .combine = 'rbind') %dopar% {
  phenotype_vars <- phenotypes_tbl_complex %>% filter(bg == "AncSR2" & specificity == REs_spec_sr2_complex[i]) %>% pull(complex) # complexes encoding specificity i
  tmp_mc <- simulate_markov_chain(states_0 = phenotype_vars,tr_mat = P_drift_sr2_ntwrk_complex,n_steps = N_MUTATIONS_B) # run markov chain from genotype set i
  tmp_pdfv <- get_PDFV_v2(tmp_mc,type="simulated mc", Bg="AncSR2",graph=net_sr2_complex_uwud,model = REs_spec_sr2_complex[i],specific = T,complex = T,pheno_tbl = phenotypes_tbl_complex) %>%
    acast(model~RE,value.var="Norm_F_prob") # compute outcome spectrum from genotype set i
  tmp_pdfv <- tmp_pdfv[,order(match(colnames(tmp_pdfv),REs[[1]]))]
}
rownames(pt_sr2_complex) <- REs_spec_sr2_complex
pt_sr2_complex <- pt_sr2_complex[, which(colSums(pt_sr2_complex) != 0)]
pt_sr2_complex <- pt_sr2_complex[,1: ncol(pt_sr2_complex)-1]
pt_sr2_complex <- t(apply(pt_sr2_complex, 1, function(x) x/sum(x)))

stopCluster(cl)
```

``` r
# Compute probabilities of conservation per phenotype
extract_prob_conservation <- function(df,RE_i){
  x <- df[rownames(df) == RE_i, colnames(df) == RE_i]
  return(x)
}

# Summarise per phenotype
# AncSR1 protein ntwrk
pt_sr1_df <- as_tibble(pt_sr1,rownames=NA) %>% rownames_to_column(var="RE") %>%
  mutate(total_prob = rowSums(across(`SRE (AA)`:`TA`))) %>%
  mutate(prob_conservation = map_dbl(.x=RE, .f=extract_prob_conservation,df=pt_sr1),
         prob_transition = total_prob - prob_conservation) %>% select(RE,prob_conservation,prob_transition) %>%
  arrange(RE = desc(prob_conservation)) %>% mutate(RE = factor(RE, levels = RE)) %>%
  pivot_longer(cols = 2:3,names_to="type",values_to="prob") %>%
  mutate(type = ifelse(type == "prob_conservation","Conservation","Transition"),
         bg="AncSR1",net="protein")

# AncSR2 protein ntwrk
pt_sr2_df <- as_tibble(pt_sr2,rownames=NA) %>% rownames_to_column(var="RE") %>%
  mutate(total_prob = rowSums(across(`SRE (AA)`:`TT`))) %>%
  mutate(prob_conservation = map_dbl(.x=RE, .f=extract_prob_conservation,df=pt_sr2),
         prob_transition = total_prob - prob_conservation) %>% select(RE,prob_conservation,prob_transition) %>%
  arrange(RE = desc(prob_conservation)) %>% mutate(RE = factor(RE, levels = RE)) %>%
  pivot_longer(cols = 2:3,names_to="type",values_to="prob") %>%
  mutate(type = ifelse(type == "prob_conservation","Conservation","Transition"),
         bg="AncSR2",net="protein")

# AncSR1 protein-DNA ntwrk
pt_sr1_complex_df <- as_tibble(pt_sr1_complex,rownames=NA) %>% rownames_to_column(var="RE") %>%
  mutate(total_prob = rowSums(across(`SRE (AA)`:`TA`))) %>%
  mutate(prob_conservation = map_dbl(.x=RE, .f=extract_prob_conservation,df=pt_sr1_complex),
         prob_transition = total_prob - prob_conservation) %>% select(RE,prob_conservation,prob_transition) %>%
  arrange(RE = desc(prob_conservation)) %>% mutate(RE = factor(RE, levels = RE)) %>%
  pivot_longer(cols = 2:3,names_to="type",values_to="prob") %>%
  mutate(type = ifelse(type == "prob_conservation","Conservation","Transition"),
         bg="AncSR1",net="protein-DNA")

# AncSR2 protein-DNA ntwrk
pt_sr2_complex_df <- as_tibble(pt_sr2_complex,rownames=NA) %>% rownames_to_column(var="RE") %>%
  mutate(total_prob = rowSums(across(`SRE (AA)`:`TT`))) %>%
  mutate(prob_conservation = map_dbl(.x=RE, .f=extract_prob_conservation,df=pt_sr2_complex),
         prob_transition = total_prob - prob_conservation) %>% select(RE,prob_conservation,prob_transition) %>%
  arrange(RE = desc(prob_conservation)) %>% mutate(RE = factor(RE, levels = RE)) %>%
  pivot_longer(cols = 2:3,names_to="type",values_to="prob") %>%
  mutate(type = ifelse(type == "prob_conservation","Conservation","Transition"),
         bg="AncSR2",net="protein-DNA")

# plot
rbind(inner_join(pt_sr1_df,pt_sr1_complex_df,by=c("RE","type")) %>% filter(type=="Conservation") %>% mutate(Bg = "AncSR1"),
             inner_join(pt_sr2_df,pt_sr2_complex_df,by=c("RE","type")) %>% filter(type=="Conservation") %>% mutate(Bg = "AncSR2")) %>%
  ggplot(aes(x=prob.x,y=prob.y)) +
  geom_point(aes(fill=RE,shape=Bg), size=4) + 
  geom_abline(slope=1,intercept = 0,linetype="dashed") + 
  theme_classic() +
  scale_shape_manual(values = c(21,24)) +
  scale_fill_manual(values = hex_RE_colors(1)) +
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) +
  labs(x=expression("Pr(conservation)"[protein]),y=expression("Pr(conservation)"[protein-DNA])) +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust = 1,vjust = 1)) +
  guides(fill = guide_legend(override.aes = list(shape=21,size=3)))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-54-1.png)

Virtually all phenotypes are more likely to remain conserved in the coevolution networks after 8 amino acid substitutions. Overall, we see that coevolution extends the influence of local bias on evolutionary outcomes and as a consequence it makes phenotypic diversification less likely to happen.

Finally, let's see what happens to the outcomes of evolution that result when trajectories begin from the ancestral genotype complex EGKA:GT.

``` r
###################################
# Phenotypic dynamics over time
###################################

# AncSR1 protein-DNA ntwrk
AA_subs <- 25
mc_iter <- round(AA_subs/AVG_AA_SUBS_PER_STEP_SR1) # Number of protein-DNA steps that results in a given number of aa subs (AA_subs)
SPEC_BINDING <- TRUE

pdfv_mc_multistep_sr1_complex <- simulate_markov_chain_multistep(REF_GENOTYPE_COMPLEX,P_drift_sr1_ntwrk_complex,mc_iter,"AncSR1",specific = SPEC_BINDING,complex = T,pheno_tbl = phenotypes_tbl_complex)
pdfv_mc_multistep_sr1_complex <- lapply(pdfv_mc_multistep_sr1_complex, function(x) x %>% mutate(model= model * AVG_AA_SUBS_PER_STEP_SR1)) # correct amino acid steps

# AncSR2 protein-DNA ntwrk
AA_subs <- 25
mc_iter <- round(AA_subs/AVG_AA_SUBS_PER_STEP_SR2) # Number of protein-DNA steps that results in a given number of aa subs (AA_subs)

pdfv_mc_multistep_sr2_complex <- simulate_markov_chain_multistep(REF_GENOTYPE_COMPLEX,P_drift_sr2_ntwrk_complex,mc_iter,"AncSR2",specific = SPEC_BINDING,complex = T,pheno_tbl = phenotypes_tbl_complex)
pdfv_mc_multistep_sr2_complex <- lapply(pdfv_mc_multistep_sr2_complex, function(x) x %>% mutate(model= model * AVG_AA_SUBS_PER_STEP_SR2)) # correct amino acid steps
```

``` r
SPEC_BINDING <- TRUE
REMOVE_PROMISCUOUS <- TRUE
NORMALIZE <- TRUE

if(!SPEC_BINDING){
  # outcome spectra at time step 0
  s0 <- data.frame(RE=REs[[1]],Norm_F_prob=c(rep(0,2),1,rep(0,13)),model=0)
} else{
  s0 <- data.frame(RE=REs[[2]],Norm_F_prob=c(rep(0,2),1,rep(0,14)),model=0)
  if(REMOVE_PROMISCUOUS){
    if(NORMALIZE){
      pdfv_mc_multistep_sr1_complex <- lapply(pdfv_mc_multistep_sr1_complex, remove_promiscuous,norm=T)
      pdfv_mc_multistep_sr2_complex <- lapply(pdfv_mc_multistep_sr2_complex, remove_promiscuous,norm=T)
    }
    else{
      pdfv_mc_multistep_sr1_complex <- lapply(pdfv_mc_multistep_sr1_complex, remove_promiscuous,norm=F)
      pdfv_mc_multistep_sr2_complex <- lapply(pdfv_mc_multistep_sr2_complex, remove_promiscuous,norm=F)
    }
    s0 <- data.frame(RE=REs[[1]],Norm_F_prob=c(rep(0,2),1,rep(0,13)),model=0)
  }
}

# plot
df_sims_sr1_complex <- do.call(rbind,pdfv_mc_multistep_sr1_complex) %>% rbind(.,s0)
df_sims_sr2_complex <- do.call(rbind,pdfv_mc_multistep_sr2_complex) %>% rbind(.,s0)

mc_iter_1 <- round(25/AVG_AA_SUBS_PER_STEP_SR1)
mc_iter_2 <- round(25/AVG_AA_SUBS_PER_STEP_SR2)

ERE_SRE_ONLY <- F # whether to color only ERE and SRE

if(ERE_SRE_ONLY){
  ere_sre_sr1 <- df_sims_sr1_complex %>% filter(RE %in% c("ERE (GT)","SRE (AA)"))
  ere_sre_sr2 <- df_sims_sr2_complex %>% filter(RE %in% c("ERE (GT)","SRE (AA)"))
  p1 <- df_sims_sr1_complex %>%
    ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
    geom_line(linewidth=1.3) + 
    geom_line(data=ere_sre,linewidth=1.8,aes(color=RE)) + 
    scale_color_manual(values = hex_RE_colors(3)) + 
    theme_classic() +
    labs(x="Substitution step",y="Probability",title = paste("AncSR1:",REF_GENOTYPE_COMPLEX),color="DNA specificity",fill="DNA specificity") +
    scale_x_continuous(breaks=seq(0,mc_iter_1,5),labels=seq(0,mc_iter_1,5)) +
    theme(axis.title = element_text(size=15),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13),
          legend.position = "none")
  
  p2 <- df_sims_sr1_complex %>%
    ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
    geom_line(linewidth=1.3) + 
    geom_line(data=ere_sre,linewidth=1.8,aes(color=RE)) + 
    scale_color_manual(values = hex_RE_colors(3)) + 
    theme_classic() +
    labs(x="Substitution step",y="Probability",title = paste("AncSR1:",REF_GENOTYPE_COMPLEX),color="DNA specificity",fill="DNA specificity") +
    scale_x_continuous(breaks=seq(0,mc_iter_2,5),labels=seq(0,mc_iter_2,5)) +
    theme(axis.title = element_text(size=15),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13)) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
} else{
  p1 <- df_sims_sr1_complex %>%
    ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
    geom_line(linewidth=1.3,aes(color=RE)) + 
    scale_color_manual(values = hex_RE_colors(1)) +
    theme_classic() +
    labs(x="Substitution step",y="Probability",title = paste("AncSR1:",REF_GENOTYPE_COMPLEX),color="DNA specificity",fill="DNA specificity") +
    scale_x_continuous(breaks=seq(0,mc_iter_1,5),labels=seq(0,mc_iter_1,5)) +
    theme(axis.title = element_text(size=15),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13),
          legend.position = "none")
  
  p2 <- df_sims_sr2_complex %>%
    ggplot(aes(x=model,y=Norm_F_prob,color=RE)) +
    geom_line(linewidth=1.3,aes(color=RE)) + 
    scale_color_manual(values = hex_RE_colors(1)) +
    theme_classic() +
    labs(x="Substitution step",y="Probability",title = paste("AncSR1:",REF_GENOTYPE_COMPLEX),color="DNA specificity",fill="DNA specificity") +
    scale_x_continuous(breaks=seq(0,mc_iter_2,5),labels=seq(0,mc_iter_2,5)) +
    theme(axis.title = element_text(size=15),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13)) +
    guides(color = guide_legend(override.aes = list(size = 5)))
}

p1 + p2
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-55-1.png)

``` r
# compare probability of ERE conservation after 3 and 10 amino acid steps between networks
knitr::kable(rbind(df_sims_sr1 %>% mutate(bg="AncSR1",ntwrk_type="protein"),
      df_sims_sr2 %>% mutate(bg="AncSR2",ntwrk_type="protein"),
      df_sims_sr1_complex %>% mutate(bg="AncSR1",ntwrk_type="protein-DNA",model=round(model)),
      df_sims_sr2_complex %>% mutate(bg="AncSR2",ntwrk_type="protein-DNA",model=round(model))) %>%
  filter(model %in% c(3,10) & RE %in% c("ERE (GT)")) %>% dplyr::rename(Probability=Norm_F_prob,step=model))
```

| RE       |  Probability|  step| bg     | ntwrk\_type |
|:---------|------------:|-----:|:-------|:------------|
| ERE (GT) |    1.0000000|     3| AncSR1 | protein     |
| ERE (GT) |    0.8424599|    10| AncSR1 | protein     |
| ERE (GT) |    0.2204272|     3| AncSR2 | protein     |
| ERE (GT) |    0.0503278|    10| AncSR2 | protein     |
| ERE (GT) |    1.0000000|     3| AncSR1 | protein-DNA |
| ERE (GT) |    0.9794327|    10| AncSR1 | protein-DNA |
| ERE (GT) |    0.9752129|     3| AncSR2 | protein-DNA |
| ERE (GT) |    0.8921665|     3| AncSR2 | protein-DNA |
| ERE (GT) |    0.2469921|    10| AncSR2 | protein-DNA |

Consistent with the previous results, we can see that ERE conservation remains a more likely outcome for a longer period of time in both GP maps when we consider coevolution. And it takes almost twice as long for SRE to become the most likely outcome in AncSR2.

### Structure of the GP maps under coevolution

Protein-DNA coevolution bias eolutionary outcomes more strongly towards conservation and delays phenotypic diversification. To understand why this happens, we will now examine the genotype networks and compare their structure to that of the protein-evolution networks.

We will compare 4 characteristics of the GP maps: 1) the number of RH neighbors per genotype, 2) the distance between phenotypes, 3) the fraction of neutral neighbors per genotype and 4) the local bias of the one-mutant neighborhood.

``` r
##########################
# Number of RH neighbors per genotype
##########################

# Detect unique protein neighbors in protein-DNA networks
detect_protein_neighbors <- function(g,graph){
  # neighbors of a given complex
  n <- neighbors(g,graph = graph)
  
  # RH genotypes of neighbor complexes and focal complex
  p_g <- str_sub(g,1,4)
  p_n <- str_sub(names(n),1,4)
  
  # number of unique RH neighbors
  r <- length(p_n[p_g != p_n])
  names(r) <- g
  return(r)
}

# compute number of neighbors per genotype
connect_sr1 <- degree(net_sr1)
connect_sr2 <- degree(net_sr2)
connect_sr1_complex <- degree(net_sr1_complex) # total number of neighbors per protein-DNA complex
connect_sr2_complex <- degree(net_sr2_complex) # total number of neighbors per protein-DNA complex
connect_sr1_complex_p <- map_int(.x=nodes_in_ntwrk_sr1_complex,.f=detect_protein_neighbors,graph=net_sr1_complex) # number of RH neighbors per protein-DNA complex
connect_sr2_complex_p <- map_int(.x=nodes_in_ntwrk_sr2_complex,.f=detect_protein_neighbors,graph=net_sr2_complex) # number of RH neighbors per protein-DNA complex

# merge data frames between networks and backgrounds
connect_df <- rbind(data.frame(bg="AncSR1",ntwrk_type="protein",neighbor_type="AA",connect=connect_sr1),
                    data.frame(bg="AncSR2",ntwrk_type="protein",neighbor_type="AA",connect=connect_sr2),
                    data.frame(bg="AncSR1",ntwrk_type="protein-DNA",neighbor_type="AA",connect=connect_sr1_complex_p),
                    data.frame(bg="AncSR2",ntwrk_type="protein-DNA",neighbor_type="AA",connect=connect_sr2_complex_p),
                    data.frame(bg="AncSR1",ntwrk_type="protein-DNA",neighbor_type="AA-DNA",connect=connect_sr1_complex),
                    data.frame(bg="AncSR2",ntwrk_type="protein-DNA",neighbor_type="AA-DNA",connect=connect_sr2_complex))

connect_df %>% filter(neighbor_type == "AA") %>%
  ggplot(aes(x = bg,y = connect)) +
  geom_violin(aes(fill=ntwrk_type),alpha=0.6) +
  geom_pointrange(aes(x=bg,group=ntwrk_type),stat = "summary",size=0.4,color="black",position = position_dodge(width=0.9)) + 
  scale_fill_manual(name="Network",values = c("white","gray30")) +
  labs(x="Background",y="Number of protein neighbors\n per genotype",color="") +
  theme_classic() +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust = 1,vjust = 1))
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-56-1.png)

``` r
knitr::kable(connect_df %>% filter(neighbor_type == "AA") %>% group_by(bg,ntwrk_type) %>%
               reframe(avg_n_RH_neighbors = mean(connect)))
```

| bg     | ntwrk\_type |  avg\_n\_RH\_neighbors|
|:-------|:------------|----------------------:|
| AncSR1 | protein     |               3.495327|
| AncSR1 | protein-DNA |               3.017857|
| AncSR2 | protein     |              10.737017|
| AncSR2 | protein-DNA |               7.720025|

``` r
# Anova
mod_connect <- lm(connect_df$connect~connect_df$bg*connect_df$ntwrk_type)
anova(mod_connect)
```

    ## Analysis of Variance Table
    ## 
    ## Response: connect_df$connect
    ##                                        Df Sum Sq Mean Sq F value    Pr(>F)    
    ## connect_df$bg                           1  11119 11118.6 611.232 < 2.2e-16 ***
    ## connect_df$ntwrk_type                   1   9643  9642.6 530.096 < 2.2e-16 ***
    ## connect_df$bg:connect_df$ntwrk_type     1    258   258.3  14.198 0.0001653 ***
    ## Residuals                           12426 226034    18.2                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

This shows that protein-DNA coevolution reduces the average number of RH neighbors relative to the protein-evolution scenario. Since genotypes are less well connected, this should affect the mutational distance between phenotypes in the network.

``` r
########################
# Average shortest path length between neutral networks
########################

# compute mutational distance betwen phenotypes
RE_combos <- expand.grid(REs[[1]],REs[[1]])

# AncSR1
pairwise_dists_sr1_complex <- RE_combos %>% mutate(links = map2_dbl(.x=as.character(Var1),.y=as.character(Var2),
                                                    pairwise_neutral_network_proximity,Bg="AncSR1",graph=net_sr1_complex,type=3,
                                                    pheno_df=phenotypes_tbl_complex,from_specific=T,to_specific=T,complex=T)) %>% acast(Var1~Var2, value.var="links")

pairwise_dists_sr1_complex <- pairwise_dists_sr1_complex[rowSums(is.na(pairwise_dists_sr1_complex)) != ncol(pairwise_dists_sr1_complex), ] # remove all NAs-rows
pairwise_dists_sr1_complex <- pairwise_dists_sr1_complex[, colSums(is.na(pairwise_dists_sr1_complex)) != nrow(pairwise_dists_sr1_complex)] # remove all NAs-cols
pairwise_dists_sr1_complex <- pairwise_dists_sr1_complex[order(match(rownames(pairwise_dists_sr1_complex),REs[[1]])),order(match(colnames(pairwise_dists_sr1_complex),REs[[1]]))] # order rows and columns

# AncSR2
pairwise_dists_sr2_complex <- RE_combos %>% mutate(links = map2_dbl(.x=as.character(Var1),.y=as.character(Var2),
                                                    pairwise_neutral_network_proximity,Bg="AncSR2",graph=net_sr2_complex,type=3,
                                                    pheno_df=phenotypes_tbl_complex,from_specific=T,to_specific=T,complex=T)) %>% acast(Var1~Var2, value.var="links")

pairwise_dists_sr2_complex <- pairwise_dists_sr2_complex[rowSums(is.na(pairwise_dists_sr2_complex)) != ncol(pairwise_dists_sr2_complex), ] # remove all NAs-rows
pairwise_dists_sr2_complex <- pairwise_dists_sr2_complex[, colSums(is.na(pairwise_dists_sr2_complex)) != nrow(pairwise_dists_sr2_complex)] # remove all NAs-cols
pairwise_dists_sr2_complex <- pairwise_dists_sr2_complex[order(match(rownames(pairwise_dists_sr2_complex),REs[[1]])),order(match(colnames(pairwise_dists_sr2_complex),REs[[1]]))] # order rows and columns
```

``` r
# correct prot-DNA steps to AA subs
pairwise_dists_sr1_complex <- pairwise_dists_sr1_complex * AVG_AA_SUBS_PER_STEP_SR1
pairwise_dists_sr2_complex <- pairwise_dists_sr2_complex * AVG_AA_SUBS_PER_STEP_SR2 

# merge data frames between networks and backgrounds
pairwise_dists_df <- rbind(inner_join(melt(pairwise_dists_sr1) %>% mutate(Pair = map2_chr(Var1, Var2, ~toString(sort(c(.x, .y)))), bg="AncSR1",net="protein") %>% distinct(Pair, .keep_all = TRUE),
                                      melt(pairwise_dists_sr1_complex) %>% mutate(Pair = map2_chr(Var1, Var2, ~toString(sort(c(.x, .y)))), bg="AncSR1",net="protein-DNA") %>% distinct(Pair, .keep_all = TRUE), by="Pair"),
                           inner_join(melt(pairwise_dists_sr2) %>% mutate(Pair = map2_chr(Var1, Var2, ~toString(sort(c(.x, .y)))), bg="AncSR2",net="protein") %>% distinct(Pair, .keep_all = TRUE),
                                      melt(pairwise_dists_sr2_complex) %>% mutate(Pair = map2_chr(Var1, Var2, ~toString(sort(c(.x, .y)))), bg="AncSR2",net="protein-DNA") %>% distinct(Pair, .keep_all = TRUE), by="Pair")) %>%
  filter(!(Var1.x == Var2.y)) # remove diagonals
 
pairwise_dists_df %>% ggplot(aes(x=value.x,y=value.y)) +
  geom_point(aes(shape=bg.x), size=2) + 
  geom_abline(slope=1,intercept = 0,linetype="dashed") + 
  theme_classic() +
  scale_shape_manual(values = c(16,17)) +
  scale_x_continuous(limits = c(3,10)) + scale_y_continuous(limits = c(3,10)) +
  labs(x=expression("Mean pairwise proximity"[protein]),y=expression("Mean pairwise proximity"[protein-DNA]),shape="Background") +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10)) +
  guides(fill = guide_legend(override.aes = list(shape=21,size=3)))
```

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-58-1.png)

As expected, coevolution also pushes phenotypes farther away in the network. Let's now explore the one-mutant neighborhoods around specific genotypes.

``` r
##########################
# Fraction of neutral neighbors
##########################

cl <- makeCluster(N_CORES,type = "FORK",outfile="")
registerDoParallel(cl)

# compute number of neutral neighbors per genotype
neighbors_sr1_complex <- foreach(i = 1:length(nodes_in_ntwrk_sr1_spec_complex), .combine = 'rbind') %dopar% {
  fraction_of_neighbors_per_genotype(from_g = nodes_in_ntwrk_sr1_spec_complex[i],Bg = "AncSR1",graph = net_sr1_complex,pheno_tbl = phenotypes_tbl_complex,complex = T)}

neighbors_sr2_complex <- foreach(i = 1:length(nodes_in_ntwrk_sr2_spec_complex), .combine = 'rbind') %dopar% {
  fraction_of_neighbors_per_genotype(from_g = nodes_in_ntwrk_sr2_spec_complex[i],Bg = "AncSR2",graph = net_sr2_complex,pheno_tbl = phenotypes_tbl_complex,complex = T)}

stopCluster(cl)
```

``` r
# merge data frames between networks and backgrounds
neutral_neighbors_df <- rbind(neighbors_sr1 %>% mutate(ntwrk_type="protein",bg="AncSR1"),
                             neighbors_sr2 %>% mutate(ntwrk_type="protein",bg="AncSR2"),
                             neighbors_sr1_complex %>% mutate(ntwrk_type="protein-DNA",bg="AncSR1"),
                             neighbors_sr2_complex %>% mutate(ntwrk_type="protein-DNA",bg="AncSR2")) %>% filter(pheno_match == TRUE)

neutral_neighbors_df %>%
  ggplot(aes(x = bg,y = prop)) +
  geom_violin(aes(fill=ntwrk_type),alpha=0.6) +
  geom_pointrange(aes(x=bg,group=ntwrk_type),stat = "summary",size=0.4,color="black",position = position_dodge(width=0.9)) + 
  scale_fill_manual(name="Network",values = c("white","gray30")) +
  labs(x="Background",y="Proportion of neutral neighbors",color="") +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25)) +
  theme_classic() +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust = 1,vjust = 1))
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-60-1.png)

``` r
knitr::kable(neutral_neighbors_df %>% group_by(bg,ntwrk_type) %>%
               reframe(avg_prop_neutral_neighbors = mean(prop)))
```

| bg     | ntwrk\_type |  avg\_prop\_neutral\_neighbors|
|:-------|:------------|------------------------------:|
| AncSR1 | protein     |                      0.7893705|
| AncSR1 | protein-DNA |                      0.8280453|
| AncSR2 | protein     |                      0.5384547|
| AncSR2 | protein-DNA |                      0.6113786|

``` r
# Anova
mod_neutral <- lm(neutral_neighbors_df$prop~neutral_neighbors_df$bg*neutral_neighbors_df$ntwrk_type)
anova(mod_neutral)
```

    ## Analysis of Variance Table
    ## 
    ## Response: neutral_neighbors_df$prop
    ##                                                           Df  Sum Sq Mean Sq
    ## neutral_neighbors_df$bg                                    1   7.117  7.1171
    ## neutral_neighbors_df$ntwrk_type                            1   3.276  3.2759
    ## neutral_neighbors_df$bg:neutral_neighbors_df$ntwrk_type    1   0.038  0.0384
    ## Residuals                                               2589 146.059  0.0564
    ##                                                          F value    Pr(>F)    
    ## neutral_neighbors_df$bg                                 126.1549 < 2.2e-16 ***
    ## neutral_neighbors_df$ntwrk_type                          58.0668 3.531e-14 ***
    ## neutral_neighbors_df$bg:neutral_neighbors_df$ntwrk_type   0.6811    0.4093    
    ## Residuals                                                                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

This shows that protein-DNA coevolution increases the average number of neutral neighbors relative to the protein-evolution scenario. Because of this, we should also expect the local bias in the one-mutant neighborhood to increase relative to the protein-evolution networks.

``` r
#####################
# Local bias in one-mutant neighborhood
#####################

cl <- makeCluster(N_CORES,type = "FORK",outfile="")
registerDoParallel(cl)

local_bias_sr1_complex <- foreach(i = 1:length(nodes_in_ntwrk_sr1_spec_complex), .combine = 'rbind') %dopar% {
  n_i <- names(neighbors(graph = net_sr1_complex,v = nodes_in_ntwrk_sr1_spec_complex[i]))
  p_n_i <- phenotypes_tbl_complex %>% ungroup() %>% filter(bg=="AncSR1") %>% filter(complex %in% n_i) %>%
    reframe(RE = REs[[1]],
            count = table(factor(specificity,levels=REs[[1]])),
            fr = count/length(n_i))
  local_b_i <- get_bias(as.vector(p_n_i$fr),input_freqs = T)
  data.frame(g=nodes_in_ntwrk_sr1_spec_complex[i],b=local_b_i)}

local_bias_sr2_complex <- foreach(i = 1:length(nodes_in_ntwrk_sr2_spec_complex), .combine = 'rbind') %dopar% {
  n_i <- names(neighbors(graph = net_sr2_complex,v = nodes_in_ntwrk_sr2_spec_complex[i]))
  p_n_i <- phenotypes_tbl_complex %>% ungroup() %>% filter(bg=="AncSR2") %>% filter(complex %in% n_i) %>%
    reframe(RE = REs[[1]],
            count = table(factor(specificity,levels=REs[[1]])),
            fr = count/length(n_i))
  local_b_i <- get_bias(as.vector(p_n_i$fr),input_freqs = T)
  data.frame(g=nodes_in_ntwrk_sr2_spec_complex[i],b=local_b_i)}

stopCluster(cl)
```

``` r
# merge data frames between networks and backgrounds
local_bias_df <- rbind(local_bias_sr1 %>% mutate(ntwrk_type="protein",bg="AncSR1"),
                       local_bias_sr2 %>% mutate(ntwrk_type="protein",bg="AncSR2"),
                       local_bias_sr1_complex %>% mutate(ntwrk_type="protein-DNA",bg="AncSR1"),
                       local_bias_sr2_complex %>% mutate(ntwrk_type="protein-DNA",bg="AncSR2"))
local_bias_df %>%
  ggplot(aes(x = bg,y = b)) +
  geom_violin(aes(fill=ntwrk_type),alpha=0.6) +
  geom_pointrange(aes(x=bg,group=ntwrk_type),stat = "summary",size=0.4,color="black",position = position_dodge(width=0.9)) + 
  scale_fill_manual(name="Network",values = c("white","gray30")) +
  labs(x="Background",y="Local bias (B)",color="") +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25)) +
  theme_classic() +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle=45,hjust = 1,vjust = 1))
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Evolutionary_modeling_GPmaps_files/figure-markdown_github/unnamed-chunk-62-1.png)

``` r
knitr::kable(local_bias_df %>% group_by(bg,ntwrk_type) %>%
               reframe(avg_local_bias_per_genotype = mean(b)))
```

| bg     | ntwrk\_type |  avg\_local\_bias\_per\_genotype|
|:-------|:------------|--------------------------------:|
| AncSR1 | protein     |                        0.9129443|
| AncSR1 | protein-DNA |                        0.9581831|
| AncSR2 | protein     |                        0.8430379|
| AncSR2 | protein-DNA |                        0.9152390|

``` r
# Anova
mod_localb <- lm(local_bias_df$b~local_bias_df$bg*local_bias_df$ntwrk_type)
anova(mod_localb)
```

    ## Analysis of Variance Table
    ## 
    ## Response: local_bias_df$b
    ##                                             Df  Sum Sq Mean Sq  F value  Pr(>F)
    ## local_bias_df$bg                             1  0.4407  0.4407  84.8423 < 2e-16
    ## local_bias_df$ntwrk_type                     1  3.4774  3.4774 669.3758 < 2e-16
    ## local_bias_df$bg:local_bias_df$ntwrk_type    1  0.0265  0.0265   5.0962 0.02406
    ## Residuals                                 2779 14.4367  0.0052                 
    ##                                              
    ## local_bias_df$bg                          ***
    ## local_bias_df$ntwrk_type                  ***
    ## local_bias_df$bg:local_bias_df$ntwrk_type *  
    ## Residuals                                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Again, as expected we see that coevolution increases the local bias in the one-mutant neighborhood. Overall, these results show that coevolution slows phenotypic diversification because it re-structures the genotype network such that evolutionary trajectories are more likely to remain in regions with the same starting specificity.

------------------------------------------------------------------------

## Exporting of genotype networks for visualization

Below is the code to export the genotype networks to a `gefx` file which can then be opened and edited in the sofware Gephi.

``` r
#######################
# networks for ploting
#######################

# protein evolution
sr1_variants_annotated <- phenotypes_tbl_prot %>% filter(bg == "AncSR1") %>% select(AA_var,specificity,max_meanF_bREs) # extract phenotypes per node
sr2_variants_annotated <- phenotypes_tbl_prot %>% filter(bg == "AncSR2") %>% select(AA_var,specificity,max_meanF_bREs) # extract phenotypes per node

net_sr1_fig <- build_genotype_network(sr1_variants_annotated,build_mat=F,adj_mat = adj_mat_sr1,type=1,node_attributes=T,cores=N_CORES) # network with annotations per node
net_sr2_fig <- build_genotype_network(sr2_variants_annotated,build_mat=F,adj_mat = adj_mat_sr2,type=1,node_attributes=T,cores=N_CORES) # network with annotations per node

#Export networks as GEFX file to open with GEPHI
library(rgexf)
output_dir = "./Networks/v5_final" # modify this

graph1 <- igraph.to.gexf(net_sr1_fig)
write.gexf(graph1$nodes[,1:2],graph1$edges[,2:3],nodesAtt=graph1$nodes[,c(1,3:4)],
           output = file.path(output_dir,"AncSR1_NET.gexf")) 

graph2 <- igraph.to.gexf(net_sr2_complex_fig)
write.gexf(graph2$nodes[,1:2],graph2$edges[,2:3],nodesAtt=graph2$nodes[,c(1,3:4)],
           output = file.path(output_dir,"AncSR2_NET.gexf"))


########


# protein-DNA coevolution
sr1_complexes_annotated <- phenotypes_tbl_complex %>% filter(bg == "AncSR1") %>% select(complex,specificity,meanF_bREs) # extract phenotypes per node
sr2_complexes_annotated <- phenotypes_tbl_complex %>% filter(bg == "AncSR2") %>% select(complex,specificity,meanF_bREs) # extract phenotypes per node


adj_mat_complex_sr1_plot <- as_adjacency_matrix(net_sr1_complex,type = "both") # adjacency matrix
links <- graph.adjacency(adj_mat_complex_sr1_plot,weighted = "1")
edges <- get.edgelist(links) # generate edgelist from matrix
net_sr1_complex_fig <- graph_from_data_frame(d=edges, vertices=sr1_complexes_annotated, directed=F) # network with annotations per node

adj_mat_complex_sr2_plot <- as_adjacency_matrix(net_sr2_complex,type = "both") # adjacency matrix
links <- graph.adjacency(adj_mat_complex_sr2_plot,weighted = "1")
edges <- get.edgelist(links) # generate edgelist from matrix
net_sr2_complex_fig <- graph_from_data_frame(d=edges, vertices=sr2_complexes_annotated, directed=F) # network with annotations per node

#Export networks as GEFX file to open with GEPHI
graph1 <- igraph.to.gexf(net_sr1_complex_fig)
write.gexf(graph1$nodes[,1:2],graph1$edges[,2:3],nodesAtt=graph1$nodes[,c(1,3:4)],
           output = file.path(output_dir,"AncSR1_NET_complex.gexf"))

graph2 <- igraph.to.gexf(net_sr2_complex_fig)
write.gexf(graph2$nodes[,1:2],graph2$edges[,2:3],nodesAtt=graph2$nodes[,c(1,3:4)],
           output = file.path(output_dir,"AncSR2_NET_complex.gexf"))
```