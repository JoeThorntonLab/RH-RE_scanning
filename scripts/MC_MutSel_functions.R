
# Functions to build, simulate and work with discrete Markov processes on RH-RE DMS data

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
REs[[2]] <- factor(c("SRE1 (AA)", "SRE2 (GA)", "ERE (GT)", "AC", 
                     "AG", "AT", "CA", "CC", 
                     "CG", "CT", "GC", "GG", 
                     "TA", "TC", "TG", "TT","Promiscuous"), 
                   levels=c("ERE (GT)", "SRE1 (AA)", "SRE2 (GA)", "AC", 
                            "AG", "AT", "CA", "CC", 
                            "CG", "CT","GC", "GG", 
                            "TA", "TC", "TG", "TT","Promiscuous"))

# Assign colors to each RE for plotting (and matching the colors in the networks)
rgb_to_hex <- function(rgb_col){
  hex <- sapply(strsplit(rgb_col, " "), function(x)
    rgb(x[1], x[2], x[3], maxColorValue=255))
  hex
}
rgb_RE_colors <- c("0 172 88", "255 63 53", "103 71 151", "37 121 77", "58 116 157", "6 183 0", "0 190 217", "0 145 145", "244 113 0",
                   "0 173 96", "0 169 254", "95 132 29", "184 78 86", "163 101 27", "255 47 136", "196 36 201")
hex_RE_colors <- rgb_to_hex(rgb_RE_colors)
names(hex_RE_colors) <- REs[[1]]

# Assign colors to DBD backgrounds
DBD_color <- function(names=T){
  colors <- c("#674797","#00AC58")
  if(names) names(colors) <- c("AncSR1","AncSR2")
  return(colors)
}

# Determine whether two amino acid sequences can be connected by a single nucleotide change given the genetic code
connect_aa_variants <- function(var1,var2){
  # Returns "1" if aa variants can be connected by single step mutation given genetic code, "0" otherwise
  var1 <- unlist(strsplit(var1,''))
  var2 <- unlist(strsplit(var2,''))
  aa_mismatches <- sum(var1 != var2)
  link <- 0
  if(aa_mismatches == 1){ #If not the same variant and not more than 2 aa diffs
    aa_diff <- which(var1 != var2) # position where variants differ
    codons_aa1 <- names(GENETIC_CODE[GENETIC_CODE==var1[aa_diff]])
    codons_aa2 <- names(GENETIC_CODE[GENETIC_CODE==var2[aa_diff]])
    # Compare all possible codon pairs: If there is at least one pair that can be interchanged by *a single* nuc change, 
    # then the link exists between AA variants
    pairs <- expand.grid(codons_aa1,codons_aa2) %>% rowwise() %>%
      mutate(nuc_mismatches = mapply(function(c1,c2) sum(c1!=c2), strsplit(as.character(Var1),''), strsplit(as.character(Var2),'')))
    if(1 %in% pairs$nuc_mismatches) link <- 1
  }
  return(link)
}

# Determine whether two amino acid sequences can be connected by a single nucleotide change given the genetic code
# The returned values can be interpreted as the relatiive propensity to go from v1-->v2 via single nucleotide mutations
# due to the genetic code.
connect_aa_variants_fraction <- function(var1,var2){
  # If aa variants can be connected by single step mutation given genetic code, returns the fraction of codons
  # of 'var1' that can access 'var2'; "0" otherwise.
  var1 <- unlist(strsplit(var1,''))
  var2 <- unlist(strsplit(var2,''))
  aa_mismatches <- sum(var1 != var2)
  link <- 0
  if(aa_mismatches == 1){ #If not the same variant and not more than 2 aa diffs
    aa_diff <- which(var1 != var2) # position where variants differ
    codons_aa1 <- names(GENETIC_CODE[GENETIC_CODE==var1[aa_diff]])
    codons_aa2 <- names(GENETIC_CODE[GENETIC_CODE==var2[aa_diff]])
    # Compare all possible codon pairs: If there is at least one pair that can be interchanged by *a single* nuc change, 
    # then the link exists between AA variants
    pairs <- expand.grid(codons_aa1,codons_aa2) %>% rowwise() %>%
      mutate(nuc_mismatches = mapply(function(c1,c2) sum(c1!=c2), strsplit(as.character(Var1),''), strsplit(as.character(Var2),'')))
    # Compute fraction of codons if at least one pair can be interconverted
    if(1 %in% pairs$nuc_mismatches) {
      link <- ungroup(pairs) %>% filter(codons_aa1 %in% Var1 & nuc_mismatches == 1) %>% with(length(unique(Var1))/length(codons_aa1))
    }
  }
  return(link)
}

# Determine whether two amino acid sequences can be connected by a single nucleotide change given the genetic code
# The returned value is the total number of mutations by which codons of var1 can access codons of var2
connect_aa_variants_count <- function(var1,var2){
  # If aa variants can be connected by single step mutation given genetic code, returns the fraction of codons
  # of 'var1' that can access 'var2'; "0" otherwise.
  var1 <- unlist(strsplit(var1,''))
  var2 <- unlist(strsplit(var2,''))
  aa_mismatches <- sum(var1 != var2)
  link <- 0
  if(aa_mismatches == 1){ #If not the same variant and not more than 2 aa diffs
    aa_diff <- which(var1 != var2) # position where variants differ
    codons_aa1 <- names(GENETIC_CODE[GENETIC_CODE==var1[aa_diff]])
    codons_aa2 <- names(GENETIC_CODE[GENETIC_CODE==var2[aa_diff]])
    # Compare all possible codon pairs: If there is at least one pair that can be interchanged by *a single* nuc change, 
    # then the link exists between AA variants
    pairs <- expand.grid(codons_aa1,codons_aa2) %>% rowwise() %>%
      mutate(nuc_mismatches = mapply(function(c1,c2) sum(c1!=c2), strsplit(as.character(Var1),''), strsplit(as.character(Var2),'')))
    # Compute number of codon paths
    n <- sum(pairs$nuc_mismatches == 1)
    # we have to multiply the number of codon paths across all invariant amino acid positions to get the total number of “backgrounds” 
    # that each amino acid mutation can occur in.
    if(length(var1 != 1)){
      bg_codons <- lapply(var1[-aa_diff],function(x) names(GENETIC_CODE)[GENETIC_CODE == x])
      nbg <- do.call(prod, lapply(bg_codons, length))
      link <- n * nbg
    }
    else link <- n
  }
  return(link)
}

# Function to calculate hamming distance between to variants
hamming_d <- function(var1,var2){
  var1 <- unlist(strsplit(var1,''))
  var2 <- unlist(strsplit(var2,''))
  d <- sum(var1 != var2)
  return(d)
}

# Function to connect protein-DNA complex variants. A link exists whether protein genotypes can be be connected by a single 
# nucleotide change given the genetic code, OR DNA variants differ by a hamming distance of one. Variants are given as 
# prot+DNA genotypes (e.g., EGKAGT)
#connect_complex <- functiona(var1,var2){
#  # split variants into protein and DNA genotypes.
#  var1_prot <- unlist(strsplit(var1,''))[1:4]
#  var1_dna <- unlist(strsplit(var1,''))[5:6]
#  var2_prot <- unlist(strsplit(var2,''))[1:4]
#  var2_dna <- unlist(strsplit(var2,''))[5:6]
#  
#  # compute links
#  link_prot <- connect_aa_variants(var1_prot,var2_prot)
#  link_dna <- hamming_d(var1_dna,var2_dna)
#  if()
#}

# Function to conect codons. Note that we are assuming no mutation bias, therefore, transitions between codons are encoded as 0 or 1.
connect_codons <- function(var1,var2){
  mismatches <- mapply(function(v1,v2) sum(v1!=v2), strsplit(var1,''), strsplit(var2,''))
  link <- 0
  if(mismatches == 1) link <- 1
  return(link)
}

# Function to find stationary distributions of amino acids from stationary distributions of codons
pi_codon_to_aa <- function(pi_codon){
  # pi_codon: vector of codon stationary distribution
  c <- data.frame(codon=names(GENETIC_CODE),aa=GENETIC_CODE) # genetic code
  pc <- data.frame(codon=names(pi_codon),freq=pi_codon) # stationary distibution of codon states
  r <- inner_join(c,pc,by="codon") %>% group_by(aa) %>% reframe(freq=sum(freq)) # sum frequencies across synonymous codons for each amino acid state
  pi_aa <- r %>% pull(freq)
  names(pi_aa) <- r$aa
  return(pi_aa)
}

# Build a mutation matrix (adjacency matrix where each entry is the mutation rate)
build_mutation_matrix <- function(nodes,type=1,cores=1){
  # nodes: a vector with vertex character names or data.frame (see argument 'node_attributes')
  # type: specify how to build the adjacency matrix. 1 - simple link through genetic code; 2 - fraction
    # of codons; 3 - mutational paths (default: 1)
  # cores: specify the number of cores to use in to use for parallel computing (default: 1)
  
  # Check for inconistent arguments
  if(is.null(nodes)){
    stop("Provide a vector of nodes to build the matrix")
  }
  if(!(type %in% seq(1,3,1))){
    stop("Invalid matrix type. Valid matrix options: 1, 2, or 3")
  }
  
  # Create adjacency matrix
  # Adjacency matrices are substantially sparse because they only store values amongst one-mutant neighbors, 
  # with the rest of cells being 0 --> Coerce outputs to to sparse matrices. Saves a significant amount of 
  # memory and speed up the processing of that data.
  # Adj. matrix --> From = rows, to = cols
  cl <- parallel::makeCluster(cores,"FORK",outfile="")
  doParallel::registerDoParallel(cl)
  
  if(type == 1){
    # Determine whether two amino acid sequences can be connected by a single nucleotide change given the genetic code
    adj_mat <- Matrix(foreach(i = 1:length(nodes), .combine = 'cbind') %dopar%
                        mapply(connect_aa_variants,nodes,nodes[i]) %>% `colnames<-`(nodes), sparse = T)
  }
  else if(type == 2){
    # fraction of codons from i that can mutate to j
    adj_mat <- Matrix(foreach(i = 1:length(nodes), .combine = 'cbind') %dopar%
                        mapply(connect_aa_variants_fraction,nodes,nodes[i]) %>% `colnames<-`(nodes), sparse = T)
  }
  else if(type == 3){
    # number of mutational paths from i to j
    adj_mat <- Matrix(foreach(i = 1:length(nodes), .combine = 'cbind') %dopar%
                        mapply(connect_aa_variants_count,nodes,nodes[i]) %>% `colnames<-`(nodes), sparse = T)
  }
  stopCluster(cl)
  
  return(adj_mat)
}

# Build a genotype network (igraph object)
build_genotype_network <- function(nodes,build_mat = TRUE,adj_mat=NULL,weighted=FALSE,node_attributes=FALSE,cores=1){
  # nodes = a vector with vertex character names or data.frame (see argument 'node_attributes')
  # build_mat = logical value to indicate whether to generate adjacency matrix (default: TRUE)
  # adj_mat = if build_mat=FALSE, then pass the corresponding adjacency matrix to generate the network (default: NULL)
  # weigthed = a string to indicate whether adjacency matrix has weighted edges: "fraction","count",FALSE (default: FALSE)
  # node_attributes = logical value to indicate whether the nodes contain annotations. If TRUE 'nodes' should be a data frame containing
    # the attrbitues per node (each column contains a node attribute). This option is mostly for plotting the network, so only undirected 
    # networks will be built. Can also work with build_mat=TRUE/FALSE.
  # cores: specify the number of cores to use in to use for parallel computing (default: 1)
  
  # Check for inconistent arguments
  if(build_mat==FALSE && is.null(adj_mat)){
    stop("Provide an adjacency matrix to build the network or set argument 'create_mat = TRUE'")
  }
  if(node_attributes && is.null(dim(nodes))){
    stop("Nodes have no attributes. Provide a data frame with annotations per node.")
  }
  
  if(build_mat==FALSE){
    # Build networks using supplied adjacency matrix
    if(weighted != FALSE){
      net <- graph_from_adjacency_matrix(adj_mat,mode="directed",weighted=TRUE)
    }
    else{
      if(node_attributes){
        links <- graph.adjacency(adj_mat,mode = "undirected")
        edges <- get.edgelist(links) # generate edgelist from matrix
        net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) # network with annotations per node
      }
      else net <- graph_from_adjacency_matrix(adj_mat,mode = "undirected")
    }
  }
  
  else if(build_mat){
    # Create adjacency matrix and build genotype network
    # Adjacency matrices are substantially sparse because they only store values amongst one-mutant neighbors, 
    # with the rest of cells being 0 --> Coerce outputs to to sparse matrices. Saves a significant amount of 
    # memory and speed up the processing of that data.
    # Adj. matrix --> From = rows, to = cols
    cl <- parallel::makeCluster(cores,"FORK",outfile="")
    doParallel::registerDoParallel(cl)
    
    if(weighted=="fraction"){
      # weights correspond to "mutational propensity" of change from var1 to var2
      adj_mat <- Matrix(foreach(i = 1:length(nodes), .combine = 'cbind') %dopar%
                          mapply(connect_aa_variants_fraction,nodes,nodes[i]) %>% `colnames<-`(nodes), sparse = T)
      
      net <- graph_from_adjacency_matrix(adj_mat,mode="directed",weighted=TRUE)
    }
    else if(weighted=="count"){
      # weights correspond to the number of mutations from var1 to var2
      adj_mat <- Matrix(foreach(i = 1:length(nodes), .combine = 'cbind') %dopar%
                          mapply(connect_aa_variants_count,nodes,nodes[i]) %>% `colnames<-`(nodes), sparse = T)
      
      net <- graph_from_adjacency_matrix(adj_mat,mode="directed",weighted=TRUE)
    }
    else{
      if(node_attributes){
        nodes_v <- nodes %>% pull(AA_var)
        adj_mat <- Matrix(foreach(i = 1:length(nodes_v), .combine = 'cbind') %dopar%
                            mapply(connect_aa_variants,nodes_v,nodes_v[i]) %>% `colnames<-`(nodes_v), sparse = T)
        links <- graph.adjacency(adj_mat,mode = "undirected")
        edges <- get.edgelist(links) # generate edgelist from matrix
        net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) # network with annotations per node
      }
      else{
        adj_mat <- Matrix(foreach(i = 1:length(nodes), .combine = 'cbind') %dopar%
                            mapply(connect_aa_variants,nodes,nodes[i]) %>% `colnames<-`(nodes), sparse = T)
        net <- graph_from_adjacency_matrix(adj_mat,mode = "undirected")
      }
    }
    stopCluster(cl)
  }
  return(net)
}

# Extract a variant's phenotypic value
get_phenotype <- function(aa_var,Bg,pheno_table,phenotype="mean"){
  # aa_var = variant or variants (can be a single string or a string vector)
  # Bg = string specifying the DBD background (= 'AncSR1' or 'AncSR2')
  # pheno_table = a data frame with phenotypic values (meanF)
  # phenotype = a string specifying which 'phenotype' to extract: 
  # 'mean' = avg meanF across bound DNA elements (default)
  # 'min' = min meanF across bound DNA elements
  # 'max' = max meanF across bound DNA elements
  p <- NA
  n <-  pheno_table %>% filter(AA_var %in% aa_var & bg == Bg) %>% pull(AA_var)
  if(phenotype == "mean"){ p <- pheno_table %>% filter(AA_var %in% aa_var & bg == Bg) %>% pull(meanF_bREs)}
  if(phenotype == "min"){ p <- pheno_table %>% filter(AA_var %in% aa_var & bg == Bg) %>% pull(min_meanF_bREs)}
  if(phenotype == "max"){ p <- pheno_table %>% filter(AA_var %in% aa_var & bg == Bg) %>% pull(max_meanF_bREs)}
  names(p) <- n
  return(p)
}

## Fitness functions: ##
# Logistic function: For scenario of directional selection and drift
fitness_logistic <- function(mF,L,k,x_o){
  # mF = meanF of allele, i.e., phenotype (dbl)
  # L = maximum fitness (dbl)
  # x_o = midpoint (dbl)
  # k = logistic grouth rate (steepness) (dbl)
  # Fitness corresponds to the exponential growth rate = N*r
  mF <- mF - AncSR1_ERE_ref # transform deltaF values to phenotypic values. dF = 0 --> meanF = AncSR1_ERE_ref
  f <- L / (1 + exp(-k*(mF-x_o)))
  f
}

# Logistic function: For scenario of directional selection and drift using a reference genotype
fitness_logistic_ref <- function(mF,L,k,x_o,mF_ref){
  # mF = meanF of allele, i.e., phenotype (dbl)
  # L = maximum fitness (dbl)
  # x_o = midpoint (dbl)
  # k = logistic grouth rate (steepness) (dbl)
  # Fitness corresponds to the exponential growth rate = N*r
  mF <- mF - AncSR1_ERE_ref # transform deltaF values to phenotypic values. dF = 0 --> meanF = AncSR1_ERE_ref
  mF_ref <- mF_ref - AncSR1_ERE_ref
  f <- ifelse(mF < mF_ref,0, L / (1 + exp(-k*(mF-x_o))))
  f
}

# Normal ('bell shaped') curve: For scenario of stabilizing selection
fitness_normal <- function(mF,u,sd,scale){
  # mF = meanF of allele, i.e., phenotype (dbl)
  # u = mean of the distribution (dbl)
  # sd = standard deviation (dbl)
  # scale = scaling factor to specify the maximum height (max fitness) (dbl)
  # Fitness corresponds to the exponential growth rate = N*r
  f <- exp(-0.5 * ((mF - u) / sd)^2) * scale
  f
}

# Normal ('bell shaped') curve: For scenario of stabilizing selection using a reference genotype
fitness_normal_ref <- function(mF,u,sd,scale,mF_ref){
  # mF = meanF of allele, i.e., phenotype (dbl)
  # u = mean of the distribution (dbl)
  # sd = standard deviation (dbl)
  # scale = scaling factor to specify the maximum height (max fitness) (dbl)
  # Fitness corresponds to the exponential growth rate = N*r
  f <- ifelse(mF < mF_ref,0, exp(-0.5 * ((mF - u) / sd)^2) * scale)
  f
}

# Step function: For scenario of Purifying selection + drift
fitness_purifying <- function(mF, mF_ref) {
  # mF = meanF of allele, i.e., phenotype (dbl)
  # mF_ref = meanF of reference genotype. Reference genotype determines the mF after which all genotypes have equivalent fitness.
  f <- ifelse(mF <= mF_ref, 0,1)
  f
}

# Calculate selection coefficient between allele i and j given a selection scenario:
get_sel_coeffs <- function(mF_i,mF_j,scenario,Ne,param){
  # mF_i = phenotype allele i (dbl)
  # mF_j = phenotype allele j (dbl)
  # scenario = a string indicating the selection scenario to model: 'drift', 'directional', 'stabilizing' (chr)
  # Ne = population size (dbl)
  # param = vector of parameters for the fitness funciton (dbl)
  
  # Use fitness functions to calculate fitness given the selection scenario
  if(scenario == "drift"){
    fit_i <- fitness_purifying(mF_i,param[1]) # N*r_i
    fit_j <- fitness_purifying(mF_j,param[1]) # N*r_j
  }
  if(scenario == "directional"){
    fit_i <- fitness_logistic(mF_i,param[1],param[2],param[3]) # N*r_i
    fit_j <- fitness_logistic(mF_j,param[1],param[2],param[3]) # N*r_j
  }
  if(scenario == "stabilizing"){
    fit_i <- fitness_normal(mF_i,param[1],param[2],param[3]) # N*r_i
    fit_j <- fitness_normal(mF_j,param[1],param[2],param[3]) # N*r_j
  }
  s_ij <- (fit_j - fit_i) / Ne # calculate selection coefficient of mutation
  return(s_ij)
}

# Kimura's equation for the probability of fixation of a new allele
pfix_kimura <- function(s_ij,Ne){
  # s_ij = selection coefficient (dbl)
  # Ne = population size (dbl)
  
  if(s_ij != 0){ p.fix <- (1-exp(-2*s_ij)) / (1-exp(-4*Ne*s_ij)) }
  else{ p.fix <- 1 / (2*Ne) } # fixation prob. under neutral case (s_ij = 0)
  
  return(p.fix)
}

# Compute transition rate R(i,j) between two variants
rate_mutational_step <- function(from_i,to_j,adj_matrix,pheno_table,scenario,Ne,fit.fn.param,Bg=NULL,mutation=FALSE,...){
  # from_i, to_j = strings specifying the name of variants
  # adj_matrix = an adjacency matrix
  # pheno_table = a data frame with phenotypic values for each variant
  # scenario: a string indicating the selection scenario to model: 'drift', 'directional', 'stabilizing' (chr)
  # Ne, fit.fn.param = population genetic parameters for the fitness functions
  # Bg = string specifying the DBD background (= 'AncSR1' or 'AncSR2')
  # mutation = logical value to indicate whether to model mutation. 
      # if TRUE, adj_matrix should contain the rates of mutation from i to j (default: FALSE)
  # ... = additional (optional) parameters 
  
  # Check for valid parameters:
    if(is.null(adj_matrix)){
    stop("Provide an adjacency matrix")
  }
  if(!(from_i %in% rownames(adj_matrix)) || !(to_j %in% colnames(adj_matrix))){
    stop("Invalid 'from_i' or 'to_j' argument: One or both variants not found in the matrix")
  }
  if(!(scenario %in% c("drift", "directional","stabilizing"))) {
    stop("Provide a valid selection scenario: 'drift', 'directional' or 'stabilizing'")
  }
  if(!(Bg %in% c("AncSR1","AncSR2")) || is.null(Bg)){
    stop("Specify DBD background: 'AncSR1' or 'AncSR2'")
  }

  r_ij <- 0 # initialize rate r_ij
  
  # If variants i and j are one-mutant neighbors, proceed to compute P(i,j)
  if(connect_aa_variants(from_i,to_j) == 1){
    # Extract phenotypes of focal variants
    mF_i = get_phenotype(from_i,Bg,pheno_table,...)
    mF_j = get_phenotype(to_j,Bg,pheno_table,...)
        
    # Calculate selection coefficient and fixation probability between i and j
    sel_coeff <- get_sel_coeffs(mF_i=mF_i,mF_j=mF_j,scenario=scenario,Ne=Ne,param=fit.fn.param)
    pfix <- pfix_kimura(s_ij=sel_coeff,Ne=Ne)
    
    # If modeling mutation, fixation probablities are weighted by the mutation rate. 
    if(mutation){
      # extract mutation rate for AA pair
      mut_rate <- adj_matrix[rownames(adj_matrix)==from_i, colnames(adj_matrix)==to_j]
      r_ij <- mut_rate * pfix
    }
    else r_ij <- pfix
    
    names(r_ij) <- to_j
  }
  return(r_ij)
}


# Compute probability P(i,j) for two variants: Prob. that *next* mutation is allele j
prob_mutational_step <- function(from_i,to_j,graph,pheno_table,scenario,Ne,fit.fn.param,Bg=NULL,mutation=FALSE,...){
  # from_i, to_j = strings specifying the name of variants
  # graph = an igraph object (the genotype network)
  # pheno_table = a data frame with phenotypic values for each variant
  # scenario: a string indicating the selection scenario to model: 'drift', 'directional', 'stabilizing' (chr)
  # Ne, fit.fn.param = population genetic parameters for the fitness functions
  # Bg = string specifying the DBD background (= 'AncSR1' or 'AncSR2')
  # mutation = logical value to indicate whether to model mutation as the relative propensities (default: FALSE)
  # ... = additional (optional) parameters 
  
  # Check for valid parameters:
  if(!(from_i %in% V(graph)$name) || !(to_j %in% V(graph)$name)){
    stop("Invalid 'from_i' or 'to_j' argument: One or both variants not found in the network")
  }
  if(!(scenario %in% c("drift", "directional","stabilizing"))) {
    stop("Provide a valid selection scenario: 'drift', 'directional' or 'stabilizing'")
  }
  if(!(Bg %in% c("AncSR1","AncSR2")) || is.null(Bg)){
    stop("Specify DBD background: 'AncSR1' or 'AncSR2'")
  }
  if(mutation && is.null(E(graph)$weight)){
    stop("Network has no weighted edges. Provide a weighted network or set 'mutation=FALSE'")
  }
  
  p_ij <- NA
  
  # If variants i and j are one-mutant neighbors, proceed to compute P(i,j)
  if(to_j %in% neighbors(graph,from_i)$name){
    # Extract phenotypes of focal variant
    mF_i = get_phenotype(from_i,Bg,pheno_table,...)
    
    #Extract local network around focal node and their phenotypes
    local_ntwrk <- neighbors(graph,from_i)$name
    pheno_neighbors <- get_phenotype(local_ntwrk,Bg,pheno_table = pheno_table,...)
    
    # Calculate selection coefficients for all k neighbors
    sel_coeffs <- map_dbl(pheno_neighbors,get_sel_coeffs,mF_i=mF_i,scenario=scenario,Ne=Ne,param=fit.fn.param)
    
    # Calculate fixation probabilities for all alternative k alleles in the local network:
    pfix_local_ntwrk <- map_dbl(sel_coeffs,pfix_kimura,Ne=Ne)
    
    # If modeling mutation, fixation probablities are weighted by the relative propensities of mutation 
    if(mutation){
      # extract weight info from local network: mutational propensities
      mut_prop <- get.data.frame(graph) %>% filter(from==from_i) %>% filter(local_ntwrk %in% to) %>% pull(weight)
      pfix_local_ntwrk <- mut_prop * pfix_local_ntwrk
    }
    
    # Calculate P(i,j): Rescale fixation probabilities
    if(length(local_ntwrk) == 1 && pfix_local_ntwrk == 0){
      p_ij <-  0 # if single neighbor is strongly deleterious (avoid 0/0 division)
    } 
    else{
      p_ij <- pfix_local_ntwrk[names(pfix_local_ntwrk)==to_j] / sum(pfix_local_ntwrk)
    }
    names(p_ij) <- to_j
  }
  return(p_ij)
}

# Build transition probability matrix (P(i,j)'s) for a DBD background and a specific selection scenario
build_transition_matrix <- function(nodes,graph,pheno_tbl,param,cores=1){
  # nodes = a vector with vertex character names
  # graph: A genotype network (igraph object)
  # pheno_tbl: a data frame containing the phenotypic information for each RH variant
  # param: A list of parameters indicating the population genetic scenario to model (in the following order):
    # Bg: string specifying the DBD background (= 'AncSR1' or 'AncSR2')
    # scenario: a string indicating the selection scenario to model: 'drift', 'directional', 'stabilizing' (chr)
    # Ne: Effective population size
    # fit.fn.param: Parameters of the fitness function
    # mutation: a logical value to indicate whether to include mutational propensities in the computation of P(i,j)
  # phenotype: phenotypic value to use from RH variants
  # cores: specify the number of cores to use in to use for parallel computing (default: 1)
  
  # Transition matrices are substatially sparse because they only store P(i,j) values amongst one-mutant neighbors, with the rest of cells
  # being NA --> Coerce outputs to to sparse matrices. (output is a dgCMatrix-like matrix (from package Matrix) but instead of zeros, NAs are 
  # not explicitly stored). Use 'dropNA2matrix(x)' to restore NA values. 
  # *CAUTION*: Zeros are represented with a very small value (.Machine$double.xmin) so they do not get dropped in the sparse representation. 
  # Actual zero ratings have a small, but non-zero value.
  
  # Begin
  print(paste("Building transition matrix for model: ",deparse(substitute(param))))
  
  # Extract population genetic and model parameters:
  Bg <- param[[1]]
  scenario <- param[[2]]
  Ne <- param[[3]]
  fit.fn.param <- param[[4]]
  mutation <- param[[5]]
  pheno <- param[[6]]
  
  
  cl <- parallel::makeCluster(cores,"FORK",outfile="")
  doParallel::registerDoParallel(cl)
  
  # Build matrix (From = rows, to = cols):
  tr_mat <- dropNA(foreach(i = 1:length(nodes), .combine = 'rbind') %dopar%
                     map_dbl(nodes,prob_mutational_step,from_i=nodes[i],
                             graph=graph,pheno_table=pheno_tbl,
                             scenario=scenario,Ne=Ne,
                             fit.fn.param=fit.fn.param,Bg=Bg,
                             mutation=mutation,phenotype=pheno) %>% 
                     `colnames<-`(nodes) %>% `rownames<-`(nodes))
  stopCluster(cl)
  
  return(tr_mat)
  print("...Done!")
}

# Build transition probability matrix (P(i,j)'s) for a DBD background and a specific selection scenario
# (uses the function 'rate_mutational_step', instead of 'prob_mutational_step')
build_transition_matrix_v2 <- function(nodes,adj_matrix,pheno_tbl,param,cores=1){
  # nodes = a vector with vertex character names
  # adj_matrix = an adjacency matrix
  # pheno_tbl: a data frame containing the phenotypic information for each RH variant
  # param: A list of parameters indicating the population genetic scenario to model (in the following order):
    # Bg: string specifying the DBD background (= 'AncSR1' or 'AncSR2')
    # scenario: a string indicating the selection scenario to model: 'drift', 'directional', 'stabilizing' (chr)
    # Ne: Effective population size
    # fit.fn.param: Parameters of the fitness function
    # mutation: a logical value to indicate whether to include mutational propensities in the computation of P(i,j)
  # phenotype: phenotypic value to use from RH variants
  # cores: specify the number of cores to use in to use for parallel computing (default: 1)
  
  # Transition matrices are substatially sparse because they only store P(i,j) values amongst one-mutant neighbors, with the rest of cells
  # being NA --> Coerce outputs to to sparse matrices. (output is a dgCMatrix-like matrix (from package Matrix) but instead of zeros, NAs are 
  # not explicitly stored). Use 'dropNA2matrix(x)' to restore NA values. 
  # *CAUTION*: Zeros are represented with a very small value (.Machine$double.xmin) so they do not get dropped in the sparse representation. 
  # Actual zero ratings have a small, but non-zero value.
  
  # Begin
  print(paste("Building transition matrix for model: ",deparse(substitute(param))))
  
  # Extract population genetic and model parameters:
  Bg <- param[[1]]
  scenario <- param[[2]]
  Ne <- param[[3]]
  fit.fn.param <- param[[4]]
  mutation <- param[[5]]
  pheno <- param[[6]]
  
  
  cl <- parallel::makeCluster(cores,"FORK",outfile="")
  doParallel::registerDoParallel(cl)
  
  # Build matrix (From = rows, to = cols):
  tr_mat <- foreach(i = 1:length(nodes), .combine = 'rbind') %dopar%
                     map_dbl(nodes,rate_mutational_step,from_i=nodes[i],
                             adj_matrix=adj_matrix,pheno_table=pheno_tbl,
                             scenario=scenario,Ne=Ne,
                             fit.fn.param=fit.fn.param,Bg=Bg,
                             mutation=mutation,phenotype=pheno) %>% 
                     `colnames<-`(nodes) %>% `rownames<-`(nodes)
  stopCluster(cl)

  # Normalize 'tr_mat' to make it a transition probability matrix (all rows sum to one)
  tr_mat <- t(apply(tr_mat, 1, function(x) x / sum(x)))
  tr_mat <- replace(tr_mat,is.nan(tr_mat),0) # replace NaN with zero (NaN are produced because variant i is disconnected from network)
  tr_mat <- as(tr_mat, "sparseMatrix")
  
  return(tr_mat)
  print("...Done!")
}


# Extract the main ('biggest') sub-network of a genotype network and modify accordingly the transition matrix.
# This function essentially returns an ergodic transition matrix which facilitiates Markov chain simulations.
extract_main_ntwrk <- function(graph,tr_mat,nodes=FALSE){
  # graph = a genotype network (igraph object)
  # tr_mat = a probability transition matrix
  # nodes = logical value to indicate whether to return the nodes of the main sub-network (filtered matrix)
  
  subnetworks <- components(graph,mode="weak")
  main_net <- which(subnetworks$csize==max(subnetworks$csize)) # index of the biggest sub-network
  genotypes_main <- names(subnetworks$membership[which(subnetworks$membership==main_net)]) # genotypes of the main sub-network
  
  if(nodes) res <- genotypes_main
  else{
    # New (ergodic) probability transition matrix
    res <- tr_mat[rownames(tr_mat) %in% genotypes_main, colnames(tr_mat) %in% genotypes_main]
  }
  return(res)
}

# Finding the stationary distribution of a discrete Markov chain
stationary_dist <- function(tr_mat){
  # tr_mat: A transition probability matrix
  
  # Check whether 'tr_mat' is an ergodic transition probability matrix
  if(!((is.matrix(tr_mat) || inherits(tr_mat,"dgCMatrix")) && nrow(tr_mat)==ncol(tr_mat) && 
       all(tr_mat>=0 & tr_mat<=1) && sum(apply(tr_mat,1,sum))==ncol(tr_mat))){
    stop("Provide a square probability transition matrix")
  }
  
  # Solve the system of linear equations of the form A*pi = b (QR decomposition), where pi is the (transposed) vector of the stationary distribution
  p <- diag(nrow(tr_mat)) - tr_mat
  A <- rbind(t(p),rep(1, ncol(tr_mat)))
  b <- c(rep(0, nrow(tr_mat)),1)
  res <- qr.solve(A, b)
  names(res) <- rownames(tr_mat)
  return(res)
}

# Run a discrete Markov chain. Returns a data frame with the distribution of states after N steps 
simulate_markov_chain <- function(states_0,tr_mat,n_steps,freqs_states_0=NULL){
  # states_0 = vector containing the states from which markov chain starts
  # tr_mat = a probability transition matrix
  # n_steps = number of (mutational) steps to run the Markov chain
  # freq_0 = optional argument to indicate the vector of initial frequencies of states_0 at n_steps=0 (default: NULL)
  
  # Check for inconsistent/valid paramenters
  # Check whether 'tr_mat' is an ergodic transition probability matrix
  if(!((is.matrix(tr_mat) || inherits(tr_mat,"dgCMatrix")) && nrow(tr_mat)==ncol(tr_mat) && 
       all(tr_mat>=0 & tr_mat<=1) && sum(apply(tr_mat,1,sum))==ncol(tr_mat))){
    stop("Provide a square probability transition matrix")
  }
  if(is.null(n_steps) || n_steps==0){
    stop("Provide a correct number of steps to run the markov chain (n_steps > 0)")
  }
  
  # Check whether any P(i,j) in the transition matrix should be zero
  # During sparse matrix transformation, zeroes were coerced to '.Machine$double.xmin' 
  if(any(tr_mat == .Machine$double.xmin)){
    tr_mat[tr_mat == .Machine$double.xmin] <- 0
  }
  
  # Possible states of the markov chain
  states <- rownames(tr_mat)
  
  if(is.null(freqs_states_0)){
    # Build vector of state frequencies at n_steps=0 (all states in 'states_0' have equal initial frequency)
    freq_0 <- rep(0,length(states)); names(freq_0) <- states
    freq_0[names(freq_0) %in% states_0] <- 1/length(states_0)
  }
  else{
    # User defined initial frequencies for states in states_0
    freq_0 <- rep(0,length(states)); names(freq_0) <- states
    freq_0[names(freq_0) %in% states_0] <- freqs_states_0
  }
  
  # operator to compute the k-th power of a sparse matrix (equivalent to `%^%` from 'expm' package)
  `%^^%` = function(x, k) Reduce(`%*%`, rep(list(x), k))
  
  # Extract the frequency distribution of states after n_steps
  freq_n <- as.vector(freq_0 %*% (tr_mat %^^% n_steps))
  names(freq_n) <- states
  return(freq_n)
}

# Compute probability distribution of functional variation (PDFV) accessible to a focal node.
get_PDFV_v2 <- function(state_freqs=NULL,type=NULL,Bg=NULL,model=NA,specific=FALSE){
  # state_freqs = a vector of state frequencies after N steps of a discrete Markov chain (output of 'simulate_markov_chain' function) (default: NULL)
  # type = a string indicating how the PDFV should be computed:
      # "network" = PDFV is proportional to the number of genotypes in the entire network (default)
      # "local network" = PDFV is proportional to the number of *accessible* genotypes
      # "simulated mc" = PDFV is computed according to the probabilities from the 'state_freqs' vector.
  # Bg = a string indicating the DBD background ("AncSR1" or "AncSR2")
  # model = a string indicating the name of model (default: NA)
  # specific = logical argument to indicate whether to compute PDFV using specific genotypes (TRUE) or including promiscuous (FALSE) (default: FALSE).
  
  # Check for inconsistent parameters
  if(!(Bg %in% c("AncSR1","AncSR2")) || is.null(Bg)){
    stop("Specify DBD background: 'AncSR1' or 'AncSR2'")
  }
  if(is.null(type) || !(type %in% c("network", "local network","simulated mc"))){
    stop("Select the type of PDFV to compute: 'network', 'local network', or 'simulated mc'")
  }
  if(is.null(state_freqs) && (type=="local network" || type=="simulated mc")){
    stop("Provide a vector of state frequencies after N steps of a discrete Markov chain")
  }
  if(type=="network" && !is.null(state_freqs)){
    message("Ignoring vector of frequencies...")
  }
  
  if(type=="network"){
    # Compute the expected PDFV as the fraction of genotypes encoding each function in the network.
    if(specific){
      # Use only specific genotypes. Make 'promiscuous genotypes' a separate 'phenotype'
      phenotypes_tbl %>% ungroup() %>% filter(bg==Bg) %>%
        reframe(RE =  REs[[2]],
                count = table(factor(specificity,levels=REs[[2]])),
                Norm_F_prob = count/sum(count),
                model=model) %>% select(RE,Norm_F_prob,model)
    }
    else{
      phenotypes_tbl %>% ungroup() %>% filter(bg==Bg) %>% 
        unnest(cols = c(bound_REs)) %>% 
        reframe(RE =  REs[[1]],
                count = table(factor(bound_REs,levels=REs[[1]])),
                Norm_F_prob = count/sum(count),
                model=model) %>% select(RE,Norm_F_prob,model)
    }
  }
  else if(type=="local network"){
    # Compute the expected PDFV as the fraction of accessible genotypes encoding each function.
    
    # Accessible RH states after N steps of a markov chain
    acc_states <- state_freqs[state_freqs>0]
    
    if(specific){
      # Use only specific genotypes. Make 'promiscuous genotypes' a separate 'phenotype'
      data.frame(AA_var=names(acc_states),bg = Bg) %>% inner_join(.,phenotypes_tbl,by=c("AA_var","bg")) %>%
        reframe(RE =  REs[[2]],
                count = table(factor(specificity,levels=REs[[2]])),
                Norm_F_prob = count/sum(count),
                model = "local GPmap") %>% select(RE,Norm_F_prob,model)
    }
    else{
      # Combine specific and promiscuous genotypes 
      data.frame(AA_var=names(acc_states),bg = Bg) %>% inner_join(.,phenotypes_tbl,by=c("AA_var","bg")) %>%
        unnest(cols = c(bound_REs)) %>% 
        reframe(RE =  REs[[1]],
                count = table(factor(bound_REs,levels=REs[[1]])),
                Norm_F_prob = count/sum(count),
                model = "local GPmap") %>% select(RE,Norm_F_prob,model)
    }
  }
  else if(type=="simulated mc"){
    # Compute the PDFV based on the Markov chain simulation
    
    # Accessible RH states after N steps of a markov chain
    acc_states <- state_freqs[state_freqs>0]
    
    if(specific){
      # Use only specific genotypes. Make 'promiscuous genotypes' a separate 'phenotype'
      data.frame(AA_var=names(acc_states),prob=acc_states,bg = Bg) %>% inner_join(.,phenotypes_tbl,by=c("AA_var","bg")) %>%
        # Add probabilitieis of trajectories ending in the DNA binding phenotype, and normalize probabilities.
        # Report final normalized probability for all 16 phenotypes.
        group_by(specificity) %>% reframe(F_prob = sum(prob)) %>% mutate(Norm_F_prob = F_prob/sum(F_prob)) %>%
        add_row(specificity=REs[[2]],Norm_F_prob=0) %>% distinct(specificity, .keep_all = TRUE) %>% 
        mutate(model = model, specificity = factor(specificity,levels(REs[[2]]))) %>% dplyr::rename(RE = specificity) %>%
        select(RE,Norm_F_prob,model)
    }
    else{
      # Combine specific and promiscuous genotypes 
      data.frame(AA_var=names(acc_states),prob=acc_states,bg = Bg) %>% inner_join(.,phenotypes_tbl,by=c("AA_var","bg")) %>%
        unnest(cols = c(bound_REs)) %>% 
        # Add probabilitieis of trajectories ending in the DNA binding phenotype, and normalize probabilities.
        # Report final normalized probability for all 16 phenotypes.
        group_by(bound_REs) %>% reframe(F_prob = sum(prob)) %>% mutate(Norm_F_prob = F_prob/sum(F_prob)) %>%
        add_row(bound_REs=REs[[1]],Norm_F_prob=0) %>% distinct(bound_REs, .keep_all = TRUE) %>% 
        mutate(model = model, bound_REs = factor(bound_REs,levels(REs[[1]]))) %>% dplyr::rename(RE = bound_REs) %>%
        select(RE,Norm_F_prob,model)
    }
  }
}

# Simulate a multi-step markov chain: Returns a list of data frames with the prob. distribution of functions at each time step
simulate_markov_chain_multistep <- function(states_0,tr_mat,n_iter,Bg=NULL,freqs_states_0=NULL,specific=FALSE){
  # states_0 = vector containing the states from which markov chain starts
  # tr_mat = a probability transition matrix
  # n_iter = number of (mutational) steps to run the Markov chain
  # Bg = a string indicating the DBD background ("AncSR1" or "AncSR2")
  # freq_0 = optional argument to indicate the vector of initial frequencies of states_0 at n_steps=0 (default: NULL)
  # specific = logical argument to indicate whether to compute PDFV using specific genotypes (TRUE) or including promiscuous (FALSE) (default: FALSE).
  
  # Check for valid paramenters
  if(!(Bg %in% c("AncSR1","AncSR2")) || is.null(Bg)){
    stop("Specify DBD background: 'AncSR1' or 'AncSR2'")
  }
  if(is.null(n_iter) || n_iter==0){
    stop("Provide a correct number of steps to run the markov chain (n_iter > 0)")
  }
  
  df <- list()
  for(i in 1:n_iter){
    # simulate Markov chain for each GP map 
    mc_tmp <- simulate_markov_chain(states_0,tr_mat,n_steps = i)
    pdfv_tmp <- get_PDFV_v2(mc_tmp,Bg = Bg,model = i,specific = specific,type="simulated mc")
    df[[i]] <- pdfv_tmp
  }
  return(df)
}

# returns a data frame indicating the number and fraction of specific and promiscuous genotypes for a mc simulation
genotype_type_mc <- function(state_freqs=NULL,Bg=NULL,model=NA){
  # state_freqs = a vector of state frequencies after N steps of a discrete Markov chain (output of 'simulate_markov_chain' function) (default: NULL)
  # Bg = a string indicating the DBD background ("AncSR1" or "AncSR2")
  # model = a string indicating the name of model (default: NA)

  # Accessible RH states after N steps of a markov chain
  acc_states <- state_freqs[state_freqs>0]
  
  data.frame(AA_var=names(acc_states),prob=acc_states,bg = Bg) %>% inner_join(.,phenotypes_tbl,by=c("AA_var","bg")) %>%
    reframe(type=c("Specific","Promiscuous"),count = table(factor(specific,levels=c("YES","NO"))), 
            fraction = count/sum(count),model=model)
}

# Function to find all the paths of a given length between two nodes
find_paths_of_length_N = function(graph, node1, node2, pathlength) {
  # graph = a network (igraph) object
  # nodeX = a string or integer indicating the node 
  # pathlength = integer indicating the hamming distance 
  
  # Initialize an empty list to store the simple paths
  SP = list()
  # If the path length is 1, check if node2 is a neighbor of node1
  if(pathlength == 1) {
    if(node2 %in% neighbors(graph, node1)) {
      # If node2 is a neighbor of node1, add the path to SP
      SP[[1]] = c(node1, node2)
    }
    # Return the list of simple paths
    return(SP)
  }
  # Find the neighbors of node2 in the graph
  Nbrs2 = neighbors(graph, node2)
  # Loop over each neighbor of node2
  for(nbr2 in Nbrs2) {
    # Recursively call the function with node1, the neighbor, and reduced path length
    ASPn2 = find_paths_of_length_N(graph, node1, nbr2, pathlength-1)
    # Loop over each simple path from node1 to the neighbor
    for(i in seq_along(ASPn2)) {
      # Check if node2 is not already in the path
      if(!(node2 %in% ASPn2[[i]])) {
        # Append the path with node2 and add it to SP
        SP = append(SP, list(c(ASPn2[[i]], node2)))
      }
    }
  }
  # Return the list of simple paths
  return(SP)
}

# Function to extract ALL paths of length N to every possible node from a focal node
extract_network_paths <- function(graph,from_i,path_length,cores=1){
  # graph = a genotype netwotk (igraph object)
  # from_i = focal genotype
  # path_length = integer with the neighborhood size
  # cores = number of cores for parallel processing (default: 1)

  # Chech whether node has single mutaant neighbors:
  if(inherits(try(neighbors(graph,from_i)), "try-error")){
    stop("Genotype has no single-mutant neighbors")
  }
  
  # graph = An igraph object (the genotype network)
  # from_i = The starting node around which the paths of length 'l' will be extracted
  # path_length = Integer giving the order of the neighborhood size = hamming distance
  # cores = number of cores for parallel processing (default: 1)
  
  cl <- parallel::makeCluster(cores,"FORK",outfile="")
  future::plan("cluster",workers=cl)
  
  # Compute a matrix containing the hamming distance between every pair of nodes
  distance_tbl = distances(graph,weights = NA)
  
  # extract variants that are at a hamming distance == 'path length' from the focal node
  v <- distance_tbl[which(rownames(distance_tbl)==from_i),]
  v <- v[v==path_length]
  
  # find all the paths of length 'path_length' between focal node and all nodes within hamming distance
  paths <- future_map(names(v),find_paths_of_length_N,graph=graph,node1=from_i,pathlength=path_length)
  paths <- unlist(paths,recursive=FALSE)
  paths <- as.data.frame(do.call(rbind, paths))
  colnames(paths) <- do.call(paste,expand.grid("step",seq(0,path_length,1))) %>% gsub(" ","",.)
  
  parallel::stopCluster(cl) # shut down cluster
  
  # replace vertex numbers with AA variant names
  for(i in 2:path_length){
    nodes <- V(graph)$name[as.integer(paths[,i])]
    paths[,i] <- nodes
  }
  
  # Return a data frame of dimensions (N_paths,path_length+1) indicating all paths from focal node to every accessible node
  return(paths)
}

# Extract the P(i,j) value from a given transition matrix
extract_step_probability <- function(var1,var2,mat){
  i <- which(rownames(mat)==var1)
  j <- which(colnames(mat)==var2)
  r <- mat[i,j]
  #Check whether P(i,j) should be zero
  if(r == .Machine$double.xmin) r <- 0
  return(r)
}

# Function to plot the PDFV as a circular ('spider') chart.
circular_PDFV_v2 <- function(data,cols = "#00AFBB",title = NULL,legend=T,fill=T,...){
  # data: a single data frame or a list of data frames to plot
  
  if(!inherits(data,"list")){
    # Fix and reorganize a single data frame
    radar_df <- data %>% spread(RE,Norm_F_prob) %>% mutate_at(vars(-model), as.numeric) %>% 
      mutate_at(vars(model),as.factor)
    max <- radar_df %>% pivot_longer(cols=2:last_col()) %>% with(max(value))
  }
  else{
    # Fix and reorganize multiple data frames
    radar_df <- do.call(rbind, data)
    radar_df <- radar_df %>% group_by(model) %>% 
      spread(RE,Norm_F_prob) %>% mutate_at(vars(-model), as.numeric) %>% 
      mutate_at(vars(model),as.factor)
    max <- radar_df %>% pivot_longer(cols=2:last_col()) %>% with(max(value))
  }
  
  # Plot
  vals <- c(0,max/2,max)
  vals <- round(vals,digits = 2)
  if(is.null(title)) title <- ""
  
  ggradar(radar_df,group.point.size = 1.7,fill = fill, group.colours= cols,
          values.radar = vals, legend.position="left",
          gridline.max.colour="black",group.line.width=0.7,axis.label.size = 3,grid.label.size=3.5,
          gridline.mid.colour="gray",grid.max = max, grid.mid = vals[2],plot.title=title,plot.legend=legend,
          legend.text.size=8,...) + 
    theme(legend.key.width = unit(0.4, 'cm'), plot.title = element_text(size=15))
}

# Compute the probability of all pairwise phenotypic transitions given a number of mutation steps.
phenotypic_transitions <- function(from=REs[[1]],to=REs[[1]],from_nodes=NULL,tr_mat,bg,n_steps,specific=F,normalize=T){
  # from, to = vectors containing the names of the DNA phenotypes to compute the trajectory (default: all REs)
  # from_nodes = alternative vector to 'from' which, instead of DNA phenotypes, includes the vector of starting genotypes.
  # tr_mat = a probability transition matrix
  # Bg = a string indicating the DBD background ("AncSR1" or "AncSR2")
  # n_steps = number of (mutational) steps to run the Markov chain
  # specific = logical argument to indicate whether to compute PDFV using specific genotypes (TRUE) or including promiscuous (FALSE) (default: FALSE).
  # normalize = logical arg. to specify whether to re-normalize the probabilities (default: TRUE)
  
  # check inconsistent/valid arguments
  if(!(bg %in% c("AncSR1","AncSR2")) || is.null(bg)){
    stop("Specify DBD background: 'AncSR1' or 'AncSR2'")
  }
  if(!is.null(from_nodes) && !is.null(from)){
    stop("Provide only one set of starting parameters; either a vector of phenotypes, or a vector of genotypes.")
  }
  if(is.null(from) && is.null(to)){
    stop("Empty 'from' and 'to' arguments. Provide the vector of phenotypes to compute the transition probabilities.")
  }
  if(is.null(tr_mat)){
    stop("Provide a probability transition matrix.")
  }
  if(is.null(n_steps) || n_steps==0){
    stop("Provide a correct number of steps to run the markov chain (n_steps > 0)")
  }
  
  if(!is.null(from_nodes)){
    # dimensions of matrix and names of columns/rows
    # extract starting phenotypes from starting nodes
    RE_from <- unique(phenotypes_tbl %>% filter(specific == "YES" & bg == bg) %>%
                        filter(AA_var %in% from_nodes) %>% pull(specificity))
    n_rows <- length(RE_from)
    
    n_cols <- length(to)
    RE_to <- REs[[1]][REs[[1]] %in% to]
  }
  else{
    # dimensions of matrix and names of columns/rows
    n_rows <- length(from)
    RE_from <- REs[[1]][REs[[1]] %in% from]
    
    n_cols <- length(to)
    RE_to <- REs[[1]][REs[[1]] %in% to]
  }
  
  # order rows and columns
  RE_from <- RE_from[order(match(RE_from,REs[[1]]))]
  RE_to <- RE_to[order(match(RE_to,REs[[1]]))]

  # create phenotypic transition matrix between all pairs of phenotypes
  pheno_transition_mat <- matrix(NA,ncol = n_cols, nrow = n_rows)
  rownames(pheno_transition_mat) <- RE_from
  colnames(pheno_transition_mat) <- RE_to
  states <- rownames(tr_mat) # states in the transition matrix
  
  for(i in RE_from){
    # extract amino acid variants from the ith neutral network 
    phenotype_vars <- phenotypes_tbl %>% filter(AA_var %in% states) %>%
      filter(bg == bg & specificity == i) %>% pull(AA_var)
    
    # if specific genotypes are not part of the main network, or no specific genotypes, continue:
    if(identical(phenotype_vars, character(0))) next
    
    # run markov chain
    mc_tmp <- simulate_markov_chain(phenotype_vars,tr_mat,n_steps = n_steps)
    pdfv_tmp <- get_PDFV_v2(mc_tmp,Bg = bg ,specific = specific,type="simulated mc") %>% filter(RE %in% RE_to)
    
    # re-normalize probabilities?
    if(normalize){
      pdfv_tmp <- pdfv_tmp %>% mutate(Norm_F_prob = Norm_F_prob /sum(Norm_F_prob))
    }
    # if excluding promiscuous, recalculate the probabilities of each transition
    if(specific && normalize){
      pdfv_tmp <- pdfv_tmp %>% filter(RE != "Promiscuous") %>% mutate(Norm_F_prob = Norm_F_prob /sum(Norm_F_prob))
    }
    
    # order REs
    x <- pdfv_tmp$Norm_F_prob; names(x) <- pdfv_tmp$RE
    x <- x[order(match(names(x),REs[[1]]))]
    pheno_transition_mat[which(rownames(pheno_transition_mat)== i),] <- x
  }
  return(pheno_transition_mat)
}

