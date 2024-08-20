#===============================================
# Evolutionary simulations on empirical GP maps
#===============================================
# Last modified: May 02, 2024 - SHA
#
# Functions to build, simulate and work with discrete Markov processes on RH-RE DMS data



# make factors with different levels for ordering of REs in plotting
REs <- list()
REs[[1]] <- factor(c("SRE (AA)", "GA", "ERE (GT)", "AC", 
                     "AG", "AT", "CA", "CC", 
                     "CG", "CT", "GC", "GG", 
                     "TA", "TC", "TG", "TT"), 
                   levels=c("ERE (GT)", "SRE (AA)", "AC", "AG", 
                            "AT", "CA", "CC", "CG", 
                            "CT", "GA", "GC", "GG", 
                            "TA", "TC", "TG", "TT"))
REs[[2]] <- factor(c("SRE (AA)", "GA", "ERE (GT)", "AC", 
                     "AG", "AT", "CA", "CC", 
                     "CG", "CT", "GC", "GG", 
                     "TA", "TC", "TG", "TT","Promiscuous"), 
                   levels=c("ERE (GT)", "SRE (AA)", "AC", "AG", 
                            "AT", "CA", "CC", "CG", 
                            "CT", "GA", "GC", "GG", 
                            "TA", "TC", "TG", "TT","Promiscuous"))

# Assign colors to each RE for plotting (and matching the colors in the networks)
hex_RE_colors <- function(type = 1){
  if(type==1){
    col <- c("#24ff24","#ffff6d","#490092","#ffb6db","#db6d00",
      "#009292","#006ddb","#594f15", "#ff6db6","#004949","#b66dff",
      "#924900","#920000","#c0c0c0","#b6dbff","#6db6ff")
    names(col) <- REs[[1]]
  }
  else if(type==2){
    col <- c("#24ff24","#ffff6d","#490092","#ffb6db","#db6d00",
      "#009292","#006ddb","#594f15", "#ff6db6","#004949","#b66dff",
      "#924900","#920000","#c0c0c0","#b6dbff","#6db6ff","white")
    names(col) <- REs[[2]]
  }
  else if(type==3){
    col <- c("#24ff24","gray60","#490092","gray60",
      "gray60","gray60","gray60","gray60",
      "gray60","gray60","gray60","gray60",
      "gray60","gray60","gray60","gray60")
    names(col) <- REs[[1]]
  }
  return(col)
}

# Assign colors to DBD backgrounds
DBD_color <- function(names=T){
  colors <- c("#674797","#00AC58")
  if(names) names(colors) <- c("AncSR1","AncSR2")
  return(colors)
}

# Determine whether two amino acid sequences can be connected by a single nucleotide change given the genetic code
connect_aa_variants <- function(var1,var2){
  # var1, var2 = amino acid variants
  
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
  # var1, var2 = amino acid variants
  
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

# Function to calculate number of nucleotide mutations between 2 amino acids
connect_single_aa_count <- function(aa1,aa2){
  # aa1, aa2 = amino acids to compare
  
  codons_aa1 <- names(GENETIC_CODE[GENETIC_CODE==aa1])
  codons_aa2 <- names(GENETIC_CODE[GENETIC_CODE==aa2])
  # Compare all possible codon pairs: If there is at least one pair that can be interchanged by *a single* nuc change, 
  # then the link exists between AAs
  pairs <- expand.grid(codons_aa1,codons_aa2) %>% rowwise() %>%
    mutate(nuc_mismatches = mapply(function(c1,c2) sum(c1!=c2), strsplit(as.character(Var1),''), strsplit(as.character(Var2),'')))
  # Compute number of codon paths
  n <- sum(pairs$nuc_mismatches == 1)
  return(n)
}

# Determine whether two amino acid sequences can be connected by a single nucleotide change given the genetic code
# The returned value is the total number of mutations by which codons of var1 can access codons of var2
connect_aa_variants_count <- function(var1,var2,syn=FALSE){
  # var1, var2 = amino acid variants
  # syn = whether synonymous mutations (between the same amino acid) are incorporated or not (default: FALSE)
  
  # If aa variants can be connected by single step mutation given genetic code, returns the number of 
  # single nuc. mutations between codons of 'var1' and 'var2'; "0" otherwise.
  
  var1 <- unlist(strsplit(var1,''))
  var2 <- unlist(strsplit(var2,''))
  
  if(length(var1) == 1 || length(var2) == 1){
    stop("For single amino acid comparison use the function `connect_single_aa_count` instead.")
  }
  
  aa_mismatches <- sum(var1 != var2)
  link <- 0
  if(aa_mismatches == 1){ #If not the same variant and not more than 2 aa diffs
    aa_diff <- which(var1 != var2) # position where variants differ
    
    # compute number of nucleotide mutations between amino acids
    n <- connect_single_aa_count(var1[aa_diff],var2[aa_diff])
    
    # we have to multiply the number of codon paths across all invariant amino acid positions to get the total number of “backgrounds” 
    # that each amino acid mutation can occur in.
    bg_codons <- map(var1[-aa_diff],function(x) names(GENETIC_CODE)[GENETIC_CODE == x])
    nbg <- do.call(prod, map(bg_codons, length))
    link <- n * nbg
  }
  else if(aa_mismatches == 0 && syn == TRUE){ #If the same amino acid AND we want to incorporate syn mutations
    # get number of synonymous mutations for each amino acid
    n_mut_aa <- map2_int(var1,var2,connect_single_aa_count)
    
    # Number of synonymous backgrounds per position
    bgs <- map(var1, function(x) names(GENETIC_CODE)[GENETIC_CODE == x])
    nbgaa <- map(bgs, length)
    npaths <- n_mut_aa * sapply(1:length(aaseq1),function(x) do.call(prod, nbgaa[-x]))
    link <- sum(npaths)
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

# Connect genotypes based on their hamming distance (link exists only if hamming distance = 1)
connect_hamming <- function(var1, var2){
  hd <- hamming_d(var1,var2)
  if(hd != 1) link = 0
  else link = hd
  return(link)
}

# Function to connect protein-DNA complex variants. A link exists whether protein genotypes can be be connected by a single 
# nucleotide change given the genetic code, OR DNA variants differ by a hamming distance of one. Variants are given as 
# prot+DNA genotypes (e.g., EGKAGT)
connect_complex <- function(var1,var2){
  # split variants into protein and DNA genotypes.
  var1_prot <- unlist(strsplit(var1,''))[1:4]
  var1_dna <- unlist(strsplit(var1,''))[5:6]
  var2_prot <- unlist(strsplit(var2,''))[1:4]
  var2_dna <- unlist(strsplit(var2,''))[5:6]
  
  link <- 0
  # Compute hamming distance for the complex
  hd_prot <- hamming_d(var1_prot,var2_prot)
  hd_dna <- hamming_d(var1_dna,var2_dna)
  
  # Compute mutation links only if complexes are *one* change appart, either in DNA or protein
  if(hd_dna + hd_prot == 1){
    if(hd_prot == 1){ # if change happened in protein, check whether variants can be mutated
      link <- connect_aa_variants_count(var1_prot,var2_prot)
    }
    else if(hd_dna == 1){
      # If change happended in DNA, find all the protein backgrounds in which mutation at the DNA can occur
      bgs <- map(var1_prot, function(x) names(GENETIC_CODE)[GENETIC_CODE == x])
      link <- do.call(prod, map(bgs, length)) * hd_dna
    }
  }
  return(link)
}

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
    # of codons; 3 - mutational paths; 4 - mutational paths between prot-DNA complexes;
    # 5 - hamming distance (default: 1)
  # cores: specify the number of cores to use in to use for parallel computing (default: 1)
  
  # Check for inconistent arguments
  if(is.null(nodes)){
    stop("Provide a vector of nodes to build the matrix")
  }
  if(!(type %in% seq(1,5,1))){
    stop("Invalid matrix type. Valid matrix options: 1, 2, 3, 4 or 5")
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
    # number of mutational paths from i to j (codon bias)
    adj_mat <- Matrix(foreach(i = 1:length(nodes), .combine = 'cbind') %dopar%
                        mapply(connect_aa_variants_count,nodes,nodes[i]) %>% `colnames<-`(nodes), sparse = T)
  }
  else if(type == 4){
    # number of mutational paths from between protein-DNA complexes
    adj_mat <- Matrix(foreach(i = 1:length(nodes), .combine = 'cbind') %dopar%
                        mapply(connect_complex,nodes,nodes[i]) %>% `colnames<-`(nodes), sparse = T)
  }
  else if(type == 5){
    # hamming distance between genotypes
    adj_mat <- Matrix(foreach(i = 1:length(nodes), .combine = 'cbind') %dopar%
                        mapply(connect_hamming,nodes,nodes[i]) %>% `colnames<-`(nodes), sparse = T)
  }
  stopCluster(cl)
  
  return(adj_mat)
}

# Build a genotype network (igraph object)
build_genotype_network <- function(nodes,build_mat = TRUE,adj_mat=NULL,type=1,node_attributes=FALSE,cores=1){
  # nodes = a vector with vertex character names or data.frame (see argument 'node_attributes')
  # build_mat = logical value to indicate whether to generate adjacency matrix (default: TRUE)
  # adj_mat = if build_mat=FALSE, then pass the corresponding adjacency matrix to generate the network (default: NULL)
  # type: specify how to build the adjacency matrix. See options for 'build_mutation_matrix' (default: 1)
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
  if(!(type %in% seq(1,5,1))){
    stop("Invalid matrix type. Valid matrix options: 1, 2, 3, 4 or 5")
  }
  
  if(build_mat==FALSE){
    # Build networks using supplied adjacency matrix
    if(type != 1 && type != 5){
      net <- graph_from_adjacency_matrix(adj_mat,mode="directed",weighted=TRUE)
    }
    else if(type == 1 || type == 5){
      if(node_attributes){
        links <- graph.adjacency(adj_mat,weighted = "1")
        edges <- get.edgelist(links) # generate edgelist from matrix
        net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) # network with annotations per node
      }
      else net <- graph_from_adjacency_matrix(adj_mat,mode = "undirected")
    }
  }
  else if(build_mat){
    if(node_attributes){
      nodes_v <- nodes %>% pull(AA_var)
      adj_mat <- build_mutation_matrix(nodes_v,type=type,cores=cores)
      links <- graph.adjacency(adj_mat,mode = "undirected")
      edges <- get.edgelist(links) # generate edgelist from matrix
      net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) # network with annotations per node
    }
    else{      
      # Build adjacency matrix and network
      adj_mat <- build_mutation_matrix(nodes,type=type,cores=cores)
      if(type == 1 || type == 5) net <- graph_from_adjacency_matrix(adj_mat,mode = "undirected")
      else if(type != 1 && type != 5) net <- graph_from_adjacency_matrix(adj_mat,weighted="1") %>% as.undirected()
    }
  }  
  return(net)
}    

# Extract a variant's phenotypic value
get_phenotype <- function(aa_var,Bg,pheno_table,phenotype="mean",complex=FALSE){
  # aa_var = variant or variants (can be a single string or a string vector)
    # Either a protein variant or a prot-DNA complex variant
  # Bg = string specifying the DBD background (= 'AncSR1' or 'AncSR2')
  # pheno_table = a data frame with phenotypic values (meanF)
  # phenotype = a string specifying which 'phenotype' to extract: 
    # 'mean' = avg meanF across bound DNA elements (default)
    # 'min' = min meanF across bound DNA elements
    # 'max' = max meanF across bound DNA elements
  # complex = whether to extract phenotypes of prot-DNA complexes (default: FALSE)
  
  # check inconsistent parameters
  if(complex && !("complex" %in% colnames(pheno_table))){
    stop("The phenotypic table provided does not include prot-DNA complexes")
  }

  p <- NA
  if(complex){
    n <- pheno_table %>% filter(complex %in% aa_var & bg == Bg) %>% pull(complex)
    p <- pheno_table %>% filter(complex %in% aa_var & bg == Bg) %>% pull(meanF_bREs)
  }
  else {
    n <-  pheno_table %>% filter(AA_var %in% aa_var & bg == Bg) %>% pull(AA_var)
    if(phenotype == "mean"){ p <- pheno_table %>% filter(AA_var %in% aa_var & bg == Bg) %>% pull(meanF_bREs)}
    if(phenotype == "min"){ p <- pheno_table %>% filter(AA_var %in% aa_var & bg == Bg) %>% pull(min_meanF_bREs)}
    if(phenotype == "max"){ p <- pheno_table %>% filter(AA_var %in% aa_var & bg == Bg) %>% pull(max_meanF_bREs)}
  }
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
  neighbors <- FALSE

  # Determine whether the genotype variants are single-mutation neighbors
  if(length(unlist(strsplit(from_i,''))) == 6){ # => protein-DNA complex
    l <- connect_complex(from_i,to_j)
    if(l > 0) neighbors <- TRUE
  }
  else{ # => protein genotype
    l <- connect_aa_variants(from_i,to_j)
    if(l == 1) neighbors <- TRUE
  }
  
  # If variants i and j are one-mutant neighbors, proceed to compute P(i,j)
  if(neighbors){
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
build_transition_matrix <- function(nodes,graph,pheno_tbl,param,complex=FALSE,cores=1){
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
  # complex: a logical value to indicate whether to extract prot-DNA complex phenotypes (default: FALSE)
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
                             mutation=mutation,phenotype=pheno,complex=complex) %>% 
                     `colnames<-`(nodes) %>% `rownames<-`(nodes))
  stopCluster(cl)
  
  return(tr_mat)
  print("...Done!")
}

# Build transition probability matrix (P(i,j)'s) for a DBD background and a specific selection scenario
# (uses the function 'rate_mutational_step', instead of 'prob_mutational_step')
build_transition_matrix_v2 <- function(nodes,adj_matrix,pheno_tbl,param,complex=FALSE,cores=1){
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
  # complex: a logical value to indicate whether to extract prot-DNA complex phenotypes (default: FALSE)
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
                             mutation=mutation,phenotype=pheno,complex=complex) %>% 
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
  
  graph = as.undirected(graph)

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
  if(!((is.matrix(tr_mat) || inherits(tr_mat,"dgCMatrix") || inherits(tr_mat,"dsCMatrix")) && 
    nrow(tr_mat)==ncol(tr_mat) && 
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
  if(!((is.matrix(tr_mat) || inherits(tr_mat,"dgCMatrix") || inherits(tr_mat,"dsCMatrix")) && 
    nrow(tr_mat)==ncol(tr_mat) && all(tr_mat>=0 & tr_mat<=1) && sum(apply(tr_mat,1,sum))==ncol(tr_mat))){
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

# Compute probability distribution of functional variation (PDFV) accessible to a focal node(s).
get_PDFV_v2 <- function(state_freqs=NULL,type=NULL,Bg=NULL,model=NA,pheno_tbl=phenotypes_tbl,specific=FALSE,complex=FALSE,graph=NULL){
  # state_freqs = a vector of state frequencies after N steps of a discrete Markov chain (output of 'simulate_markov_chain' function) (default: NULL)
  # type = a string indicating how the PDFV should be computed:
      # "network" = PDFV is proportional to the number of genotypes in the entire network
      # "main network" = PDFV is proportional to the number of genotypes in the main component of the network
      # "simulated mc" = PDFV is computed according to the probabilities from the 'state_freqs' vector.
  # Bg = a string indicating the DBD background ("AncSR1" or "AncSR2")
  # model = a string indicating the name of model (default: NA)
  # pheno_tbl = table with annotated phenotypes per variant
  # specific = logical argument to indicate whether to compute PDFV using specific genotypes (TRUE) or including promiscuous (FALSE) (default: FALSE).
  # complex = a logical argument to indicate whether PDFV is to be computed for prot-DNA complexes
  # graph = a genotype network (igraph object)

  # Check for inconsistent parameters
  if(!(Bg %in% c("AncSR1","AncSR2")) || is.null(Bg)){
    stop("Specify DBD background: 'AncSR1' or 'AncSR2'")
  }
  if(is.null(type) || !(type %in% c("production", "main network","simulated mc"))){
    stop("Select the type of PDFV to compute: 'production', 'main network', or 'simulated mc'")
  }
  if(is.null(state_freqs) && (type=="simulated mc")){
    stop("Provide a vector of state frequencies after N steps of a discrete Markov chain")
  }
  if(type=="production" && !is.null(state_freqs)){
    message("Ignoring vector of frequencies...")
  }
  if(complex & !("complex" %in% colnames(pheno_tbl))){
    stop("The phenotypic table does not contain prot-DNA complexes")
  }
  if(type == "main network" && is.null(graph)){
    stop("Need to provide a graph to extract the main sub-component")
  }
  
  if(type=="production"){
    # Compute the expected PDFV as the fraction of protein genotypes encoding each phenotype.
    if(specific){
      # Use only specific genotypes. Make 'promiscuous genotypes' a separate 'phenotype'
      pheno_tbl %>% ungroup() %>% filter(bg==Bg) %>%
        reframe(RE =  REs[[2]],
                count = table(factor(specificity,levels=REs[[2]])),
                Norm_F_prob = count/sum(count),
                model=model) %>% select(RE,Norm_F_prob,model)
    }
    else{
      pheno_tbl %>% ungroup() %>% filter(bg==Bg) %>% 
        unnest(cols = c(bound_REs)) %>% 
        reframe(RE =  REs[[1]],
                count = table(factor(bound_REs,levels=REs[[1]])),
                Norm_F_prob = count/sum(count),
                model=model) %>% select(RE,Norm_F_prob,model)
    }
  }
  else if(type=="main network"){
    # Compute the expected PDFV as the fraction of genotypes encoding each function in the main network
    
    # Extract genotypes from main network
    net_vars <- extract_main_ntwrk(graph=graph,tr_mat=NULL,nodes=TRUE)
    
    if(complex){
      # Compute PDFV from prot-DNA complexes
      if(specific){
        # Use only specific genotypes. Make 'promiscuous genotypes' a separate 'phenotype'
        pheno_tbl %>% ungroup() %>% filter(bg == Bg & complex %in% net_vars) %>% 
        reframe(RE =  REs[[2]],
                count = table(factor(specificity,levels=REs[[2]])),
                Norm_F_prob = count/sum(count),
                model=model) %>% select(RE,Norm_F_prob,model)
      }
      else{
      # Combine specific and promiscuous genotypes 
        pheno_tbl %>% ungroup() %>% filter(bg == Bg & complex %in% net_vars) %>% 
        unnest(cols = c(bound_REs)) %>% 
        reframe(RE =  REs[[1]],
                count = table(factor(bound_REs,levels=REs[[1]])),
                Norm_F_prob = count/sum(count),
                model=model) %>% select(RE,Norm_F_prob,model)
      }
    }
    else{
      # Compute PDFV from protein variants
      if(specific){
        # Use only specific genotypes. Make 'promiscuous genotypes' a separate 'phenotype'
        pheno_tbl %>% ungroup() %>% filter(bg == Bg & AA_var %in% net_vars) %>% 
        reframe(RE =  REs[[2]],
                count = table(factor(specificity,levels=REs[[2]])),
                Norm_F_prob = count/sum(count),
                model=model) %>% select(RE,Norm_F_prob,model)
      }
      else{
      # Combine specific and promiscuous genotypes 
        pheno_tbl %>% ungroup() %>% filter(bg == Bg & AA_var %in% net_vars) %>% 
        unnest(cols = c(bound_REs)) %>% 
        reframe(RE =  REs[[1]],
                count = table(factor(bound_REs,levels=REs[[1]])),
                Norm_F_prob = count/sum(count),
                model=model) %>% select(RE,Norm_F_prob,model)
      }
    }
  }
  else if(type=="simulated mc"){
    # Compute the PDFV based on the Markov chain simulation
    
    # Accessible RH states after N steps of a markov chain
    acc_states <- state_freqs[state_freqs>0]
    
    if(complex){
      # Compute PDFV from prot-DNA complexes
      if(specific){
        # Use only specific genotypes. Make 'promiscuous genotypes' a separate 'phenotype'
        data.frame(complex=names(acc_states),prob=acc_states,bg = Bg) %>% inner_join(.,pheno_tbl,by=c("complex","bg")) %>%
          # Add probabilitieis of trajectories ending in the DNA binding phenotype, and normalize probabilities.
          # Report final normalized probability for all 16 phenotypes.
          group_by(specificity) %>% reframe(F_prob = sum(prob)) %>% mutate(Norm_F_prob = F_prob/sum(F_prob)) %>%
          add_row(specificity=REs[[2]],Norm_F_prob=0) %>% distinct(specificity, .keep_all = TRUE) %>% 
          mutate(model = model, specificity = factor(specificity,levels(REs[[2]]))) %>% dplyr::rename(RE = specificity) %>%
          select(RE,Norm_F_prob,model)
      }
      else{
        # Combine specific and promiscuous genotypes 
        data.frame(complex=names(acc_states),prob=acc_states,bg = Bg) %>% inner_join(.,pheno_tbl,by=c("complex","bg")) %>%
          unnest(cols = c(bound_REs)) %>% 
          # Add probabilitieis of trajectories ending in the DNA binding phenotype, and normalize probabilities.
          # Report final normalized probability for all 16 phenotypes.
          group_by(bound_REs) %>% reframe(F_prob = sum(prob)) %>% mutate(Norm_F_prob = F_prob/sum(F_prob)) %>%
          add_row(bound_REs=REs[[1]],Norm_F_prob=0) %>% distinct(bound_REs, .keep_all = TRUE) %>% 
          mutate(model = model, bound_REs = factor(bound_REs,levels(REs[[1]]))) %>% dplyr::rename(RE = bound_REs) %>%
          select(RE,Norm_F_prob,model)
      }
    }
    else{
      if(specific){
        # Use only specific genotypes. Make 'promiscuous genotypes' a separate 'phenotype'
        data.frame(AA_var=names(acc_states),prob=acc_states,bg = Bg) %>% inner_join(.,pheno_tbl,by=c("AA_var","bg")) %>%
          # Add probabilitieis of trajectories ending in the DNA binding phenotype, and normalize probabilities.
          # Report final normalized probability for all 16 phenotypes.
          group_by(specificity) %>% reframe(F_prob = sum(prob)) %>% mutate(Norm_F_prob = F_prob/sum(F_prob)) %>%
          add_row(specificity=REs[[2]],Norm_F_prob=0) %>% distinct(specificity, .keep_all = TRUE) %>% 
          mutate(model = model, specificity = factor(specificity,levels(REs[[2]]))) %>% dplyr::rename(RE = specificity) %>%
          select(RE,Norm_F_prob,model)
      }
      else{
        # Combine specific and promiscuous genotypes 
        data.frame(AA_var=names(acc_states),prob=acc_states,bg = Bg) %>% inner_join(.,pheno_tbl,by=c("AA_var","bg")) %>%
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
}

# Simulate a multi-step markov chain: Returns a list of data frames with the prob. distribution of functions at each time step
simulate_markov_chain_multistep <- function(states_0,tr_mat,n_iter,Bg=NULL,freqs_states_0=NULL,specific=FALSE,complex=FALSE,pheno_tbl=phenotypes_tbl){
  # states_0 = vector containing the states from which markov chain starts
  # tr_mat = a probability transition matrix
  # n_iter = number of (mutational) steps to run the Markov chain
  # Bg = a string indicating the DBD background ("AncSR1" or "AncSR2")
  # freq_0 = optional argument to indicate the vector of initial frequencies of states_0 at n_steps=0 (default: NULL)
  # specific = logical argument to indicate whether to compute PDFV using specific genotypes (TRUE) or including promiscuous (FALSE) (default: FALSE).
  # complex = whether to compute PDFV for protein-DNA complexes (default: FALSE)

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
    pdfv_tmp <- get_PDFV_v2(mc_tmp,Bg = Bg,model = i,specific = specific,type="simulated mc",complex=complex,pheno_tbl = pheno_tbl)
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
circular_PDFV_v2 <- function(data,cols = "#00AFBB",title = NULL,legend=T,fill=T,point_size=1.7,line_width=0.7,leg_position="left",...){
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
      mutate_at(vars(model),as.factor) %>% mutate(model = factor(model,levels = unique(.$model),ordered=T))
    max <- radar_df %>% pivot_longer(cols=2:last_col()) %>% with(max(value))
  }
  
  # Plot
  vals <- c(0,max/2,max)
  vals <- round(vals,digits = 2)
  if(is.null(title)) title <- ""
  
  ggradar(radar_df,group.point.size = point_size,fill = fill, group.colours= cols,
          values.radar = vals, legend.position=leg_position,
          gridline.max.colour="black",group.line.width=line_width,axis.label.size = 3,grid.label.size=3.5,
          gridline.mid.colour="gray",grid.max = max, grid.mid = vals[2],plot.title=title,plot.legend=legend,
          legend.text.size=8,...) + 
    theme(legend.key.width = unit(0.4, 'cm'), plot.title = element_text(size=15)) +
    guides(color = guide_legend(override.aes = list(size = 5)))
}

# Function to plot the PDFV as a barplot.
bars_PDFV <- function(data,fill_by_RE = FALSE, fill_cols = "#00AFBB",col="black",title = NULL,legend=T,arrange_by=NULL,hex_id=NULL){
  # data: a single data frame or a list of data frames to plot
  
  if(is.null(title)) title <- ""
  
  if(inherits(data,"list")){
    data <- do.call(rbind, data)
  }
  if(fill_by_RE){
    if(!is.null(arrange_by)){
      # create arranged levels based on selected model label
      ord_levs <- as.character(data %>% filter(model == arrange_by) %>% arrange(desc(Norm_F_prob)) %>% pull(RE))
      
      p <- data %>% mutate(RE = factor(RE,levels = ord_levs),
                           model = factor(model,levels = unique(sort(.$model)))) %>% 
        ggplot(aes(x=RE,y=Norm_F_prob,fill=RE)) + 
        geom_bar(stat="identity",position="dodge",color=col) +
        scale_fill_manual(values = hex_RE_colors(hex_id))
    }
    else{
      p <- data %>% mutate(RE = factor(RE,levels = REs[[1]]),
                           model = factor(model,levels = unique(sort(.$model)))) %>% 
        ggplot(aes(x=RE,y=Norm_F_prob,fill=RE)) + 
        geom_bar(stat="identity",position="dodge",color=col) +
        scale_fill_manual(values = hex_RE_colors(hex_id))
    }
  }
  else{
    if(!is.null(arrange_by)){
      # create arranged levels based on selected model label
      ord_levs <- as.character(data %>% filter(model == arrange_by) %>% arrange(desc(Norm_F_prob)) %>% pull(RE))
      
      p <- data %>% mutate(RE = factor(RE,levels = ord_levs),
                           model = factor(model,levels = unique(sort(.$model)))) %>% 
        ggplot(aes(x=RE,y=Norm_F_prob,fill=model)) + 
        geom_bar(stat="identity",position="dodge",color=col) +
        scale_fill_manual(values = fill_cols)
    }
    else{
      p <- data %>% mutate(RE = factor(RE,levels = REs[[1]]),
                           model = factor(model,levels = unique(sort(.$model)))) %>% 
        ggplot(aes(x=RE,y=Norm_F_prob,fill=model)) + 
        geom_bar(stat="identity",position="dodge",color=col) +
        scale_fill_manual(values = fill_cols)
    }
  }
  
  p <- p + theme_classic() + ggtitle(title) +
    labs(x="DNA element",y="Probability") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "none") + 
    facet_grid(rows = vars(model))
  
  if(legend == FALSE) p <- p + theme(legend.position = "none")
  else if(legend == "bottom") p <- p + theme(legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 5)))
  
  return(p)
}

# Compute the proximity between two neutral networks
pairwise_neutral_network_proximity <- function(n_ntwrk1,n_ntwrk2,Bg,graph,type=1,pheno_df=phenotypes_tbl,complex=FALSE,from_specific=T,to_specific=T){
  # n_ntwrk1,n_ntwrk2 = strings specifying the names of RE functions
  # Bg = string specifying the DBD background (= 'AncSR1' or 'AncSR2')
  # graph = a genotype network (igraph object)
  # type = proximity calculation:
    # 1 - fraction overlap due to promiscuous genotypes; 2 - number of direct links betwen neutral networks; 3 - average shortest path length between neutral networks
  # pheno_df = df with phenotypic annotations per genotype
  # complex = whether to use protein-DNA complex genoytpes (default: FALSE)
  # from_specific, to_specific = whether to select specific genotypes 
  
  # Check inconsistent parameters
  if(!(type %in% c(1,2,3))){
    stop("Provide a valid value for type: 1, 2 or 3")
  }
  if(!(Bg %in% c("AncSR1","AncSR2")) || is.null(Bg)){
    stop("Specify DBD background: 'AncSR1' or 'AncSR2'")
  }
  if(complex && type == 1){
    stop("protein-DNA complex networks don't have explicit promiscuous genotypes")
  }

  prox <- NA

  if(type == 1){
    # Compute fraction overlap due to promiscuous genotypes: AnB/AuB

    # Extract genotypes from each neutral network: all binders
    all_vars_ntwrk1 <- unique(pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == n_ntwrk1 & bg == Bg) %>% pull(AA_var))
    all_vars_ntwrk2 <- unique(pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == n_ntwrk2 & bg == Bg) %>% pull(AA_var))

    # Filter variants to retain those present in main component genotype network
    net_vars <- extract_main_ntwrk(graph=graph,tr_mat=NULL,nodes=TRUE)
    all_vars_ntwrk1 <- all_vars_ntwrk1[all_vars_ntwrk1 %in% net_vars]
    all_vars_ntwrk2 <- all_vars_ntwrk2[all_vars_ntwrk2 %in% net_vars]

    # Check that neutral networks are encoded by at least one genotype 
    if(!identical(all_vars_ntwrk1, character(0)) && !identical(all_vars_ntwrk2, character(0))){
      # Compute fraction overlap due to promiscuous genotypes: AnB/AuB
      U <- union(all_vars_ntwrk1,all_vars_ntwrk2)
      I <- intersect(all_vars_ntwrk1,all_vars_ntwrk2)
      prox <- length(I) / length(U)
    }
  }
  else if(type == 2){
    # Compute number of direct links between specific genotypes

    # Extract genotypes from each neutral network: specific binders
    if(complex){
      if(from_specific) vars_ntwrk1 <- pheno_df %>% filter(specificity == n_ntwrk1 & bg == Bg) %>% pull(complex) # specific
      else vars_ntwrk1 <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == n_ntwrk1 & bg == Bg) %>% pull(complex) # binding

      if(to_specific) vars_ntwrk2 <- pheno_df %>% filter(specificity == n_ntwrk2 & bg == Bg) %>% pull(complex) # specific
      else vars_ntwrk2 <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == n_ntwrk2 & bg == Bg) %>% pull(complex) # binding
    }
    else{
      if(from_specific) vars_ntwrk1 <- pheno_df %>% filter(specificity == n_ntwrk1 & bg == Bg) %>% pull(AA_var) # specific
      else vars_ntwrk1 <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == n_ntwrk1 & bg == Bg) %>% pull(AA_var) # binding
      
      if(to_specific) vars_ntwrk2 <- pheno_df %>% filter(specificity == n_ntwrk2 & bg == Bg) %>% pull(AA_var) # specific
      else vars_ntwrk2 <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == n_ntwrk2 & bg == Bg) %>% pull(AA_var) # binding
    }

    # Filter variants to retain those present in main component genotype network
    net_vars <- extract_main_ntwrk(graph=graph,tr_mat=NULL,nodes=TRUE)
    vars_ntwrk1 <- vars_ntwrk1[vars_ntwrk1 %in% net_vars]
    vars_ntwrk2 <- vars_ntwrk2[vars_ntwrk2 %in% net_vars]

    # Check that neutral networks are encoded by at least one genotype 
    if(!identical(vars_ntwrk1, character(0)) && !identical(vars_ntwrk2, character(0))){
      # Compute all direct links
      tmp1 <- as_data_frame(graph) %>% filter(from %in% vars_ntwrk1 & to %in% vars_ntwrk2)
      tmp2 <- as_data_frame(graph) %>% filter(from %in% vars_ntwrk2 & to %in% vars_ntwrk1)

      prox <- dim(tmp1)[1] + dim(tmp2)[1]
    }
  }
  else if(type == 3){
    # Compute average shortest path length between neutral networks

    # Extract genotypes from each neutral network: specific binders
    if(complex){
      if(from_specific) vars_ntwrk1 <- pheno_df %>% filter(specificity == n_ntwrk1 & bg == Bg) %>% pull(complex) # specific
      else vars_ntwrk1 <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == n_ntwrk1 & bg == Bg) %>% pull(complex) # binding
      
      if(to_specific) vars_ntwrk2 <- pheno_df %>% filter(specificity == n_ntwrk2 & bg == Bg) %>% pull(complex) # specific
      else vars_ntwrk2 <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == n_ntwrk2 & bg == Bg) %>% pull(complex) # binding
    }
    else{
      if(from_specific) vars_ntwrk1 <- pheno_df %>% filter(specificity == n_ntwrk1 & bg == Bg) %>% pull(AA_var) # specific
      else vars_ntwrk1 <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == n_ntwrk1 & bg == Bg) %>% pull(AA_var) # binding
      
      if(to_specific) vars_ntwrk2 <- pheno_df %>% filter(specificity == n_ntwrk2 & bg == Bg) %>% pull(AA_var) # specific
      else vars_ntwrk2 <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == n_ntwrk2 & bg == Bg) %>% pull(AA_var) # binding
    }

    # Filter variants to retain those present in main component genotype network
    net_vars <- extract_main_ntwrk(graph=graph,tr_mat=NULL,nodes=TRUE)
    vars_ntwrk1 <- vars_ntwrk1[vars_ntwrk1 %in% net_vars]
    vars_ntwrk2 <- vars_ntwrk2[vars_ntwrk2 %in% net_vars]

    # Check that neutral networks are encoded by at least one genotype 
    if(!identical(vars_ntwrk1, character(0)) && !identical(vars_ntwrk2, character(0))){
      # Compute shortest path lengths between every pair og genotypes in the network
      d_table <- distances(graph)
      d_table <- d_table[rownames(d_table) %in% vars_ntwrk1,colnames(d_table) %in% vars_ntwrk2]
      prox <- mean(d_table)
    }
  }
  return(prox) 
}

# Compute the fraction of neighbor nodes with same or different phenotype for a focal node
fraction_of_neighbors_per_genotype <- function(from_g,graph,Bg,pheno_tbl=phenotypes_tbl_prot,complex=FALSE){
  # from_g = focal genotype
  # Bg = string specifying the DBD background (= 'AncSR1' or 'AncSR2')
  # graph = a genotype network (igraph object)
  # pheno_tbl = df with phenotypic annotations per genotype
  # complex = whether to use protein-DNA complex genoytpes (default: FALSE)
  
  # Check inconsistent parameters
  if(!(Bg %in% c("AncSR1","AncSR2")) || is.null(Bg)){
    stop("Specify DBD background: 'AncSR1' or 'AncSR2'")
  }
  
  match_phenotypes <- function(g_i,g_j,complex){
    if(complex){
      p_i <- pheno_tbl %>% filter(complex == g_i & bg == Bg) %>% pull(specificity)
      p_j <- pheno_tbl %>% filter(complex == g_j & bg == Bg) %>% pull(specificity)
    }
    else{
      p_i <- pheno_tbl %>% filter(AA_var == g_i & bg == Bg) %>% pull(specificity)
      p_j <- pheno_tbl %>% filter(AA_var == g_j & bg == Bg) %>% pull(specificity)
    }
    return(p_i == p_j)
  }
  
  # extract neighbors of focal genotype and their phenotypes
  neighbors_df <- get.data.frame(as.directed(graph)) %>% filter(from == from_g) %>% mutate(pheno_match = map2_lgl(.x=from,.y=to,.f=match_phenotypes,complex=complex)) %>% 
    group_by(pheno_match) %>% reframe(count = n()) %>% ungroup() %>% mutate(total = sum(count)) %>%
    group_by(pheno_match) %>% mutate(prop = count / total,g=from_g)
  
  return(neighbors_df)
}

# Simulate genotype networks
simulate_GPmap <- function(graph,type=1,which="net",cores=1,seed = NULL,n_sample = NULL,genotypes=NULL){
  # graph: reference genotype network (igraph object) - for 'type' options 1,2 and 3
  # type: an integer indicating the type of simulation to produce
    # 1 - network built from hamming distance (not accounting for genetic code)
    # 2 - same number of nodes and edges but randomized connections between nodes
    # 3 - maximally connected network (same number of nodes)
    # 4 - network from random sample of genotypes 
    # 5 - network from a user-specified sample of genotypes
  # which: specify whether to return the genotype network ("net"), the adjacency matrix ("mat") or both ("both").(default: "net")
  # cores: specify the number of cores to use in to use for parallel computing (default: 1)
  # seed: integer (seed) to make random networks reproducible (default: NULL)
  # n_sample: Number of genotypes to randomly sample (default: NULL)
  # genotypes: User-specified sample of genotypes (default: NULL)
  
  # check inconsistent parameters
  if(!(which %in% c("net","mat","both"))){
    stop("Provide a valid output ('which') paramter: 'net','mat' or 'both'")
  }
  if(type %in% c(1,2,3) && is.null(graph)){
    stop("Must provide a reference genotype network to do the simulation")
  }
  if(type %in% c(4) && is.null(n_sample)){
    stop("Must provide a value > 0 for 'n_sample' argument")
  }
  if(type == 4 && !is.null(graph)){
    message("Ignoring genotype network...")
  }
  if(type == 5 && is.null(genotypes)){
    stop("Must provide a vector of genotypes")
  }

  # keep the number of nodes and edges
  if(type %in% c(1,2,3)){
    n_nodes <- length(V(graph))
    n_edges <- length(E(graph))
  }
  
  res <- NULL
  
  # simulate genotype network maintaining node identity from original network
  if(type == 1){
    # network built from hamming distance (not accounting for genetic code)
    adj_mat <- build_mutation_matrix(nodes = names(V(graph)),type = 5,cores = cores)
    tmp_net <- build_genotype_network(nodes = names(V(graph)),build_mat = FALSE,adj_mat = adj_mat,type = 5,cores = cores)
  }
  if(type == 2){
    # same number of edges but randomized connections between nodes
    if(!is.null(seed)){
      set.seed(seed)
    }
    tmp_net <- sample_gnm(n=n_nodes, m=n_edges)
    tmp_net <- set.vertex.attribute(tmp_net, "name", value=names(V(graph))) # rename vertices
    adj_mat <- as_adjacency_matrix(tmp_net,type = "both") # adjacency matrix
  }
  else if(type == 3){
    # maximally connected network
    tmp_net <- make_full_graph(n_nodes,loops = FALSE)
    tmp_net <- set.vertex.attribute(tmp_net, "name", value=names(V(graph))) # rename vertices
    adj_mat <- as_adjacency_matrix(tmp_net,type = "both") # adjacency matrix
  }
  else if(type == 4){
    # make network from random sample of all possible genotypes
    set.seed(seed)
    possible_AAvar <- do.call(paste,expand.grid(AA_STANDARD,AA_STANDARD,AA_STANDARD,AA_STANDARD)) %>% gsub(" ","",.)
    g <- sample(possible_AAvar,n_sample,replace=FALSE)
    adj_mat <- build_mutation_matrix(nodes = g,type = 3,cores = cores)
    tmp_net <- build_genotype_network(nodes = g,build_mat = FALSE,adj_mat = adj_mat,type = 3,cores = cores)
  }
  else if(type == 5){
    # network from a user-specified sample of genotypes
    g <- genotypes
    adj_mat <- build_mutation_matrix(nodes = g,type = 3,cores = cores)
    tmp_net <- build_genotype_network(nodes = g,build_mat = FALSE,adj_mat = adj_mat,type = 3,cores = cores)
  }
  
  if(which == "net") res <- tmp_net
  else if(which == "mat") res <- adj_mat
  else if(which == "both") res <- list(tmp_net,adj_mat)
  return(res)
}

# Function to check whether phenotype is in encoded by at least one variant
is.encoded <- function(RE,graph,pheno_df,Bg,complex=FALSE){
  # extract binders
  if(complex) vars_RE<- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == RE & bg == Bg) %>% pull(complex)
  else vars_RE <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == RE & bg == Bg) %>% pull(AA_var)

  return(!identical(vars_RE, character(0)))
}

# Function to check whether phenotype is in main network
is.in.ntwrk <- function(RE,graph,pheno_df,Bg,specific=FALSE,complex=FALSE){
  # genotypes in main component genotype network
  net_vars <- extract_main_ntwrk(graph=graph,tr_mat=NULL,nodes=TRUE)

  # extract binders
  if(complex){
    if(specific){
      vars_RE <- pheno_df %>% filter(specificity == RE & bg == Bg) %>% pull(complex)
    }
    else{
      vars_RE <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == RE & bg == Bg) %>% pull(complex)
    }
  }
  else{
    if(specific){
      vars_RE <- pheno_df %>% filter(specificity == RE & bg == Bg) %>% pull(AA_var)
    }
    else vars_RE <- pheno_df %>% unnest(bound_REs) %>% filter(bound_REs == RE & bg == Bg) %>% pull(AA_var)
  }

  #filter variants to those present in the main network
  vars_RE <- vars_RE[vars_RE %in% net_vars]

  return(!identical(vars_RE, character(0)))
}

# Function to check whether phenotype is in bound specifically (for protein netwrok)
is.bound.specific <- function(RE,graph,pheno_df,Bg,complex=FALSE){
  # extract specific binders
  if(complex) vars_RE <- pheno_df %>% filter(specificity == RE & bg == Bg) %>% pull(complex)
  else vars_RE <- pheno_df %>% filter(specificity == RE & bg == Bg) %>% pull(AA_var)

  return(!identical(vars_RE, character(0)))
}

# Function to remove promiscuous from PDFV table and renormalize probabilities
remove_promiscuous <- function(pdfv_tbl,norm=F){
  if(norm) pdfv_tbl %>% filter(RE != "Promiscuous") %>% mutate(Norm_F_prob = Norm_F_prob /sum(Norm_F_prob))
  else pdfv_tbl %>% filter(RE != "Promiscuous")
}

# Function to get number of genotypes and phenotypes accessed during evolutionary exploration
network_exploration <- function(n_steps,tr_mat,Bg,from=REF_GENOTYPE,pheno_tbl,complex=F,get_pheno_id=F){
  # n_steps = number of iterations to run mc chain
  # tr_mat = transition matrix
  # Bg = ancestral background ("AncSR1" or "AncSR2")
  # from = character string or vector of characters indicating the initial states of the mc chain (default: REF_GENOTYPE)
  # pheno_tbl = phenotypic table with annotations per genotype
  # complex = whether to compute analysis for protein-DNA complexes
  # get_pheno_id = whether to return the IDs of bound phenotypes (T) or just counts (F). (default: F)

  mc_tmp <- simulate_markov_chain(states_0 = from,tr_mat = tr_mat,n_steps = n_steps)
  
  # Get number of accessible genotypes
  n_genotypes <- length(mc_tmp[mc_tmp > 0])
  
  # Get number of accessible phenotypes
  pdfv_tmp_all <- get_PDFV_v2(mc_tmp,type = "simulated mc",Bg = Bg,model = n_steps,pheno_tbl = pheno_tbl,complex = complex,specific = F)
  pdfv_tmp_spec <- get_PDFV_v2(mc_tmp,type = "simulated mc",Bg = Bg,model = n_steps,pheno_tbl = pheno_tbl,complex = complex,specific = T)
  
  n_bound <- pdfv_tmp_all %>% filter(Norm_F_prob > 0) %>% summarize(c = n()) %>% pull(c)
  n_bound_spec <- remove_promiscuous(pdfv_tmp_spec) %>% filter(Norm_F_prob > 0) %>% summarize(c = n()) %>% pull(c)
  
  df <- data.frame(Bg=Bg,mut_step = n_steps, n_genotypes = n_genotypes, n_bound = n_bound, n_bound_spec = n_bound_spec)
  
  # Get identity of bound phenotypes
  if(get_pheno_id){
    id_bound <- pdfv_tmp_all %>% filter(Norm_F_prob > 0) %>% pull(RE)
    id_bound_spec <- pdfv_tmp_spec %>% filter(RE != "Promiscuous") %>% filter(Norm_F_prob > 0) %>% pull(RE)
  
    df <- tibble(Bg=Bg,mut_step = n_steps, n_genotypes = n_genotypes, n_bound = n_bound, n_bound_spec = n_bound_spec, 
                  id_bound = list(id_bound), id_bound_spec = list(id_bound_spec))
  }
  return(df)
}


# Fraction of neighbor nodes that result in amino acid change (useful only for prot-DNA networks)
frac_aa_subs_per_g <- function(graph,nodes_in_ntwrk){
  # graph = a genotype network (igraph object)
  # nodes_in_ntwrk = nodes in the main component of the network

  # neighbors per focal genotype
  n <-  map(.x=nodes_in_ntwrk,.f=neighbors,graph=graph)
  n <- lapply(n,names)
  
  # protein genotypes of neighbor complexes and focal genotypes
  p_g <- map(.x=nodes_in_ntwrk,.f=str_sub,start=1,end=4)
  p_n <- lapply(n,str_sub,start=1,end=4)
  
  # fraction of unique protein neighbors
  fr <- map2_dbl(.x=p_g,.y=p_n,.f=function(x,y) length(y[x != y])/length(y))

  names(fr) <- nodes_in_ntwrk
  return(fr)
}

# Compute the average number of AA substitutions per step (matrix multiplication) for protein-DNA simulations:
# For every substitution step that includes protein OR DNA substitutions, compute the expected number of subtitutions
# that are due to AA subs only 
compute_aa_subs_from_protDNA_steps <- function(states_0,n_steps,tr_mat,fr_aa_subs){
  # states_0 = vector containing the states from which markov chain starts
  # n_steps = number of (mutational) steps to run the Markov chain
  # tr_mat = a probability transition matrix
  # fr_aa_neighbors = fraction of subs that result in amino acid changes (per genotype)

  # Check for inconsistent/valid paramenters
  # Check whether 'tr_mat' is an ergodic transition probability matrix
  if(!((is.matrix(tr_mat) || inherits(tr_mat,"dgCMatrix") || inherits(tr_mat,"dsCMatrix")) && 
    nrow(tr_mat)==ncol(tr_mat) && all(tr_mat>=0 & tr_mat<=1) && sum(apply(tr_mat,1,sum))==ncol(tr_mat))){
    stop("Provide a square probability transition matrix")
  }
  if(is.null(n_steps) || n_steps==0){
    stop("Provide a correct number of steps to run the markov chain (n_steps > 0)")
  }

  # Compute the average number of amino acid substitutions per prot-DNA step in the markov chain
  state_probs <- simulate_markov_chain(states_0=states_0,n_steps=n_steps,tr_mat=tr_mat)
  avg_aa_subs <- sum(state_probs*fr_aa_subs)
  return(avg_aa_subs)
}


# Function to get the amount of bias (and direction) from a GP map:
# Either due to production, equilibrium freqs., or around the local neighborhood of a variant (or set of variants)
get_bias <- function(from,Bg,tr_mat=NULL,n_steps=NULL,complex=F,pheno_tbl=NULL,main_direction=FALSE,input_freqs=FALSE,specific=T){
  # from = Either "production", "equilibrium", or a specific genotype
  # bg = a string indicating the DBD background ("AncSR1" or "AncSR2")
  # tr_mat = a transition matrix
  # pheno_tbl = table with annotated phenotypes per variant
  # complex = a logical argument to indicate whether PDFV is to be computed for prot-DNA complexes
  # main_direction = whether to get the main direction of the bias
  # input_freqs = whether to simulate a mc or use input frequency distribution
  
  # Check for inconsistent parameters
  if(input_freqs == F && (!(Bg %in% c("AncSR1","AncSR2")) || is.null(Bg))){
    stop("Specify DBD background: 'AncSR1' or 'AncSR2'")
  }
  if(input_freqs == T && !is.null(n_steps)){
    message("Ignoring number of steps...")
  }
  if(input_freqs && (is.null(from))){
    stop("Provide a data frame or vector with the frequency distribution")
  }
  if(input_freqs == F && is.null(from)){
    stop("Provide a starting genotype or set of genotypes to compute the bias")
  }
  
  if(input_freqs){
    if(is.null(dim(from))){
      # vector of frequencies
      f <- from[from > 0]
      dir <- from[from == max(from)]
    }
    else {
      f <- from %>% remove_promiscuous(norm = T) %>% filter(Norm_F_prob > 0) %>% pull(Norm_F_prob)
      dir <- from %>% remove_promiscuous(norm = T) %>% filter(Norm_F_prob==max(Norm_F_prob)) %>% pull(RE)
    }
  }
  else{
    x <- simulate_markov_chain(states_0 = from,tr_mat = tr_mat,n_steps = n_steps)
    f <- get_PDFV_v2(x,type = "simulated mc",Bg = Bg,pheno_tbl = pheno_tbl,specific = specific,model="NA",complex = complex)
    dir <- f %>% remove_promiscuous(norm = T) %>% filter(Norm_F_prob==max(Norm_F_prob)) %>% pull(RE)
    f <- f %>% remove_promiscuous(norm = T) %>% filter(Norm_F_prob > 0) %>% pull(Norm_F_prob)
  }
  
  # Compute bias
  B <- 1-(-sum(f * log(f,base=16)))
  if(main_direction){
    r <- list(B,dir)
    return(r)
  }
  else return(B)
}


# Function to compute the Kullback-Leibler divergence metric (relative entropy). Measures the "distance" between
# a distribution and a reference distribution. It is asymmetric (depends on which distribution is the reference, Qx)
KL_div <- function(Px,Qx){
  if(length(Px) != 16){
    extra <- 16 - length(Px)
    Px <- c(Px,rep(0,extra))
  }

  if(0 %in% Px | 0 %in% Qx){
    Px_tmp <- Px + .Machine$double.xmin
    Px <- Px_tmp / sum(Px_tmp)
    Qx_tmp <- Qx + .Machine$double.xmin
    Qx <- Qx_tmp / sum(Qx_tmp)
  }
  sum(Px * log(Px/Qx)) / log(16) # normalized so that KL ranges from 0 to 1
}

# Function to compute the Population Stability Index (PSI), a symmetric derivation of the KL divergence
# metric. Captures the total difference in information between two distributions
PSI_div <- function(Px,Qx){
  # PSI(Px,Qx) = KL(P||Q) + KL(Q||P)
  if(0 %in% Px | 0 %in% Qx){
    Px_tmp <- Px + .Machine$double.xmin
    Px <- Px_tmp / sum(Px_tmp)
    Qx_tmp <- Qx + .Machine$double.xmin
    Qx <- Qx_tmp / sum(Qx_tmp)
  }
  sum((Px-Qx) * log(Px/Qx))
}

# function to compute Shannon's entropy per site for a sequence alignment
shannon_entropy <- function(AA_vars,n_sites=4,which="bits"){
  # vector of amino acid variants
  # n_sites = number of sites in the alignment
  # which = entropy (H per site) - bits (height of state X at position i) - freq (prop of every AA state per site)
  
  vars_df <- as.list(AA_vars)
  vars_df <- lapply(vars_df,strsplit,split = '')
  vars_df <- matrix(unlist(vars_df,recursive = T), ncol = n_sites, byrow = TRUE)
  
  # Compute Entropy (H) per site (from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC332411/pdf/nar00204-0153.pdf)
  entropy_per_site <- c()
  bits_per_state_per_site <-  matrix(NA,ncol = length(AA_STANDARD),nrow = n_sites)
  freq_states_per_site <- matrix(NA,ncol = length(AA_STANDARD),nrow = n_sites)
  colnames(freq_states_per_site) <- colnames(bits_per_state_per_site) <- AA_STANDARD
  rownames(freq_states_per_site) <- rownames(bits_per_state_per_site) <- 1:n_sites
  
  for(i in 1:n_sites){
    # Entropy
    freqs_s_i <- table(vars_df[,i])
    prop_s_i <- freqs_s_i / dim(vars_df)[1]
    H_i <- -sum(prop_s_i * log(prop_s_i,base=20))
    entropy_per_site <- c(entropy_per_site,H_i)
    
    # frequencies for all 20 states per site
    freqs <- table(factor(vars_df[,i],levels = AA_STANDARD)) / dim(vars_df)[1]
    freq_states_per_site[i,] <- freqs
    
    # Bits per state per site
    total_info_s_i <- log2(20) - (-sum(prop_s_i * log2(prop_s_i)))
    bits_per_state_site_i <- freqs * total_info_s_i
    bits_per_state_per_site[i,] <- bits_per_state_site_i
  }
  if(which == "bits"){
    return(bits_per_state_per_site)
  }
  else if(which == "freq"){
    return(freq_states_per_site)
  }
  else if(which == "entropy"){
    return(entropy_per_site)
  }
}

# Returns a matrix (or plot) with the marginal frequencies for each possible state*site combination, 
# and the observed co-frequency (analagous to a variance-covariance matrix)
freq_cofreq_mat <- function(df,plot=T){
  # df = data frame with production spectra
  
  df.heatmap <- remove_promiscuous(df,norm=T) %>% 
    mutate(RE = case_when(RE == "SRE (AA)" ~ "AA",RE == "ERE (GT)" ~ "GT",T ~ RE),
           RE_1 = str_sub(RE,1,1),
           RE_2 = str_sub(RE,2,2)) %>%
    mutate_at(vars(Norm_F_prob),list(as.numeric)) 
  
  
  if(plot==T){
    h_total <- df.heatmap %>%
      group_by(RE_2) %>%
      reframe(Norm_F_prob = sum(Norm_F_prob)) %>%
      mutate(RE_1 = "Total")
    
    v_total <- df.heatmap %>%
      group_by(RE_1) %>%
      reframe(Norm_F_prob = sum(Norm_F_prob)) %>%
      mutate(RE_2 = "Total")
    
    midpoint_marginal <- mean(c(v_total$Norm_F_prob,h_total$Norm_F_prob))
    midpoint_covar <- mean(df.heatmap$Norm_F_prob)
    
    ggplot(as.data.frame(df.heatmap), aes(x = RE_1, y = RE_2)) +
      geom_tile(aes(fill = Norm_F_prob)) +
      geom_text(size = 3, aes(label = round(Norm_F_prob,3))) + #displays cell values
      scale_fill_gradient2(high = "purple", #colors
                           mid = "lightblue",
                           low = "white",
                           midpoint = midpoint_covar) +
      geom_point(data = h_total, 
                 aes(color = Norm_F_prob), 
                 size = 10, 
                 shape = 19) +
      geom_point(data = v_total, 
                 aes(color = Norm_F_prob), 
                 size = 10, 
                 shape = 19) +
      scale_color_gradient2(high = "red", #colors
                            mid = "orange",
                            low = "white",
                            midpoint = midpoint_marginal) +
      geom_text(data = h_total, size = 3, aes(label = round(Norm_F_prob,3))) +
      geom_text(data = v_total, size = 3, aes(label = round(Norm_F_prob,3))) +
      theme(panel.grid.major.x=element_blank(), #no gridlines
            panel.grid.minor.x=element_blank(), 
            panel.grid.major.y=element_blank(), 
            panel.grid.minor.y=element_blank(),
            panel.background=element_rect(fill="white"),
            axis.text.x = element_text(angle=0, hjust = 0.5,vjust=0.5, size = 10,face = NULL),
            axis.text.y = element_text(size = 10,face = NULL),
            legend.title=element_text(face="bold", size=8))  + 
      labs(x="site 1",y="site 2",fill="Co-frequency",color="Marginal freq")
  }
  else{
    df.heatmap <- df.heatmap %>% acast(RE_1 ~ RE_2,value.var="Norm_F_prob")
    total_s1 <- apply(df.heatmap, 1, sum)
    total_s2 <- apply(df.heatmap, 2, sum); total_s2 <- c(total_s2,NA)
    df.heatmap <- cbind(df.heatmap,total_s1)
    df.heatmap <- rbind(df.heatmap,total_s2)
    df.heatmap
  }
}

