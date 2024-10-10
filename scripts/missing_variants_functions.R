################################################################################
# functions for fitting reference-free mutation effects model for RH-RE DMS data
################################################################################

# Build a covariate matrix of RH/RE mutation effects of the specified order.
# RH.order is the maximum order of RH effects.
# RE.order is the maximum order of RE effects.
# RH.RE.order is either a 2-column matrix of RHxRE interactions to include, 
# where the first column is the order of RH effects and the second column is the 
# corresponding order of RE effects in the interaction, or 0 if no RHxRE 
# interactions.
# Note that the model will always include lower order "component" effects up to
# the specified order.
# If reffree = FALSE, then a single amino acid/RE state is chosen as the 
# reference state.
# Assume that RH variants are ordered alphabetically in the data.
# cores indicates the number of cores to use for parallel computing.
# its indicates the number of iterations to use for serialization.
make_cov_matrix <- function(data, RH.order = 1, RE.order = 0, RH.RE.order = 0, 
                            sparse = TRUE, reffree = TRUE, cores = 1, its = 1) {
  
  # create data frame of AAs at each RH site for all variants in the data
  statemat <- str_split(as.character(data$AA_var), pattern = "", 
                        simplify = TRUE)
  statemat <- as.data.frame(statemat, stringsAsFactors = TRUE)
  colnames(statemat) <- c("AA1", "AA2", "AA3", "AA4")
  
  if(RE.order > 0) {
    # create data frame of nts at each RE site for all variants in the data
    REstatemat <- gsub("[ES]RE.+\\(|\\)", "", data$RE)
    REstatemat <- str_split(REstatemat, pattern = "", simplify = TRUE)
    REstatemat <- as.data.frame(REstatemat, stringsAsFactors = TRUE)
    colnames(REstatemat) <- c("RE1", "RE2")
    
    statemat <- cbind(statemat, REstatemat)
  }
  
  if(reffree) {
    # modify contrasts for each site so that they include all states in the data
    for(i in 1:ncol(statemat)) {
      attributes(statemat[,i])$contrasts <- contrasts(statemat[,i], 
                                                      contrasts = FALSE)
    }
  }
  
  # construct formula for covariate matrix
  fmla <- "~ 0 + (AA1 + AA2 + AA3 + AA4)"
  if(RH.order > 1) fmla <- paste0(fmla, "^", RH.order)
  if(RE.order > 0) fmla <- paste0(fmla, " + (RE1 + RE2)")
  if(RE.order > 1) fmla <- paste0(fmla, "^", RE.order)
  if(!identical(RH.RE.order, 0)) {
    for(i in 1:nrow(RH.RE.order)) {
      fmla <- paste0(fmla,
                     " + (AA1 + AA2 + AA3 + AA4)",
                     ifelse(RH.RE.order[i,1] > 1, paste0("^", RH.RE.order[i,1]), ""),
                     " * (RE1 + RE2)",
                     ifelse(RH.RE.order[i,2] > 1, paste0("^", RH.RE.order[i,2]), ""))
    }
  }
  fmla <- as.formula(fmla)
  
  # construct covariate matrix
  cov_mat <- list()
  size1 <- nrow(statemat)
  chunksize1 <- ceiling(size1 / its)
  
  # serial processing for matrix construction
  for(i in 1:its) {
    statemat_chunk1 <- statemat[((i-1)*chunksize1+1):min(i*chunksize1, size1),]
    
    # parallel processing for matrix construction
    cl <- parallel::makeCluster(cores, "FORK", outfile="")
    doParallel::registerDoParallel(cl)
    # construct covariate matrix
    cov_mat[[i]] <- foreach(j = 1:cores, .combine = "rbind") %dopar% {
      # split statemat into chunks
      size2 <- nrow(statemat_chunk1)
      chunksize2 <- ceiling(size2 / cores)
      statemat_chunk2 <- statemat_chunk1[((j-1)*chunksize2+1):min(j*chunksize2, size2),]
      
      # construct covariate matrix for chunk
      model.Matrix(fmla, statemat_chunk2, sparse = sparse)
    }
    stopCluster(cl)
  }
  
  cov_mat <- do.call("rbind", cov_mat)
  return(cov_mat)
}


# objective function for nonlinear least squares with logistic nonspecific
# epistasis model
obj_LU <- function(par, X, y) {
  beta <- par[1:(length(par) - 2)]
  L <- par[length(par) - 1]
  U <- par[length(par)]
  
  return(sum((y - L - (U - L) / (1 + exp(-X %*% beta)))^2))
}


# gradient function for NLS with logistic nonspecific epistasis model, 
# inferring lower and upper bound
gr_LU <- function(par, X, y) {
  beta <- par[1:(length(par) - 2)]
  L <- par[length(par) - 1]
  U <- par[length(par)]
  
  z0 <- as.vector(exp(-X %*% beta))
  z1 <- (U - L) / (1 + z0)
  
  # partial derivatives for beta
  partial.beta <- (-2 * (y - L - z1) * z1^2 / (U - L) * z0) %*% X
  partial.beta <- as.vector(partial.beta)
  
  # partial derivative for L
  partial.L <- sum(-2 * (y - L - z1) * z0 / (1 + z0))
  
  # partial derivative for U
  partial.U <- sum(-2 * (y - L - z1) * 1 / (1 + z0))
  
  return(c(partial.beta, partial.L, partial.U))
}


# Fit an unregularized mutation effects model with logistic nonspecific epistasis.
# y is a vector of observed phenotypes.
# cov_mat is the covariance matrix of mutation effects contributing to the
# phenotype of each variant.
# method is the algorithm used by the optim function for model fitting
#
# Outputs the fitted model from optim; the first k-2 arguments of the par
# attribute are the betas; the k-1th element is L and the kth element is U,
# where k is the length of par
fit_LU <- function(y, cov_mat, method = "L-BFGS-B", maxit = 100) {
  
  # starting values for parameters
  start <- c(rep(0, ncol(cov_mat)), min(y), max(y))
  
  if(method == "L-BFGS-B") {
    # lower and upper bounds for parameters
    lower <- c(rep(-Inf, ncol(cov_mat)), rep(min(y), 2))
    upper <- c(rep(Inf, ncol(cov_mat)), rep(max(y), 2))
    
    model <- optim(start, obj_LU, gr_LU, method = method, X = cov_mat, y = y,
                   lower = lower, upper = upper, control = list(maxit = maxit))
  } 
  
  else if(method == "BFGS" || method == "CG") {
    model <- optim(start, obj_LU, gr_LU, method = method, X = cov_mat, y = y,
                   control = list(maxit = maxit))
  } 
  
  else if(method == "Nelder-Mead" || method == "SANN") {
    model <- optim(start, obj_LU, method = method, X = cov_mat, y = y, 
                   control = list(maxit = maxit))
  }
  
  return(model)
}

# get genetic scores from data fitted using fit_LU
# fit is the output from fit_LU
# its is the number of times to loop through the data
# cores is the number of CPUs to use for parallel processing
get_gs_LU_fit <- function(cov_mat, fit, its, cores) {
  gs <- numeric(length = nrow(cov_mat))
  for(i in 1:its) {
    size1 <- nrow(cov_mat)
    chunksize1 <- ceiling(size1 / its)
    chunk1 <- cov_mat[((i-1)*chunksize1+1):min(i*chunksize1, size1),]
    
    # parallel processing
    cl <- parallel::makeCluster(cores, "FORK", outfile="")
    doParallel::registerDoParallel(cl)
    chunk1_gs <- foreach(j = 1:cores, .combine = "c") %dopar% {
      # split covariate matrix into chunks
      size2 <- nrow(chunk1)
      chunksize2 <- ceiling(size2 / cores)
      chunk2 <- chunk1[((j-1)*chunksize2+1):min(j*chunksize2, size2),]
      
      # compute genetic scores for chunk
      as.numeric(chunk2 %*% fit$par[1:(length(fit$par)-2)])
    }
    stopCluster(cl)
    
    gs[((i-1)*chunksize1+1):min(i*chunksize1, size1)] <- chunk1_gs
  }
  return(gs)
}

#####################################
# functions for fitting with glmnet #
#####################################

# generalized logistic function with lower and upper bounds L and U
# this is the inverse link function; this is the inverse link function for
# glmnet
logistic <- function(x, L, U) {
  return(L + (U - L) / (1 + exp(-x)))
}

# generalized logit link function (inverse of the generalized logistic function);
# this is the link function for glmnet
logit <- function(x, L, U) {
  l <- log((x - L) / (U - x))
  # define outputs when x is outside the domain of the function 
  # so that glmnet will converge
  l[x <= L] <- -Inf
  l[x >= U] <- Inf
  return(l)
}

# gradient of generalized logistic function with respect to x
logistic.gr <- function(x, L, U) {
  d <- ((U - L) * exp(-x)) / (1 + exp(-x))^2
  d[is.nan(d)] <- 0  # NaNs occur if x is very small; derivative should be almost 0
  return(d)
}

# Build a generalized logit link model.
glogit <- function(L, U) {
  linkfun <- function(mu) logit(mu, L, U)
  linkinv <- function(eta) logistic(eta, L, U)
  mu.eta <- function(eta) logistic.gr(eta, L, U)
  valideta <- function(eta) TRUE  # whether a linear predictor eta is within the domain of linkinv
  link <- paste0("glogit(L=", L, ",U=", U, ")")  # name of the link function
  return(structure(list(linkfun = linkfun, linkinv = linkinv,
                        mu.eta = mu.eta, valideta = valideta, name = link),
                   class = "link-glm"))
}


##########################################################
# functions for glmnet v4.1.6 with large sparse matrices #
##########################################################


# Helper function to compute the weighted mean and sd
# If x is sparse, then use the sparseMatrixStats package to compute the
# weighted sd. This gives different values than the original implementation,
# so use with caution if standardizing the matrix.
weighted_mean_sd_large <- function(x, weights=rep(1,nrow(x))){
  xm <- drop(t(weights/sum(weights))%*%x)
  if(!inherits(x, "sparseMatrix"))
    xv <- drop(t(weights/sum(weights))%*%scale(x,xm,FALSE)^2)
  else {
    x <- as(x, "dgCMatrix")
    xv <- sparseMatrixStats::colWeightedVars(x, weights)
  }
  
  xv[xv < 10*.Machine$double.eps] <- 0
  list(mean = xm, sd = sqrt(xv))
}


# glmnet.path from glmnet v4.1.6 that can handle large sparse matrices

# can also be fed parameters for a manual warmstart using init; this should
# be a list with elements a0 and beta giving the intercept and remaining 
# parameters

# set standardize to FALSE (this implementation uses the sparseMatrixStats 
# package to calculated weighted sds, which is different from how they are
# calculated in the original glmnet package; this only affects the fit when
# standardize = TRUE)
glmnet.path.large <- function(x, y, init = NULL, weights = NULL, 
                              lambda = NULL, nlambda = 100,
                              lambda.min.ratio = ifelse(nobs < nvars, 0.01, 0.0001),
                              alpha = 1.0, offset = NULL, family = gaussian(),
                              standardize = FALSE, intercept = TRUE, 
                              thresh = 1e-10, maxit = 100000,
                              penalty.factor = rep(1.0, nvars), 
                              exclude = integer(0), lower.limits = -Inf,
                              upper.limits = Inf, trace.it = 0) {
  
  ### Check on family argument
  if(is.function(family)) family = family()
  if(!inherits(family, "family"))
    stop("Invalid family argument; must be either character, 
         function or family object")
  
  ### Prepare all the generic arguments
  if (alpha > 1) {
    warning("alpha > 1; set to 1")
    alpha = 1
  } else if (alpha < 0) {
    warning("alpha < 0; set to 0")
    alpha = 0
  }
  alpha = as.double(alpha)
  
  this.call <- match.call()
  
  np = dim(x)
  if(is.null(np) || (np[2] <= 1)) 
    stop("x should be a matrix with 2 or more columns")
  nobs = as.integer(np[1]); nvars = as.integer(np[2])
  
  # get feature variable names
  vnames = colnames(x)
  if(is.null(vnames)) vnames = paste("V", seq(nvars), sep = "")
  
  # check weights
  if(is.null(weights)) weights = rep(1,nobs)
  else if (length(weights) != nobs)
    stop(paste("Number of elements in weights (",length(weights),
               ") not equal to the number of rows of x (",nobs,")",sep=""))
  weights <- as.double(weights)
  
  # initialize from family function. Makes y a vector in case of binomial, and 
  # possibly changes weights. Expects nobs to be defined, and creates n and 
  # mustart (neither used here).  Some cases expect to see things, so we set it 
  # up just to make it work
  etastart = 0; mustart = NULL; start = NULL
  eval(family$initialize)
  
  ## Just in case this was not done in initialize()
  y <- drop(y)  # we don't like matrix responses
  
  is.offset <- !(is.null(offset))
  if (is.offset == FALSE) {
    offset <- as.double(y * 0) #keeps the shape of y
  }
  # infinite penalty factor vars are excluded
  if(any(penalty.factor == Inf)) {
    exclude = c(exclude, seq(nvars)[penalty.factor == Inf])
    exclude = sort(unique(exclude))
  }
  
  # Compute weighted mean and variance of columns of x, sensitive to sparse 
  # matrix needed to detect constant columns below, and later if standardization
  meansd <- weighted_mean_sd_large(x, weights)
  
  ## look for constant variables, and if any, then add to exclude
  
  const_vars <- meansd$sd == 0
  nzvar <- setdiff(which(!const_vars), exclude)
  # if all the non-excluded variables have zero variance, throw error
  if (length(nzvar) == 0) stop("All used predictors have zero variance")
  
  ## if any constant vars, add to exclude
  if(any(const_vars)) {
    exclude <- sort(unique(c(which(const_vars), exclude)))
    meansd$sd[const_vars] <- 1.0 ## we divide later, and do not want bad numbers
  }
  if(length(exclude) > 0) {
    jd = match(exclude, seq(nvars), 0)
    if(!all(jd > 0)) stop ("Some excluded variables out of range")
    penalty.factor[jd] = 1 # ow can change lambda sequence
  }
  # check and standardize penalty factors (to sum to nvars)
  vp = pmax(0, penalty.factor)
  if (max(vp) <= 0) stop("All penalty factors are <= 0")
  vp = as.double(vp * nvars / sum(vp))
  
  
  ### check on limits
  control <- glmnet.control()
  if (thresh >= control$epsnr)
    warning("thresh should be smaller than glmnet.control()$epsnr",
            call. = FALSE)
  
  if(any(lower.limits > 0)) { stop("Lower limits should be non-positive") }
  if(any(upper.limits < 0)) { stop("Upper limits should be non-negative") }
  lower.limits[lower.limits == -Inf] = -control$big
  upper.limits[upper.limits == Inf] = control$big
  if (length(lower.limits) < nvars) {
    if(length(lower.limits) == 1) lower.limits = rep(lower.limits, nvars) else
      stop("Require length 1 or nvars lower.limits")
  } else lower.limits = lower.limits[seq(nvars)]
  if (length(upper.limits) < nvars) {
    if(length(upper.limits) == 1) upper.limits = rep(upper.limits, nvars) else
      stop("Require length 1 or nvars upper.limits")
  } else upper.limits = upper.limits[seq(nvars)]
  
  if (any(lower.limits == 0) || any(upper.limits == 0)) {
    ###Bounds of zero can mess with the lambda sequence and fdev;
    ###ie nothing happens and if fdev is not zero, the path can stop
    fdev <- glmnet.control()$fdev
    if(fdev!= 0) {
      glmnet.control(fdev = 0)
      on.exit(glmnet.control(fdev = fdev))
    }
  }
  ### end check on limits
  ### end preparation of generic arguments
  
  # standardize x if necessary
  if (intercept) {
    xm <- meansd$mean
  } else {
    xm <- rep(0.0, times = nvars)
  }
  if (standardize) {
    xs <- meansd$sd
  } else {
    xs <- rep(1.0, times = nvars)
  }
  if (!inherits(x, "sparseMatrix")) {
    x <- scale(x, xm, xs)
  } else {
    attr(x, "xm") <- xm
    attr(x, "xs") <- xs
  }
  lower.limits <- lower.limits * xs
  upper.limits <- upper.limits * xs
  
  # get null deviance and lambda max
  start_val <- glmnet:::get_start(x, y, weights, family, intercept, is.offset,
                                  offset, exclude, vp, alpha)
  
  # work out lambda values
  nlam = as.integer(nlambda)
  user_lambda = FALSE   # did user provide their own lambda values?
  if (is.null(lambda)) {
    if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
    
    # compute lambda max: to add code here
    lambda_max <- start_val$lambda_max
    
    # compute lambda sequence
    ulam <- exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
                    length.out = nlam))
  } else {  # user provided lambda values
    user_lambda = TRUE
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam = as.double(rev(sort(lambda)))
    nlam = as.integer(length(lambda))
  }
  
  # start progress bar
  if (trace.it == 1) pb <- utils::txtProgressBar(min = 0, max = nlam, style = 3)
  
  a0 <- rep(NA, length = nlam)
  beta <- matrix(0, nrow = nvars, ncol = nlam)
  dev.ratio <- rep(NA, length = nlam)
  # fit <- NULL
  fit <- init
  mnl <- min(nlam, control$mnlam)
  for (k in 1:nlam) {
    # get the correct lambda value to fit
    if (k > 1) {
      cur_lambda <- ulam[k]
    } else {
      cur_lambda <- ifelse(user_lambda, ulam[k], control$big)
    }
    
    if (trace.it == 2) cat("Fitting lambda index", k, ":", ulam[k], fill = TRUE)
      
    fit <- glmnet:::glmnet.fit(x, y, weights / sum(weights), cur_lambda, 
                               alpha = alpha, offset = offset, family = family, 
                               intercept = intercept, thresh = thresh,
                               maxit = maxit, penalty.factor = vp, 
                               exclude = exclude, lower.limits = lower.limits, 
                               upper.limits = upper.limits, warm = fit, 
                               from.glmnet.path = TRUE, save.fit = TRUE,
                               trace.it = trace.it)
    if (trace.it == 1) utils::setTxtProgressBar(pb, k)
    # if error code non-zero, a non-fatal error must have occurred
    # print warning, ignore this lambda value and return result
    # for all previous lambda values
    if (fit$jerr != 0) {
      errmsg <- glmnet:::jerr.glmnetfit(fit$jerr, maxit, k)
      warning(errmsg$msg, call. = FALSE)
      k <- k - 1
      break
    }
    
    a0[k] <- fit$a0
    beta[, k] <- as.matrix(fit$beta)
    dev.ratio[k] <- fit$dev.ratio
    
    # early stopping if dev.ratio almost 1 or no improvement
    if (k >= mnl && user_lambda == FALSE) {
      if (dev.ratio[k] > control$devmax) break
      else if (k > 1) {
        if (family$family == "gaussian") {
          if (dev.ratio[k] - dev.ratio[k-1] < control$fdev * dev.ratio[k])
            break
        } else if (family$family == "poisson") {
          if (dev.ratio[k] - dev.ratio[k - mnl + 1] <
              10 * control$fdev * dev.ratio[k])
            break
        } else if (dev.ratio[k] - dev.ratio[k-1] < control$fdev) break
      }
    }
  }
  if (trace.it == 1) {
    utils::setTxtProgressBar(pb, nlam)
    cat("", fill = TRUE)
  }
  
  # truncate a0, beta, dev.ratio, lambda if necessary
  if (k < nlam) {
    a0 <- a0[1:k]
    beta <- beta[, 1:k, drop = FALSE]
    dev.ratio <- dev.ratio[1:k]
    ulam <- ulam[1:k]
  }
  
  # return coefficients to original scale (because of x standardization)
  beta <- beta / xs
  a0 <- a0 - colSums(beta * xm)
  
  # output
  stepnames <- paste0("s", 0:(length(ulam) - 1))
  out <- list()
  out$a0 <- a0
  names(out$a0) <- stepnames
  out$beta <- Matrix::Matrix(beta, sparse = TRUE,
                             dimnames = list(vnames, stepnames))
  out$df <- colSums(abs(beta) > 0)
  out$dim <- dim(beta)
  out$lambda <- ulam
  out$dev.ratio <- dev.ratio
  out$nulldev <- start_val$nulldev
  out$npasses <- fit$npasses
  out$jerr <- fit$jerr
  out$offset <- is.offset
  out$call <- this.call
  out$family <- family
  out$nobs <- nobs
  class(out) <- c("glmnetfit", "glmnet")
  
  return(out)
}


# glmnet for large matrices

# can also be fed parameters for a manual warmstart using init; this should
# be a list with elements a0 and beta giving the intercept and remaining 
# parameters (both numeric vectors, with a0 length 1)

# defaults to standardize = FALSE

glmnet.large <- function(x, y, family = c("gaussian", "binomial", "poisson",
                                          "multinomial", "cox", "mgaussian"),
                         init = NULL, weights = NULL, offset = NULL, 
                         alpha = 1.0, nlambda = 100, 
                         lambda.min.ratio = ifelse(nobs<nvars, 1e-2, 1e-4),
                         lambda = NULL, standardize = FALSE, intercept = TRUE,
                         thresh = 1e-7, dfmax = nvars+1, 
                         pmax = min(dfmax*2+20, nvars), exclude = NULL, 
                         penalty.factor=rep(1, nvars), lower.limits = -Inf,
                         upper.limits = Inf, maxit = 100000,
                         type.gaussian = ifelse(nvars<500, "covariance", "naive"),
                         type.logistic = c("Newton", "modified.Newton"),
                         standardize.response = FALSE,
                         type.multinomial = c("ungrouped", "grouped"), 
                         relax = FALSE, trace.it = 0, ...) {
  
  this.call = match.call()
  ### Need to do this first so defaults in call can be satisfied
  np = dim(x)
  ##check dims
  if(is.null(np) | (np[2] <= 1)) 
    stop("x should be a matrix with 2 or more columns")
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])
  ##check for NAs
  if(any(is.na(x)))
    stop("x has missing values; consider using makeX() to impute them")
  if(is.null(weights)) weights = rep(1, nobs)
  else if(length(weights) != nobs)
    stop(paste("number of elements in weights (", length(weights), 
               ") not equal to the number of rows of x (", nobs, ")", sep = ""))
  if(is.function(exclude))
    exclude <- check.exclude(exclude(x = x, y = y, weights = weights), nvars)
  if (length(penalty.factor) != nvars)
    stop("the length of penalty.factor does not match the number of variables")
  
  ### See whether its a call to glmnet or to glmnet.path, based on family arg
  if(!is.character(family)) {
    ## new.call=this.call
    ## new.call[[1]]=as.name("glmnet.path")
    ## fit=eval(new.call, parent.frame())
    fit = glmnet.path.large(x, y, init, weights, lambda, nlambda, 
                            lambda.min.ratio, alpha, offset, family, 
                            standardize, intercept, thresh = thresh, maxit, 
                            penalty.factor, exclude, lower.limits, upper.limits, 
                            trace.it = trace.it)
    fit$call = this.call
  } else {
    family = match.arg(family)
    if (family == "cox" && use.cox.path(x, y)) {
      # we should call the new cox.path()
      fit <- cox.path(x, y, weights, offset, alpha, nlambda, lambda.min.ratio,
                      lambda, standardize, thresh, exclude,penalty.factor,
                      lower.limits, upper.limits, maxit, trace.it, ...)
      fit$call <- this.call
    } else {
      ### Must have been a call to old glmnet
      ### Prepare all the generic arguments, then hand off to family functions
      if(alpha > 1){
        warning("alpha >1; set to 1")
        alpha = 1
      }
      if(alpha < 0) {
        warning("alpha<0; set to 0")
        alpha = 0
      }
      alpha = as.double(alpha)
      nlam = as.integer(nlambda)
      y = drop(y) # we dont like matrix responses unless we need them
      dimy = dim(y)
      nrowy = ifelse(is.null(dimy), length(y), dimy[1])
      if(nrowy != nobs)
        stop(paste("number of observations in y (", nrowy, 
                   ") not equal to the number of rows of x (", nobs, ")", 
                   sep=""))
      vnames = colnames(x)
      if(is.null(vnames)) vnames = paste("V", seq(nvars), sep="")
      ne = as.integer(dfmax)
      nx = as.integer(pmax)
      if(is.null(exclude)) exclude = integer(0)
      if(any(penalty.factor == Inf)) {
        exclude = c(exclude, seq(nvars)[penalty.factor == Inf])
        exclude = sort(unique(exclude))
      }
      if(length(exclude) > 0) {
        jd = match(exclude, seq(nvars), 0)
        if(!all(jd > 0)) stop("Some excluded variables out of range")
        penalty.factor[jd] = 1 #ow can change lambda sequence
        jd = as.integer(c(length(jd), jd))
      } else jd = as.integer(0)
      vp = as.double(penalty.factor)
      internal.parms = glmnet.control()
      if(internal.parms$itrace) trace.it = 1
      else{
        if(trace.it) {
          glmnet.control(itrace = 1)
          on.exit(glmnet.control(itrace = 0))
        }
      }
      ###check on limits
      if(any(lower.limits > 0)) {stop("Lower limits should be non-positive")}
      if(any(upper.limits < 0)) {stop("Upper limits should be non-negative")}
      lower.limits[lower.limits == -Inf] = -internal.parms$big
      upper.limits[upper.limits == Inf] = internal.parms$big
      if(length(lower.limits) < nvars) {
        if(length(lower.limits) == 1) lower.limits = rep(lower.limits, nvars)
        else stop("Require length 1 or nvars lower.limits")
      }
      else lower.limits=lower.limits[seq(nvars)]
      if(length(upper.limits) < nvars){
        if(length(upper.limits) == 1) upper.limits = rep(upper.limits, nvars)
        else stop("Require length 1 or nvars upper.limits")
      }
      else upper.limits = upper.limits[seq(nvars)]
      cl = rbind(lower.limits, upper.limits)
      if(any(cl == 0)) {
        ###Bounds of zero can mess with the lambda sequence and fdev; ie nothing 
        ###happens and if fdev is not zero, the path can stop
        fdev = glmnet.control()$fdev
        if(fdev != 0) {
          glmnet.control(fdev = 0)
          on.exit(glmnet.control(fdev = fdev))
        }
      }
      storage.mode(cl) = "double"
      ### end check on limits
      
      isd = as.integer(standardize)
      intr = as.integer(intercept)
      if(!missing(intercept) && family == "cox") 
        warning("Cox model has no intercept")
      jsd = as.integer(standardize.response)
      thresh = as.double(thresh)
      if(is.null(lambda)) {
        if(lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
        flmin = as.double(lambda.min.ratio)
        ulam = double(1)
      }
      else{
        flmin = as.double(1)
        if(any(lambda < 0)) stop("lambdas should be non-negative")
        ulam = as.double(rev(sort(lambda)))
        nlam = as.integer(length(lambda))
      }
      is.sparse = FALSE
      ix = jx = NULL
      if(inherits(x, "sparseMatrix")) {##Sparse case
        is.sparse = TRUE
        x = as(x, "CsparseMatrix")
        x = as(x, "dMatrix")
        #        x=as(x,"generalMatrix")
        ix = as.integer(x@p + 1)
        jx = as.integer(x@i + 1)
        
        # TODO: changed everything except cox to C++ implementation.
        # C++ version takes xd as the dgCMatrix itself.
        # Fortran requires xd to be the compressed data vector.
        if (family != "cox") {
          xd <- x
        } else {
          xd <- x@x
        }
        
      } else if (!inherits(x, "matrix")) {
        xd <- data.matrix(x)
      } else {
        xd <- x
      }
      # TODO: only coerce if xd is not sparse
      if(!inherits(xd, "sparseMatrix")) {
        storage.mode(xd) <- "double"
      }
      if (trace.it) {
        if (relax) cat("Training Fit\n")
        pb  <- createPB(min = 0, max = nlam, initial = 0, style = 3)
      } else {
        pb <- NULL # dummy initialize (won't be used, but still need to pass)
      }
      kopt=switch(match.arg(type.logistic),
                  "Newton" = 0,#This means to use the exact Hessian
                  "modified.Newton" = 1 # Use the upper bound
      )
      if(family == "multinomial") {
        type.multinomial = match.arg(type.multinomial)
        if(type.multinomial == "grouped") kopt=2 #overrules previous kopt
      }
      kopt = as.integer(kopt)
      
      fit = switch(family,
                   "gaussian" = elnet(xd, is.sparse, y, weights, offset, 
                                      type.gaussian, alpha, nobs, nvars, jd, vp,
                                      cl, ne, nx, nlam, flmin, ulam, thresh, 
                                      isd, intr, vnames, maxit, pb),
                   "poisson" = fishnet(xd, is.sparse, y, weights, offset, alpha,
                                       nobs, nvars, jd, vp, cl, ne, nx, nlam, 
                                       flmin, ulam, thresh, isd, intr, vnames,
                                       maxit, pb),
                   "binomial" = lognet(xd, is.sparse, ix, jx, y, weights, 
                                       offset, alpha, nobs, nvars, jd, vp, cl, 
                                       ne, nx, nlam, flmin, ulam, thresh, isd,
                                       intr, vnames, maxit, kopt, family, pb),
                   "multinomial" = lognet(xd, is.sparse, ix, jx, y, weights,
                                          offset, alpha, nobs, nvars, jd, vp, 
                                          cl, ne, nx, nlam, flmin, ulam, thresh,
                                          isd, intr, vnames, maxit, kopt, 
                                          family, pb),
                   "cox" = coxnet(xd, is.sparse, ix, jx, y, weights, offset, 
                                  alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam,
                                  flmin, ulam, thresh, isd, vnames, maxit),
                   "mgaussian" = mrelnet(xd, is.sparse, y, weights, offset, 
                                         alpha, nobs, nvars, jd, vp, cl, ne, nx,
                                         nlam, flmin, ulam, thresh, isd, jsd, 
                                         intr, vnames, maxit, pb)
      )
      if (trace.it) {
        utils::setTxtProgressBar(pb, nlam)
        close(pb)
      }
      if(is.null(lambda)) fit$lambda = fix.lam(fit$lambda) ##first lambda is infinity; changed to entry point
      fit$call = this.call
      fit$nobs = nobs
      class(fit) = c(class(fit), "glmnet")
    }
  }
  
  if(relax)
    relax.glmnet(fit, x = x, y = y, weights = weights, offset = offset,
                 lower.limits = lower.limits, upper.limits = upper.limits,
                 penalty.factor = penalty.factor, check.args=FALSE,...)
  else
    fit
}

