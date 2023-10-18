# this script defines functions for calculating mean GFP fluorescence for 
# isogenic yeast strains measured using flow cytometry

# install and load packages
packages <- c("flowCore")
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages[installed_packages], library, character.only = TRUE))


# functions

# Extracts log-normalized GFP values from FCS files.
# mCherryfluor and GFPfluor provide strings with which to identify channels for
# the corresponding fluorescent proteins. If mCherryfluor is NULL, then mCherry-
# positive cells are not filtered out. If not NULL, then mCherry-positive cells
# will be identified and filtered out.
# If norm_method = 1, then fluorescence is normalized by cell area (FSC_A). If
# norm_method = 1.5, then fluorescence is normalized by cell volume (FCS_A^1.5).
extract_GFP <- function(fcs_dir, mCherryfluor = "PEDazzle594", 
                        GFPfluor = 'FITC', norm_method = 1.5) {
  
  # read in FCS data
  fcs <- flowCore::read.FCS(fcs_dir)
  data <- flowCore::exprs(fcs)
  
  # extract GFP fluorescence
  GFP <- data[, grepl(GFPfluor, colnames(data)) & grepl('-A', colnames(data))]
  FSC_A <- data[, 'FSC-A']
  
  if(!is.null(mCherryfluor)) {
    # extract mCherry fluorescence
    mCherry <- data[, grepl(mCherryfluor, colnames(data)) & 
                      grepl('-A', colnames(data))]
    
    # remove cells that are at the upper and lower bounds of the FSC-A range
    filter_out <- (FSC_A == 2 ^ 18 - 1) | (GFP <= 0) | (mCherry <= 0)
    mCherry <- mCherry[!filter_out]
    GFP <- GFP[!filter_out]
    FSC_A <- FSC_A[!filter_out]
    
    # fitting a mixture of two logistic distributions to mCherry
    norm_mCherry <- log(mCherry / FSC_A, 10)
    param <- infer(norm_mCherry)
    
    # identifying the threshold of mCherry over which cells with such values are 
    # more than 99.9% likely to be active
    f <- function(x) {
      abs(dlogis(x, param[1], param[2], TRUE) - 
            dlogis(x, param[3], param[4], TRUE) - log(0.999))
    }
    
    threshold <- optimize(f, interval = c(param[3], param[1]))$minimum[1]
    
    # normalize GFP by FSC_A^norm_method
    norm_GFP <- log(GFP / FSC_A ^ norm_method, 10)
    
    return(norm_GFP[norm_mCherry > threshold])
    
  } else {
    # remove cells that are at the upper and lower bounds of the FSC-A range
    filter_out <- (FSC_A == 2 ^ 18 - 1) | (GFP <= 0)
    GFP <- GFP[!filter_out]
    FSC_A <- FSC_A[!filter_out]
    
    # normalize GFP by FSC_A^norm_method
    norm_GFP <- log(GFP / FSC_A^norm_method, 10)
    
    return(norm_GFP)
  }
}


# calculates the log-likelihood of a mixture of two logistic distributions given 
# data and distribution parameters
calc_ll <- function(x, distr_param, null_distr_param, lambda) {
  # checking parameter range:
  # if mixture value is non-sensical
  if(lambda < 0 || lambda > 1) return(NaN) 
  
  # if dispersion value is negative
  if(distr_param[2] < 0 || null_distr_param[2] < 0) return(NaN) 
  
  # calculating log-likelihood
  sum(log(lambda * dlogis(x, distr_param[1], distr_param[2]) +
            (1 - lambda) * dlogis(x, null_distr_param[1], null_distr_param[2])))
}


# infers the mixture two-logistic distribution that best fits the data
infer <- function(x) {
  
  # setting initial parameter values for optimization
  
  param_init <- rep(0.5, 5)
  
  param_init[1] <- mean(x[x > quantile(x, 0.75, na.rm=TRUE)])
  param_init[2] <- sd(x[x >quantile(x, 0.75, na.rm=TRUE)]) * sqrt(3) / pi
  param_init[3] <- mean(x[x <quantile(x, 0.25, na.rm=TRUE)])
  param_init[4] <- sd(x[x < quantile(x, 0.25, na.rm=TRUE)]) * sqrt(3) / pi
  
  # inference
  
  f <- function(arg) calc_ll(x, arg[1:2], arg[3:4], arg[5])
  param_inf <- optim(param_init, f, control = list(fnscale = -1))
  res <- c(param_inf$par, param_inf$value) # the last value is the total log-likelihood
  
  # ordering the parameters
  
  if(res[3] > res[1]) res <- c(res[3:4], res[1:2], res[5:6])
  
  res
}
