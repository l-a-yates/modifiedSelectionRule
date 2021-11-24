#------------------------------------------------------------
#
# Goby survival models
# First example from "Parsimonious Model Selection using Information Theory: a Modified Selection Rule"
# Authors: Yates, L.A, Richards, S.A, Brook, B.W.
# 
# This script contains function definitions to:
#   - Specify and fit Goby survival models, following Bolker(2008);
#   - Compute model scores using leave-one-out cross-validation;
# 
# Date: Original version Dec 2019; this version Sept 2020
# Data: This data is due originally to Wilson(2004)
#
# Please see manuscript for further details and references
#------------------------------------------------------------

# N.B. maxCores variable is set globally in main script

# The following three functions are adapted from Chapter 8, Bolker (2008)
# See also https://ms.mcmaster.ca/~bolker/emdbook/chap8.R (chunks 46 - 60)

# computes log-Likelihood for full model.
# called from GS.fit.all() which supplies the data within its environment
nll.xqdi = function(lscale0,lscale.q,lscale.d,lscale.i,
                       lscale.x2,lscale.x3,lscale.x4,lscale.x5,
                       lshape) {
  lscalediff = c(0,lscale.x2,lscale.x3,lscale.x4,lscale.x5)
  scale=exp(lscale0+lscalediff[exper]+
              lscale.q*qual+(lscale.d+lscale.i*qual)*density)
  shape = exp(lshape)
  suppressWarnings(-sum(log(pweibull(day2,shape,scale)-
             pweibull(day1,shape,scale))))
}

# computes predictive, negative log-likelihood for a fitted model (mle.fit) and new data (data)
# used for cross validation
nll.pred = function(mle.fit,data) {
  with(data, with(as.list(coef(mle.fit)),{
        lscalediff = c(0,lscale.x2,lscale.x3,lscale.x4,lscale.x5)
        scale=exp(lscale0+lscalediff[exper]+lscale.q*qual+(lscale.d+lscale.i*qual)*density)
        shape = exp(lshape)
        -sum(log(pweibull(day2,shape,scale)- pweibull(day1,shape,scale)))
      }
      ))
}

# specifies and fits all ten models to the supplied data
# input: data
# output: list of fitted objects (mle2)
fit.goby <- function(data){
  startvals.GS = list(lscale0=log(mean((data$d1+data$d2)/2)),
                      lscale.x2=0,lscale.x3=0,lscale.x4=0,lscale.x5=0,
                      lscale.q=0,lscale.d=0,lscale.i=0,lshape=0)
  GS <- list()
  GS$xqdi = mle2(nll.xqdi,startvals.GS, data = data)
  GS$qdi <-  mle2(nll.xqdi,startvals.GS,
                  fixed=list(lscale.x2=0,lscale.x3=0,lscale.x4=0,lscale.x5=0), data = data)
  GS$qd <-  mle2(nll.xqdi,startvals.GS, 
                 fixed=list(lscale.i=0,lscale.x2=0,lscale.x3=0,lscale.x4=0,lscale.x5=0), data = data)
  GS$xq <-  mle2(nll.xqdi,startvals.GS,
                 fixed=list(lscale.i=0,lscale.d=0), data = data)
  GS$d <-  mle2(nll.xqdi,startvals.GS,
                fixed=list(lscale.i=0,lscale.x2=0,lscale.x3=0,lscale.x4=0, lscale.x5=0,lscale.q=0), data = data)
  GS$q <-  mle2(nll.xqdi,startvals.GS,
                fixed=list(lscale.i=0, lscale.x2=0,lscale.x3=0,lscale.x4=0,lscale.x5=0,lscale.d=0), data = data)
  GS$zero <-  mle2(nll.xqdi,startvals.GS,
                   fixed=list(lscale.i=0,lscale.x2=0,lscale.x3=0,lscale.x4=0,lscale.x5=0,lscale.d=0,lscale.q=0), data = data)
  GS$xd <-  mle2(nll.xqdi,startvals.GS, fixed=list(lscale.i=0,lscale.q=0), data = data)
  GS$x <-  mle2(nll.xqdi,startvals.GS, fixed=list(lscale.i=0,lscale.q=0,lscale.d=0), data = data)
  GS$xqd <-  mle2(nll.xqdi,startvals.GS,fixed=list(lscale.i=0), data = data)
  return(GS)
}

# computes the leave-one-out cross-validated deviance for each data point for each model
# input: data
# output: tibble of CV values with one row per datum and one column per model
GS.loo <- function(data){
  pbmclapply(1:nrow(data), function(i){
    fit.train <- fit.goby(data[-i,])
    sapply(fit.train, function(m) nll.pred(m, data[i,]))*2
  }, mc.cores = MAX_CORES) %>% sapply(identity) %>% t %>%  as_tibble
 }

# calculates the deviance (-2 x log density) for non-parametric bootstrap samples (all models)
# input: data; number of bootstrap samples
# output: a tibble of deviance values with one row per sample and one column per model
GS.boot <- function(data, nboot){
  pbmclapply(1:nboot, function(i){
    fitted <- F
    while (!fitted){
      GS.boot = try(fit.goby(sample_n(data,nrow(data), replace = T)))
      fitted = sum(sapply(GS.boot, is.character))==0
    }
    sapply(GS.boot, function(x) logLik(x)*-2)
  }, mc.cores = MAX_CORES) %>% sapply(identity) %>% t %>%  as_tibble
}
