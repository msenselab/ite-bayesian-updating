#' refine Fredrik's code for factorial modeling and computational efficiency. 
#' Date: 3rd August, 2017 by Strongway

library(tidyverse)
library(parallel)
library(doParallel)
#' Approximate pdf for the Wiener diffusion model 
#' 
#' (from Lee, Fuss & Navarro, 2007)
#' @param rt the reaction time, 
#' @param theta is the threshold, 
#' @param mu is the drift rate, 
#' @param delta delta is the non-decision time (or ter)
ddiffusion <- function(rt,theta,delta,mu) { 
  1/sqrt(2*pi)*(2*theta+mu*(rt-delta))/(2*(rt-delta)^1.5)*exp(-(2*theta-mu*(rt-delta))^2/(2*(rt-delta))) }


#' Logarithm of the DDM PDF
#'
#' optimized compared to using log(ddiffusion) by cancelling the exponential against the logarithm
#' @param rt the reaction time, 
#' @param theta is the threshold, 
#' @param mu is the drift rate, 
#' @param delta delta is the non-decision time (or ter)
log_ddiffusion <- function(rt,theta,delta,mu) {  log((2*theta+mu*(rt-delta))/(2*(rt-delta)^1.5)) - log(sqrt(2*pi)) - (2*theta-mu*(rt-delta))^2/(2*(rt-delta)) }

#' Recinormal PDF
#' 
#' @param rt the reaction time
#' @param mu the reciprocal mu
#' @param sig the reciprocal sig
drecinorm <- function(rt, mu, sig){
  1/(rt*rt*sqrt(2*pi*sig*sig))*exp(-(mu*rt-1)^2/2/sig^2/rt^2)
}

#' so the Bayesian updating framework can be separated as follows:
#' 1. Calculate trial-wise prior and posterior updates
#'    two methods: simulation, and analytical approach with hyper parameter integration 
#'    This includes memory component (forgetting): $u_i = (1-m) u_0 + m * u_{i-1}$ (similar to Kalman filter)
#' 2. Estimate starting point $S_0$ based one logPrior
#' 3. Estimate likelihood of all trials for given parameters and models

#' prior updates with parameter integration
#' 
#' Updating priors based on Beta distribution
#' When m = 1, keep all updates (full memory). Otherwise, partial leakage
#' In this model, we assume participants has some knowlege of the mean of the prior (i.e., mu)
#' But it is somehow forgetting or partially integrated. 
#' Participants still use some original prior (or partially integrated). Mathematically:
#' mu_update = (1-m) * mu_0 + m * mu
#' Another parameter - sample size v remains. Thus, Beta distribution parameters alpha and beta
#' can be expressed as : alpha = mu * v, beta = (1-u)v
#' It turns out that two beta parameters are updated as follows:
#' alpha_update = (1-m)* alpha_0 + m * alpha
#' beta_update = (1-m)*beta_0 + m * beta
#' 
#' The above approach is very similar to the simulation approach where two prior distributions are fused. 
#' It is a normal approximation of fusion of two beta distribution. 
#' @param targets An array indicate target or non-target. logical FALSE, TRUE
#' @param a Beta distribution parameter a
#' @param b Beta distribution parameter b
#' @param m memory component. When m = 0, no updating, i.e., complete ignore prior history. 
updatePriorA <- function(targets, a = 1, b = 1, m = 0.5){
  a0 = a
  b0 = b
  # with memory approach
  mu = rep(0,length(targets))
  mu[1] = a/(a+b) # initial prior: equal (a/(a+b))
  for(idx in 2: (length(targets))) {
    # update parameters after trial n
    if (targets[idx-1]){
      a = a + 1
    } else {
      b = b + 1
    }
    mu_cur = a/(a+b) # current prior
    # partial reset (forgetting)
    mu[idx] = (1-m)*mu[1] + mu_cur*m
    
    # similarly, update a and b 
    a = (1-m)*a0 + m*a
    b = (1-m)*b0 + m*b
  }
  return(mu)
}

#' prior updates - simulation approach
#' 
#' Updating priors based on Beta distribution
#' @param targets An array indicate target or non-target. It must be 0,1 or FALSE, TRUE
#' @param a Beta distribution parameter a
#' @param b Beta distribution parameter b
#' @param m memory component. When m = 0, no updating, i.e., complete ignore prior history. 
updatePrior <- function(targets, a = 1, b = 1, m = 0.5){
  x<-seq(0,1,0.001)
  p0<-dbeta(x,a,b) # assume a and b are the same initial
  p0 <- ifelse(p0==Inf, 20, p0) # replace Inf with 20 (this happens for Jeffrey prior)
  pd <- p0
  mu = rep(a/(a+b),length(targets)+1)
  idx = 2
  for(target in targets) {
    # Update the distribution based on the Bernoulli likelihood 
    pd <- pd*x*target + (1-x)*pd*(1-target)
    # Implement with 'forgetting'
    pd <- (1-m)*p0 + m * pd
    # Normalize the probability distribution
    pd<-pd/mean(pd) 
    mu[idx] <- mean(pd*x)
    idx  = idx + 1
  }
  # return priors
  return(mu[1:length(targets)])
}

#' update of the drift rate
#' 
#' Updating drift rate so that it is reduced on dimension switch trials, this version has a memory of one trial back and doesn't
#' care whether there has been many repeats/switches in a row. 
#' @param intertrial An array that indicate dimension repeat (1), switch (0) or target absent (NA)
#' @param sc Scaling parameter, should be <= 1: the drift rate is scaled down by this amount after a dimension switch
#' When m = 1, keep all updates (full memory). Otherwise, partial leakage
updateRate <- function(intertrial, sc) {
  scale <- intertrial + (1-intertrial)*sc
  scale[is.na(scale)] <- 1
  return(scale)
}

#' update of the drift rate - with memory
#' 
#' Updating drift rate so that it is reduced on dimension switch trials, this version has a memory of more than one trial back
#' the scaling factor is reduced after a switch and increased after a repeat, but with exponential discounting of old trials 
#' @param intertrial An array that indicate dimension repeat (1), switch (0) or target absent (NA)
#' @param mem memory parameter, m=0 means no updating, m=1 means old trials matter as much as recent trials
#' @param delta determines the amount that the scaling factor is changed by on each trial
updateRateMem <- function(intertrial, mem, delta) {
  intertrial <- 2*(intertrial-0.5) # rescale intertrial so that a switch is represented by -1 instead of 0
  intertrial[is.na(intertrial)] <- 0 # Target absent trials are neither a dimension switch nor a repeat
  s = rep(1, length(intertrial))
  for (i in 2:length(intertrial)) {
    s[i] = max(0.001, s[i-1] + intertrial[i]*delta) # make sure it is positive
    s[i] = s[i]*mem + (1-mem)
  }
  s[intertrial==0] <- 1 # Dimension switch costs don't affect target absent trials
  return(s)
}

#' update of the drift rate - dimension weighting
#' 
#' Updating drift rate so that it is reduced on dimension switch trials, this version attempts to do something closer to the spirit
#' of dimension weighting by updating based on the dimension itself rather than whether there has been a repeat or switch: 
#' ds reflects how much "weight" has been shifted from orientation to color, consequently the rate is scaled by s+ds on a color trial
#' and by s-ds on an orientation trial.
#' 
#' @param dim An array that indicate the target dimension color, orientation or absent (NA)
#' @param mem memory parameter, m=0 means no updating, m=1 means old trials matter as much as recent trials
#' @param delta determines the amount that the scaling factor is changed by on each trial
updateRateWeight <- function(dim, mem, delta) {
  dim[dim>0] <- 2*(dim[dim>0]-1.5) # rescale intertrial so that color orientation and absent are -1, 1 and 0
  ds = rep(0, length(dim))
  s = rep(1, length(dim))
  for (i in 2:length(dim)) {
    ds[i] = ds[i-1] + dim[i-1]*delta
    ds[i] = ds[i]*mem 
    s[i]<-s[i]+dim[i]*ds[i]
    s[i]=max(0,s[i])
  }
  return(s)
}


#' Estimate loglikelihood
#'
#' Estimate log likelihood for a given model (theta, mu, sig)
#' @param data input data, which should have one column rt, one column of updated prior, called beta, and one column target
#' @param theta decision boundary (theta1, theta2)
#' @param mu drift rate / ergodic rate (mu1, mu2)
#' @param sig sigma of the drift rate (sig1, sig2)
#' @param ter nondecision time (assume it is independent from response/target condition)
#' @param later_diffusion a flag for using LATER or Diffusion model as RT distribution
logll <- function(data, theta=4, mu=5,  sig=1, ter  = 0, later_diffusion =1) {
  
  # data <-   data %>% ungroup(.) %>% filter(!(error | outlier)) %>%
  #   mutate( #beta = max(0.001, min(0.999, beta)), # constrain beta in 0.01 - 0.99  # min/max not working in mutate!!!!
  #     s0 = log(beta/(1-beta)) * ((targets>0) *2 -1),
  #     mu = mu*scale,  sig = sig,  theta = theta,
  #     rs = 1/(rt - ter), # add non-decision time for LATER model
  #     delta =  theta-s0) 
  # LATER model  
  rt = data$rt
  scale = data$scale
  delta = theta - data$s0
  if(later_diffusion==1) {
    # use recinormal pdf (i.e., rt pdf, comparible to DDM)
    nll = -sum(log(drecinorm(rt-ter, mu*scale/delta, sig/delta)))
    # nll = -sum(log(dnorm(1/(rt - ter), mu*scale/delta, sig/delta)))
  } else {
    # put back scaling parameter sig back, given that delta is not scaled
    nll = -sum(log_ddiffusion(rt, delta/sig, ter, mu/sig*scale))
  }
  if (nll == Inf | nll == -Inf) # avoid optim error with L-BFGS-B method
    nll = 1e10
  return(nll)
#  if(later_diffusion==1) {
#    data %>% mutate(delta = theta - s0, ll =  log(dnorm(1/(rt - ter), mu*scale/delta, sig/delta))) %>%
#      summarise(logll = -sum(ll))
#  } else {
#    data %>% mutate(delta = theta - s0, ll = log_ddiffusion(rt, delta, ter, mu*scale)) %>%
#      summarise(logll = -sum(ll))    
#  }

}

# Old version used in first round of optimization
findParameters_V1 <- function(par1, data, fixed, op = TRUE) {
  tryCatch({
    ll<-0 # log likelihood
    # recombine to-be-fitted parameterd and fixed parameters
    x = fixed
    x[is.na(fixed)] = par1
    # parameters for outer fit: m = x1, beta_a = x2, scale x3, um = x4, u_delta = x5, ter = x6,
    # flags: later_diffusion = x7, s0_update = x8, u_update = x9
    # 
    # update starting point
    switch(x[8],
           data$beta <- 0.5, # no update 1
           data$beta <- updatePrior(data$targets > 0, a = x[2], b=x[2], m = x[1]), # simulate
           data$beta <- updatePriorA(data$targets > 0, a = x[2], b=x[2], m = x[1]) # analytical
    )
    # update the rate
    switch(x[9], # 1: no update, 2: step update, 3: drift
           data$scale <- 1,
           data$scale <- updateRate(data$inttrial, x[3]),
           data$scale <- updateRateMem(data$inttrial, x[4], x[5]) 
    )
    
    # later or diffusion
    # calculate log likelihood based how many diffusion /later processes 
    # the targets has implicit indication of process (e.g., targets contains 0,1,2 indicates 3 processes)
    # each find optimal parameters here (inner optimization)
    # Inner parameters: theta=x[1], mu=x[2],  sig=x[3], ter  = x[4] , later_diffusion =x[5]
    pars = x # parameters
    inner_par0 = c(theta = 4, mu = 5, sig = 1, ter = x[6], later_diffusion = fixed[7])
    inner_fixed = c(NA, NA, NA, x[6], fixed[7]) 
    
    data <- data %>% ungroup(.) %>% filter(!(error | outlier)) %>%
      mutate( s0 = log(beta/(1-beta)) * ((targets>0) *2 -1)) 
    # add contraints for parameters - theta, mu, sig, ter
    # not working at the moment: return L-BFGS-B needs finite values of 'fn'
    #i_lowers = c(max(abs(data$s0)), 0.001, 0.001) #theta must be greater than s0
    #i_uppers = c(50, 50, 50) 
    ui1 = rbind(c(1,0,0), c(0,1,0),c(0,0,1))
    ci1 = c(max(abs(data$s0)), 0.001, 0.001)
    ui2 = -ui1
    ci2 = c(-50, -50, -50)
    ui = rbind(ui1, ui2)
    ci = c(ci1, ci2)
    
    for (itarget in sort(unique(data$targets))){
      subdata <- data %>% filter(targets == itarget)
      inner_par1 <- inner_par0[is.na(inner_fixed)] # find to-be-fixed parameters
      par <- constrOptim(inner_par1,  findInnerParameters, NULL,
                         ui = ui, ci = ci, 
                         data = subdata, fixed = inner_fixed)
      #      par <- optim(inner_par1,  findInnerParameters, data = subdata, fixed = inner_fixed)
      #      par <- optim(inner_par1,  findInnerParameters, data = subdata, fixed = inner_fixed,
      #                   lower = i_lowers, upper = i_uppers, method = 'BFGS')
      ll <- ll + par$value
      pars = c(pars, par$par)
    }
    pars = c(pars, nll = ll)
    if (op)
      return(ll)
    else
      return(pars)
  },
  error = function(e){
    print(e)
    stop(e)
  })
}

#' for given parameters calculate negative log likelihood
#' 
#' @param x to-be-optimzed parameter
#' @param data data frame contain RT data
#' @param fixed pass fixed parameters here. Same length as x. NA for to-be-optimzed parameter
#' @param op logical flag for optimization or returning parameters
#' @return negloglikelihood return negative log likelihood
findParameters <- function(par1, data, fixed, op = TRUE) {
  tryCatch({
    ll<-0 # log likelihood
    # recombine to-be-fitted parameterd and fixed parameters
    x = fixed
    x[is.na(fixed)] = par1
    data$beta <- 0.5
    data$scale <-1

    # Response based updating 
    switch(x['resp_update'], # If 0 do nothing
           data$beta <- updatePrior(data$targets > 0, a = x["beta_resp"], b=x["beta_resp"], m = 1), # 1: S0 updating with full memory
           data$beta <- updatePrior(data$targets > 0, a = x["beta_resp"], b=x["beta_resp"], m = x["mem_resp"]), # 2: S0 updating with forgetting
           data$scale <- updateRate(data$inttrial_resp, sc=x["scale_resp"]), # 3: single trial back switch cost on DR
           data$scale <- updateRateMem(data$inttrial_resp, mem=x["u_mem_resp"], delta=x["u_delta_resp"]), # 4: switch cost on DR with longer memory
           data$scale <- updateRateWeight(as.numeric(data$targets>0)+1, mem=x["u_mem_resp"], delta=x["u_delta_resp"]) # 5: response based DR weighting 
    )
    # Dimension based updating
    switch(x['dim_update'], # If 0 do nothing
           {
             dims <- c(1, 2)
             for(dim in dims) {
               curdim <- as.numeric(data$dimension)==dim # 1 if current trial has the dimension dim, 0 otherwise
               beta_update = updatePrior(curdim, a = x["beta_dim"], b = x["beta_dim"], m = 1) - 0.5 # how beta differs from "neutral"
               data[curdim,"beta"] <- data[curdim,"beta"] + beta_update[curdim] # Use these changes to beta only on trials where the dimension is dim
             }
           }, # 1: S0 updating with full memory
           {
             dims <- c(1, 2)
             for(dim in dims) {
               curdim <- as.numeric(data$dimension)==dim # 1 if current trial has the dimension dim, 0 otherwise
               beta_update = updatePrior(curdim, a = x["beta_dim"], b = x["beta_dim"], m = x["mem_dim"]) - 0.5 # how beta differs from "neutral"
               data[curdim,"beta"] <- data[curdim,"beta"] + beta_update[curdim] # Use these changes to beta only on trials where the dimension is dim
             }
           }, # 2: S0 updating with forgetting
           data$scale <- data$scale*updateRate(data$inttrial_dim, sc=x["scale_dim"]), # 3: single trial back switch cost on DR
           data$scale <- data$scale*updateRateMem(data$inttrial_dim, mem=x["u_mem_dim"], delta=x["u_delta_dim"]), # 4: switch cost on DR with longer memory
           data$scale <- data$scale*updateRateWeight(as.numeric(data$dimension), mem=x["u_mem_dim"], delta=x["u_delta_dim"]) # 5: response based DR weighting 
    )
    
    # later or diffusion
    # calculate log likelihood based how many diffusion /later processes 
    # the targets has implicit indication of process (e.g., targets contains 0,1,2 indicates 3 processes)
    # each find optimal parameters here (inner optimization)
    # Inner parameters: theta=x[1], mu=x[2],  sig=x[3], ter  = x[4] , later_diffusion =x[5]
    pars = x # parameters
    # Ran exps 1-2 with these parameters:
    #inner_par0 = c(theta = 4, mu = 5, sig = 1, ter = x['ter'], later_diffusion = fixed['later_diffusion']) 
    # Ran exp 3 with these parameters because fit (subject "fra") failed with theta = 4, mu = 5:    
    inner_par0 = c(theta = 3, mu = 6, sig = 1, ter = x['ter'], later_diffusion = fixed['later_diffusion'])
    inner_fixed = c(NA, NA, NA, x['ter'], fixed['later_diffusion']) 
    
    data <- data %>% ungroup(.) %>% filter(!(error | outlier)) %>%
      mutate( s0 = log(beta/(1-beta)) * ((targets>0) *2 -1)) 
    # add contraints for parameters - theta, mu, sig, ter
    # not working at the moment: return L-BFGS-B needs finite values of 'fn'
    #i_lowers = c(max(abs(data$s0)), 0.001, 0.001) #theta must be greater than s0
    #i_uppers = c(50, 50, 50) 
    ui1 = rbind(c(1,0,0), c(0,1,0),c(0,0,1))
    ci1 = c(max(abs(data$s0)), 0.001, 0.001)
    ui2 = -ui1
    ci2 = c(-50, -50, -50)
    ui = rbind(ui1, ui2)
    ci = c(ci1, ci2)
    
    if(sum(is.na(data$s0)>0))
    {
      ll <- 1e10
      return(ll)
    } else if(max(abs(data$s0)>inner_par0[1])) {
      ll <- 1e10
      return(ll)
    }
    
    for (itarget in sort(unique(data$targets))){
      subdata <- data %>% filter(targets == itarget)
      inner_par1 <- inner_par0[is.na(inner_fixed)] # find to-be-fixed parameters
      par <- constrOptim(inner_par1,  findInnerParameters, NULL,
                         ui = ui, ci = ci, 
                         data = subdata, fixed = inner_fixed)
#      par <- optim(inner_par1,  findInnerParameters, data = subdata, fixed = inner_fixed)
#      par <- optim(inner_par1,  findInnerParameters, data = subdata, fixed = inner_fixed,
#                   lower = i_lowers, upper = i_uppers, method = 'BFGS')
      ll <- ll + par$value
      pars = c(pars, par$par)
    }
    pars = c(pars, nll = ll)
    if (op)
      return(ll)
    else
      return(pars)
  },
  error = function(e){
    print(e)
    stop(e)
  })
}

# Generates a sequence of starting points and  drift rates
genSeq <- function(x, data) {
  
  data$beta <- 0.5
  data$scale <- 1
  
  switch(x$resp_update, # If 0 do nothing
         data$beta <- updatePrior(data$targets>0, a = x[["beta_resp"]], b=x[["beta_resp"]], m = 1), # 1: S0 updating with full memory
         data$beta <- updatePrior(data$targets>0, a = x[["beta_resp"]], b=x[["beta_resp"]], m = x[["mem_resp"]]), # 2: S0 updating with forgetting
         data$scale <- updateRate(data$inttrial_resp, sc=x[["scale_resp"]]), # 3: single trial back switch cost on DR
         data$scale <- updateRateMem(data$inttrial_resp, mem=x[["u_mem_resp"]], delta=x[["u_delta_resp"]]), # 4: switch cost on DR with longer memory
         data$scale <- updateRateWeight(as.numeric(data$target!="Absent")+1, mem=x[["u_mem_resp"]], delta=x[["u_delta_resp"]]) # 5: response based DR weighting 
  )
  # Dimension based updating
  switch(x$dim_update, # If 0 do nothing
         {
           dims <- c(1, 2)
           for(dim in dims) {
             curdim <- as.numeric(data$dimension)==dim # 1 if current trial has the dimension dim, 0 otherwise
             beta_update = updatePrior(curdim, a = x[["beta_dim"]], b = x[["beta_dim"]], m = 1) - 0.5 # how beta differs from "neutral"
             data[curdim,"beta"] <- data[curdim,"beta"] + beta_update[curdim] # Use these changes to beta only on trials were the dimension is dim
           }
         }, # 1: S0 updating with full memory
         {
           dims <- c(1, 2)
           for(dim in dims) {
             curdim <- as.numeric(data$dimension)==dim # 1 if current trial has the dimension dim, 0 otherwise
             beta_update = updatePrior(curdim, a = x[["beta_dim"]], b = x[["beta_dim"]], m = x[["mem_dim"]]) - 0.5 # how beta differs from "neutral"
             data[curdim,"beta"] <- data[curdim,"beta"] + beta_update[curdim] # Use these changes to beta only on trials were the dimension is dim
           }
         }, # 2: S0 updating with forgetting
         data$scale <- data$scale*updateRate(data$inttrial_dim, sc=x[["scale_dim"]]), # 3: single trial back switch cost on DR
         data$scale <- data$scale*updateRateMem(data$inttrial_dim, mem=x[["u_mem_dim"]], delta=x[["u_delta_dim"]]), # 4: switch cost on DR with longer memory
         data$scale <- data$scale*updateRateWeight(as.numeric(data$dimension), mem=x[["u_mem_dim"]], delta=x[["u_delta_dim"]]) # 5: response based DR weighting 
  )
  
  data$rate <- data$scale  
  data$theta <- 0 
  data$sigma <- 0
  dims <- unique(data$dimension)
  if (length(dims)==3) {
    data$dimension <- factor(data$dimension, levels=0:2, labels=c("Absent", "Color", "Orientation"))
    dims <- unique(data$dimension)
    for (dim in dims) {
      data[data$dimension==dim,]$rate <- switch(dim, "Absent"= x$mu1, "Color" = x$mu2,  "Orientation" = x$mu3)*data[data$dimension==dim,]$rate
      data[data$dimension==dim,]$theta <- switch(dim, "Absent"= x$theta1, "Color" = x$theta2,  "Orientation" = x$theta3)
      data[data$dimension==dim,]$sigma <- switch(dim, "Absent"= x$sig1, "Color" = x$sig2,  "Orientation" = x$sig3)
    }
  } else {
    for (dim in dims) {
      data[data$dimension==dim,]$rate <- switch(dim, "Color" = x$mu2,  "Orientation" = x$mu1)*data[data$dimension==dim,]$rate
      data[data$dimension==dim,]$theta <- switch(dim, "Color" = x$theta2,  "Orientation" = x$theta1)
      data[data$dimension==dim,]$sigma <- switch(dim, "Color" = x$sig2,  "Orientation" = x$sig1)
    }
  }
  data<-mutate(data, s0=log(beta/(1-beta)) * ((targets>0) *2 -1)) 
  
  t <- seq(0, 3, 0.01)
  delta <- data$theta - data$s0
  sigma <- data$sigma
  mu <- data$rate 
  predRT <- rep(0, nrow(data))
  for(n in 1:nrow(data)) {
    if(x$later_diffusion=="LATER") {
      dist <- drecinorm(t-x$ter, mu[n]/delta[n], sigma[n]/delta[n])
    } else {
      dist <- ddiffusion(t, delta[n]/sigma[n], x$ter, mu[n]/sigma[n])
    }
    predRT[n] <- mean(t*dist, na.rm=TRUE)/mean(dist, na.rm=TRUE)
  }
  data$predRT <- predRT
  
  return(data)
}

# return inner LATER/DIFFUSION parameters
findInnerParameters <- function(par, data, fixed){
  x = fixed
  x[is.na(fixed)] = par
  logll(data, theta=x[1], mu=x[2],  sig=x[3], ter  = x[4] , later_diffusion =x[5]) 
}

#' fit parameters with given models (old version used in the first round of optimization, for the current version see fitPar below)
#' 
#' This function will find optimal parameters for a given model and give subject data
#' @param sub subject rt data frame
#' @param modelname model names. 
#' model name rules: Before '_' indicates models, after '_' indicates parameters
#' a parameter followed by a number (eg T0) means that the parameter is fixed to have that value
#' a parameter followed by V means the parameter is free to vary in the optimization
#' e.g., LS_MVB1T0R2 means later model using simulation, memory parameter is free (MV) but prior is fixed to beta(1,1) (B1) 
#' and non-decision time is fixed to 0 (T0), rate is updated using the step model (R2)
#' Note that non-updating models (LN or DN) do not use the M and B parameters 
fitParV1 <- function(sub, modelname = 'LS_MVBVT0R1') {
  # parameters for outer fit: m = x1, beta_a = x2, scale x3, um = x4, u_delta = x5, ter = x6,
  # flags: later_diffusion = x7, s0_update = x8 (1: no update, 2: simulation, 3: analytical), 
  # u_update = x9 (1: no update, 2: step, 3: drift)
  par0 <- c(m = 0.5, beta_a = 2, u_scale = 1, u_m =0.5, u_delta = 0.1, ter = 0.01, 
            later_diffusion = 1, s0_update = 2, u_update = 1)
  switch(substr(modelname, 1, 2),
         'LS' = {# LATER SIM
           fixed = c(NA, NA, NA, NA, NA, NA, 1, 2, NA)
         },
         'DS' = {# Diffusion SIM
           fixed = c(NA, NA, NA, NA, NA, NA, 0, 2, NA)
         },
         'LA' = {# LATER analytical
           fixed = c(NA, NA, NA, NA, NA, NA, 1, 3, NA)
         },
         'DA' = {# Diffusion analytical
           fixed = c(NA, NA, NA, NA, NA, NA, 0, 3, NA)
         },
         'LN' = {# LATER noupdate
           fixed = c(0, 1, NA, NA, NA, NA, 1, 1, NA)
         },
         'DN' = {# Diffusion noupdate
           fixed = c(0, 1, NA, NA, NA, NA, 0, 1, NA)
         }
  )
  if (substr(modelname,2,2)!="N") { # fixing or varying updating parameters only makes sense if updating is turned on
    switch(substr(modelname, 4, 7),
           'M1B1' = { # No forgetting, prior fixed to beta(1,1)
             fixed[1:2] <- c(1, 1)  
           },
           'M1BV' = { # No forgetting, prior free to vary
             fixed[1:2] <- c(1, NA)   
           },
           'MVB1' = { # memory parameter free to vary, prior fixed to beta(1,1)
             fixed[1:2] <- c(NA, 1)   
           },
           'MVBV' = { # memory parameter and prior both free to vary
             fixed[1:2] <- c(NA, NA)   
           }
    )
  }
  N = nchar(modelname)
  switch(substr(modelname, N-3, N),
         'T0R1' = { # No non-decision time, no drift rate updating
           fixed[c(3:6,9)] <- c(1, 0, 0, 0, 1)
         }, 
         'T0R2' = { #  No non-decision time, drift rate updated based on rep/switch from one trial back
           fixed[c(3:6,9)] <- c(NA, 0, 0, 0, 2) # step update only (on parameter)
         },
         'T0R3' = { # No non-decision time, drift rate updating with exponential discounting of older trials
           fixed[c(3:6,9)] <- c(0, NA, NA, 0, 3) # drift rate
         },
         'TVR1' = { # decision time free to vary, no drift rate updating
           fixed[c(3:6,9)] <- c(1, 0, 0, NA, 1)
         },
         'TVR2' = { # decision time free to vary, drift rate updated based on rep/switch from one trial back
           fixed[c(3:6,9)] <- c(NA, 0, 0, NA, 2) # only one parameter 
         },
         'TVR3' = { # decision time free to vary, drift rate updating with exponential discounting of older trials
           fixed[c(3:6,9)] <- c(0, NA, NA, NA, 3) # drift rate
         }
  )
  par1 <- par0[is.na(fixed)] # find to-be-fixed parameters
#  t0 = proc.time()
  # add constraints on parameters (mem, beta, u_scale, u_m, u_delta, non_decision)
  lowers = c(0, 0.5, 0, 0, 0, 0)
  lowers = lowers[is.na(fixed[1:6])] # contraints for to-be-fixed params
  uppers = c(1, 100, 1, 1, 1, min(sub$rt-0.01)) # fixed u_m <=1, u_delta avoid Inf
  uppers = uppers[is.na(fixed[1:6])]
  par <- optim(par1,  findParameters, data = sub, fixed = fixed,
               lower = lowers, upper = uppers, method = 'L-BFGS-B')
#  proc.time() - t0
  # get all parameters
  fixed_par = fixed
  fixed_par[is.na(fixed_par)] = par$par # all fixed
  parameters <- findParameters(fixed_par, sub, fixed_par, op = FALSE)
  N <- nrow(sub)
  npar <- length(par1) + 3*length(unique(sub$targets))
  nll <- parameters[length(parameters)]
  parameters <- c(parameters, AIC=2*nll + 2*npar)
  parameters <- c(parameters, BIC=2*nll + log(N)*npar)
  parameters <- c(parameters, sub = sub$sub[1], model = modelname)
  return(parameters)
}

#' fit parameters with given models
#' 
#' This function will find optimal parameters for a given model and give subject data
#' @param sub subject rt data frame
#' @param modelname model names. 
#' model name rules: 1 ("L" or "D"): First letter indicates which model is used L = LATER model, D = DDM
#' 2-3 ("TX"): indicates whether non-decision time is fixed to 0 (T0) or a free parameter (TV)
#' 4-5 ("RX"): indicates which type of updating is done based on the response history (R0 = no updating, R1-R2: starting point updating, 
#' R3-R5: drift rate updating)
#' 6-7 ("DX"): indicates which type of updating is done based on the dimension history (D0 = no updating, D1-D2: starting point updating,
#' D3-D5: drift rate updating)
fitPar <- function(sub, modelname = 'LT0R0D0') {
  # parameters for outer fit: mem, beta = parameters of prior updating (memory and beta prior),  
  # scale,  u_mem, u_delta = parameters of rate updating,  ter = non-decision time,
  # flags: later_diffusion, resp_update (response related updating, 0: no update, 1-2: prior updating, 3-5 rate updating), 
  # dim_update (dimension related updating, 0: no update, 1-2: prior updating, 3-5 rate updating)
  par0 <- c(mem_resp = 0.5, beta_resp = 2, mem_dim = 0.5, beta_dim = 2, scale_resp = 1, u_mem_resp =0.5, u_delta_resp = 0.1,
            scale_dim = 1, u_mem_dim =0.5, u_delta_dim = 0.1, ter = 0.01, 
            later_diffusion = 1, resp_update = 2, u_update = 1)
  fixed = c(mem_resp=NA, beta_resp=NA, mem_dim=NA, beta_dim=NA, scale_resp=NA, u_mem_resp=NA,  u_delta_resp=NA, 
            scale_dim=NA, u_mem_dim=NA,  u_delta_dim=NA, ter=NA, later_diffusion=1, resp_update=0, dim_update=0)
  # Later model or DDM
  switch(substr(modelname, 1, 1),
         'L' = {# LATER 
           fixed['later_diffusion'] = 1
         },
         'D' = {# Diffusion 
           fixed['later_diffusion'] = 0
         }
  )
  # Non decision time or not
  switch(substr(modelname, 2, 3),
         'T0' = { # No non-decision time
           fixed['ter'] <- 0
         }, 
         'TV' = { # non-decision time free to vary
           fixed['ter'] <- NA
         }
  )
  # Response history based updating: 0 = no updating, 1 = Bayesian S0 updating with full memory, 2 = Bayesian S0 updating with forgetting,
  # 3 = single trial back switch cost on DR, 4 = switch cost on DR with longer memory, 5 = response based DR weighting 
  switch(substr(modelname, 4, 5),
         'R0' = { # No updating based on response history
           fixed[c('mem_resp','beta_resp','scale_resp','u_mem_resp','u_delta_resp','resp_update')] <- c(0, 1, 1, 0, 0, 0)
         }, 
         'R1' = { # Bayesian S0 updating with full memory
           fixed[c('mem_resp','beta_resp','scale_resp','u_mem_resp','u_delta_resp','resp_update')] <- c(1, NA, 1, 0, 0, 1)
         },
         'R2' = { # Bayesian S0 updating with forgetting
           fixed[c('mem_resp','beta_resp','scale_resp','u_mem_resp','u_delta_resp','resp_update')] <- c(NA, NA, 1, 0, 0, 2)
         }, 
         'R3' = { # single trial back switch cost on DR
           fixed[c('mem_resp','beta_resp','scale_resp','u_mem_resp','u_delta_resp','resp_update')] <- c(0, 1, NA, 0, 0, 3)
         },
         'R4' = { # switch cost on DR with longer memory
           fixed[c('mem_resp','beta_resp','scale_resp','u_mem_resp','u_delta_resp','resp_update')] <- c(0, 1, 1, NA, NA, 4)
         },
         'R5' = { # response based DR weighting 
           fixed[c('mem_resp','beta_resp','scale_resp','u_mem_resp','u_delta_resp','resp_update')] <- c(0, 1, 1, NA, NA, 5)
         }
  )
  # Dimension history based updating: 0 = no updating, 1 = Bayesian S0 updating with full memory, 2 = Bayesian S0 updating with forgetting,
  # 3 = single trial back switch cost on DR, 4 = switch cost on DR with longer memory, 5 = dimension based DR weighting 
  switch(substr(modelname, 6, 7),
         'D0' = { #  no updating
           fixed[c('mem_dim','beta_dim','scale_dim','u_mem_dim','u_delta_dim','dim_update')] <- c(0, 1, 1, 0, 0, 0)
         }, 
         'D1' = { # Bayesian S0 updating with full memory
           fixed[c('mem_dim','beta_dim','scale_dim','u_mem_dim','u_delta_dim','dim_update')] <- c(1, NA, 1, 0, 0, 1)
         }, 
         'D2' = { # Bayesian S0 updating with forgetting
           fixed[c('mem_dim','beta_dim','scale_dim','u_mem_dim','u_delta_dim','dim_update')] <- c(NA, NA, 1, 0, 0, 2)
         }, 
         'D3' = { # single trial back switch cost on DR
           fixed[c('mem_dim','beta_dim','scale_dim','u_mem_dim','u_delta_dim','dim_update')] <- c(0, 1, NA, 0, 0, 3)
         }, 
         'D4' = { # switch cost on DR with longer memory
           fixed[c('mem_dim','beta_dim','scale_dim','u_mem_dim','u_delta_dim','dim_update')] <- c(0, 1, 1, NA, NA, 4)
         }, 
         'D5' = { # dimension based DR weighting 
           fixed[c('mem_dim','beta_dim','scale_dim','u_mem_dim','u_delta_dim','dim_update')] <- c(0, 1, 1, NA, NA, 5)
         }
  )
  
  par1 <- par0[is.na(fixed)] # find to-be-fixed parameters
  #  t0 = proc.time()
  
  # add constraints on parameters (mem, beta, u_scale, u_mem, u_delta, non_decision)
  lowers = c(mem_resp=0, beta_resp=0.5, mem_dim=0, beta_dim=0.5, scale_resp=0, u_mem_resp=0,  u_delta_resp=0,
             scale_dim=0, u_mem_dim=0,  u_delta_dim=0, ter=0)
  lowers = lowers[is.na(fixed[1:11])] # contraints for to-be-fixed params
  uppers = c(mem_resp=1, beta_resp=100, mem_dim=1, beta_dim=100, scale_resp=1, u_mem_resp=1, u_delta_resp=1,
             scale_dim=1, u_mem_dim=1, u_delta_dim=1, ter=min(filter(sub, outlier==FALSE)$rt-0.01)) # fixed u_m <=1, u_delta avoid Inf
  uppers = uppers[is.na(fixed[1:11])]

  par <- optim(par1,  findParameters, data = sub, fixed = fixed,
               lower = lowers, upper = uppers, method = 'L-BFGS-B')
  #  proc.time() - t0
  # get all parameters
  fixed_par = fixed
  fixed_par[is.na(fixed_par)] = par$par # all fixed
  parameters <- findParameters(fixed_par, sub, fixed_par, op = FALSE)
  N <- nrow(sub)
  npar <- length(par1) + 3*length(unique(sub$targets))
  nll <- parameters[length(parameters)]
  parameters <- c(parameters, AIC=2*nll + 2*npar)
  parameters <- c(parameters, BIC=2*nll + log(N)*npar)
  parameters <- c(parameters, sub = sub$sub[1], model = modelname)
  return(parameters)
}

#' Parallel compute parameters
#' 
#' Parallel compute parameters for each subject each model
#' @param lsubs list of subject data, which should include 'targets' and 'inttrial'
#' @param modelnames a list of modelnames (see optimization function)
#' @return paras return a list of parameters
parEstimate <- function(lsubs, modelnames){
  no_cores <- detectCores() - 1
  print(no_cores)
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {library(dplyr) })
  clusterExport(cl, c('findParameters','findInnerParameters', 'updateRate','updateRateMem', 'updateRateWeight',
                      'fitPar','drecinorm','updatePriorA','updatePrior','logll','log_ddiffusion'))
  
  t0 = proc.time()
  paras <- clusterMap(cl, fitPar, lsubs, modelnames)
  stopCluster(cl)
  print(proc.time()-t0)
  return(paras)
}
