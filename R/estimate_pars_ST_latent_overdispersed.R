if(!dir.exists(here::here("output"))) {
  dir.create(here::here("output"))
}
library(raster)
library(geosphere)
library(tidyverse)
library(sp)
library(rgdal)
library(rdist)
library(rgeos)
library(ggplot2)
library(invgamma)
library(mvnfast)
library(splines)
library(pgdraw)
library(fields)
library(geoR)
library(dplyr)
library(data.table)
# install.packages("/home/adawson/Documents/projects/pgR", repos=NULL, type="source")
library(pgR)

version='4.0'

#### DATA PREP ####
y <- readRDS(here::here('data', paste0('pollen_dat_', version, '.RDS')))
taxa.keep <- readRDS(here::here('data', paste0('taxa_', version, '.RDS')))
locs <- readRDS(here::here('data', paste0('pollen_locs_', version, '.RDS')))
rescale <- 1e3

#### RUNNING THE MODEL & SAVING OUTPUT####
# scale locations so distance values aren't too large (from 1m to 100km units)
locs_scaled <- locs/rescale
N_locs_dat = nrow(locs_scaled)
X <- matrix(rep(1, N_locs_dat), N_locs_dat, 1)

params <- default_params()
params$n_adapt <- 2000
params$n_mcmc <- 5000
params$n_message <- 100
params$n_thin <- 5
priors <- default_priors_pg_stlm(y, X, corr_fun = "matern")
# > priors
# $mu_beta
# [1] 0
# $Sigma_beta
# [,1]
# [1,]    1
# $alpha_tau
# [1] 0.1
# $beta_tau
# [1] 0.1
# $mean_nu
# [1] -1
# $sd_nu
# [1] 1
# $mean_range
# [1] 0
# $sd_range
# [1] 10

# mu_sigma_prop    = 1
priors$mean_nu     = -0.7  # change from -0.9 to -0.7 on 6/29/21
priors$sd_nu       = sqrt(0.005)
priors$mean_range  = 7  # change from 4.6 to 7 on 6/29/21
priors$sd_range    = sqrt(0.2)
priors$alpha_tau   = 2  # changed from default to 2 (checked by AS/AD 10/29/21)
priors$beta_tau    = 1  # changed from default to 1 (checked by AS/AD 10/29/21)
priors$mu_beta     = 0  #current default in pgR (changed by AS/AD 10/29/21)
priors$Sigma_beta  = 10 #current default in pgR (changed by AS/AD 10/29/21)

J <- ncol(y)
theta_mean <- c(priors$mean_range, priors$mean_nu)
theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)

inits <- list(
  beta  = t(mvnfast::rmvn(J-1, priors$mu_beta, priors$Sigma_beta)),
  tau2  = rgamma(J-1, priors$alpha_tau, priors$beta_tau),
  theta = mvnfast::rmvn(J-1, theta_mean, theta_var),
  rho = 0.8  # changed from runif(1, 0, 1) to 0.8 on 6/29
)

# XXX: need this to work around undefined variable in pgSPLM
d <- ncol(y)
Y = y
X = as.matrix(X)
locs = as.matrix(locs_scaled)
n_cores = 12L
n_chain = 1

for (n in 1:dim(Y)[1]){
  for (t in 1:dim(Y)[3]){
    if (!any(is.na(Y[n,,t]))){
      if (sum(Y[n,,t])==0){
        #print(Y[n,,t])
        Y[n,,t] = rep(NA, dim(Y)[2])
      }
    }
  }
}

# code to run matern model
if (!file.exists( here::here('output', paste0('polya-gamma-posts_', version, '_latent_overdispersed.RDS')))) {
  
  out <- pg_stlm_latent_overdispersed(Y = Y,
                 X = X,
                 locs = locs,
                 params,
                 priors,
                 n_cores = n_cores,
                 n_chain = n_chain,
                 corr_fun = "matern",
                 shared_covariance_params = FALSE,
                 inits = inits)
  saveRDS(out, here::here('output', paste0('polya-gamma-posts_', version, '_latent_overdispersed.RDS')),
          compress = FALSE)
  
#  pushoverr::pushover(message = "Finished fitting latent overdispersed Matern model")
}

dat <- list(y = y,
            X = X, 
            locs = locs_scaled,
            rescale = rescale,
            taxa.keep = taxa.keep)
if (!file.exists( here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))) {
  saveRDS(dat, here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))
}


