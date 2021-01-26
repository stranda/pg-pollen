
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

version='Dec15'

#### DATA PREP ####
y <- readRDS(paste0('data/', 'pollen_dat_', version, '.RDS'))
taxa.keep <- readRDS(paste0('data/', 'pollen_taxa_', version, '.RDS'))
locs <- readRDS(paste0('data/', 'pollen_locs_', version, '.RDS'))
rescale <- 1e3

#### RUNNING THE MODEL & SAVING OUTPUT####
# scale locations so distance values aren't too large (from 1m to 100km units)
locs_scaled <- locs/rescale
N_locs_dat = nrow(locs_scaled)
X <- matrix(rep(1, N_locs_dat),N_locs_dat, 1)

params <- default_params()
params$n_adapt <- 1000
params$n_mcmc <- 1800
params$n_message <- 100
params$n_thin <- 1
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
priors$mean_nu     = -0.9
priors$sd_nu       = sqrt(0.005)
priors$mean_range  = 4.6
priors$sd_range    = sqrt(0.2)
# priors$alpha_tau  = 0.5
# priors$beta_tau   = 10

J <- ncol(y)
theta_mean <- c(priors$mean_range, priors$mean_nu)
theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)

inits <- list(
  beta  = t(mvnfast::rmvn(J-1, priors$mu_beta, priors$Sigma_beta)),
  tau2  = rgamma(J-1, priors$alpha_tau, priors$beta_tau),
  theta = mvnfast::rmvn(J-1, theta_mean, theta_var),
  rho = runif(1, 0, 1)
)

# XXX: need this to work around undefined variable in pgSPLM
d <- ncol(y)
Y = y
X = as.matrix(X)
locs = as.matrix(locs_scaled)
n_cores = 20L
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
out <- pgSTLM(Y = Y,
              X = X,
              locs = locs,
              params,
              priors,
              n_cores = n_cores,
              n_chain = n_chain,
              corr_fun = "matern",
              shared_covariance_params = FALSE,
              inits = inits,
              verbose=TRUE)
saveRDS(out, paste0('output/', 'polya-gamma-posts_', version, '.RDS'))

dat <- list(y = y,
            X = X, 
            locs = locs_scaled,
            rescale = rescale,
            taxa.keep = taxa.keep)
saveRDS(dat, paste0('output/', 'polya-gamma-dat_', version,'.RDS'))
