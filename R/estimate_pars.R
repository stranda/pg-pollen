
require(raster)
require(geosphere)
require(tidyverse)
require(sp)
require(rgdal)
require(rdist)
require(rgeos)
require(ggplot2)
require(invgamma)

require(mvnfast)
require(splines)
require(pgdraw)
require(fields)
require(geoR)
library(dplyr)

library(data.table)

library(pgR)

#
# after loading pgR, make sure to source "default-inits.R" and "pgSPLM.R"
#

run='pgSPLM_test'

#### DATA PREP ####
# dato <- readRDS("../data/pollen_data.RData")
# dat <- readRDS("../data/6taxa_dat.RData")
dat <- readRDS("../data/12taxa_dat.RDS")

locs_grid = readRDS('../polya-gamma/data/grid.RDS')

##### WRITING THE MODEL #####
K <- 200
locs = as.data.frame(dat[,c('x', 'y')])
rescale=1e3

D_inter   <- cdist(as.matrix(locs_grid/rescale), as.matrix(locs/rescale)) # N_locs x N_cores
any(D_inter == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
# D_inter <- ifelse(D_inter == 0, 0.007, D_inter)  # remove them
closest.cell = apply(D_inter, 2, function(x) which.min(x))

dat.mod = data.frame(closest=closest.cell, dat)

dat.merged = as.data.table(dat.mod[, !(colnames(dat.mod) %in% c('x', 'y'))])[, lapply(.SD, sum, na.rm=TRUE), 
                                                                             by = list(closest)]
locs = locs_grid[dat.merged$closest,]

# taxa.keep =  as.vector(colnames(dat)[!(colnames(dat) %in% c('x', 'y', 'closest'))])
taxa.keep =  as.vector(colnames(dat.merged)[!(colnames(dat.merged) %in% c('closest'))])
y = as.data.frame(dat.merged[,..taxa.keep])

# ## calculate the Matern correlation using parameters theta on the log scale
# correlation_function <- function(D, theta) {
#   geoR::matern(D, exp(theta[1]), exp(theta[2]))
# }

#### RUNNING THE MODEL & SAVING OUTPUT####
locs_scaled = locs/rescale

X <- matrix(rep(1, nrow(y)), nrow(y), 1)

params <- default_params()
params$n_adapt <- 100
params$n_mcmc <- 100
params$n_message <- 50
params$n_thin <- 1
priors <- default_priors_pgSPLM(y, X)
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

# mu_sigma_prop = 1
# mean_nu       = -0.9
# sd_nu         = 0.005
# mean_range    = 4.6
# sd_range      = 0.2
# alpha_tau     = 1/2
# beta_tau      = 10

priors$mean_nu   = -0.9
priors$sd_nu      = 0.005
priors$mean_range = 4.6
priors$sd_range  = 0.2
#priors$alpha_tau  = 0.5
#priors$beta_tau   = 10

inits <- default_inits_pgSPLM(y, 
                              X, 
                              priors, 
                              corr_fun = "matern", 
                              shared_theta=FALSE, 
                              shared_tau=TRUE)

#inits <- NULL

# XXX: need this to work around undefined variable in pgSPLM
d <- ncol(y)

Y = as.matrix(y)
X = as.matrix(X)
locs = as.matrix(locs_scaled)
n_cores = 8L
shared_covariance_params = FALSE
shared_variance_params   = TRUE

# out <- pgSPLMAD(Y = Y, 
#               X = X, 
#               locs = locs, 
#               params, 
#               priors, 
#               n_cores = n_cores)

# out <- pgSPLM(Y = Y, 
#               X = X, 
#               locs = locs, 
#               params, 
#               priors, 
#               n_cores = n_cores,
#               shared_covariance_params = shared_covariance_params)
# 
# out <- pgSPLMAD2(Y = Y, 
#               X = X, 
#               locs = locs, 
#               params, 
#               priors, 
#               n_cores = 8L,
#               shared_covariance_params = shared_covariance_params,
#               shared_variance_params = shared_variance_params,
#               verbose=TRUE)

out <- pgSPLM(Y = Y, 
              X = X, 
              locs = locs, 
              params, 
              priors, 
              n_cores = 6L,
              corr_fun = "matern",
              shared_theta = FALSE,
              shared_tau = TRUE,
              verbose=TRUE)

# out <- pgSPLMAD3(Y = Y, 
#                  X = X, 
#                  locs = locs, 
#                  params, 
#                  priors, 
#                  n_cores = 8L,
#                  shared_covariance_params = shared_covariance_params,
#                  shared_variance_params = shared_variance_params,
#                  verbose=TRUE)

# out <- mcmc_fix(y,
#                locs_scaled,
#                K = 100,
#                message = 100,
#                mean_nu = mean_nu,
#                sd_nu = sd_nu,
#                mean_range = mean_range,
#                sd_range = sd_range)
saveRDS(out, paste0('polya-gamma-posts_pgR_', run,'.RDS'))

dat = list(y=y,
           locs=locs_scaled,
           rescale=rescale)
saveRDS(dat, paste0('polya-gamma-dat_pgR_', run,'.RDS'))