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

run='pgR-ST-short'

#### DATA PREP ####
dat  = readRDS("data/pollen_dat_vtime_1.0.RDS")
taxa.keep = taxa.keep = c('Acer', 'Alnus','Betula', 'Cyperaceae', 'Fagus', 'Ostrya.Carpinus',
                          'Picea', 'Pinus', 'Quercus', 'Tsuga', 'Ulmus', 'Other')#c('Acer','Betula','Ostrya.Carpinus','Ulmus')

locs = readRDS('data/pollen_locs_vtime_1.0.RDS')[,1:2]

locs_grid = readRDS('data/grid.RDS')

##### WRITING THE MODEL #####
K <- 200
# locs <- as.data.frame(dat[,c('x', 'y')])
rescale <- 1e3

# find cells within prediction grid that are closest to pollen core sites, then 
# replace core site location with location of closest grid cell
D_inter <- cdist(as.matrix(locs_grid/rescale), as.matrix(locs/rescale)) # N_locs x N_cores
any(D_inter == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
# D_inter <- ifelse(D_inter == 0, 0.007, D_inter)  # if so, remove them
closest.cell <- apply(D_inter, 2, function(x) which.min(x))

dat.merged = list()
# for (i in 1:length(dat)){
# dat.mod <- data.frame(closest=closest.cell, dat)
# dat.merged <- as.data.table(dat.mod[, !(colnames(dat.mod) %in% c('x', 'y'))])[, lapply(.SD, sum, na.rm=TRUE), 
#                                                                              by <- list(closest)]
# }
# 
# locs <- locs_grid[dat.merged$by,]

dat.merged = dat
# taxa.keep <- as.vector(colnames(dat.merged)[!(colnames(dat.merged) %in% c('by'))])
# y <- as.data.frame(dat.merged[,..taxa.keep])
y <- dat.merged

#### RUNNING THE MODEL & SAVING OUTPUT####
# scale locations so distance values aren't too large (from 1m to 100km units)
locs_scaled <- locs/rescale
N_locs_dat = nrow(locs_scaled)
X <- matrix(rep(1, N_locs_dat),N_locs_dat, 1)

params <- default_params()
params$n_adapt <- 500
params$n_mcmc <- 1000
params$n_message <- 50
params$n_thin <- 1
priors <- default_priors_pgSTLM(y, X, corr_fun = "matern")
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

# inits <- default_inits_pgSPLM(y, 
#                               X, 
#                               priors, 
#                               corr_fun = "matern", 
#                               shared_covariance_params = FALSE)

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
# Y = as.matrix(y)
Y = y
X = as.matrix(X)
locs = as.matrix(locs_scaled)
n_cores = 1L

source("/home/adawson/Documents/projects/pgR/R/helper-functions.R")
source("/home/adawson/Documents/projects/pgR/R/check-input-spatial.R")
source("/home/adawson/Documents/projects/pgR/R/check-params.R")
source("/home/adawson/Documents/projects/pgR/R/pgSTLM.R")
source("/home/adawson/Documents/projects/pgR/R/update-tuning.R")
out <- pgSTLM(Y = Y,
              X = X,
              locs = locs,
              params,
              priors,
              n_cores = n_cores,
              corr_fun = "matern",
              shared_covariance_params = FALSE,
              inits = inits,
              verbose=TRUE)

dir.create(file.path('output'), showWarnings = FALSE)
saveRDS(out, paste0('output/polya-gamma-posts_pgR_', run,'.RDS'))

dat <- list(y = y,
            X = X, 
           locs = locs_scaled,
           rescale = rescale,
           taxa.keep = taxa.keep)
saveRDS(dat, paste0('output/polya-gamma-dat_pgR_', run,'.RDS'))
