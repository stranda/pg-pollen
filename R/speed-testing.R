# Testing model fit speed as a function of number of cores ---------------------
if(!dir.exists(here::here("results"))) {
    dir.create(here::here("results"))
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

version='3.1'  # v. 3.1 uses same time bins as ABC, ENM, and uses updated priors

#### DATA PREP ####
y <- readRDS(here::here('data', paste0('paleo_pollen_dat_', version, '.RDS')))
taxa.keep <- readRDS(here::here('data', paste0('pollen_taxa_', version, '.RDS')))
locs <- readRDS(here::here('data', paste0('paleo_pollen_locs_', version, '.RDS')))
rescale <- 1e3

#### RUNNING THE MODEL & SAVING OUTPUT####
# scale locations so distance values aren't too large (from 1m to 100km units)
locs_scaled <- locs/rescale
N_locs_dat = nrow(locs_scaled)
X <- matrix(rep(1, N_locs_dat),N_locs_dat, 1)

params <- default_params()
params$n_adapt <- 5
params$n_mcmc <- 5
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
priors$mean_nu     = -0.7  # change from -0.9 to -0.7 on 6/29/21
priors$sd_nu       = sqrt(0.005)
priors$mean_range  = 7  # change from 4.6 to 7 on 6/29/21
priors$sd_range    = sqrt(0.2)
priors$alpha_tau   = 2  # changed from default to 2
priors$beta_tau    = 1  # changed from default to 1
priors$mu_beta     = -5  # changed from default to -5 on 6/29/21
priors$Sigma_beta  = 0.5  # changed from default to 0.5 on 6/29/21

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

# convert Y into an integer matrix
# mode(Y) <- "integer"

# Matern -----------------------------------------------------------------------
if (file.exists(here::here("results", "matern-timings.Rdata"))) {
    load(here::here("results", "matern-timings.Rdata"))
} else {
    
  # dat_matern <- data.frame(n_cores = c(1L),
  #                          time = rep(0, 1))
  dat_matern <- data.frame(n_cores = c(1L, 10L, 20L, 30L, 40L, 50L, 60L),
                             time = rep(0, 7))
    for (i in dat_matern$n_cores) {
        message("Fitting matern model with ", i, " cores")
        n_cores = i
        
        # code to run matern model
        start <- Sys.time()
        pg_stlm(Y = Y,
                X = X,
                locs = locs,
                params,
                priors,
                n_cores = n_cores,
                n_chain = n_chain,
                corr_fun = "matern",
                shared_covariance_params = FALSE,
                inits = inits)
        finish <- Sys.time()

        dat_matern$time[dat_matern$n_cores == i] <- finish - start
        
    }
    save(dat_matern, file = here::here("results", "matern-timings.Rdata"))
    pushoverr::pushover(message = "Finished testing Matern model timings")
}



# BayesMRA ---------------------------------------------------------------------
# use the default model priors just to get the code working
priors <- default_priors_pg_stlm(Y, X)

if (file.exists(here::here("results", "mra-timings.Rdata"))) {
    load(here::here("results", "mra-timings.Rdata"))
} else {
    
  # dat_mra <- data.frame(n_cores = c(1L),
  #                       time = rep(0, 1))
  dat_mra <- expand_grid(
    M = 2:4, 
    n_coarse_grid = c(5, 10, 15, 20, 25)) %>%
    mutate(time = 0)
    for (i in 1:nrow(dat_mra)) {
        
        message("Fitting mra model with ", dat_mra$M[i], " resolutions and ", dat_mra$n_coarse_grid[i], " coarse grid points")
      
        # code to run mra model
        start <- Sys.time()
        pg_stlm_mra(Y = Y,
                    X = X,
                    locs = locs,
                    params,
                    priors,
                    n_cores = 1,
                    n_chain = n_chain,
                    M = dat_mra$M[i], 
                    n_coarse_grid = dat_mra$n_coarse_grid[i])
        finish <- Sys.time()
        
        dat_mra$time[i] <- difftime(finish, start, units = "secs")
        
    }
    
    save(dat_mra, file = here::here("results", "mra-timings.Rdata"))
    pushoverr::pushover(message = "Finished testing MRA model timings")
}

ggsave(
  dat_mra %>%
    mutate(M = as.factor(M)) %>%
    # rbind(mutate(dat_matern, model = "Matern"),
    #     mutate(dat_mra, model = "MRA"),
    #     mutate(dat_matern_new, model = "Matern approx"),
    #     mutate(dat_mra_new, model = "MRA approx")) %>%
    ggplot(aes(x = n_coarse_grid, y = time, group = M, color = M)) +
    scale_color_viridis_d(end = 0.7) +
    geom_line() +
    ylab("time for 10 iterations") +
    xlab("number of coarse grid functions"),
  file = here::here("figures", "MRA-timings.png"),
  width = 16,
  height = 9)

dat_mra$n_basis <- 0
for (i in 1:nrow(dat_mra)) {
  dat_mra$n_basis[i] <- sum(BayesMRA::mra_wendland_2d(locs, M = dat_mra$M[i], n_coarse_grid = dat_mra$n_coarse_grid[i])$n_dims)
}

ggsave(
  dat_mra %>%
    mutate(M = as.factor(M)) %>%
    # rbind(mutate(dat_matern, model = "Matern"),
    #     mutate(dat_mra, model = "MRA"),
    #     mutate(dat_matern_new, model = "Matern approx"),
    #     mutate(dat_mra_new, model = "MRA approx")) %>%
    ggplot(aes(x = n_basis, y = time, group = M, color = M)) +
    geom_line() +
    scale_color_viridis_d(end = 0.7) +
    ylab("time for 10 iterations") +
    xlab("number of basis functions"),
  file = here::here("figures", "MRA-timings-basis.png"),
  width = 16,
  height = 9)

ggsave(
  dat_mra %>%
    mutate(M = as.factor(M)) %>%
    # rbind(mutate(dat_matern, model = "Matern"),
    #     mutate(dat_mra, model = "MRA"),
    #     mutate(dat_matern_new, model = "Matern approx"),
    #     mutate(dat_mra_new, model = "MRA approx")) %>%
    ggplot(aes(x = n_coarse_grid, y = n_basis, group = M, color = M)) +
    geom_line() +
    scale_color_viridis_d(end = 0.7) +
    ylab("total number of basis functions") +
    xlab("number of coarse grid functions"),
  file = here::here("figures", "MRA-timings-basis-resolution.png"),
  width = 16,
  height = 9)

# ggsave(
#   rbind(mutate(dat_matern, model = "Matern"),
#         mutate(dat_mra, model = "MRA")) %>%
#     # rbind(mutate(dat_matern, model = "Matern"),
#     #     mutate(dat_mra, model = "MRA"),
#     #     mutate(dat_matern_new, model = "Matern approx"),
#     #     mutate(dat_mra_new, model = "MRA approx")) %>%
#     # ggplot(aes(x = n_cores, y = time, color = model)) +
#     # geom_line() +
#     ggplot(aes(x = model, y = time, color = model)) +
#     geom_col() +
#     ylab("time for 10 iterations") +
#     # xlab("number of cores"),
#     xlab("model"),
#   file = here::here("figures", "model-timings.png"),
#   width = 16,
#   height = 9)

