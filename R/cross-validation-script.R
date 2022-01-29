# TODO ------------------------------------------------------------------------
# 1. Write code to validate MCMC convergence -- plot and save trace plots
# 2. Write code to evaluate cross-validation predictive log-likelihood
# 3. Write code to calculate loo log-likelihood for each model

# generate the cross-validation datasets
if(!dir.exists(here::here("output"))) {
    dir.create(here::here("output"))
}
if(!dir.exists(here::here("output", "cross-validate"))) {
    dir.create(here::here("output", "cross-validate"))
}


if(!dir.exists(here::here("figures", "cross-validation"))) {
    dir.create(here::here("figures", "cross-validation"))
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
library(xtable)
# install.packages("/home/adawson/Documents/projects/pgR", repos=NULL, type="source")
library(pgR)

source(here::here("R", "xtable_print_bold.R"))

version='4.1'  # v. 4.1 uses same time bins as ABC, ENM, and uses updated priors

#### DATA PREP ----------------------------------------------------------------
y <- readRDS(here::here('data', paste0('pollen_dat_', version, '.RDS')))
taxa.keep <- readRDS(here::here('data', paste0('taxa_', version, '.RDS')))
locs <- readRDS(here::here('data', paste0('pollen_locs_', version, '.RDS')))
rescale <- 1e3

locs_scaled <- locs/rescale
N_locs_dat = nrow(locs_scaled)
X <- matrix(rep(1, N_locs_dat),N_locs_dat, 1)

# setup cross-validation ------------------------------------------------------

set.seed(2021)
K <- 5
total_obs <- dim(y)[1] * dim(y)[3]
# make sure this is an exact integer (which it is in our case -- could be more general but is not needed)
n_cv <- total_obs / K

if (file.exists(here::here("output", "cross-validate", "cv_idx.RData"))) {
    load(here::here("output", "cross-validate", "cv_idx.RData"))
} else {
    cv_idx <- matrix(sample(1:total_obs), n_cv, K)
    obs_idx <- expand_grid(row_id = 1:dim(y)[1], col_id = 1:dim(y)[3])
    save(cv_idx, obs_idx, file = here::here("output", "cross-validate", "cv_idx.RData"))
}

# fit the Matern model for each of the K folds --------------------------------

for (k in 1:K) {
    # set up the cross-validation
    y_train <- y
    y_test <- array(NA, dim = dim(y))
    obs_fold <- obs_idx[cv_idx[, k], ]
    for (i in 1:n_cv) {
        y_train[obs_fold$row_id[i], , obs_fold$col_id[i]] <- NA
        y_test[obs_fold$row_id[i], , obs_fold$col_id[i]] <- y[obs_fold$row_id[i], , obs_fold$col_id[i]] 
    }
    # fit the Matern model on the fold
    
    params <- default_params()
    params$n_adapt <- 2000
    params$n_mcmc <- 5000
    params$n_message <- 100
    params$n_thin <- 5
    priors <- default_priors_pg_stlm(y_train, X, corr_fun = "matern")
    
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
    
    # Y = y_train
    X = as.matrix(X)
    locs = as.matrix(locs_scaled)
    n_cores = 12L
    n_chain = 1
    
    if (!file.exists( here::here('output', 'cross-validate' , paste0('matern-cv-', k, '.RData')))) {
        
        out <- pg_stlm(Y = y_train,
                       X = X,
                       locs = locs,
                       params,
                       priors,
                       n_cores = n_cores,
                       n_chain = n_chain,
                       corr_fun = "matern",
                       shared_covariance_params = FALSE,
                       inits = inits)
        save(out, y_train, y_test, obs_fold, 
             file = here::here('output', 'cross-validate' , paste0('matern-cv-', k, '.RData')),
             compress = FALSE)
        
###        pushoverr::pushover(message = paste0("Finished fitting Matern model for fold ", k, " out of ", K))
        print(paste0("Finished fitting Matern model for fold ", k, " out of ", K))
    }
    
    if (!file.exists( here::here('output', 'cross-validate' , paste0('overdispersed-cv-', k, '.RData')))) {
        
        out <- pg_stlm_overdispersed(Y = y_train,
                                     X = X,
                                     locs = locs,
                                     params,
                                     priors,
                                     n_cores = n_cores,
                                     n_chain = n_chain,
                                     corr_fun = "matern",
                                     shared_covariance_params = FALSE,
                                     inits = inits)
        save(out, y_train, y_test, obs_fold, 
             file = here::here('output', 'cross-validate' , paste0('overdispersed-cv-', k, '.RData')),
             compress = FALSE)
        
###        pushoverr::pushover(message = paste0("Finished fitting Overdispersed model for fold ", k, " out of ", K))
print(paste0("Finished fitting Overdispersed model for fold ", k, " out of ", K))
    }
    
    if (!file.exists( here::here('output', 'cross-validate' , paste0('latent-cv-', k, '.RData')))) {
        
        out <- pg_stlm_latent_overdispersed(Y = y_train,
                                            X = X,
                                            locs = locs,
                                            params,
                                            priors,
                                            n_cores = n_cores,
                                            n_chain = n_chain,
                                            corr_fun = "matern",
                                            shared_covariance_params = FALSE,
                                            inits = inits)
        save(out, y_train, y_test, obs_fold, 
             file = here::here('output', 'cross-validate' , paste0('latent-cv-', k, '.RData')),
             compress = FALSE)
        
###        pushoverr::pushover(message = paste0("Finished fitting Latent model for fold ", k, " out of ", K))
        print(paste0("Finished fitting Latent model for fold ", k, " out of ", K))
    }
    
    if (!file.exists( here::here('output', 'cross-validate' , paste0('mra-cv-', k, '.RData')))) {
        
        
        inits <- list(
            beta  = t(mvnfast::rmvn(J-1, priors$mu_beta, priors$Sigma_beta)),
            # tau2  = rgamma(J-1, priors$alpha_tau, priors$beta_tau),
            theta = mvnfast::rmvn(J-1, theta_mean, theta_var),
            rho = rep(0.8, J-1)  # changed from runif(1, 0, 1) to 0.8 on 6/29
        )
        
        out <- pg_stlm_mra(Y = y_train,
                       X = X,
                       locs = locs,
                       params,
                       priors,
                       n_cores = n_cores,
                       M             = 3,
                       n_coarse_grid = 10,
                       n_chain = n_chain,
                       inits = inits)
        save(out, y_train, y_test, obs_fold, 
             file = here::here('output', 'cross-validate' , paste0('mra-cv-', k, '.RData')),
             compress = FALSE)
        
###        pushoverr::pushover(message = paste0("Finished fitting MRA model for fold ", k, " out of ", K))
        print(paste0("Finished fitting MRA model for fold ", k, " out of ", K))
    }
}



# Check that the models have converged ----

for (k in 1:K) {
    # set up the cross-validation
    y_train <- y
    y_test <- array(NA, dim = dim(y))
    obs_fold <- obs_idx[cv_idx[, k], ]
    row_id <- c()
    col_id <- c()
    for (i in 1:n_cv) {
        y_train[obs_fold$row_id[i], , obs_fold$col_id[i]] <- NA
        y_test[obs_fold$row_id[i], , obs_fold$col_id[i]] <- y[obs_fold$row_id[i], , obs_fold$col_id[i]] 
        if (!any(is.na(y_test[obs_fold$row_id[i], , obs_fold$col_id[i]]))) {
            row_id <- c(row_id, obs_fold$row_id[i])
            col_id <- c(col_id, obs_fold$col_id[i])
        }
    }
    
    
    # load the Matern model
    X = as.matrix(X)
    locs = as.matrix(locs_scaled)
    
    source("R/plot-trace.R")
    source("R/plot-trace_latent.R")
    
    # generate trace plots for matern model
    load( here::here('output', 'cross-validate' , paste0('matern-cv-', k, '.RData')))
    plot_trace(out, base_size = 7, file = here::here("figures", "cross-validation", paste0("matern-trace-plots-", k, ".png")))
    plot_trace_latent
    
    # generate trace plots for overdispesed model
    load( here::here('output', 'cross-validate' , paste0('overdispersed-cv-', k, '.RData')))
    plot_trace(out, base_size = 7, file = here::here("figures", "cross-validation", paste0("overdispersed-trace-plots-", k, ".png")))    
    
    # generate trace plots for latent model
    load( here::here('output', 'cross-validate' , paste0('latent-cv-', k, '.RData')))
    plot_trace(out, base_size = 7, file = here::here("figures", "cross-validation", paste0("latent-trace-plots-", k, ".png")))    
    
    
    
    
    
}



# calculate the cross-validation predictive log-likelihood --------------------

for (k in 1:K) {
    # set up the cross-validation
    y_train <- y
    y_test <- array(NA, dim = dim(y))
    obs_fold <- obs_idx[cv_idx[, k], ]
    row_id <- c()
    col_id <- c()
    for (i in 1:n_cv) {
        y_train[obs_fold$row_id[i], , obs_fold$col_id[i]] <- NA
        y_test[obs_fold$row_id[i], , obs_fold$col_id[i]] <- y[obs_fold$row_id[i], , obs_fold$col_id[i]] 
        if (!any(is.na(y_test[obs_fold$row_id[i], , obs_fold$col_id[i]]))) {
            row_id <- c(row_id, obs_fold$row_id[i])
            col_id <- c(col_id, obs_fold$col_id[i])
        }
    }
    
    
    # load the Matern model
    X = as.matrix(X)
    locs = as.matrix(locs_scaled)
    
    if (!file.exists( here::here('output', 'cross-validate' , paste0('predictive-ll-matern-cv-', k, '.RData')))) {
        load( here::here('output', 'cross-validate' , paste0('matern-cv-', k, '.RData')))
        
        ll_matern_raw <- calc_ll_pg_stlm(y_test, X, out)
        # drop the NAs and add in row/column id
        ll_matern <- matrix(NA, dim(ll_matern_raw$ll)[1], length(row_id))
        for (i in 1:length(row_id)) {
            ll_matern[, i] <- ll_matern_raw$ll[, row_id[i], col_id[i]]
        }
        
        save(ll_matern_raw, ll_matern,
             row_id, col_id,
             file = here::here('output', 'cross-validate' , paste0('predictive-ll-matern-cv-', k, '.RData')))
    }
    
    if (!file.exists( here::here('output', 'cross-validate' , paste0('predictive-ll-overdispersed-cv-', k, '.RData')))) {
        load( here::here('output', 'cross-validate' , paste0('overdispersed-cv-', k, '.RData')))
        
        ll_overdispersed_raw <- calc_ll_pg_stlm(y_test, X, out)
        # drop the NAs and add in row/column id
        ll_overdispersed <- matrix(NA, dim(ll_overdispersed_raw$ll)[1], length(row_id))
        for (i in 1:length(row_id)) {
            ll_overdispersed[, i] <- ll_overdispersed_raw$ll[, row_id[i], col_id[i]]
        }
        
        save(ll_overdispersed_raw, ll_overdispersed, 
             row_id, col_id,
             file = here::here('output', 'cross-validate' , paste0('predictive-ll-overdispersed-cv-', k, '.RData')))
    }
    
    if (!file.exists( here::here('output', 'cross-validate' , paste0('predictive-ll-latent-cv-', k, '.RData')))) {
        load( here::here('output', 'cross-validate' , paste0('latent-cv-', k, '.RData')))
        
        ll_latent_raw <- calc_ll_pg_stlm(y_test, X, out)
        # drop the NAs and add in row/column id
        ll_latent <- matrix(NA, dim(ll_latent_raw$ll)[1], length(row_id))
        for (i in 1:length(row_id)) {
            ll_latent[, i] <- ll_latent_raw$ll[, row_id[i], col_id[i]]
        }

        save(ll_latent_raw, ll_latent, 
             row_id, col_id,
             file = here::here('output', 'cross-validate' , paste0('predictive-ll-latent-cv-', k, '.RData')))
    }
}

# process the cross-validation predictive likelihoods -------------------------
if (file.exists( here::here('output', 'cross-validate', 'predictive-ll.RData'))) {
    load( here::here('output', 'cross-validate' , 'predictive-ll.RData'))
} else {    
    dat_ll <- data.frame()
    for (k in 1:K) {
        # extract the Matern model predictive log-likelihood
        load(here::here('output', 'cross-validate' , paste0('predictive-ll-matern-cv-', k, '.RData')))
        dimnames(ll_matern) <- list(iteration = 1:nrow(ll_matern),
                                    obs_idx = 1:ncol(ll_matern))
        dat_matern <- as.data.frame.table(ll_matern, responseName = 'll') %>%
            mutate(model = "matern", 
                   fold = k, 
                   iteration = as.numeric(iteration), 
                   obs_idx = as.numeric(obs_idx)) %>% 
            left_join(
                data.frame(
                    obs_idx = 1:ncol(ll_matern),
                    row_id = row_id,
                    col_id = col_id), by = "obs_idx")

        
        # extract the overdispersed model predictive log-likelihood
        load(here::here('output', 'cross-validate' , paste0('predictive-ll-overdispersed-cv-', k, '.RData')))
        dimnames(ll_overdispersed) <- list(iteration = 1:nrow(ll_overdispersed),
                                           obs_idx = 1:ncol(ll_overdispersed))
        dat_overdispersed <- as.data.frame.table(ll_overdispersed, responseName = 'll') %>%
            mutate(model = "overdispersed", 
                   fold = k, 
                   iteration = as.numeric(iteration), 
                   obs_idx = as.numeric(obs_idx)) %>% 
            left_join(
                data.frame(
                    obs_idx = 1:ncol(ll_overdispersed),
                    row_id = row_id,
                    col_id = col_id), by = "obs_idx")
        
        # extract the latent model predictive log-likelihood    
        load(here::here('output', 'cross-validate' , paste0('predictive-ll-latent-cv-', k, '.RData')))
        dimnames(ll_latent) <- list(iteration = 1:nrow(ll_latent),
                                    obs_idx = 1:ncol(ll_latent))
        dat_latent <- as.data.frame.table(ll_latent, responseName = 'll') %>%
            mutate(model = "latent", 
                   fold = k, 
                   iteration = as.numeric(iteration), 
                   obs_idx = as.numeric(obs_idx))%>% 
            left_join(
                data.frame(
                    obs_idx = 1:ncol(ll_latent),
                    row_id = row_id,
                    col_id = col_id), by = "obs_idx")
        
        dat_ll <- rbind(dat_ll, do.call("rbind", list(dat_matern, dat_overdispersed, dat_latent)))
    }
    save(dat_ll, file = here::here('output', 'cross-validate' , 'predictive-ll.RData'))
}


dat_ll %>%
    group_by(model) %>%
    summarize(mean_ll = mean(ll), 
              median_ll = median(ll),
              sd_ll = sd(ll), 
              n = n(),
              ci_lower = mean_ll - 1.96 * sd_ll / sqrt(n/1000),
              ci_upper = mean_ll + 1.96 * sd_ll / sqrt(n/1000)) %>%
    # reorder rows
    slice(match(c("matern", "overdispersed", "latent"), model)) %>%
    select(-n) %>%
    xtable(caption = "Cross-validation results") %>%
    printbold(each = "column", max = c(NA, TRUE, TRUE, rep(NA, 3)), type = "latex",
              file = here::here("tables", "cross-validation-ll.tex"), 
              include.rownames = FALSE)


dat_ll %>%
    group_by(model, fold) %>%
    summarize(mean_ll_fold = mean(ll), 
              median_ll_fold = median(ll),
              sd_ll_fold = sd(ll)) %>%
    ungroup() %>%
    group_by(model) %>%
    summarize(mean_ll = mean(mean_ll_fold),
              sd_ll = sd(mean_ll_fold),
              n = n(), 
              ci_lower = mean_ll + qt(0.025, df = n-1) * sd_ll / sqrt(n),
              ci_upper = mean_ll + qt(0.975, df = n-1) * sd_ll / sqrt(n)) %>%
    ungroup() %>%
    select(-n) %>%
    # reorder rows
    slice(match(c("matern", "overdispersed", "latent"), model)) %>%
    xtable(caption = "Cross-validation results") %>%
    printbold(each = "column", max = c(NA, TRUE, rep(NA, 3)), type = "latex",
              file = here::here("tables", "cross-validation-fold.tex"), 
              include.rownames = FALSE)



p0 <- dat_ll %>%
    group_by(model, fold) %>%
    summarize(mean_ll = mean(ll)) %>%
    ungroup() %>%
    ggplot(aes(x = model, y = mean_ll)) +
    geom_boxplot() +
    geom_point(aes(color = factor(fold)), position = position_dodge(width = 0.1))

p1 <- dat_ll %>%
    group_by(obs_idx, fold, model) %>%
    summarize(mean_ll = mean(ll)) %>%
    pivot_wider(names_from = c(model), values_from = mean_ll) %>%
    ggplot(aes(x = latent, y = matern)) +
    geom_point(alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, color = "red", alpha = 0.5)


p2 <- dat_ll %>%
    group_by(obs_idx, fold, model) %>%
    summarize(mean_ll = mean(ll)) %>%
    pivot_wider(names_from = c(model), values_from = mean_ll) %>%
    ggplot(aes(x = latent, y = overdispersed)) +
    geom_point(alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, color = "red", alpha = 0.5)


p3 <- dat_ll %>%
    group_by(obs_idx, fold, model) %>%
    summarize(mean_ll = mean(ll)) %>%
    pivot_wider(names_from = c(model), values_from = mean_ll) %>%
    ggplot(aes(x = matern, y = overdispersed)) +
    geom_point(alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, color = "red", alpha = 0.5)

p0 / (p1 + p2 + p3)



# find outlying observations
n_outlier <- 100
inliers <- dat_ll %>%
    group_by(model, row_id, col_id) %>%
    summarize(mean_ll = mean(ll)) %>%
    arrange(desc(mean_ll)) %>%
    head(n = n_outlier)

outliers <- dat_ll %>%
    group_by(model, row_id, col_id) %>%
    summarize(mean_ll = mean(ll)) %>%
    arrange(mean_ll) %>%
    head(n = n_outlier)

locs_inliers <- locs[inliers$row_id, ]
locs_outliers <- locs[outliers$row_id, ]

# doesn't seem to be much of a pttern here
plot(locs)
points(locs_inliers, col = "blue")
points(locs_outliers, col = "red")

data.frame(time = c(inliers$col_id, outliers$col_id), 
           observation = rep(c("inlier", "outlier"), each = n_outlier)) %>%
ggplot(aes(x = time, group = observation, fill = observation)) +
    geom_histogram(position = "identity", alpha = 0.5)


y_inlier <- matrix(NA, n_outlier, dim(y)[2])
y_outlier <- matrix(NA, n_outlier, dim(y)[2])
for (i in 1:n_outlier) {
    y_inlier[i, ] <- y[inliers$row_id[i], , inliers$col_id[i]]
    y_outlier[i, ] <- y[outliers$row_id[i], , outliers$col_id[i]]
}
dimnames(y_inlier) <- list(observation = 1:n_outlier, species = 1:dim(y)[2])
dimnames(y_outlier) <- list(observation = 1:n_outlier, species = 1:dim(y)[2])

p_inlier <- as.data.frame.table(y_inlier, responseName = "y") %>%
    ggplot(aes(x = observation, y = y, fill = species)) +
    geom_bar(position = "fill", stat = "identity") +
    ggtitle("inliers")
p_outlier <- as.data.frame.table(y_outlier, responseName = "y") %>%
    ggplot(aes(x = observation, y = y, fill = species)) +
    geom_bar(position = "fill", stat = "identity") +
    ggtitle("outliers")


p_inlier / p_outlier

hist(rowSums(y_outlier), breaks = 100)
hist(rowSums(y_inlier))
hist(apply(y, c(1, 3), sum), breaks = 100)

# TODO ----
# plot held-out observations vs. posterior predictive - DONE
# investigate outlier log-likelhoods - DONE
# change in binning over time - Alissa will do


# predicted vs. fitted values

for (k in 1:K) {
    # set up the cross-validation
    y_train <- y
    y_test <- array(NA, dim = dim(y))
    obs_fold <- obs_idx[cv_idx[, k], ]
    row_id <- c()
    col_id <- c()
    for (i in 1:n_cv) {
        y_train[obs_fold$row_id[i], , obs_fold$col_id[i]] <- NA
        y_test[obs_fold$row_id[i], , obs_fold$col_id[i]] <- y[obs_fold$row_id[i], , obs_fold$col_id[i]] 
        if (!any(is.na(y_test[obs_fold$row_id[i], , obs_fold$col_id[i]]))) {
            row_id <- c(row_id, obs_fold$row_id[i])
            col_id <- c(col_id, obs_fold$col_id[i])
        }
    }
    
    
    # load the Matern model
    X = as.matrix(X)
    locs = as.matrix(locs_scaled)
    
    if (!file.exists( here::here('output', 'cross-validate' , paste0('pis-matern-cv-', k, '.RData')))) {
        load( here::here('output', 'cross-validate' , paste0('matern-cv-', k, '.RData')))
        
        pi_matern <- array(NA, dim = c(dim(out$pi)[1], dim(out$pi)[3], length(row_id)))
        for (i in 1:length(row_id)) {
            pi_matern[, , i] <- out$pi[, row_id[i], , col_id[i]]
        }
        
        dimnames(pi_matern) <- list(
            iteration = 1:dim(pi_matern)[1],
            species = 1:dim(pi_matern)[2],
            obs_idx = 1:dim(pi_matern)[3])
        dat_pi_matern <- as.data.frame.table(pi_matern, responseName = "pi") %>%
            mutate(iteration = as.numeric(iteration), obs_idx = as.numeric(obs_idx)) %>%
            left_join(
                data.frame(
                    obs_idx = 1:dim(pi_matern)[3],
                    row_id = row_id,
                    col_id = col_id), by = "obs_idx") %>% 
            select(-obs_idx) %>%
            mutate(model = "matern")
        
        
        save(dat_pi_matern,
             row_id, col_id,
             file = here::here('output', 'cross-validate' , paste0('pis-matern-cv-', k, '.RData')))
    }
    
    if (!file.exists( here::here('output', 'cross-validate' , paste0('pis-overdispersed-cv-', k, '.RData')))) {
        load( here::here('output', 'cross-validate' , paste0('overdispersed-cv-', k, '.RData')))
        
        pi_overdispersed <- array(NA, dim = c(dim(out$pi)[1], dim(out$pi)[3], length(row_id)))
        for (i in 1:length(row_id)) {
            pi_overdispersed[, , i] <- out$pi[, row_id[i], , col_id[i]]
        }
        
        dimnames(pi_overdispersed) <- list(
            iteration = 1:dim(pi_overdispersed)[1],
            species = 1:dim(pi_overdispersed)[2],
            obs_idx = 1:dim(pi_overdispersed)[3])
        dat_pi_overdispersed <- as.data.frame.table(pi_overdispersed, responseName = "pi") %>%
            mutate(iteration = as.numeric(iteration), obs_idx = as.numeric(obs_idx)) %>%
            left_join(
                data.frame(
                    obs_idx = 1:dim(pi_overdispersed)[3],
                    row_id = row_id,
                    col_id = col_id), by = "obs_idx") %>% 
            select(-obs_idx) %>%
            mutate(model = "overdispersed")
        
        save(dat_pi_overdispersed,
             row_id, col_id,
             file = here::here('output', 'cross-validate' , paste0('pis-overdispersed-cv-', k, '.RData')))
    }
    
    if (!file.exists( here::here('output', 'cross-validate' , paste0('pis-latent-cv-', k, '.RData')))) {
        load( here::here('output', 'cross-validate' , paste0('latent-cv-', k, '.RData')))
        
        pi_latent <- array(NA, dim = c(dim(out$pi)[1], dim(out$pi)[3], length(row_id)))
        for (i in 1:length(row_id)) {
            pi_latent[, , i] <- out$pi[, row_id[i], , col_id[i]]
        }
        
        dimnames(pi_latent) <- list(
            iteration = 1:dim(pi_latent)[1],
            species = 1:dim(pi_latent)[2],
            obs_idx = 1:dim(pi_latent)[3])
        dat_pi_latent <- as.data.frame.table(pi_latent, responseName = "pi") %>%
            mutate(iteration = as.numeric(iteration), obs_idx = as.numeric(obs_idx)) %>%
            left_join(
                data.frame(
                    obs_idx = 1:dim(pi_latent)[3],
                    row_id = row_id,
                    col_id = col_id), by = "obs_idx") %>% 
            select(-obs_idx) %>%
            mutate(model = "latent")
        
        save(dat_pi_latent,
             row_id, col_id,
             file = here::here('output', 'cross-validate' , paste0('pis-latent-cv-', k, '.RData')))
    }
}


if (file.exists( here::here('output', 'cross-validate', 'pis-all.RData'))) {
    load( here::here('output', 'cross-validate' , 'pis-all.RData'))
} else {    
    dat_pis <- data.frame()
    for (k in 1:K) {
        # extract the Matern model predictive log-likelihood
        load(here::here('output', 'cross-validate' , paste0('pis-matern-cv-', k, '.RData')))

        # extract the overdispersed model predictive log-likelihood
        load(here::here('output', 'cross-validate' , paste0('pis-overdispersed-cv-', k, '.RData')))

        # extract the latent model predictive log-likelihood    
        load(here::here('output', 'cross-validate' , paste0('pis-latent-cv-', k, '.RData')))
        
        dat_pis <- rbind(dat_pis, do.call("rbind", list(dat_pi_matern, dat_pi_overdispersed, dat_pi_latent)))
    }
    save(dat_pis, file = here::here('output', 'cross-validate' , 'pis-all.RData'))
}

y_prop <- y
dimnames(y_prop) <- list(row_id = 1:dim(y)[1], species = 1:dim(y)[2], col_id = 1:dim(y)[3])
for (i in 1:dim(y_prop)[1]) {
    for (tt in 1:dim(y_prop)[3]) {
        y_prop[i, , tt] <- y_prop[i, , tt] / sum(y_prop[i, , tt])
    }
    
}
dat_y <- as.data.frame.table(y_prop, responseName = "y") %>%
    mutate(row_id = as.numeric(row_id), col_id = as.numeric(col_id))

dat_plot <- dat_pis %>%
    group_by(row_id, col_id, species, model) %>%
    summarize(pi_mean = mean(pi),
            pi_lower = quantile(pi, prob = 0.025), 
            pi_upper = quantile(pi, prob = 0.975)) %>%
    ungroup() %>%
    left_join(dat_y, by = c("row_id", "col_id", "species"))

p_fitted <- ggplot(dat_plot, aes(x = pi_mean, y = y, color = model)) +
    geom_point(alpha = 0.5) +
    facet_grid(col_id ~ species)
p_fitted


# plot fitted vs. observed for each model ----

# change the time numbering and the species
time_bins <- seq(-285, 21000, by=990)[-1]
time_vec <- paste(time_bins[1:21], "ybp")
names(time_vec) <- paste(1:21)
species_names <- readRDS(here::here('data', 'species-names.RDS'))
species_names <- str_replace(species_names, "\\.", " \n ")
names(species_names) <- paste(1:13)
base_size <- 18



calibration_matern <- dat_plot %>%
    filter(model == "matern") %>%
    ggplot(aes(x = y, y = pi_mean)) +
    geom_point(alpha = 0.5) +
    facet_grid(col_id ~ species, 
               labeller = labeller(.rows = time_vec,
                                   .cols = species_names)) +
    geom_linerange(aes(ymin = pi_lower, ymax = pi_upper)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Matern model") +
    xlab("Observed proportions") +
    ylab("Estimated proportions") +
    theme_bw(base_size = base_size) +
    scale_x_continuous(breaks = c(0.25, 0.75)) +
    scale_y_continuous(breaks = c(0.25, 0.75)) +
    theme(strip.text.y.right = element_text(angle = 0),
          strip.text.x.top = element_text(size = 0.56 * base_size),
          axis.text.y = element_text(size = 0.5 * base_size))


ggsave(calibration_matern, 
       file = here::here("figures", "matern-calibration.png"),
       width = 16, height = 9, device = "png", units = "in")

calibration_overdispersed <- dat_plot %>%
    filter(model == "overdispersed") %>%
    ggplot(aes(x = y, y = pi_mean)) +
    geom_point(alpha = 0.5) +
    facet_grid(col_id ~ species, 
               labeller = labeller(.rows = time_vec,
                                   .cols = species_names)) +
    geom_linerange(aes(ymin = pi_lower, ymax = pi_upper)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Overdispersed model") +
    xlab("Observed proportions") +
    ylab("Estimated proportions") +
    theme_bw(base_size = base_size) +
    scale_x_continuous(breaks = c(0.25, 0.75)) +
    scale_y_continuous(breaks = c(0.25, 0.75)) +  
    theme(strip.text.y.right = element_text(angle = 0),
          strip.text.x.top = element_text(size = 0.56 * base_size),
          axis.text.y = element_text(size = 0.5 * base_size))


ggsave(calibration_overdispersed, 
       file = here::here("figures", "overdispersed-calibration.png"),
       width = 16, height = 9, device = "png", units = "in")


calibration_latent <- dat_plot %>%
    filter(model == "latent") %>%
    ggplot(aes(x = y, y = pi_mean)) +
    geom_point(alpha = 0.5) +
    facet_grid(col_id ~ species, 
               labeller = labeller(.rows = time_vec,
                                   .cols = species_names)) +
    geom_linerange(aes(ymin = pi_lower, ymax = pi_upper)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Latent model") +
    xlab("Latent proportions") +
    ylab("Estimated proportions") +
    theme_bw(base_size = base_size) +
    scale_x_continuous(breaks = c(0.25, 0.75)) +
    scale_y_continuous(breaks = c(0.25, 0.75)) +
    theme(strip.text.y.right = element_text(angle = 0),
          strip.text.x.top = element_text(size = 0.56 * base_size),
          axis.text.y = element_text(size = 0.5 * base_size))

# calibration_latent

ggsave(calibration_latent, 
       file = here::here("figures", "latent-calibration.png"),
       width = 16, height = 9, device = "png", units = "in")

dat_all <- dat_ll %>%
    group_by(row_id, col_id, model, fold) %>%
    summarize(mean_ll = mean(ll)) %>%
    ungroup() %>%
    left_join(dat_plot, by = c("row_id", "col_id", "model")) %>%
    mutate(type = case_when(mean_ll < -2500 ~ "outlier", TRUE ~ "inlier"))

p_all <- ggplot(dat_all, aes(x = pi_mean, y = y, color = type)) +
    geom_point(aes(alpha = type)) +
    scale_alpha_discrete(range = c(0.125, 0.25)) +
    facet_grid(col_id ~ species)
p_all

ggsave(
    filename = here::here("figures", "cross-validation-outliers.png"),
    p_all, width = 16, height = 9, units = "in")

# subset for the inliers -- need to better subset these
p_inlier <- dat_plot %>%
    filter(
        (row_id == inliers$row_id[1] & col_id == inliers$col_id[1]) |
            (row_id == inliers$row_id[2] & col_id == inliers$col_id[2]) |
            (row_id == inliers$row_id[3] & col_id == inliers$col_id[3]) |
            (row_id == inliers$row_id[4] & col_id == inliers$col_id[4]) |
            (row_id == inliers$row_id[5] & col_id == inliers$col_id[5]) |
            (row_id == inliers$row_id[6] & col_id == inliers$col_id[6]) |
            (row_id == inliers$row_id[7] & col_id == inliers$col_id[7]) |
            (row_id == inliers$row_id[8] & col_id == inliers$col_id[8]) |
            (row_id == inliers$row_id[9] & col_id == inliers$col_id[9]) |
            (row_id == inliers$row_id[10] & col_id == inliers$col_id[10]) ) %>%
    ggplot(aes(x = y, y = pi_mean, color = model)) +
    geom_point(alpha = 0.5) +
    facet_grid(col_id ~ species) +
    geom_linerange(aes(ymin = pi_lower, ymax = pi_upper), alpha = 0.25) +
    geom_abline(slope = 1, intercept = 0, color = "red")


# subset for the outliers
p_outlier <- dat_plot %>%
    filter(
        (row_id == outliers$row_id[1] & col_id == outliers$col_id[1]) |
            (row_id == outliers$row_id[2] & col_id == outliers$col_id[2]) |
            (row_id == outliers$row_id[3] & col_id == outliers$col_id[3]) |
            (row_id == outliers$row_id[4] & col_id == outliers$col_id[4]) |
            (row_id == outliers$row_id[5] & col_id == outliers$col_id[5]) |
            (row_id == outliers$row_id[6] & col_id == outliers$col_id[6]) |
            (row_id == outliers$row_id[7] & col_id == outliers$col_id[7]) |
            (row_id == outliers$row_id[8] & col_id == outliers$col_id[8]) |
            (row_id == outliers$row_id[9] & col_id == outliers$col_id[9]) |
            (row_id == outliers$row_id[10] & col_id == outliers$col_id[10]) ) %>%
    ggplot(aes(x = y, y = pi_mean, color = model)) +
    geom_point(alpha = 0.5) +
    facet_grid(col_id ~ species) +
    geom_linerange(aes(ymin = pi_lower, ymax = pi_upper), alpha = 0.25) +
    geom_abline(slope = 1, intercept = 0, color = "red")

p_inlier / p_outlier
    

# compare models without the "outliers"

p_no_outlier <- dat_ll %>%
    group_by(model, fold, row_id, col_id) %>%
    summarize(mean_ll = mean(ll)) %>%
    mutate(type = case_when(mean_ll < -750 ~ "outlier", TRUE ~ "inlier")) %>%
    filter(type == "inlier") %>%
    ungroup() %>%
    ggplot(aes(x = model, y = mean_ll)) +
    geom_boxplot() +
    geom_point(aes(color = factor(fold)), position = position_dodge(width = 0.1))
p_no_outlier
