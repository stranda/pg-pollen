library(pgR)
library(ggplot2)
library(fields)
library(rgdal)
library(loo)
library(tidyverse)

version <- '3.1'
dat <- readRDS(here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))

# check the Matern model

out_matern <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '.RDS')))

ll_matern_raw <- calc_ll_pg_stlm(dat$y, dat$X, out_matern)
# reduce memory usage if needed
# rm(out_matern)

# check for convergence of log-likelihood
layout(matrix(1:4, 2, 2))
matplot(ll_matern_raw$ll[, , 1], type = 'l')
matplot(ll_matern_raw$ll[, , 2], type = 'l')
matplot(ll_matern_raw$ll[, , 3], type = 'l')
matplot(ll_matern_raw$ll[, , 4], type = 'l')

# convert ll to a matrix for loo and WAIC (need to figure out missing values)
ll_matern <- t(apply(ll_matern_raw$ll, 1, c))
# drop the missing values
na_idx <- rep(NA, ncol(ll_matern))
for (i in 1:ncol(ll_matern)) {
    na_idx[i] <- any(is.na(ll_matern[, i]))
}
ll_matern <- ll_matern[, !na_idx]


# check the overdispersed Matern model

out_matern_overdispersed <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_overdispersed.RDS')))

ll_matern_overdispersed_raw <- calc_ll_pg_stlm(dat$y, dat$X, out_matern_overdispersed)
# reduce memory usage if needed
# rm(out_matern)

# check for convergence of log-likelihood
layout(matrix(1:4, 2, 2))
matplot(ll_matern_overdispersed_raw$ll[, , 1], type = 'l')
matplot(ll_matern_overdispersed_raw$ll[, , 2], type = 'l')
matplot(ll_matern_overdispersed_raw$ll[, , 3], type = 'l')
matplot(ll_matern_overdispersed_raw$ll[, , 4], type = 'l')

# convert ll to a matrix for loo and WAIC (need to figure out missing values)
ll_matern_overdispersed <- t(apply(ll_matern_overdispersed_raw$ll, 1, c))
# drop the missing values
na_idx <- rep(NA, ncol(ll_matern_overdispersed))
for (i in 1:ncol(ll_matern_overdispersed)) {
    na_idx[i] <- any(is.na(ll_matern_overdispersed[, i]))
}
ll_matern_overdispersed <- ll_matern_overdispersed[, !na_idx]


# # check the MRA model
# 
# out_mra <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_MRA.RDS')))
# 
# # temporary class fix -- fix the calc_ll_pg_stlm function
# # class(out_mra) <- "pg_stlm"
# 
# ll_mra_raw <- calc_ll_pg_stlm(dat$y, dat$X, out_mra)
# 
# # reduce memory usage if needed
# # rm(out_mra)
# 
# # check for convergence of log-likelihood
# layout(matrix(1:4, 2, 2))
# matplot(ll_mra_raw$ll[, , 1], type = 'l')
# matplot(ll_mra_raw$ll[, , 2], type = 'l')
# matplot(ll_mra_raw$ll[, , 3], type = 'l')
# matplot(ll_mra_raw$ll[, , 4], type = 'l')
# 
# # convert ll to a matrix for loo and WAIC (need to figure out missing values)
# ll_mra <- t(apply(ll_mra_raw$ll, 1, c))
# # drop the missing values
# na_idx <- rep(NA, ncol(ll_mra))
# for (i in 1:ncol(ll_mra)) {
#     na_idx[i] <- any(is.na(ll_mra[, i]))
# }
# ll_mra <- ll_mra[, !na_idx]





# WAIC using LaplacesDemon package
WAIC_matern <- LaplacesDemon::WAIC(ll_matern)
WAIC_matern_overdispersed <- LaplacesDemon::WAIC(ll_matern_overdispersed)
# WAIC_mra <- LaplacesDemon::WAIC(ll_mra)
WAIC_matern$WAIC
WAIC_matern_overdispersed$WAIC
# WAIC_mra$WAIC

if (WAIC_matern$WAIC < WAIC_matern_overdispersed$WAIC) {
    message("Matern Model Fits Better")
} else {
    message("Overdispersed Matern Model Fits Better")
}
# if (WAIC_matern$WAIC < WAIC_mra$WAIC) {
#     message("Matern Model ")
# }


# I think the warnings from loo come from the autocorrelation
# I'm not really sure how to address this going forward but maybe the ideas
# presented here are a good start... https://mc-stan.org/loo/articles/loo2-non-factorizable.html

# WAIC using loo package
loo_waic_matern <- waic(ll_matern)
loo_waic_matern_oversdispersed <- waic(ll_matern_overdispersed)
# loo_waic_mra <- waic(ll_mra)
# loo_compare(loo_waic_matern, loo_waic_mra)
loo_compare(loo_waic_matern, loo_waic_matern_oversdispersed)

r_eff_matern <- relative_eff(ll_matern, chain_id = rep(1, nrow(ll_matern)))
r_eff_matern_overdispersed <- relative_eff(ll_matern_overdispersed, chain_id = rep(1, nrow(ll_matern_overdispersed)))
# r_eff_mra <- relative_eff(ll_mra, chain_id = rep(1, nrow(ll_mra)))
loo_matern <- loo(ll_matern, cores = 32L, r_eff = r_eff_matern)
loo_matern_overdispersed <- loo(ll_matern_overdispersed, cores = 32L, r_eff = r_eff_matern_overdispersed)
# loo_mra <- loo(ll_mra, cores = 32L, r_eff = r_eff_mra)
comp <- loo_compare(loo_matern, loo_matern_overdispersed)
# comp <- loo_compare(loo_matern, loo_mra)

# the loo diagnostics don't look good -- maybe if we increase the sample size?
print(comp)

layout(matrix(1:2, 2, 1))
plot(loo_matern)
plot(loo_matern_overdispersed)
# plot(loo_mra)

dat1 <- data.frame(k = loo_matern$diagnostics$pareto_k, model = "matern") %>%
    mutate(sample_id = 1:length(loo_matern$diagnostics$pareto_k))
dat2 <- data.frame(k = loo_matern_overdispersed$diagnostics$pareto_k, model = "overdispersed") %>%
    mutate(sample_id = 1:length(loo_matern_overdispersed$diagnostics$pareto_k))

# on average, the Matern k parameter is larger
rbind(dat1, dat2) %>%
    pivot_wider(names_from = model, values_from = k) %>%
    mutate(k_diff = matern - overdispersed) %>%
    summarize(mean_diff = mean(k_diff))

rbind(dat1, dat2) %>%
    pivot_wider(names_from = model, values_from = k) %>%
    mutate(k_diff = matern - overdispersed) %>%
    ggplot(aes(x = k_diff)) +
    geom_histogram()

rbind(dat1, dat2) %>%
    ggplot(aes(x = k, fill = model)) +
    geom_density(alpha = 0.5)
