library(pgR)
library(ggplot2)
library(fields)
library(rgdal)
library(loo)
library(tidyverse)

version <- '4.0'
dat <- readRDS(here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))

# check the Matern model -------------------------------------------------------

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


# check the overdispersed Matern model -----------------------------------------

out_matern_overdispersed <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_overdispersed.RDS')))

ll_matern_overdispersed_raw <- calc_ll_pg_stlm(dat$y, dat$X, out_matern_overdispersed)
# reduce memory usage if needed
# rm(out_matern_overdispersed)

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

# check the latent overdispersed Matern model ----------------------------------

out_matern_latent <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_latent_overdispersed.RDS')))

ll_matern_latent_raw <- calc_ll_pg_stlm(dat$y, dat$X, out_matern_latent)
# reduce memory usage if needed
# rm(out_matern_latent)

# check for convergence of log-likelihood
layout(matrix(1:4, 2, 2))
matplot(ll_matern_latent_raw$ll[, , 1], type = 'l')
matplot(ll_matern_latent_raw$ll[, , 2], type = 'l')
matplot(ll_matern_latent_raw$ll[, , 3], type = 'l')
matplot(ll_matern_latent_raw$ll[, , 4], type = 'l')

# convert ll to a matrix for loo and WAIC (need to figure out missing values)
ll_matern_latent <- t(apply(ll_matern_latent_raw$ll, 1, c))
# drop the missing values
na_idx <- rep(NA, ncol(ll_matern_latent))
for (i in 1:ncol(ll_matern_latent)) {
    na_idx[i] <- any(is.na(ll_matern_latent[, i]))
}
ll_matern_latent <- ll_matern_latent[, !na_idx]


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
WAIC_matern_latent <- LaplacesDemon::WAIC(ll_matern_latent)
# WAIC_mra <- LaplacesDemon::WAIC(ll_mra)
WAIC_matern$WAIC
WAIC_matern_overdispersed$WAIC
WAIC_matern_latent$WAIC

dat_scores <- data.frame(
    WAIC = c(WAIC_matern$WAIC, WAIC_matern_overdispersed$WAIC, WAIC_matern_latent$WAIC),
    lppd = c(WAIC_matern$lppd, WAIC_matern_overdispersed$lppd, WAIC_matern_latent$lppd),
    pWAIC = c(WAIC_matern$pWAIC, WAIC_matern_overdispersed$pWAIC, WAIC_matern_latent$pWAIC),
    pWAIC1 = c(WAIC_matern$pWAIC1, WAIC_matern_overdispersed$pWAIC1, WAIC_matern_latent$pWAIC1),
    model = factor(c("matern", "overdispersed", "latent"), levels = c("matern", "overdispersed", "latent"))) 

dat_scores %>%
    pivot_longer(cols = 1:4, names_to = "score") %>%
    ggplot(aes(x = model, y = value)) +
    geom_col() +
    facet_wrap(~ score, scales = "free")

cut_off <- 200000
p_WAIC <- dat_scores %>%
    mutate(WAIC = WAIC - cut_off) %>%
    ggplot(aes(x = model, y = WAIC)) +
    geom_text(aes(label = round(WAIC + cut_off, digits = 2), vjust = 0)) +
    scale_y_continuous(labels = function(x) x + cut_off)  +
    geom_col() +
    ggtitle("WAIC") +
    theme_bw(base_size = 22)

ggsave(p_WAIC, 
       file = here::here("figures", "WAIC-scores.png"),
       width = 16, height = 9, device = "png", units = "in")

cut_off <- -6000
dat_scores %>%
    mutate(lppd = lppd - cut_off) %>%
    ggplot(aes(x = model, y = lppd)) +
    geom_text(aes(label = round(lppd + cut_off, digits = 2), vjust = 1)) +
    scale_y_continuous(labels = function(x) x + cut_off)  +
    geom_col() +
    ggtitle("lppd")


library(dplyr)
library(ggplot2)

cut.off = 500                                            # (= 1 in your use case)

diamonds %>%
    filter(clarity %in% c("SI1", "VS2")) %>%
    count(cut, clarity) %>%
    mutate(n = n - cut.off) %>%                            # subtract cut.off from y values
    ggplot(aes(x = cut, y = n, fill = cut)) +
    geom_col() +
    geom_text(aes(label = n + cut.off,                     # label original values (optional)
                  vjust = ifelse(n > 0, 0, 1))) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(labels = function(x) x + cut.off) + # add cut.off to label values
    facet_grid(clarity ~ .)

# WAIC_mra$WAIC



if (WAIC_matern$WAIC < WAIC_matern_overdispersed$WAIC & WAIC_matern$WAIC < WAIC_matern_latent$WAIC) {
    message("Matern Model Fits Better")
} else if (WAIC_matern_overdispersed$WAIC < WAIC_matern_latent$WAIC) {
    message("Overdispersed Matern Model Fits Better")
} else {
    message("Latent Matern Model Fits Better")
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
loo_waic_matern_latent <- waic(ll_matern_latent)
# loo_waic_mra <- waic(ll_mra)
# loo_compare(loo_waic_matern, loo_waic_mra)
loo_compare(loo_waic_matern, loo_waic_matern_oversdispersed, loo_waic_matern_latent)

r_eff_matern <- relative_eff(ll_matern, chain_id = rep(1, nrow(ll_matern)))
r_eff_matern_overdispersed <- relative_eff(ll_matern_overdispersed, chain_id = rep(1, nrow(ll_matern_overdispersed)))
r_eff_matern_latent <- relative_eff(ll_matern_latent, chain_id = rep(1, nrow(ll_matern_latent)))
# r_eff_mra <- relative_eff(ll_mra, chain_id = rep(1, nrow(ll_mra)))
loo_matern <- loo(ll_matern, cores = 32L, r_eff = r_eff_matern)
loo_matern_overdispersed <- loo(ll_matern_overdispersed, cores = 32L, r_eff = r_eff_matern_overdispersed)
loo_matern_latent <- loo(ll_matern_latent, cores = 32L, r_eff = r_eff_matern_latent)
# loo_mra <- loo(ll_mra, cores = 32L, r_eff = r_eff_mra)
comp <- loo_compare(loo_matern, loo_matern_overdispersed, loo_matern_latent)
# comp <- loo_compare(loo_matern, loo_mra)

# the loo diagnostics don't look good -- maybe if we increase the sample size?
print(comp)

layout(matrix(1:3, 3, 1))
plot(loo_matern)
plot(loo_matern_overdispersed)
plot(loo_matern_latent)
# plot(loo_mra)

dat1 <- data.frame(k = loo_matern$diagnostics$pareto_k, model = "matern") %>%
    mutate(sample_id = 1:length(loo_matern$diagnostics$pareto_k))
dat2 <- data.frame(k = loo_matern_overdispersed$diagnostics$pareto_k, model = "overdispersed") %>%
    mutate(sample_id = 1:length(loo_matern_overdispersed$diagnostics$pareto_k))
dat3 <- data.frame(k = loo_matern_latent$diagnostics$pareto_k, model = "latent") %>%
    mutate(sample_id = 1:length(loo_matern_latent$diagnostics$pareto_k))

# on average, the Matern k parameter is larger
rbind(dat1, dat2, dat3) %>%
    pivot_wider(names_from = model, values_from = k) %>%
    mutate(k_diff_1 = matern - overdispersed,
           k_diff_2 = matern - latent,
           k_diff_3 = overdispersed - latent) %>%
    summarize(mean_diff_1 = mean(k_diff_1),
              mean_diff_2 = mean(k_diff_2),
              mean_diff_3 = mean(k_diff_3))

rbind(dat1, dat2, dat3) %>%
    pivot_wider(names_from = model, values_from = k) %>%
    mutate(k_diff = matern - overdispersed) %>%
    ggplot(aes(x = k_diff)) +
    geom_histogram()

rbind(dat1, dat2, dat3) %>%
    ggplot(aes(x = k, fill = model)) +
    geom_density(alpha = 0.5)
