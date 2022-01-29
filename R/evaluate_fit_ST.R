library(pgR)
library(ggplot2)
library(fields)
library(rgdal)
library(loo)
library(tidyverse)

source(here::here("R", "xtable_print_bold.R"))

version <- '3.1'
dat <- readRDS(here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))

# check the Matern model -------------------------------------------------------

if (file.exists(here::here("output", paste0('polya-gamma-ll_', version,'.RData')))) {
    load(here::here("output", paste0('polya-gamma-ll_', version,'.RData')))
} else {
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
    save(ll_matern_raw, ll_matern, file = here::here("output", paste0('polya-gamma-ll_', version,'.RData')))
}

# check the overdispersed Matern model -----------------------------------------

if (file.exists(here::here("output", paste0('polya-gamma-ll_', version,'_overdispersed.RData')))) {
    load(here::here("output", paste0('polya-gamma-ll_', version,'_overdispersed.RData')))
} else {
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

    save(ll_matern_overdispersed_raw, ll_matern_overdispersed, file = here::here("output", paste0('polya-gamma-ll_', version,'_overdispersed.RData')))
}

# check the latent overdispersed Matern model ----------------------------------
if (file.exists(here::here("output", paste0('polya-gamma-ll_', version,'_latent_overdispersed.RData')))) {
    load(here::here("output", paste0('polya-gamma-ll_', version,'_latent_overdispersed.RData')))
} else {
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
    save(ll_matern_latent_raw, ll_matern_latent, file = here::here("output", paste0('polya-gamma-ll_', version,'_latent_overdispersed.RData')))
}

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

dat_scores %>%
    # reorder rows
    slice(match(c("matern", "overdispersed", "latent"), model)) %>%
    # reorder columns
    select(model, everything()) %>%
    xtable(caption = "WAIC results") %>%
    printbold(each = "column", max = c(NA, FALSE, rep(NA, 3)), type = "latex",
              file = here::here("tables", "WAIC.tex"), 
              include.rownames = FALSE)



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


# Calculate loo-ll ------------------------------------------------------------
if (file.exists(here::here("output", paste0('polya-gamma-ll_loo_', version,'.RData')))) {
    load(here::here("output", paste0('polya-gamma-ll_loo_', version,'.RData')))
} else {
    out_matern <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '.RDS')))
    
    # start time 4:23PM 2021-08-27
    ll_matern_loo <- calc_ll_pg_stlm_loo(dat$y, dat$X, dat$locs, out_matern, n_message = 1)
    # finished time about 2:40PM 2021-08-28
    
    save(ll_matern_loo, file = here::here("output", paste0('polya-gamma-ll_loo_', version,'.RData')))
    pushoverr::pushover(message = "Finished ll-loo Matern model")
}

# sanity check of results



dat_ll <- data.frame(ll = c(apply(ll_matern_loo$ll, c(2, 3), mean)),
                     ll_loo = c(apply(ll_matern_loo$ll_loo, c(2, 3), mean)))

ggplot(dat_ll, aes(x = ll, y = ll_loo)) +
    geom_point(alpha = 0.5)
    # dat <- data.frame(ll = c(ll_matern_loo$ll), ll_loo = c(ll_matern_loo$ll_loo))
# ggplot(dat, aes(x = ll, y = ll_loo)) +
#     geom_hex()

dat_pi <- data.frame(pi = c(apply(ll_matern_loo$pi, c(2, 3, 4), mean)),
                     pi_loo = c(apply(ll_matern_loo$pi_loo, c(2, 3, 4), mean)))

ggplot(dat_pi, aes(x = pi, y = pi_loo)) +
    geom_point(alpha = 0.5)



# dat <- data.frame(pi = c(ll_matern_loo$pi), pi_loo = c(ll_matern_loo$pi_loo))
# ggplot(dat, aes(x = pi, y = pi_loo)) +
#     geom_hex()



if (file.exists(here::here("output", paste0('polya-gamma-ll_loo_', version,'_overdispersed.RData')))) {
    load(here::here("output", paste0('polya-gamma-ll_loo_', version,'_overdispersed.RData')))
} else {
    out_overdispersed <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_overdispersed.RDS')))
    
    
    ll_overdispersed_loo <- calc_ll_pg_stlm_loo(dat$y, dat$X, dat$locs, out_matern, n_message = 1)
    
    save(ll_overdispersed_loo, file = here::here("output", paste0('polya-gamma-ll_loo_', version,'_overdispersed.RData'))) 
    pushoverr::pushover(message = "Finished ll-loo overdispersed model")
}

# sanity check of results
dat <- data.frame(ll = c(ll_overdispersed_loo$ll), ll_loo = c(ll_overdispersed_loo$ll_loo))
ggplot(dat, aes(x = ll, y = ll_loo)) +
    geom_hex()

dat_pi <- data.frame(pi = c(apply(ll_overdispersed_loo$pi, c(2, 3, 4), mean)),
                     pi_loo = c(apply(ll_overdispersed_loo$pi_loo, c(2, 3, 4), mean)))

ggplot(dat_pi, aes(x = pi, y = pi_loo)) +
    geom_point(alpha = 0.5)

# dat <- data.frame(pi = c(ll_overdispersed_loo$pi), pi_loo = c(ll_overdispersed_loo$pi_loo))
# ggplot(dat, aes(x = pi, y = pi_loo)) +
#     geom_hex()

if (file.exists(here::here("output", paste0('polya-gamma-ll_loo_', version,'_latent.RData')))) {
    load(here::here("output", paste0('polya-gamma-ll_loo_', version,'_latent.RData')))
} else {
    out_latent <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_latent_overdispersed.RDS')))
    
    # start time 4:23PM 2021-08-27
    ll_latent_loo <- calc_ll_pg_stlm_loo(dat$y, dat$X, dat$locs, out_latent, n_message = 1)
    # finished time about 2:40PM 2021-08-28
    
    save(ll_latent_loo, file = here::here("output", paste0('polya-gamma-ll_loo_', version,'_latent.RData')))
    pushoverr::pushover(message = "Finished ll-loo latent model")
}

# sanity check of results
dat <- data.frame(ll = c(ll_latent_loo$ll), ll_loo = c(ll_latent_loo$ll_loo))
ggplot(dat, aes(x = ll, y = ll_loo)) +
    geom_hex()

dat_pi <- data.frame(pi = c(apply(ll_latent_loo$pi, c(2, 3, 4), mean)),
                     pi_loo = c(apply(ll_latent_loo$pi_loo, c(2, 3, 4), mean)))

ggplot(dat_pi, aes(x = pi, y = pi_loo)) +
    geom_point(alpha = 0.5)

# dat <- data.frame(pi = c(ll_latent_loo$pi), pi_loo = c(ll_latent_loo$pi_loo))
# ggplot(dat, aes(x = pi, y = pi_loo)) +
#     geom_hex()


# Try loo with loo-likelihood

# convert ll to a matrix for loo and WAIC (need to figure out missing values)
ll_matern_loo_mat <- t(apply(ll_matern_loo$ll_loo, 1, c))
# drop the missing values
na_idx <- rep(NA, ncol(ll_matern_loo_mat))
for (i in 1:ncol(ll_matern_loo_mat)) {
    na_idx[i] <- any(is.na(ll_matern_loo_mat[, i]))
}
ll_matern_loo_mat <- ll_matern_loo_mat[, !na_idx]


ll_overdispersed_loo_mat <- t(apply(ll_overdispersed_loo$ll_loo, 1, c))
# drop the missing values
na_idx <- rep(NA, ncol(ll_overdispersed_loo_mat))
for (i in 1:ncol(ll_overdispersed_loo_mat)) {
    na_idx[i] <- any(is.na(ll_overdispersed_loo_mat[, i]))
}
ll_overdispersed_loo_mat <- ll_overdispersed_loo_mat[, !na_idx]


ll_latent_loo_mat <- t(apply(ll_latent_loo$ll_loo, 1, c))
# drop the missing values
na_idx <- rep(NA, ncol(ll_latent_loo_mat))
for (i in 1:ncol(ll_latent_loo_mat)) {
    na_idx[i] <- any(is.na(ll_latent_loo_mat[, i]))
}
ll_latent_loo_mat <- ll_latent_loo_mat[, !na_idx]

# I think the warnings from loo come from the autocorrelation
# I'm not really sure how to address this going forward but maybe the ideas
# presented here are a good start... https://mc-stan.org/loo/articles/loo2-non-factorizable.html

# WAIC using loo package
loo_waic_matern <- waic(ll_matern_loo_mat)
loo_waic_matern_oversdispersed <- waic(ll_overdispersed_loo_mat)
loo_waic_matern_latent <- waic(ll_latent_loo_mat)
# loo_waic_mra <- waic(ll_mra)
# loo_compare(loo_waic_matern, loo_waic_mra)
loo_compare(loo_waic_matern, loo_waic_matern_oversdispersed, loo_waic_matern_latent)

r_eff_matern <- relative_eff(ll_matern_loo_mat, chain_id = rep(1, nrow(ll_matern_loo_mat)))
r_eff_matern_overdispersed <- relative_eff(ll_overdispersed_loo_mat, chain_id = rep(1, nrow(ll_overdispersed_loo_mat)))
r_eff_matern_latent <- relative_eff(ll_latent_loo_mat, chain_id = rep(1, nrow(ll_latent_loo_mat)))
# r_eff_mra <- relative_eff(ll_mra, chain_id = rep(1, nrow(ll_mra)))
loo_matern <- loo(ll_matern_loo_mat, cores = 32L, r_eff = r_eff_matern)
loo_matern_overdispersed <- loo(ll_overdispersed_loo_mat, cores = 32L, r_eff = r_eff_matern_overdispersed)
loo_matern_latent <- loo(ll_latent_loo_mat, cores = 32L, r_eff = r_eff_matern_latent)
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
