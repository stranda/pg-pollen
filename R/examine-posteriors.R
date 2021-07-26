library(tidyverse)
library(pgR)
library(latex2exp)


version <- "3.1"

out <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '.RDS')))
out_mra <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_MRA.RDS')))
dat <- readRDS(here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))

source("R/plot-trace.R")
plot_trace(out, base_size = 7, file = here::here("figures", "matern", "trace-plots.png"))
plot_trace(out_mra, base_size = 7, file = here::here("figures", "MRA", "trace-plots.png"))

matplot(out$beta[, 1, ], type = 'l')
matplot(out_mra$beta[, 1, ], type = 'l')

layout(matrix(1:4, 2, 2))
matplot(exp(out$theta[, , 1]) / 1e3, type = 'l')
matplot(exp(out$theta[, , 2]), type = 'l')


matplot(out$tau2, type = 'l')
layout(matrix(1:3, 1, 3))
matplot(out_mra$tau2[, 1, ], type = 'l')
matplot(out_mra$tau2[, 2, ], type = 'l')
matplot(out_mra$tau2[, 3, ], type = 'l')


matplot(out_mra$sigma2, type = 'l')

hist(out$eta)

## Calculate Mi
N <- dim(dat$y)[1]
J <- dim(dat$y)[2]
n_time <- dim(dat$y)[3]
Y <- dat$y
missing_idx <- matrix(FALSE, N, n_time)
for (i in 1:N) {
    for (tt in 1:n_time) {
        missing_idx[i, tt] <- any(is.na(Y[i, , tt]))
    }
}

Mi    <- array(0, dim = c(N, J - 1, n_time))
kappa <- array(0, dim = c(N, J - 1, n_time))
for (tt in 1:n_time) {
    Mi[, , tt] <- calc_Mi(Y[, , tt])
    kappa[, , tt] <- calc_kappa(Y[, , tt], Mi[, , tt])
}



plot(apply(out$eta, c(2, 3, 4), mean), kappa)
plot(apply(out_mra$eta, c(2, 3, 4), mean), kappa)

plot(apply(out$eta, c(2, 3, 4), mean), apply(out$omega, c(2, 3, 4), mean))
plot(apply(out_mra$eta, c(2, 3, 4), mean), apply(out_mra$omega, c(2, 3, 4), mean))
plot(apply(out$omega, c(2, 3, 4), mean), kappa)
plot(apply(out_mra$omega, c(2, 3, 4), mean), kappa)

cor(apply(out$omega, c(2, 3, 4), mean), kappa)
cor(apply(out_mra$omega, c(2, 3, 4), mean), kappa)


plot(apply(out$eta, c(2, 3, 4), mean), kappa)
plot(apply(out_mra$eta, c(2, 3, 4), mean), kappa)



# check the predicted probabilities between the different models

preds <- readRDS(here::here("output", paste0('polya-gamma-predictions_', version, '.RDS')))
preds_mra <- readRDS(here::here("output", paste0('polya-gamma-predictions_', version, '_MRA.RDS')))

pi_mean <- apply(preds$pi, c(2, 3, 4), mean)
pi_mra_mean <- apply(preds_mra$pi, c(2, 3, 4), mean)
dimnames(pi_mean) <- list(
    location = 1:dim(pi_mean)[1],
    species = 1:dim(pi_mean)[2],
    time = 1:dim(pi_mean)[3]
)
dimnames(pi_mra_mean) <- list(
    location = 1:dim(pi_mra_mean)[1],
    species = 1:dim(pi_mra_mean)[2],
    time = 1:dim(pi_mra_mean)[3]
)

dat_pi <- as.data.frame.table(pi_mean, responseName = "pi") %>%
    mutate(model = "Matern", location = as.numeric(location))
dat_pi_mra <- as.data.frame.table(pi_mra_mean, responseName = "pi") %>%
    mutate(model = "MRA", location = as.numeric(location))

y_prop = dat$y
for(i in 1:dim(y_prop)[1]) {
    for (tt in 1:dim(y_prop)[3]) {
        y_prop[i, , tt] <- y_prop[i, , tt] / sum(y_prop[i, , tt])
    }
}
dimnames(y_prop) <- list(
    location = 1:dim(y_prop)[1],
    species = 1:dim(y_prop)[2],
    time = 1:dim(y_prop)[3]
)
dat_obs <- as.data.frame.table(y_prop, responseName = "pi") %>%
    mutate(model = "observed", location = as.numeric(location))

rbind(dat_pi, dat_pi_mra, dat_obs) %>% 
    pivot_wider(names_from = model, values_from = pi) %>%
    ggplot(aes(x = Matern, y = MRA)) +
    geom_hex() +
    facet_grid(species ~ time)


rbind(dat_pi, dat_pi_mra, dat_obs) %>% 
    pivot_wider(names_from = model, values_from = pi) %>%
    ggplot(aes(x = Matern, y = observed)) +
    geom_hex() +
    facet_grid(species ~ time)

rbind(dat_pi, dat_pi_mra, dat_obs) %>% 
    pivot_wider(names_from = model, values_from = pi) %>%
    ggplot(aes(x = MRA, y = observed)) +
    geom_hex() +
    facet_grid(species ~ time)

plot(pi_mean, pi_mra_mean)
plot(pi_mean[, 1, 2], pi_mra_mean[, 1, 2])

