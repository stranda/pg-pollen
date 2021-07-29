library(tidyverse)
library(pgR)
library(latex2exp)


version <- "3.1"

out <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '.RDS')))
out_overdispersed <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_overdispersed.RDS')))
out_latent <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_latent_overdispersed.RDS')))
# out_mra <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_MRA.RDS')))
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

# preds <- readRDS(here::here("output", paste0('polya-gamma-predictions_', version, '.RDS')))
# preds_mra <- readRDS(here::here("output", paste0('polya-gamma-predictions_', version, '_MRA.RDS')))

pi_mean <- apply(out$pi, c(2, 3, 4), mean)
pi_mean_overdispersed <- apply(out_overdispersed$pi, c(2, 3, 4), mean)
pi_mean_latent <- apply(out_latent$pi, c(2, 3, 4), mean)

# latent values only -- no overdispersion
# X <- matrix(rep(1, nrow(dat$y)), nrow(dat$y), 1)
# 
# pi_latent <- array(0, dim = c(nrow(out_latent$beta), dim(dat$y)))
# for (k in 1:nrow(out_latent$beta)) {
#     for (tt in 1:dim(dat$y)[3]) {
#         pi_latent[k, , , tt] <- eta_to_pi(X %*% out_latent$beta[k, , ] + out_latent$psi[k, , , tt])
#     }
# }
# pi_mean_latent <- apply(pi_latent, c(2, 3, 4), mean)

# pi_mra_mean <- apply(out_mra$pi, c(2, 3, 4), mean)
dimnames(pi_mean) <- list(
    location = 1:dim(pi_mean)[1],
    species = 1:dim(pi_mean)[2],
    time = 1:dim(pi_mean)[3]
)
dimnames(pi_mean_overdispersed) <- list(
    location = 1:dim(pi_mean_overdispersed)[1],
    species = 1:dim(pi_mean_overdispersed)[2],
    time = 1:dim(pi_mean_overdispersed)[3]
)
dimnames(pi_mean_latent) <- list(
    location = 1:dim(pi_mean_latent)[1],
    species = 1:dim(pi_mean_latent)[2],
    time = 1:dim(pi_mean_latent)[3]
)
# dimnames(pi_mra_mean) <- list(
#     location = 1:dim(pi_mra_mean)[1],
#     species = 1:dim(pi_mra_mean)[2],
#     time = 1:dim(pi_mra_mean)[3]
# )

dat_pi <- as.data.frame.table(pi_mean, responseName = "pi") %>%
    mutate(model = "Matern", location = as.numeric(location))
dat_pi_overdispersed <- as.data.frame.table(pi_mean_overdispersed, responseName = "pi") %>%
    mutate(model = "Overdispersed", location = as.numeric(location))
dat_pi_latent <- as.data.frame.table(pi_mean_latent, responseName = "pi") %>%
    mutate(model = "Latent", location = as.numeric(location))
# dat_pi_mra <- as.data.frame.table(pi_mra_mean, responseName = "pi") %>%
#     mutate(model = "MRA", location = as.numeric(location))

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

# check each of the models to the observed proportions
p1 <- rbind(dat_pi, dat_obs) %>%
    pivot_wider(names_from = model, values_from = pi) %>%
    ggplot(aes(x = Matern, y = observed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red")  + 
    facet_grid(species ~ time)
p2 <- rbind(dat_pi_overdispersed, dat_obs) %>%
    pivot_wider(names_from = model, values_from = pi) %>%
    ggplot(aes(x = Overdispersed, y = observed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red")  + 
    facet_grid(species ~ time)
p3 <- rbind(dat_pi_latent, dat_obs) %>%
    pivot_wider(names_from = model, values_from = pi) %>%
    ggplot(aes(x = Latent, y = observed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red")  + 
    facet_grid(species ~ time)

p1 / p2 / p3


rbind(dat_pi, dat_pi_overdispersed) %>%
    pivot_wider(names_from = model, values_from = pi) %>%
    ggplot(aes(x = Matern, y = Overdispersed)) +
    geom_point() +
    facet_grid(species ~ time)
rbind(dat_pi_latent, dat_pi_latent, dat_obs) %>%
    pivot_wider(names_from = model, values_from = pi) %>%
    ggplot(aes(x = Matern, y = Latent)) +
    geom_hex() +
    facet_grid(species ~ time)

# rbind(dat_pi, dat_pi_mra, dat_obs) %>% 
#     rbind(dat_pi, dat_pi_mra, dat_obs) %>% 
#     pivot_wider(names_from = model, values_from = pi) %>%
#     ggplot(aes(x = Matern, y = MRA)) +
#     geom_hex() +
#     facet_grid(species ~ time)


rbind(dat_pi, dat_obs) %>% 
    pivot_wider(names_from = model, values_from = pi) %>%
    ggplot(aes(x = Matern, y = observed)) +
    geom_point() +
    facet_grid(species ~ time)

# rbind(dat_pi_mra, dat_obs) %>% 
#     pivot_wider(names_from = model, values_from = pi) %>%
#     ggplot(aes(x = MRA, y = observed)) +
#     geom_hex() +
#     facet_grid(species ~ time)

plot(pi_mean, pi_mra_mean)
plot(pi_mean[, 1, 2], pi_mra_mean[, 1, 2])








# Plot the full fitted vs. observed proportions
pis_matern <- out$pi
pis_overdispersed <- out_overdispersed$pi
pis_latent <- out_latent$pi


dimnames(pis_matern) <- list(
    iteration = 1:dim(pis_matern)[1],
    location = 1:dim(pis_matern)[2],
    species = 1:dim(pis_matern)[3],
    time = 1:dim(pis_matern)[4])

dimnames(pis_overdispersed) <- list(
    iteration = 1:dim(pis_overdispersed)[1],
    location = 1:dim(pis_overdispersed)[2],
    species = 1:dim(pis_overdispersed)[3],
    time = 1:dim(pis_overdispersed)[4])

dimnames(pis_latent) <- list(
    iteration = 1:dim(pis_latent)[1],
    location = 1:dim(pis_latent)[2],
    species = 1:dim(pis_latent)[3],
    time = 1:dim(pis_matern)[4])

dat_matern <- as.data.frame.table(pis_matern, responseName = "pi") %>%
    mutate(model = "matern", 
           iteration = as.integer(iteration),
           location = as.integer(location),  
           species = as.integer(species), 
           time = as.integer(time))
dat_overdispersed <- as.data.frame.table(pis_overdispersed, responseName = "pi") %>%
    mutate(model = "overdispersed", 
           iteration = as.integer(iteration),
           location = as.integer(location),  
           species = as.integer(species), 
           time = as.integer(time))
dat_latent <- as.data.frame.table(pis_latent, responseName = "pi") %>%
    mutate(model = "latent", 
           iteration = as.integer(iteration),
           location = as.integer(location),  
           species = as.integer(species), 
           time = as.integer(time))


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
    mutate(model = "observed", 
           location = as.numeric(location), 
           species = as.integer(species), 
           time = as.integer(time))

# check each of the models to the observed proportions
calibration_matern <- dat_matern %>%
    group_by(location, species, time) %>%
    summarize(pi_mean = mean(pi),
              pi_lower  = quantile(pi, prob = 0.025),
              pi_upper  = quantile(pi, prob = 0.975)) %>%
    ungroup() %>%
    left_join(dat_obs) %>%
    ggplot(aes(x = pi, y = pi_mean)) +
    geom_point(alpha = 0.5) +
    facet_grid(species ~ time) +
    geom_errorbar(aes(ymin = pi_lower, ymax = pi_upper)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Matern model") +
    xlab("Observed proportions") +
    ylab("Estimated proportions") +
    theme_bw()

calibration_overdispersed <- dat_overdispersed %>%
    group_by(location, species, time) %>%
    summarize(pi_mean = mean(pi),
              pi_lower  = quantile(pi, prob = 0.025),
              pi_upper  = quantile(pi, prob = 0.975)) %>%
    ungroup() %>%
    left_join(dat_obs) %>%
    ggplot(aes(x = pi, y = pi_mean)) +
    geom_point(alpha = 0.5) +
    facet_grid(species ~ time) +
    geom_errorbar(aes(ymin = pi_lower, ymax = pi_upper)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Overdispersed model") +
    xlab("Observed proportions") +
    ylab("Estimated proportions") +
    theme_bw()

calibration_latent <- dat_latent %>%
    group_by(location, species, time) %>%
    summarize(pi_mean = mean(pi),
              pi_lower  = quantile(pi, prob = 0.025),
              pi_upper  = quantile(pi, prob = 0.975)) %>%
    ungroup() %>%
    left_join(dat_obs) %>%
    ggplot(aes(x = pi, y = pi_mean)) +
    geom_point(alpha = 0.5) +
    facet_grid(species ~ time) +
    geom_errorbar(aes(ymin = pi_lower, ymax = pi_upper)) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    ggtitle("Latent model") +
    xlab("Observed proportions") +
    ylab("Estimated proportions") +
    theme_bw()
    

ggsave(calibration_matern, 
       file = here::here("figures", "matern-calibration.png"),
       width = 16, height = 9, device = "png", units = "in")
ggsave(calibration_overdispersed, 
       file = here::here("figures", "overdispersed-calibration.png"),
       width = 16, height = 9, device = "png", units = "in")
ggsave(calibration_latent, 
       file = here::here("figures", "latent-calibration.png"),
       width = 16, height = 9, device = "png", units = "in")

