library(tidyverse)
library(pgR)
library(latex2exp)
library(Cairo)
source("R/plot-trace.R")
source("R/plot_trace_latent.R")

version <- "4.0"

dat <- readRDS(here::here('output.first', paste0('polya-gamma-dat_', version,'.RDS')))
out_latent <- readRDS(here::here('output.first', paste0('polya-gamma-posts_', version, '_latent_overdispersed.RDS')))
plot_trace_latent(out_latent, base_size = 7, file = "figures/trace-plots_latent.pdf")

out=out_latent

CairoPDF("figures/eta_distribution.pdf")
hist(out_latent$eta)
dev.off()



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


CairoPDF("figures/params_vs_kappa.pdf")
plot(apply(out$eta, c(2, 3, 4), mean), kappa)

plot(apply(out$eta, c(2, 3, 4), mean), apply(out$omega, c(2, 3, 4), mean))

plot(apply(out$omega, c(2, 3, 4), mean), kappa)

cor(apply(out$omega, c(2, 3, 4), mean), kappa)

plot(apply(out$eta, c(2, 3, 4), mean), kappa)

dev.off()



# check the predicted probabilities between the different models

# preds <- readRDS(here::here("output", paste0('polya-gamma-predictions_', version, '.RDS')))
# preds_mra <- readRDS(here::here("output", paste0('polya-gamma-predictions_', version, '_MRA.RDS')))

pi_mean <- apply(out$pi, c(2, 3, 4), mean)

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
# dimnames(pi_mra_mean) <- list(
#     location = 1:dim(pi_mra_mean)[1],
#     species = 1:dim(pi_mra_mean)[2],
#     time = 1:dim(pi_mra_mean)[3]
# )

dat_pi <- as.data.frame.table(pi_mean, responseName = "pi") %>%
    mutate(model = "Matern", location = as.numeric(location))

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

CairoPDF("figures/model_vs_obs_prop.pdf")
print(p1)
dev.off()


# Plot the full fitted vs. observed proportions
pis_matern <- out$pi


dimnames(pis_matern) <- list(
    iteration = 1:dim(pis_matern)[1],
    location = 1:dim(pis_matern)[2],
    species = 1:dim(pis_matern)[3],
    time = 1:dim(pis_matern)[4])


dat_matern <- as.data.frame.table(pis_matern, responseName = "pi") %>%
    mutate(model = "matern", 
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


ggsave(calibration_matern, 
       file = "figures/calibration.pdf",
       width = 16, height = 9, device = cairo_pdf, units = "in")

