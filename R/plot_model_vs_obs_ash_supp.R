library(tidyverse)
library(pgR)
library(latex2exp)
library(Cairo)
source("R/plot-trace.R")


version <- "4.1"

dat <- readRDS(here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))
out <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_overdispersed.RDS')))

out=out


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



pi_mean <- apply(out$pi, c(2, 3, 4), mean)

dimnames(pi_mean) <- list(
    location = 1:dim(pi_mean)[1],
    species = 1:dim(pi_mean)[2],
    time = 1:dim(pi_mean)[3]
)

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

# specify time bins; use median of time bins for subsequent work
timev <- c(-70, seq(705, by = 990, length.out = 22)) * -1
time <- data.frame(YBP = timev[-1])
time$time <- rownames(time)


# check each of the models to the observed proportions
p1 <- rbind(dat_pi, dat_obs) %>%
    pivot_wider(names_from = model, values_from = pi) %>%
    filter(species==8, time %in% c(1,9,14,22)) %>%
    left_join(time) %>% 
    ggplot(aes(x = Matern, y = observed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red")  + 
    facet_grid(~YBP, labeller=label_both) + xlab("Model prediction") + ylab("Observed pollen abundance") +
    ggtitle("Model vs observed comparison for Fraxinus")

CairoPDF(paste0("figures/model_vs_obs_prop_FRAX",version,".pdf"))
print(p1)
dev.off()


CairoPNG(paste0("figures/model_vs_obs_prop_FRAX",version,".png"))
print(p1)
dev.off()



