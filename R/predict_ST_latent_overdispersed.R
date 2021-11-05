library(pgR)
library(ggplot2)
library(fields)
library(rgdal)

version <- '4.0'

out <- readRDS(here::here('output.first', paste0('polya-gamma-posts_', version, '_latent_overdispersed.RDS')))
dat <- readRDS(here::here('output.first', paste0('polya-gamma-dat_', version,'.RDS')))

# note that locations were scaled to fit the model
# unscaling to think in meters, then will rescale again before prediction
rescale = dat$rescale
locs_pollen <- dat$locs*rescale 
names(locs_pollen) <- c("x", "y")

taxa.keep = dat$taxa.keep
y = dat$y
X = dat$X
N_cores = nrow(locs_pollen)

#### DISTANCE MATRICES ####
D_pollen <- fields::rdist(locs_pollen/rescale)# N_cores x N_cores
any(D_pollen == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
D_pollen <- ifelse(D_pollen == 0, 0.007, D_pollen)  # remove them
diag(D_pollen) <- 0

locs_grid = readRDS(here::here('data', paste0('grid_', version, '.RDS')))
X_pred <- matrix(rep(1, nrow(locs_grid)), nrow(locs_grid), 1)
locs = locs_pollen/rescale
locs_pred = locs_grid/rescale

#### MAKE PREDICTIONS ####
# class(out) <- "pg_stlm"

if (!file.exists(here::here("output.first", paste0('polya-gamma-predictions_', version, '_latent_overdispersed.RDS')))) {
  preds = predict_pg_stlm_latent_overdispersed(
    out,
    X,
    X_pred,
    locs = locs,
    locs_pred = locs_pred,
    corr_fun = "matern",
    shared_covariance_params = FALSE,
    progress = TRUE, 
    verbose = TRUE
  )
  saveRDS(preds, here::here("output", paste0('polya-gamma-predictions_', version, '_latent_overdispersed.RDS')),
          compress = FALSE)
  #pushoverr::pushover(message = "Finished predicting latent overdispersed Matern model")
}

