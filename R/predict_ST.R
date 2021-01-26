library(pgR)
library(ggplot2)
library(fields)
library(rgdal)

version <- 'Dec15'

out <- readRDS(paste0('output/polya-gamma-posts_', version, '.RDS'))
dat <- readRDS(paste0('output/polya-gamma-dat_', version,'.RDS'))

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

locs_grid = readRDS(paste0('data/grid_', version, '.RDS'))
X_pred <- matrix(rep(1, nrow(locs_grid)), nrow(locs_grid), 1)
locs = locs_pollen/rescale
locs_pred = locs_grid/rescale

#### MAKE PREDICTIONS ####
class(out) <- "pg_stlm"
preds = predict_pg_stlm(
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
saveRDS(preds, paste0('output/preds_', version, '.RDS'))

