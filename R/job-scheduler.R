library(here)
library(rstudioapi)

# Model speed comparison
jobRunScript(here::here("R", "speed-testing.R"))

# Build prediction grid
# this file needs to be resolved

# Matern model
jobRunScript(here::here("R", "estimate_pars_ST.R"))
jobRunScript(here::here("R", "predict_ST.R"))
jobRunScript(here::here("R", "generate-maps_ST.R"))

# Overdispersed Matern model
jobRunScript(here::here("R", "estimate_pars_ST_overdispersed.R"))
jobRunScript(here::here("R", "predict_ST_overdispersed.R"))
jobRunScript(here::here("R", "generate-maps_ST_overdispersed.R"))

# Latent Overdispersed Matern model
jobRunScript(here::here("R", "estimate_pars_ST_latent_overdispersed.R"))
jobRunScript(here::here("R", "predict_ST_latent_overdispersed.R"))
jobRunScript(here::here("R", "generate-maps_ST_latent_overdispersed.R"))

# MRA model -- currently fits slow and has prediction (indexing?) errors
# jobRunScript(here::here("R", "estimate_pars_ST_MRA.R"))
# jobRunScript(here::here("R", "predict_ST_MRA.R"))
# jobRunScript(here::here("R", "generate-maps_ST_MRA.R"))

