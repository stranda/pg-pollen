# prediction maps

if(!dir.exists(here::here("figures"))) {
    dir.create(here::here("figures"))
}
if(!dir.exists(here::here("figures", "matern"))) {
    dir.create(here::here("figures", "matern"))
}

library(tidyverse)
library(patchwork)

version <- '3.1'

# load the species names
species_names <- readRDS(here::here('data', 'species-names.RDS'))

# load the preds data
preds <- readRDS(here::here("output", paste0('polya-gamma-predictions_', version, '.RDS')))

# load the preds grid
locs_grid <- readRDS(here::here('data', paste0('grid_', version, '.RDS')))

pi_mean <- apply(preds$pi, c(2, 3, 4), mean)
pi_sd <- apply(preds$pi, c(2, 3, 4), sd)

# modify this based on Alissa's input
dimnames(pi_mean) <- list(
    location = 1:dim(pi_mean)[1],
    species = species_names,
    time = 1:dim(pi_mean)[3]
)


dat_pi_mean <- as.data.frame.table(pi_mean, responseName = "pi_mean") %>%
    mutate(location = as.numeric(location), time = as.numeric(time)) %>%
    left_join(locs_grid %>%
    mutate(location = 1:dim(locs_grid)[1]))

dimnames(pi_sd) <- list(
    location = 1:dim(pi_sd)[1],
    species = species_names,
    time = 1:dim(pi_sd)[3]
)

dat_pi_sd <- as.data.frame.table(pi_sd, responseName = "pi_sd") %>%
    mutate(location = as.numeric(location), time = as.numeric(time)) %>%
    left_join(locs_grid %>%
                  mutate(location = 1:dim(locs_grid)[1]))

# generate the plots

for (species_to_plot in unique(dat_pi_mean$species)) {
    
    if (!file.exists(here::here("figures", "matern", paste0("predictions-", species_to_plot, ".png")))) {
        p_mean <- dat_pi_mean %>%
            # filter(species %in% species_to_plot) %>%
            filter(species == species_to_plot) %>%
            ggplot(aes(x = x, y = y, fill = pi_mean)) +
            geom_raster() +
            facet_wrap(~ time) +
            scale_fill_viridis_c()
        
        p_sd <- dat_pi_sd %>%
            # filter(species %in% species_to_plot) %>%
            filter(species == species_to_plot) %>%
            ggplot(aes(x = x, y = y, fill = pi_sd)) +
            geom_raster() +
            facet_wrap(~ time) +
            scale_fill_viridis_c()
        ggsave(p_mean + p_sd, 
               file = here::here("figures", "matern", paste0("predictions-", species_to_plot, ".png")), 
               height = 9,
               width = 16)
    }
}
