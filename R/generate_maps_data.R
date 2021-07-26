# data maps

if(!dir.exists(here::here("figures"))) {
    dir.create(here::here("figures"))
}
if(!dir.exists(here::here("figures", "data"))) {
    dir.create(here::here("figures", "data"))
}

library(tidyverse)
library(patchwork)

# load the data

version <- '3.1'

# load the observation data
dat <- readRDS(here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))

# load the species names
species_names <- readRDS(here::here('data', 'species-names.RDS'))

ys <- dat$y
for (tt in 1:dim(ys)[3]) {
    na_idx <- which(!is.na(ys[, , tt]), arr.ind = TRUE)
    ys[, , tt][na_idx] <- counts_to_proportions(matrix(ys[, , tt][na_idx], ncol = ncol(ys[, , tt])))
}
dimnames(ys) <- list(
    location = 1:dim(ys)[1],
    species  = species_names,
    time     = 1:dim(ys)[3])

locss <- dat$locs
dimnames(locss) <- list(    
    location = 1:dim(ys)[1],
    variable = c("lon", "lat"))
                            
dat_pollen <- as.data.frame.table(ys, responseName = "y") %>%
    left_join(as.data.frame(dat$locs) %>%
                  transmute(lon = x, lat = y) %>%
                  mutate(location = factor(1:n())))

for (species_to_plot in unique(dat_pollen$species)) {
    p_pollen <- dat_pollen %>% 
        filter(species == species_to_plot) %>%
        drop_na() %>%
        ggplot(aes(x = lon, y = lat, fill = y, color = y)) +
        geom_point() +
        scale_fill_viridis_c(end = 0.8) + 
        scale_color_viridis_c(end = 0.8) + 
        facet_wrap( ~ time)
    ggsave(p_pollen, 
           file = here::here("figures", "data", paste0("observations-", species_to_plot, ".png")), 
           height = 9,
           width = 16)
}
