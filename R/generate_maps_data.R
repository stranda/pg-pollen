# data maps

if(!dir.exists(here::here("figures"))) {
    dir.create(here::here("figures"))
}
if(!dir.exists(here::here("figures", "data"))) {
    dir.create(here::here("figures", "data"))
}

library(tidyverse)
library(patchwork)
require(rasterVis)
require(fields)
require(rgdal)
require(raster)
require(enmSdm)
require(rgeos)


# load the data

version <- '3.1'

# load the observation data
dat <- readRDS(here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))

# load the species names
species_names <- readRDS(here::here('data', 'species-names.RDS'))

# load in the predition grid (to make all plots on the same scale)
locs_grid <- readRDS(here::here('data', paste0('grid_', version, '.RDS')))


ys <- dat$y
for (tt in 1:dim(ys)[3]) {
    na_idx <- which(!is.na(ys[, , tt]), arr.ind = TRUE)
    ys[, , tt][na_idx] <- counts_to_proportions(matrix(ys[, , tt][na_idx], ncol = ncol(ys[, , tt])))
}
dimnames(ys) <- list(
    location = 1:dim(ys)[1],
    species  = species_names,
    time     = 1:dim(ys)[3])

locss <- dat$locs * 1000
dimnames(locss) <- list(    
    location = 1:dim(ys)[1],
    variable = c("lon", "lat"))
                            
dat_pollen <- as.data.frame.table(ys, responseName = "y") %>%
    left_join(as.data.frame(dat$locs * 1000) %>%
                  transmute(lon = x, lat = y) %>%
                  mutate(location = factor(1:n())))


#### READ MAP DATA ####
# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
# stack <- stack(here::here('data', 'map-data', 'study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'))
# names(stack) <- 1:701
# sub <- seq(1, 701, by = 33)
# stack_sub <- subset(stack, subset = paste0('X', sub))  # only want mask every 990 years
# stack_sub[stack_sub > 0.6] <- NA
# proj <- proj4string(stack)
proj <- "+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

na_shp <- readOGR(here::here("data", "map-data", "NA_States_Provinces_Albers.shp"), layer = "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR(here::here("data", "map-data", "Great_Lakes.shp"), "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj)

# getting bounding box to specify spatial domain when plotting
bbox_tran <- function(x, coord_formula = '~ x + y', from, to) {
    sp::coordinates(x) <- formula(coord_formula)
    sp::proj4string(x) <- sp::CRS(from)
    bbox <- as.vector(sp::bbox(sp::spTransform(x, CRSobj = sp::CRS(to))))
    return(bbox)
}

grid_box <- bbox_tran(locs_grid, '~ x + y',
                      proj,
                      proj)
xlim = c(grid_box[1], grid_box[3])
ylim = c(grid_box[2], grid_box[4])



dat_cont_shp <- fortify(cont_shp)
dat_lake_shp <- fortify(lake_shp)
# glimpse(dat_ar_shp)


# plot the precincts
map <- ggplot() +
    geom_polygon(data = dat_cont_shp,
                 aes(x = long, y = lat, group = group), fill = NA,
                 color = 'black', size = .2) +
    geom_polygon(data = dat_lake_shp,
                 aes(x = long, y = lat, group = group), fill = "blue",
                 color = 'black', size = .2)
map


# change the time numbering
time_bins <- seq(-285, 21000, by=990)[-1]

time_vec <- paste(time_bins[1:21], "ybp")
names(time_vec) <- paste(1:21)
time_vec


# convert to propotions
dat_pollen %>%
    group_by(location, time) %>%
    mutate(p = y / sum(y)) %>%
    ungroup


# generate the plots
base_size <- 22


for (species_to_plot in unique(dat_pollen$species)) {
    p_pollen <- dat_pollen %>% 
        group_by(location, time) %>%
        mutate(p = y / sum(y)) %>%
        ungroup %>%
        filter(species == species_to_plot) %>%
        drop_na() %>%
        ggplot(aes(x = lon, y = lat, fill = p, color = p)) +
        geom_point() +
        geom_polygon(data = dat_cont_shp,
                     aes(x = long, y = lat, group = group), fill = NA,
                     color = 'black', size = .2) +
        geom_polygon(data = dat_lake_shp,
                     aes(x = long, y = lat, group = group), fill = "black",
                     color = 'black', size = .2) +
        facet_wrap(~ time, nrow = 3, labeller = as_labeller(time_vec)) +
        scale_fill_viridis_c(end = 0.8, name = "p") + 
        scale_color_viridis_c(end = 0.8, name = "p") + 
        ggtitle(paste0("Observed proportion ", species_to_plot)) +
        theme_bw(base_size = base_size) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        xlab("") +
        ylab("") +         
        coord_cartesian(xlim = xlim, ylim = ylim)
    
    p_pollen
    
    
    ggsave(p_pollen, 
           file = paste0("~/pg-pollen/figures/data/observations-", species_to_plot, ".png") ,
               # here::here("figures", "data", paste0("observations-", species_to_plot, ".png")), 
           height = 9,
           width = 16)
}
