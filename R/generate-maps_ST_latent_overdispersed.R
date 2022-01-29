# prediction maps

if(!dir.exists(here::here("figures"))) {
    dir.create(here::here("figures"))
}
if(!dir.exists(here::here("figures", "matern_latent_overdispersed"))) {
    dir.create(here::here("figures", "matern_latent_overdispersed"))
}

library(tidyverse)
library(patchwork)
require(rasterVis)
require(fields)
require(rgdal)
require(raster)
require(enmSdm)
require(rgeos)


version <- '3.1'

# load the species names
species_names <- readRDS(here::here('data', 'species-names.RDS'))

# load the preds data
preds <- readRDS(here::here("output", paste0('polya-gamma-predictions_', version, '_latent_overdispersed.RDS')))

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
# dat_pi_mean <- dat_pi_mean %>%
#     mutate(time = recode(time, 
#                          "1" = time_bins[1],
#                          "2" = time_bins[2],
#                          "3" = time_bins[3],
#                          "4" = time_bins[4],
#                          "5" = time_bins[5],
#                          "6" = time_bins[6],
#                          "7" = time_bins[7],
#                          "8" = time_bins[8],
#                          "9" = time_bins[9],
#                          "10" = time_bins[10],
#                          "11" = time_bins[11],
#                          "12" = time_bins[12],
#                          "13" = time_bins[13],
#                          "14" = time_bins[14],
#                          "15" = time_bins[15],
#                          "16" = time_bins[16],
#                          "17" = time_bins[17],
#                          "18" = time_bins[18],
#                          "19" = time_bins[19],
#                          "20" = time_bins[20],
#                          "21" = time_bins[21]))
# dat_pi_sd <- dat_pi_sd %>%
#     mutate(time = recode(time, 
#                             "1" = time_bins[1],
#                             "2" = time_bins[2],
#                             "3" = time_bins[3],
#                             "4" = time_bins[4],
#                             "5" = time_bins[5],
#                             "6" = time_bins[6],
#                             "7" = time_bins[7],
#                             "8" = time_bins[8],
#                             "9" = time_bins[9],
#                             "10" = time_bins[10],
#                             "11" = time_bins[11],
#                             "12" = time_bins[12],
#                             "13" = time_bins[13],
#                             "14" = time_bins[14],
#                             "15" = time_bins[15],
#                             "16" = time_bins[16],
#                             "17" = time_bins[17],
#                             "18" = time_bins[18],
#                             "19" = time_bins[19],
#                             "20" = time_bins[20],
#                             "21" = time_bins[21]))

# time_vec <- paste(1:21)
# names(time_vec) <- paste(time_bins[1:21], "ybp")
time_vec <- paste(time_bins[1:21], "ybp")
names(time_vec) <- paste(1:21)
time_vec

# as.character(1:21))
# time_labeller <- list(
#     paste(time_bins[1], "ybp") = "1",
#     paste(time_bins[2], "ybp") = "2", 
#     paste(time_bins[3], "ybp") = "3",
#     paste(time_bins[4], "ybp") = "4",
#     paste(time_bins[5], "ybp") = "5",
#     paste(time_bins[6], "ybp") = "6",  
#     paste(time_bins[7], "ybp") = "7",  
#     paste(time_bins[8], "ybp") = "8",  
#     paste(time_bins[9], "ybp") = "9",  
#     paste(time_bins[10], "ybp") = "10", 
#     paste(time_bins[11], "ybp") = "11",  
#     paste(time_bins[12], "ybp") = "12",  
#     paste(time_bins[13], "ybp") = "13",  
#     paste(time_bins[14], "ybp") = "14",  
#     paste(time_bins[15], "ybp") = "15",  
#     paste(time_bins[16], "ybp") = "16",  
#     paste(time_bins[17], "ybp") = "17",  
#     paste(time_bins[18], "ybp") = "18",  
#     paste(time_bins[19], "ybp") = "19",  
#     paste(time_bins[20], "ybp") = "20",  
#     paste(time_bins[21], "ybp") = "21")  

# generate the plots
base_size <- 14






# generate the plots

for (species_to_plot in unique(dat_pi_mean$species)) {
    
    if (!file.exists(here::here("figures", "matern_latent_overdispersed", paste0("predictions-", species_to_plot, "_latent_overdispersed.png")))) {
        p_mean <- dat_pi_mean %>%
            # filter(species %in% species_to_plot) %>%
            filter(species == species_to_plot) %>%
            ggplot(aes(x = x, y = y, fill = pi_mean)) +
            geom_raster() +
            geom_polygon(data = dat_cont_shp,
                         aes(x = long, y = lat, group = group), fill = NA,
                         color = 'black', size = .2) +
            geom_polygon(data = dat_lake_shp,
                         aes(x = long, y = lat, group = group), fill = "black",
                         color = 'black', size = .2) +
            facet_wrap(~ time, nrow = 3, labeller = as_labeller(time_vec)) +
            scale_fill_viridis_c() +
            ggtitle(paste0("Posterior proportion mean ", species_to_plot)) +
            theme_bw(base_size = base_size) +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y=element_blank()) +
            xlab("") +
            ylab("") +         
            coord_cartesian(xlim = xlim, ylim = ylim)
        
        p_sd <- dat_pi_sd %>%
            # filter(species %in% species_to_plot) %>%
            filter(species == species_to_plot) %>%
            ggplot(aes(x = x, y = y, fill = pi_sd)) +
            geom_raster() +
            geom_polygon(data = dat_cont_shp,
                         aes(x = long, y = lat, group = group), fill = NA,
                         color = 'black', size = .2) +
            geom_polygon(data = dat_lake_shp,
                         aes(x = long, y = lat, group = group), fill = "black",
                         color = 'black', size = .2) +
            facet_wrap(~ time, nrow = 3, labeller = as_labeller(time_vec)) +
            scale_fill_viridis_c() +
            ggtitle(paste0("Posterior proportion sd ", species_to_plot)) +
            theme_bw(base_size = base_size) +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y=element_blank()) +
            xlab("") +
            ylab("") +
            coord_cartesian(xlim = xlim, ylim = ylim)
        scale_fill_viridis_c()
        
        
        ggsave(p_mean / p_sd, 
               file = here::here("figures", "matern_latent_overdispersed", paste0("predictions-", species_to_plot, "_latent_overdispersed.png")), 
               height = 9,
               width = 16)
    }
}
