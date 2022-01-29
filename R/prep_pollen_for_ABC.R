# setwd('C:/Users/abrow/Documents/pg-pollen')
require(raster)
require(dplyr)
require(ggplot2)
require(tidyr)
require(enmSdm)
require(ggforce)
require(gridExtra)
require(holoSimCell)
require(rgdal)

# read in predictions, convert to list of lists (times within iterations)
preds <- readRDS('output/preds_frax_n200_v4.1.RDS')
locs <- readRDS('data/grid_4.1.RDS')
n_iter <- dim(preds)[1]
n_times <- dim(preds)[3]

preds_list <- list()
for(i in 1:n_iter){
  preds_list[[i]] <- cbind(locs, data.frame(preds[i,,]))
  preds_list[[i]] <- pivot_longer(preds_list[[i]], cols = 3:ncol(preds_list[[i]]),
                                  names_to = 'time', values_to = 'pred')
  preds_list[[i]]$time <- as.integer(substr(preds_list[[i]]$time, start = 2, 
                                            stop = nchar(preds_list[[i]]$time)))
  preds_list[[i]] <- split(preds_list[[i]], f = preds_list[[i]]$time)
  preds_list[[i]] <- lapply(preds_list[[i]], function(x) { x['time'] <- NULL; x })
}


# convert dataframes to rasters
proj <- '+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs'
raster_list <- preds_list
for(i in 1:n_iter){
  for(j in 1:n_times){
    raster_list[[i]][[j]] <- rasterFromXYZ(raster_list[[i]][[j]], crs = CRS(proj))
  }
}


# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
stack[stack == 1] <- NA
saveRDS(stack, 'data/map-data/study_region_mask_glacierNA.RDS')
stack <- readRDS('data/map-data/study_region_mask_glacierNA.RDS')
names(stack) <- 700:0

# convert rasterstack masks into list of rasters, reverse order of items to correspond with pollen time
mask_list <- unstack(stack)
n_times <- length(mask_list)
mask_names <- NA
for(i in 1:n_times){
  mask_names[i] <- names(mask_list[[i]])
  mask_names[i] <- as.integer(substr(mask_names[i], start = 2, stop = nchar(mask_names[i]))) * 30
}
names(mask_list) <- mask_names
mask_list <- rev(mask_list)
mask_names <- rev(mask_names)
interp_times <- as.numeric(mask_names)

# create vector of pollen times
time <- c(-70, seq(705, by = 990, length.out = 22))
time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
time$time_mid <- (time$timeFrom + time$timeTo) / 2  # to get midpoint of time bins
times <- time$time_mid

# interpolate pollen rasters to temporal resolution of 30yrs
# convert list (iterations) of lists (times) to list of rasterstacks (stacked by time)
rasterstack_list <- list()
for(i in 1:n_iter){
  rasterstack_list[[i]] <- stack(raster_list[[i]])
}

# interpolate for each iteration
interp_list <- list()
for(i in 1:n_iter){
  interp_list[[i]] <- interpolateRasters(rasterstack_list[[i]], interpFrom = times, interpTo = interp_times)
  interp_list[[i]] <- unstack(interp_list[[i]])
}
saveRDS(interp_list, 'output/interp_preds_n200_v4.1.RDS')


# mask out ice/water
masked_interp_list <- interp_list
for(i in 1:n_iter){
  for(j in 1:n_times){
  masked_interp_list[[i]][[j]] <- raster::resample(masked_interp_list[[i]][[j]], y = mask_list[[j]])
  masked_interp_list[[i]][[j]] <- mask(masked_interp_list[[i]][[j]], mask = mask_list[[j]])
  }
}
saveRDS(masked_interp_list, 'output/interp_masked_preds_n200_v4.1.RDS')


# weight relative abundance by proportion of cell covered by land (for cells with partial glacial coverage)
mask_list_revalue <- mask_list
for(i in 1:n_times){
  values(mask_list_revalue[[i]]) <- 1 - values(mask_list_revalue[[i]])
}

masked_interp_revalue <- masked_interp_list
for(i in 1:n_iter){
  for(j in 1:n_times){
    masked_interp_revalue[[i]][[j]] <- overlay(masked_interp_list[[i]][[j]], mask_list_revalue[[j]],
                                             fun = function(x, y){(x * y)})
  }
}

# convert masked, interpolated, revalued rasters to list of rasterstacks (each stack = 1 iteration)
# first, reverse order of rasters within each iteration
for(i in 1:n_iter){
  masked_interp_revalue[[i]] <- rev(masked_interp_revalue[[i]])
  masked_interp_revalue[[i]] <- stack(masked_interp_revalue[[i]])
}
saveRDS(masked_interp_revalue, 'output/preds_for_ABC_n200_v4.1.RDS')

# first 11 time slices contain no maps - fill those in with most recent map
interp_revalue_fill <- masked_interp_revalue
for(i in 1:n_iter){
  interp_revalue_fill[[i]] <- unstack(interp_revalue_fill[[i]])
  interp_revalue_fill[[i]][691:701] <- interp_revalue_fill[[i]][[690]]
  interp_revalue_fill[[i]] <- stack(interp_revalue_fill[[i]])
}
saveRDS(interp_revalue_fill, 'output/preds_for_ABC_n200_v4.1_final.RDS')



################################################
## CALCULATE AND PLOT REFUGE LOCATIONS
################################################

# read in 'simulated' dataframe
sim <- raster('data/map-data/study_region_resampled_to_genetic_demographic_simulation_resolution.tif')
proj <- proj4string(sim)

# read in fraxinus predictions for each of the 50 iterations
# subset to only LGM predictions
grid <- readRDS('data/grid_4.1.RDS')
frax <- readRDS('output/preds_frax_n200_v4.1.RDS')
frax_lgm <- frax[,,22]

# convert data array to list of dataframes (1 dataframe for each iteration)
n_iter <- dim(frax_lgm)[1]
lgm_list <- list()
for(i in 1:n_iter){
  lgm_list[[i]] <- data.frame(cbind(grid, frax_lgm[i,]))
  names(lgm_list[[i]]) <- c('x','y','pred')
}

# convert dfs to rasters
lgm_raster <- list()
for(i in 1:n_iter){
  lgm_raster[[i]] <- rasterFromXYZ(lgm_list[[i]], crs = CRS(proj))
}

# read in masking rasters, extract mask at LGM
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
LGM <- stack$study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.1  # LGM
LGM[LGM == 1] <- NA  # convert glacial values to NAs

# mask out glacial coverage and lakes/ocean
lgm_resamp <- list()
for(i in 1:n_iter){
  lgm_resamp[[i]] <- raster::resample(x = lgm_raster[[i]], y = LGM)
  lgm_resamp[[i]] <- mask(x = lgm_resamp[[i]], mask = LGM)
}

# weight relative abundance by proportion of cell covered by land (for cells with partial glacial coverage)
# partial glacial coverage is indicated as proportion of cell covered by glacier
# first, convert this value to proportion of cell covered by land
lgm_revalue <- LGM
values(lgm_revalue) <- 1 - values(lgm_revalue)

# multiply relative abundance predictions by proportion of cell covered by land
lgm_preds_revalue <- lgm_resamp
for(i in 1:n_iter){
  lgm_preds_revalue[[i]] <- overlay(lgm_preds_revalue[[i]], lgm_revalue,
                                    fun = function(x, y){(x * y)})
}

# identify refuge locations
# can specify the threshold you want, in units of relative abundance
refuge_list <- list()
for(i in 1:n_iter){
  refuge_list[[i]] <- assignRefugiaFromAbundanceRaster(abund = lgm_preds_revalue[[i]],
                                                       sim = sim, threshold = 0.015)
}

# check to make sure all iterations have been assigned a refuge
test <- rep(NA, n_iter)
for(i in 1:n_iter){
  test[i] <- refuge_list[[i]]$meanRefugeAbund
}
anyNA(test)  # should be FALSE

# save as list of rasters
# THIS WILL BE IN REFUGIA INPUT FOR POLLEN-INTEGRATED MODEL
refuge_for_allan <- list()
for(i in 1:n_iter){
  refuge_for_allan[[i]] <- refuge_list[[i]][["simulationScale"]]@layers[[2]]
}
saveRDS(refuge_for_allan, 'output/refuge_rasters_n200_V4.1.RDS')


# PLOT 50 REFUGIA ON SINGLE PDF
# convert refugia rasters to dataframes for plotting in ggplot2
refuge_dfs <- list()
for(i in 1:n_iter){
  refuge_dfs[[i]] <- raster::as.data.frame(refuge_list[[i]]$simulationScale@layers[[2]], xy = TRUE)
  refuge_dfs[[i]]$iter <- i
}
refuge_df <- bind_rows(refuge_dfs)

# north america shapefiles
na_shp <- readOGR("data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("data/map-data/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj)

ylim <- c(min(refuge_df$y), max(refuge_df$y))
xlim <- c(min(refuge_df$x), max(refuge_df$x))

# create discrete color scale for plotting
# make sure 0s aren't color-coded, since they aren't refugia
breaks_vec <- c(0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.45)
breaks <- cut(refuge_df$refugiaAbund, breaks = breaks_vec, include.lowest = TRUE, labels = FALSE)
refuge_df$breaks <- breaks
breaks_vec <- c('0.025 - 0.05', '0.05 - 0.075', '0.075 - 0.1', '0.1 - 0.15', '0.15 - 0.2', 
                '0.2 - 0.3', '0.3 - 0.45')
my_colors <- RColorBrewer::brewer.pal(9, 'YlOrRd')[3:9]

# plot across multiple PDF pages
plot_list <- list()
for(i in 1:5){
  plot_list[[i]] <- ggplot(refuge_df) + 
    geom_tile(aes(x = x, y = y, fill = as.factor(breaks))) + 
    scale_fill_manual(values = my_colors, labels = breaks_vec, na.value = 'white') +
    scale_y_continuous(limits = ylim) + 
    scale_x_continuous(limits = xlim) +
    geom_path(data = cont_shp, aes(x = long, y = lat, group = group), alpha = 0.5) +
    geom_path(data = lake_shp, aes(x = long, y = lat, group = group), alpha = 0.5) +
    labs(fill = 'Relative\nabund') +
    facet_wrap_paginate(~ iter, nrow = 5, ncol = 2, page = i) +
    theme_void() +
    theme(strip.text = element_blank(),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10)) +
    coord_equal()
}
plots <- marrangeGrob(plot_list, nrow = 1, ncol = 1)  # takes several minutes to run
ggsave('figures/frax_refugia_n50_v4.0.pdf', plots, 
       height = 15, width = 8.5, units = 'in')
