# setwd('C:/Users/abrow/Documents/pg-pollen')
require(raster)
require(dplyr)
require(ggplot2)
require(tidyr)
require(enmSdm)

# read in predictions, convert to list of lists (times within iterations)
preds <- readRDS('output/frax_preds_n50_3.2_latent_overdispersed.RDS')
locs <- readRDS('data/grid_3.1.RDS')
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
stack <- readRDS('data/map-data/study_region_mask_glacierNA.RDS')
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
stack[stack == 1] <- NA
names(stack) <- 700:0
# saveRDS(stack, 'data/map-data/study_region_mask_glacierNA.RDS')
# stack <- readRDS('data/map-data/study_region_mask_glacierNA.RDS')

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

# create vector of pollen times
time <- c(-70, seq(705, by = 990, length.out = 22))
time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
time$time_mid <- (time$timeFrom + time$timeTo) / 2  # to get midpoint of time bins
times <- time$time_mid

# create vector of ice mask times
interp_times <- as.numeric(mask_names)

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
saveRDS(interp_list, 'output/interp_preds_n50_latent_overdispersed_v3.2.RDS')


# mask out ice/water
masked_interp_list <- interp_list
for(i in 1:n_iter){
  for(j in 1:n_times){
  masked_interp_list[[i]][[j]] <- raster::resample(masked_interp_list[[i]][[j]], y = mask_list[[j]])
  masked_interp_list[[i]][[j]] <- mask(masked_interp_list[[i]][[j]], mask = mask_list[[j]])
  }
}
saveRDS(masked_interp_list, 'output/interp_masked_preds_n50_latent_overdispersed_v3.2.RDS')



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
saveRDS(masked_interp_revalue, 'output/preds_for_ABC_n50_latent_overdispersed_v3.2.RDS')

# first 11 time slices contain no maps - fill those in with most recent map
interp_revalue_fill <- masked_interp_revalue
for(i in 1:n_iter){
  interp_revalue_fill[[i]] <- unstack(interp_revalue_fill[[i]])
  interp_revalue_fill[[i]][691:701] <- interp_revalue_fill[[i]][[690]]
}
saveRDS(interp_revalue_fill, 'output/preds_for_ABC_n50_latent_overdispersed_v3.2_final.RDS')


