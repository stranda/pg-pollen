
require(tidyr)
require(ggplot2)
require(rasterVis)
require(fields)
require(rgdal)
require(raster)
require(enmSdm)
require(rgeos)
require(sp)
require(dplyr)

setwd('C:/Users/abrow/Documents/pg-pollen')
version <- 'Dec15'

# READ IN DATA
# read in raw data
dat <- readRDS(paste0('output/polya-gamma-dat_', version, '.RDS'))
taxa <- dat$taxa.keep
frax_num <- which(taxa == 'Fraxinus')

# read in prediction grid locations
locs_grid <- readRDS('data/grid_Oct12.RDS')

# read in predictions
preds <- readRDS('output/Dec2020-inter_frax_pis_small.RDS')

# specify projection
proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
+lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"


# Get raster masks and spatial domain you want to use
# (Adam Smith's .tif from: NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
names(stack) <- 1:701
sub <- seq(0, 701, by = 33)
stack_sub <- subset(stack, subset = paste0('X', sub))  # only want mask every 990 years
proj <- proj4string(stack)

# reproject coordinates to match with ENM and genetic rasters
sp::coordinates(locs_grid) <- ~x + y
sp::proj4string(locs_grid) <- CRS(proj_out)
locs_grid <- sp::spTransform(locs_grid, CRS(proj))
locs_grid <- data.frame(locs_grid)

# convert prediction array into a list of raster objects
ntimes <- dim(preds)[3]

preds_list <- list()
for(i in 1:dim(preds)[1]){
  preds_list[[i]] <- preds[i,,]
}

preds_list <- lapply(preds_list, as.data.frame)
preds_list <- lapply(preds_list, function(x) cbind(locs_grid, x))

pivot_fxn <- function(x) pivot_longer(x, cols = 3:ncol(x), 
                                      names_to = 'time',
                                      values_to = 'pred')
preds_list <- lapply(preds_list, pivot_fxn)

for(i in 1:length(preds_list)){
  preds_list[[i]]$time <- as.integer(substr(preds_list[[i]]$time, start = 2, 
                                  stop = nchar(preds_list[[i]]$time)))
}

preds_list2 <- preds_list
for(i in 1:length(preds_list)){
  preds_list2[[i]] <- split(preds_list2[[i]], f = preds_list2[[i]]$time)
}

# remove time column from all dataframes - I think this will make it easier
# to deal with converting dfs to rasters - should only have x, y, pred columns
for(i in 1:length(preds_list2)){
  for(j in 1:length(ntimes)){
    preds_list2[[i]] <- lapply(preds_list2[[i]], 
                                    function(x) x[!(names(x) %in% c('time'))])
  }
}

# convert dfs to rasters
rast_list <- preds_list2
for(i in 1:length(rast_list)){
  for(j in 1:ntimes){
    rast_list[[i]][[j]] <- rasterFromXYZ(rast_list[[i]][[j]])
    proj4string(rast_list[[i]][[j]]) <- proj
  }
}

# saveRDS(rast_list, 'output/frax_preds_rasters_by_iteration.RDS')
# rast_list <- readRDS('output/frax_preds_rasters_by_iteration.RDS')

\

# resample rasters so they're at the correct res and extent
# Code to use:
# resampledPollen <- resample(pollenRaster, template)
# resampledPollen <- calc(resampledPollen, fun=function(x) ifelse(x > 1, 1, x))
# resampledPollen <- calc(resampledPollen, fun=function(x) ifelse(x < 0, 0, x))
# template <- raster('data/map-data/study_region_resampled_to_genetic_demographic_simulation_resolution.tif')
# for(i in 1:length(rast_list)){
#   for(j in 1:ntimes){
#     rast_list[[i]][[j]] <- resample(rast_list[[i]][[j]], template)
#   }
# }



# FOR NOW, TO FIX EXTENT ISSUE, CROP BOTTOM OF POLLEN RASTERS TO SAME YMIN
# VALUE AS ADAM'S STACK
# BUT FIRST NEED TO RESAMPLE, OTHERWISE CROPPING WON'T FIX EXTENT DIFFERENCES

for(i in 1:length(rast_list)){
  for(j in 1:ntimes){
    rast_list[[i]][[j]] <- resample(rast_list[[i]][[j]], stack_sub)
  }
}

# e <- extent(rast_list[[1]][[1]]@extent@xmin, 
#             rast_list[[1]][[1]]@extent@xmax, 
#             stack_sub@extent@ymin, 
#             rast_list[[1]][[1]]@extent@ymax)

# rast_list_cr <- rast_list
# for(i in 1:100){
#   for(j in 1:21){
#     rast_list_cr[[i]][[j]] <- crop(rast_list[[i]][[j]], e)
#   }
# }

# stack_sub_cr <- crop(stack_sub, rast_list_cr[[1]][[1]])


# mask out ocean, lakes, glaciers
unstack <- unstack(stack_sub)
mask_list <- rast_list
for(i in 1:length(mask_list)){
  for(j in 1:ntimes){
    mask_list[[i]][[j]] <- mask(mask_list[[i]][[j]], 
                                mask = unstack[[j]])
  }
}
# STARTED AT 12:33 - 12:47

# ERROR: Error in compareRaster(x, mask) : different extent
# my rasters have different extent than Adam's rasterstack
# 
# might need to crop the larger one to the extent of the other, if the other is
# actually nested
# croppedRast <- raster::crop(oldRast, study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation)


# convert list of rasters into raster stack
mask_stack <- list()
for(i in 1:length(mask_list)){
  mask_stack[[i]] <- raster::stack(mask_list[[i]])
}

# saveRDS(mask_stack, 'output/frax_preds_rasterstack_by_iteration.RDS')

# iterate bioticVelocity function through each iteration
# sub <- seq(0, 701, by = 33)
times <- sub[-1] * 30

bv_list <- list()
for(i in 1:100){
  bv_list[[i]] <- bioticVelocity(mask_stack[[i]],
                                 times = times,
                                 onlyInSharedCells = TRUE)
}
# 12:50 - 12:58

# plot BVs with uncertainty
require(dplyr)
bv <- bind_rows(bv_list)
bv$timeTo <- as.factor(bv$timeTo)

par(mfrow = c(5,5))
for(i in 1:20){
  h <- hist(bv[levels(bv$timeTo)[i],]$centroidVelocity)
  print(h)
}

ggplot(bv, aes(x = timeTo, y = centroidVelocity)) +
  geom_boxplot()



##################################
## CALCULATING REFUGE LOCATIONS ##
##################################
require(holoSimCell)

# read in 'simulated' dataframe
sim <- raster('data/map-data/study_region_resampled_to_genetic_demographic_simulation_resolution.tif')
proj <- proj4string(sim)

# read in fraxinus predictions at the LGM (use mean of iterations)
grid <- readRDS('data/grid_Dec15.RDS')
frax_preds <- readRDS('output/preds_frax_LGM_mean.RDS')
frax_df <- data.frame(cbind(frax_preds, grid))  # cbind means with locs

# convert df to raster
coordinates(frax_df) <- ~x + y
gridded(frax_df) <- TRUE
frax_raster <- raster(frax_df)
proj4string(frax_raster) <- proj

# mask out glacial coverage and lakes/ocean
# read in masking rasters, extract mask at LGM
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
LGM <- stack$study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.1  # LGM
LGM[LGM > 0.6] <- NA  # convert glacial values to NAs
LGM_resamp <- raster::resample(x = LGM, y = frax_raster)  # convert LGM to same res/extent as frax
frax_masked <- mask(x = frax_raster, mask = LGM_resamp)


# refuge locations
refuge <- assign_refugia_from_abundance_raster(abund = frax_masked, 
                                     sim = sim, quant = 0.9)
plot(refuge)

# convert refuge locations to dataframe so you can plot on a map
refuge_df <- as.data.frame(refuge, xy = TRUE)
refuge_df$id <- as.integer(refuge_df$id)
refuge_df$id_fac <- as.factor(refuge_df$id)

hist(refuge_df$abundance, 
     main = 'Pollen abund preds on refugia\n0.9 quantile threshold', 
     xlab = '',
     breaks = 10)

ylim <- c(frax_masked@extent@ymin, frax_masked@extent@ymax)
xlim <- c(frax_masked@extent@xmin, frax_masked@extent@xmax)

# north america shapefiles
na_shp <- readOGR("data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("data/map-data/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj)

# if you want to plot raw data on top of refugia, run this code
# READ RAW DATA
version <- 'Dec15'
dat <- readRDS(paste0('output/polya-gamma-dat_', version, '.RDS'))
taxa <- dat$taxa.keep
frax_num <- which(taxa == 'Fraxinus')

frax <- data.frame(dat$y[, frax_num,])
frax <- cbind(dat$locs, frax)
frax$x <- frax$x * 1000
frax$y <- frax$y * 1000

frax <- pivot_longer(data = frax, cols = X1:X22, 
                     names_to = 'time', values_to = 'abund')
frax$time <- substr(frax$time, start = 2, stop = nchar(frax$time))
frax$time <- as.integer(frax$time)



ggplot(data = refuge_df) +
  geom_tile(aes(x = x, y = y, fill = id_fac)) +
  scale_fill_discrete(na.value = 'white') +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), alpha = 0.5) +
  geom_path(data = lake_shp, aes(x = long, y = lat, group = group), alpha = 0.5) +
  scale_y_continuous(limits = ylim) +
  scale_x_continuous(limits = xlim) +
  geom_point(data = frax[frax$time == 21 & !is.na(frax$abund),], aes(x = x, y = y),
             pch = 21, color = 'black', fill = 'black', alpha = 0.7, size = 3) +
  labs(fill = 'Refuge ID', title = 'Threshold: 0.8 quantile \n6 refugia') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        # legend.title = element_text(size = 16),
        # legend.text = element_text(size = 14),
        legend.position = 'none') +
  coord_equal()

