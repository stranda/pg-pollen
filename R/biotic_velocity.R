
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
require(holoSimCell)

setwd('C:/Users/abrow/Documents/pg-pollen')
version <- 'Dec15'

#### READ IN DATA ####
# read in raw data
dat <- readRDS(paste0('output/polya-gamma-dat_', version, '.RDS'))
taxa <- dat$taxa.keep
# frax_num <- which(taxa == 'Fraxinus')

# read in prediction grid locations
locs_grid <- readRDS(paste0('data/grid_', version, '.RDS'))

# read in predictions from all simulations (ie, not summarized as mean/median)
preds <- readRDS(paste0('output/preds_', version, '_all_taxa_n20.RDS'))


# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
names(stack) <- 1:701
sub <- seq(1, 701, by = 33)
stack_sub <- subset(stack, subset = paste0('X', sub))  # only want mask every 990 years
stack_sub[stack_sub > 0.6] <- NA
proj <- proj4string(stack)


#### DATA PREPARATION ####
# convert rasterstack masks into list of rasters, reverse order of items to correspond with pollen time
mask_list <- unstack(stack_sub)
mask_names <- NA
for(i in 1:n_times){
  mask_names[i] <- names(mask_list[[i]])
  mask_names[i] <- as.integer(substr(mask_names[i], start = 2, stop = nchar(mask_names[i]))) * 30
}
mask_names <- rev(mask_names)
names(mask_list) <- mask_names
mask_list <- rev(mask_list)

# convert prediction array to list of rasterstacks
n_iter <- dim(preds)[1]
n_taxa <- dim(preds)[3]
n_times <- dim(preds)[4]

# this takes about a minute to run with dimensions [20, 7221, 14, 22]
preds_list <- rep(list(list()), n_taxa) 
for(i in 1:n_taxa){
  for(j in 1:n_iter){
    preds_list[[i]][[j]] <- cbind(locs_grid, data.frame(preds[j,,i,]))
    preds_list[[i]][[j]] <- pivot_longer(preds_list[[i]][[j]], cols = 3:ncol(preds_list[[i]][[j]]),
                                         names_to = 'time', values_to = 'pred')
    preds_list[[i]][[j]]$time <- as.integer(substr(preds_list[[i]][[j]]$time, start = 2, 
                                                 stop = nchar(preds_list[[i]][[j]]$time)))
    preds_list[[i]][[j]] <- split(preds_list[[i]][[j]], f = preds_list[[i]][[j]]$time)
    preds_list[[i]][[j]] <- lapply(preds_list[[i]][[j]], function(x) { x['time'] <- NULL; x })
  }
}

# takes about 20 minutes to run with dimensions [20, 7221, 14, 22]
rasterstack_list <- preds_list
for(i in 1:n_taxa){
  for(j in 1:n_iter){
    for(k in 1:n_times){
      rasterstack_list[[i]][[j]][[k]] <- rasterFromXYZ(rasterstack_list[[i]][[j]][[k]])
      proj4string(rasterstack_list[[i]][[j]][[k]]) <- proj
      rasterstack_list[[i]][[j]][[k]] <- resample(rasterstack_list[[i]][[j]][[k]], y = mask_list[[k]])
      rasterstack_list[[i]][[j]][[k]] <- mask(rasterstack_list[[i]][[j]][[k]], mask = mask_list[[k]])
    }
  }
}

# convert list of rasters into raster stack
for(i in 1:n_taxa){
  for(j in 1:n_iter){
    rasterstack_list[[i]][[j]] <- raster::stack(rasterstack_list[[i]][[j]])
  }
}
saveRDS(rasterstack_list, 'output/preds_rasterstack_list_for_bv_calc.RDS')


##############################################
## CALCULATE AND PLOT BIOTIC VELOCITIES
##############################################
# iterate bioticVelocity function through each iteration
# sub <- seq(1, 701, by = 33)
times <- rev(-(sub * 30))

bv_list <- replicate(n_taxa, list(replicate(n_iter, list())))
bv_list <- replicate(n_taxa, list(replicate(n_iter, data.frame())))

for(i in 1:n_taxa){
  for(j in 1:n_iter){
    bv_list[[i]][[j]] <- bioticVelocity(rasterstack_list[[i]][[j]],
                                 times = times,
                                 onlyInSharedCells = TRUE)
  }
}
# 12:50 - 12:58


# LEFT OFF HERE
# STUCK ON ABOVE LOOP - GIVES ERROR: 
# Error in array(x, c(length(x), 1L), if (!is.null(names(x))) list(names(x),  : 
# length of 'dimnames' [1] not equal to array extent


#### PLOT BVS WITH UNCERTAINTY ####

bv_list2 <- bv_list
for(i in 1:n_taxa){
  bv_list2[[i]] <- bind_rows(bv_list[[i]])
}

bv$timeTo <- as.factor(bv$timeTo)

par(mfrow = c(5,5))
for(i in 1:20){
  h <- hist(bv[levels(bv$timeTo)[i],]$centroidVelocity)
  print(h)
}

ggplot(bv, aes(x = timeTo, y = centroidVelocity)) +
  geom_boxplot()



##################################
## CALCULATE REFUGE LOCATIONS
##################################

# read in 'simulated' dataframe
sim <- raster('data/map-data/study_region_resampled_to_genetic_demographic_simulation_resolution.tif')
proj <- proj4string(sim)

# read in fraxinus predictions at the LGM (use mean of iterations)
grid <- readRDS('data/grid_3.0.RDS')
frax_preds <- readRDS('output/preds_paleo_mean_frax.RDS')
frax_df <- data.frame(cbind(frax_preds, grid))  # cbind means with locs
frax_df <- frax_df[,c('X42','x','y')]  # use only LGM data
names(frax_df) <- c('pred','x','y')

# convert df to raster
coordinates(frax_df) <- ~x + y
gridded(frax_df) <- TRUE
frax_raster <- raster(frax_df)
proj4string(frax_raster) <- proj

# mask out glacial coverage and lakes/ocean
# read in masking rasters, extract mask at LGM
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
LGM <- stack$study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.1  # LGM
LGM[LGM == 1] <- NA  # convert glacial values to NAs
LGM_resamp <- raster::resample(x = LGM, y = frax_raster)  # convert LGM to same res/extent as frax
frax_masked <- mask(x = frax_raster, mask = LGM_resamp)
saveRDS(frax_masked, 'output/LGM_frax_masked_for_refuge.RDS')

# refuge locations
refuge <- assignRefugiaFromAbundanceRaster(abund = frax_masked, 
                                     sim = sim, threshold = 0.02)
plot(refuge$simulationScale)

# convert refuge locations to dataframe so you can plot on a map
refuge_df <- as.data.frame(refuge$simulationScale, xy = TRUE)
refuge_df$refugiaId <- as.integer(refuge_df$refugiaId)
refuge_df$refugiaId_fac <- as.factor(refuge_df$refugiaId)

hist(refuge_df$refugiaAbund, 
     main = 'Pollen abund preds on refugia\n0.02 rel abund threshold', 
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
  geom_tile(aes(x = x, y = y, fill = refugiaId_fac)) +
  scale_fill_discrete(na.value = 'white') +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), alpha = 0.5) +
  geom_path(data = lake_shp, aes(x = long, y = lat, group = group), alpha = 0.5) +
  scale_y_continuous(limits = ylim) +
  scale_x_continuous(limits = xlim) +
  # geom_point(data = frax[frax$time == 21 & !is.na(frax$abund),], aes(x = x, y = y),
  #            pch = 21, color = 'black', fill = 'black', alpha = 0.7, size = 3) +
  labs(fill = 'Refuge ID', title = 'Threshold: 0.02 rel abund \n1 refuge') +
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





#### INTERPOLATING TO 30-YEAR TIME RESOLUTION
# require(enmSdm)
# require(raster)

paleo <- readRDS('output/preds_paleo_mean_frax.RDS')
modern <- readRDS('output/preds_modern_mean_frax.RDS')
# taxa <- readRDS('data/pollen_taxa_3.0.RDS')
locs <- readRDS('data/grid_3.0.RDS')
proj <- '+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs'

# CONVERT PREDICTION ARRAYS TO RASTERSTACK
n_times <- dim(modern)[2]
modern_list <- list()
for(i in 1:n_times){
  modern_list[[i]] <- as.data.frame(modern[,i])
  names(modern_list[[i]]) <- 'pred'
  modern_list[[i]] <- cbind(locs, modern_list[[i]])
  modern_list[[i]] <- rasterFromXYZ(modern_list[[i]], crs = CRS(proj))
}

n_times <- dim(paleo)[2]
paleo_list <- list()
for(i in 1:n_times){
  paleo_list[[i]] <- as.data.frame(paleo[,i])
  names(paleo_list[[i]]) <- 'pred'
  paleo_list[[i]] <- cbind(locs, paleo_list[[i]])
  paleo_list[[i]] <- rasterFromXYZ(paleo_list[[i]], crs = CRS(proj))
}

all_list <- c(modern_list, paleo_list)

# MASK OUT SEAS, LAKES, ICE
# rasterstack masks are arranged LGM (item 1) to modern (item 701)
# I think modern = 2000 A.D.
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
stack_list <- unstack(stack)
stack_list <- rev(stack_list)

stack_list_c <- stack_list
for(i in 1:length(stack_list_c)){
  stack_list_c[[i]][stack_list_c[[i]] == 1] <- NA
}

# Specify pollen times (46 total; 4 modern and 42 paleo)
times <- c(seq(0, 150, by = 50), seq(from = 650, by = 500, length.out = n_times))

# create vector of ice mask times
adam_times <- seq(0, 700, 1) * 30

# interpolate pollen rasters to temporal resolution of 30yrs
all_stack <- stack(all_list)
interp <- interpolateRasters(all_stack, interpFrom = times, interpTo = adam_times)
interp_list <- unstack(interp)

# mask out ice/water
n_times <- length(stack_list_c)
mask_list <- interp_list
for(i in 1:n_times){
  mask_list[[i]] <- raster::resample(mask_list[[i]], y = stack_list_c[[i]])
  mask_list[[i]] <- mask(mask_list[[i]], mask = stack_list_c[[i]])
}

# GET BIOTIC VELOCITY FROM INTERPOLATED/MASKED RASTERS
mask_list <- rev(mask_list)
mask_stack <- stack(mask_list)

bv_times <- rev(adam_times * -1)
bvs_interp <- bioticVelocity(mask_stack, times = bv_times, onlyInSharedCells = TRUE)
saveRDS(bvs_interp, 'output/interp_bvs.RDS')

minor_breaks <- rev(times$from) * -1
breaks <- minor_breaks[seq(1,42, by = 2)]

ggplot(data = bvs_interp, aes(x = timeFrom, y = centroidVelocity)) +
  geom_line() + 
  scale_x_continuous(minor_breaks = minor_breaks, breaks = breaks) +
  labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))
ggsave('figures/bvs_after_30yr_interp.png')




# OLD CODE

# CALCULATE BV USING FRAX PREDS
bv_stack <- stack(frax_mask)
sub <- seq(0, 701, by = 33)
times <- sub[-1] * 30
bv <- bioticVelocity(bv_stack, times = times) # , onlyInSharedCells = TRUE)
saveRDS(bv, 'output/frax_bvs.RDS')

# PLOTTING BVS OVER TIME
par(mfrow = c(3,1))
par(mar = c(4, 5, 2, 2))
plot(bv$timeTo, bv$centroidVelocity, type = 'b', cex = 2, lwd = 2,
     ylab = 'Centroid velocity (m/yr)', xlab = '', main = 'Fraxinus pollen biotic velocities',
     xlim = rev(range(bv$timeTo)), cex.lab = 2, cex.axis = 2, cex.main = 2)
grid(col = 'gray10')
plot(bv$timeTo, bv$nsQuantVelocity_quant0p95, type = 'b', cex = 2, lwd = 2,
     ylab = '95% quantile velocity', xlab = '', main = '',
     xlim = rev(range(bv$timeTo)), cex.lab = 2, cex.axis = 2)
grid(col = 'gray10')
abline(h = 0, col = 'black', lwd = 2)
plot(bv$timeTo, bv$nsQuantVelocity_quant0p05, type = 'b', cex = 2, lwd = 2,
     ylab = '5% quantile velocity', xlab = 'Years before present', main = '',
     xlim = rev(range(bv$timeTo)), cex.lab = 2, cex.axis = 2)
grid(col = 'gray10')
abline(h = 0, col = 'black', lwd = 2)
# dev.off()

# PLOT CENTROIDS OF FRAX RANGE OVER TIME ON MAP
# use biotic velocity dataframe (bv)
bvs <- readRDS('output/no-interp_bvs.RDS')
ggplot(data = bvs) +
  geom_point(aes(x = centroidLong, y = centroidLat, fill = timeTo),
             size = 4, alpha = 0.9, pch = 21) +
  # geom_text(aes(x = centroidLong, y = centroidLat, 
  #               label = ifelse(timeTo > -7000, as.character(timeTo),''))) +
  geom_label_repel(aes(x = centroidLong, y = centroidLat, 
                       label = ifelse(timeTo > -3000, as.character(timeTo),'')),
                   box.padding   = 0.5,
                   point.padding = 0.5,
                   segment.color = 'red',
                   max.overlaps = 30) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(-700000,1500000)) +
  scale_x_continuous(limits = c(-500000,2500000)) +
  labs(title = 'Centroid BVs (no interpolation)') +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16)) +
  coord_equal()

bvs_interp <- readRDS('output/interp_bvs.RDS')
ggplot(data = bvs_interp) +
  geom_point(aes(x = centroidLong, y = centroidLat, fill = timeTo),
             size = 4, alpha = 0.9, pch = 21) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(-700000,1500000)) +
  scale_x_continuous(limits = c(-500000,2500000)) +
  labs(title = 'Centroid BVs (no interpolation)') +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16)) +
  coord_equal()


# INSTALLING MOST RECENT VERSION OF ENMSDM
# download zip file from repo: https://github.com/adamlilith/enmSdm/tree/master/zipTarFiles
# might also need to update omnibus and satisfactory, also on adam's github page
install.packages('C:/Users/abrow/Desktop/postdoc/enmSdm_0.5.2.9.zip', lib='C:/Users/abrow/Documents/R/win-library/4.0',repos = NULL)

install.packages('C:/Users/abrow/Desktop/postdoc/omnibus_0.3.4.7.zip', lib='C:/Users/abrow/Documents/R/win-library/4.0',repos = NULL)
install.packages('C:/Users/abrow/Desktop/postdoc/statisfactory_0.3.4.zip', lib='C:/Users/abrow/Documents/R/win-library/4.0',repos = NULL)
