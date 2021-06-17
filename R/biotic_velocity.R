
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
version <- '3.0'

#### READ IN FRAXINUS DATA ####
# prediction grid, paleo/modern predictions
locs_grid <- readRDS(paste0('data/grid_', version, '.RDS'))
preds_paleo <- readRDS('output/preds_paleo_n50_3.0.RDS')
preds_mod <- readRDS('output/preds_modern_n50_3.0.RDS')
preds <- array(c(preds_mod, preds_paleo), dim = c(50, 3951, 46))

modern_bins <- seq(-50, 149, by = 50)
paleo_bins <- seq(150, by = 500, length.out = 43)
time <- c(modern_bins, paleo_bins)
time <- data.frame(timeFrom = time[-47], timeTo = time[-1])
time$time_mid <- (time$timeFrom + time$timeTo) / 2  # to get midpoint of time bins
time$time_mid <- time$time_mid + 25  # to align with 0 start date


# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
proj <- proj4string(stack)
stack_times <- rev(seq(0, by = 30, length.out = 701))
time$closest <- NA
for(i in 1:nrow(time)){
  time$closest[i] <- which.min(abs(stack_times - time$time_mid[i]))
}
# stack_times_sub <- stack_times[time$closest]
# time$closest_time <- stack_times_sub
stack_list <- unstack(stack)
stack_sub <- stack_list[time$closest]

for(i in 1:length(stack_sub)){
  stack_sub[[i]][stack_sub[[i]] == 1] <- NA  # convert areas completely covered by ice to NAs
}



#### DATA PREPARATION ####
# convert prediction array to list of rasterstacks
n_iter <- dim(preds)[1]
n_times <- dim(preds)[3]

preds_list <- list() 
for(i in 1:n_iter){
    preds_list[[i]] <- cbind(locs_grid, data.frame(preds[i,,]))
    preds_list[[i]] <- pivot_longer(preds_list[[i]], cols = 3:ncol(preds_list[[i]]),
                                         names_to = 'time', values_to = 'pred')
    preds_list[[i]]$time <- as.integer(substr(preds_list[[i]]$time, start = 2, 
                                                 stop = nchar(preds_list[[i]]$time)))
    preds_list[[i]] <- split(preds_list[[i]], f = preds_list[[i]]$time)
    preds_list[[i]] <- lapply(preds_list[[i]], function(x) { x['time'] <- NULL; x })
}


# takes about 2 minutes to run with dimensions [50, 3951, 46]
rasterstack_list <- preds_list
for(i in 1:n_iter){
  for(j in 1:n_times){
      rasterstack_list[[i]][[j]] <- rasterFromXYZ(rasterstack_list[[i]][[j]])
      proj4string(rasterstack_list[[i]][[j]]) <- proj
      rasterstack_list[[i]][[j]] <- resample(rasterstack_list[[i]][[j]], y = stack_sub[[j]])
      rasterstack_list[[i]][[j]] <- mask(rasterstack_list[[i]][[j]], mask = stack_sub[[j]])
  }
}

# convert list of rasters into raster stack
for(i in 1:n_iter){
    rasterstack_list[[i]] <- raster::stack(rasterstack_list[[i]])
}
saveRDS(rasterstack_list, 'output/rasterstacks_for_bv_calc_n50.RDS')


##############################################
## CALCULATE AND PLOT BIOTIC VELOCITIES
##############################################
# iterate bioticVelocity function through each iteration (takes ~6 minutes)
bv_list <- list()
for(i in 1:n_iter){
    bv_list[[i]] <- bioticVelocity(rasterstack_list[[i]],
                                 times = time$time_mid,
                                 onlyInSharedCells = TRUE)
}

for(i in 1:n_iter){
  bv_list[[i]]$iter <- i
}
bvs <- bind_rows(bv_list)
saveRDS(bvs, 'output/no_interp_bvs_n50.RDS')

#### PLOT BVS WITH UNCERTAINTY ####
bvs$median_time <- as.factor((bvs$timeFrom + bvs$timeTo) / 2)

ggplot(bvs, aes(x = median_time, y = centroidVelocity)) +
  geom_boxplot() + 
  theme_classic()




##################################
## CALCULATE REFUGE LOCATIONS
##################################

# read in 'simulated' dataframe
sim <- raster('data/map-data/study_region_resampled_to_genetic_demographic_simulation_resolution.tif')
proj <- proj4string(sim)

# read in fraxinus predictions at the LGM for each of the 50 iterations
grid <- readRDS('data/grid_3.0.RDS')
frax_lgm <- readRDS('output/preds_paleo_n50_3.0.RDS')
frax_lgm <- frax_lgm[,,42]

n_iter <- 50
lgm_list <- list()
for(i in 1:n_iter){
  lgm_list[[i]] <- data.frame(cbind(grid, frax_lgm[i,]))
  names(lgm_list[[i]]) <- c('x','y','pred')
}

# convert df to raster
lgm_raster <- list()
for(i in 1:n_iter){
  lgm_raster[[i]] <- rasterFromXYZ(lgm_list[[i]], crs = CRS(proj))
}


# mask out glacial coverage and lakes/ocean
# read in masking rasters, extract mask at LGM
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
LGM <- stack$study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.1  # LGM
LGM[LGM == 1] <- NA  # convert glacial values to NAs

lgm_resamp <- list()
for(i in 1:n_iter){
  lgm_resamp[[i]] <- raster::resample(x = lgm_raster[[i]], y = LGM)
  lgm_resamp[[i]] <- mask(x = lgm_resamp[[i]], mask = LGM)
}
saveRDS(lgm_resamp, 'output/lgm_frax_masked_for_refuge_n50.RDS')

# refuge locations
refuge_list <- list()
for(i in 1:n_iter){
  refuge_list[[i]] <- assignRefugiaFromAbundanceRaster(abund = lgm_resamp[[i]],
                                                       sim = sim, threshold = 0.01)
}

# check to make sure all iterations have been assigned a refuge
test <- rep(NA, n_iter)
for(i in 1:n_iter){
  test[i] <- refuge_list[[i]]$meanRefugeAbund
}
anyNA(test)  # should be FALSE

# save as list of rasters
refuge_for_allan <- list()
for(i in 1:n_iter){
  refuge_for_allan[[i]] <- refuge_list[[i]][["simulationScale"]]@layers[[2]]
}
saveRDS(refuge_for_allan, 'output/refuge_rasters_n50.RDS')

# check out all refugia plots
for(i in 1:n_iter){
  jpeg(paste('figures/refugia_gif/iter', i, '.jpeg', sep = ''))
  plot(refuge_list[[i]]$simulationScale@layers[[2]])
  dev.off()
}
# make gif of plots
require(magick)
frames <- paste('figures/refugia_gif/iter', 1:50, '.jpeg', sep = '')
m <- image_read(frames)
m <- image_animate(m, fps = 1)
image_write(m, "figures/refugia_gif/refugia_n50.gif")


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
  labs(fill = 'Refuge ID', title = 'Threshold: 0.04 rel abund \n2 refugia') +
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
ggsave('figures/frax_refuge_0.04threshold_map.png')




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
# saveRDS(all_list, 'output/preds_mean_frax_all_times_rasters.RDS')

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
# saveRDS(interp, 'output/interp_preds_mean_frax_rasterstack.RDS')

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





# PLOT BV FOR ALL 5 METHODS
require(gridExtra)
# read bv data
setwd('C:/Users/abrow/Documents/green_ash')
abc_enm <- read.csv('BV table ENM nnet_0.1.csv', stringsAsFactors = FALSE)
abc <- read.csv('BV table NAIVE nnet_0.1.csv', stringsAsFactors = FALSE)
abc_pollen <- read.csv('BV_table_POLLEN_nnet_0.1.csv', stringsAsFactors = FALSE)
load('biotic_velocities.rda')
enm <- velocities

setwd('C:/Users/abrow/Documents/pg-pollen')
pollen <- readRDS('output/no_interp_bvs_n50.RDS')

# calculate midpoint of time range for each dataset (if it doesn't already exist)
abc$timeFrom <- as.numeric(sub('-.*', '', abc$time))
abc$timeTo <- as.numeric(sub('.*-', '', abc$time))
abc$median_time <- (abc$timeFrom + abc$timeTo) / 2

abc_enm$timeFrom <- as.numeric(sub('-.*', '', abc_enm$time))
abc_enm$timeTo <- as.numeric(sub('.*-', '', abc_enm$time))
abc_enm$median_time <- (abc_enm$timeFrom + abc_enm$timeTo) / 2

abc_pollen$timeFrom <- as.numeric(sub('-.*', '', abc_pollen$time))
abc_pollen$timeTo <- as.numeric(sub('.*-', '', abc_pollen$time))
abc_pollen$median_time <- (abc_pollen$timeFrom + abc_pollen$timeTo) / 2

enm$median_time <- (abs(enm$timeFrom) + abs(enm$timeTo)) / 2

# plot
# pollen
# for now, use 21 times most closely associated with 990-yr intervals
levels_use <- seq(5, 46, by = 2)
times_use <- levels(bvs$median_time)[levels_use]
bvs_sub <- bvs %>% filter(median_time %in% times_use)
bvs_sub$median_time <- droplevels(bvs_sub$median_time)

plot_pollen <- ggplot(bvs_sub, aes(x = median_time, y = centroidVelocity)) +
  geom_boxplot() + 
  ylab('\n ') +
  geom_text(x = 3, y = 1850, label = 'Pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 12))

# enm
enm_sub <- enm %>% filter(timeSpan == 990)
enm_sub$median_time <- as.factor(enm_sub$median_time)

plot_enm <- ggplot(enm_sub, aes(x = median_time, y = centroidVelocity)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(0, 400)) +
  ylab('\n ') +
  geom_text(x = 3, y = 375, label = 'ENM', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# abc naive
abc$median_timeF <- as.factor(abc$median_time)

plot_abc <- ggplot(abc, aes(x = median_time, y = BVcent))+ 
  geom_boxplot(middle = abc$BVcent, 
               lower = abc$BVcent0p025,
               upper = abc$BVcent0p975)+
  scale_y_continuous(limits = c(0, 400)) +
  ylab('Centroid BV (m/yr)\n') +
  geom_text(x = 3, y = 375, label = 'Naive ABC', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# abc-enm
abc_enm$median_time <- as.factor(abc_enm$median_time)

plot_abc_enm <- ggplot(abc_enm, aes(x = median_time, y = BVcent))+ 
  geom_boxplot(middle = abc_enm$BVcent, 
               lower = abc_enm$BVcent0p025,
               upper = abc_enm$BVcent0p975)+
  scale_y_continuous(limits = c(0, 400)) +
  ylab('\n ') +
  xlab('Years before present') +
  geom_text(x = 3, y = 375, label = 'ABC-ENM', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# abc-pollen
abc_pollen$median_time <- as.factor(abc_pollen$median_time)

plot_abc_pollen <- ggplot(abc_pollen, aes(x = median_time, y = BVcent))+ 
  geom_boxplot(middle = abc_pollen$BVcent, 
               lower = abc_pollen$BVcent0p025,
               upper = abc_pollen$BVcent0p975)+
  scale_y_continuous(limits = c(0, 400)) +
  ylab('\n ') +
  geom_text(x = 3, y = 375, label = 'ABC-pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# grid.arrange(plot_pollen, plot_enm, plot_abc, plot_abc_pollen, plot_abc_enm, ncol = 1)
g <- arrangeGrob(plot_pollen, plot_enm, plot_abc, plot_abc_pollen, plot_abc_enm, 
                 ncol = 1, heights = c(1,1,1,1,1.35))
ggsave(file = 'figures/BVs_all_methods.png', plot = g, height = 10, width = 7, units = 'in')


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



# OLD CODE, BUT MIGHT USE AGAIN - FOR DEALING WITH 4D ARRAY OF PREDS USING ALL TAXA
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
