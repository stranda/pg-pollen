
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





# PLOT CENTROID BV FOR ALL 5 METHODS
require(gridExtra)
# read bv data
setwd('C:/Users/abrow/Documents/green_ash')
abc_enm <- read.csv('BV table ENM nnet_0.1.csv', stringsAsFactors = FALSE)
abc <- read.csv('BV table NAIVE nnet_0.1.csv', stringsAsFactors = FALSE)
abc_pollen <- read.csv('BV_table_POLLEN_nnet_0.1.csv', stringsAsFactors = FALSE)
load('enm_biotic_velocities.rda')
enm <- velocities
pollen <- readRDS('pollen_no_interp_bvs_n50.RDS')

# pollen
# for now, use 21 times most closely associated with 990-yr intervals
levels_use <- seq(5, 46, by = 2)
times_use <- levels(pollen$median_time)[levels_use]
pollen_sub <- pollen %>% filter(median_time %in% times_use)
pollen_sub$median_time <- droplevels(pollen_sub$median_time)

# summarize across 50 iterations to get 97.5 and 2.5 quantiles
pollen_summ <- pollen_sub %>%
  group_by(median_time) %>%
  summarize(BVcent = quantile(centroidVelocity, 0.5),
            BVcent0p025 = quantile(centroidVelocity, 0.025),
            BVcent0p975 = quantile(centroidVelocity, 0.975),
            BVN = quantile(nCentroidVelocity, 0.5),
            BVN0p025 = quantile(nCentroidVelocity, 0.025),
            BVN0p975 = quantile(nCentroidVelocity, 0.975),
            BVS = quantile(sCentroidVelocity, 0.5),
            BVS0p025 = quantile(sCentroidVelocity, 0.025),
            BVS0p975 = quantile(sCentroidVelocity, 0.975)
            )
saveRDS(pollen_summ,'pollen_bvs_quants_for_plotting.RDS')

pollen <- readRDS('pollen_bvs_quants_for_plotting.RDS')
pollen <- pollen %>% arrange(desc(median_time))  # arrange so time goes from LGM to modern (top to bottom)

plot_pollen <- ggplot(pollen, aes(x = median_time, y = BVcent)) +
  geom_boxplot(middle = pollen$BVcent,
               lower = pollen$BVcent0p025,
               upper = pollen$BVcent0p975) + 
  scale_y_continuous(limits = c(0,1500)) +
  ylab('\n ') +
  geom_text(x = 3, y = 1450, label = 'Pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 12))

# enm
enm_sub <- enm %>% filter(timeSpan == 990, onlyInSharedCells == TRUE, onlyInContinuouslyExposedLand == FALSE)
enm_sub$median_time <- as.factor(enm_sub$median_time)

# summarize across 24 models to get 97.5 and 2.5 quantiles
enm_summ <- enm_sub %>% 
  group_by(median_time) %>% 
  summarize(BVcent = quantile(centroidVelocity, 0.5),
            BVcent0p025 = quantile(centroidVelocity, 0.025),
            BVcent0p975 = quantile(centroidVelocity, 0.975), 
            BVN = quantile(nCentroidVelocity, 0.5),
            BVN0p025 = quantile(nCentroidVelocity, 0.025),
            BVN0p975 = quantile(nCentroidVelocity, 0.975),
            BVS = quantile(sCentroidVelocity, 0.5),
            BVS0p025 = quantile(sCentroidVelocity, 0.025),
            BVS0p975 = quantile(sCentroidVelocity, 0.975)
            )
saveRDS(enm_summ,'enm_bvs_quants_for_plotting.RDS')

enm <- readRDS('enm_bvs_quants_for_plotting.RDS')
enm <- enm %>% arrange(desc(median_time))

plot_enm <- ggplot(enm, aes(x = median_time, y = BVcent)) + 
  geom_boxplot(middle = enm$BVcent,
               lower = enm$BVcent0p025,
               upper = enm$BVcent0p975) + 
  scale_y_continuous(limits = c(0, 400)) +
  ylab('\n ') +
  geom_text(x = 3, y = 375, label = 'ENM', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# abc naive
abc$median_timeF <- as.factor(abc$median_time)

plot_abc <- ggplot(abc, aes(x = time, y = BVcent))+ 
  geom_boxplot(middle = abc$BVcent, 
               lower = abc$BVcent0p025,
               upper = abc$BVcent0p975)+
  scale_y_continuous(limits = c(0, 400)) +
  ylab('Centroid BV (m/yr)\n') +
  geom_text(x = 3, y = 375, label = 'Naive ABC', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# abc-enm
abc_enm$median_time <- as.factor(abc_enm$median_time)

plot_abc_enm <- ggplot(abc_enm, aes(x = time, y = BVcent))+ 
  geom_boxplot(middle = abc_enm$BVcent, 
               lower = abc_enm$BVcent0p025,
               upper = abc_enm$BVcent0p975)+
  scale_x_discrete(limits = c("21000-20010","20010-19020","19020-18030","18030-17040",
                              "17040-16050","16050-15060","15060-14070","14070-13080",
                              "13080-12090","12090-11100","11100-10110","10110-9120",
                              "9120-8130","8130-7140","7140-6150","6150-5160",
                              "5160-4170","4170-3180","3180-2190","2190-1200","1200-210"))+
  scale_y_continuous(limits = c(0, 400)) +
  ylab('\n ') +
  xlab('Years before present') +
  geom_text(x = 3, y = 375, label = 'ABC-ENM', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# abc-pollen
abc_pollen$median_time <- as.factor(abc_pollen$median_time)

plot_abc_pollen <- ggplot(abc_pollen, aes(x = time, y = BVcent))+ 
  geom_boxplot(middle = abc_pollen$BVcent, 
               lower = abc_pollen$BVcent0p025,
               upper = abc_pollen$BVcent0p975)+
  scale_y_continuous(limits = c(0, 400)) +
  ylab('\n ') +
  geom_text(x = 3, y = 375, label = 'ABC-pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# grid.arrange(plot_pollen, plot_enm, plot_abc, plot_abc_pollen, plot_abc_enm, ncol = 1)
g <- arrangeGrob(plot_pollen, plot_enm, plot_abc, plot_abc_pollen, plot_abc_enm, 
                 ncol = 1, heights = c(1,1,1,1,1.6))
ggsave(file = 'figures/BVs_all_methods.png', plot = g, height = 10, width = 7, units = 'in')




# PLOT NORTH CENTROID BVS FOR ALL 5 METHODS
plot_pollen <- ggplot(pollen, aes(x = median_time, y = BVN)) +
  geom_boxplot(middle = pollen$BVN,
               lower = pollen$BVN0p025,
               upper = pollen$BVN0p975) + 
  scale_y_continuous(limits = c(-200,1700)) +
  ylab('\n\n ') +
  geom_text(x = 3, y = 1600, label = 'Pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12))

plot_enm <- ggplot(enm, aes(x = median_time, y = BVN)) + 
  geom_boxplot(middle = enm$BVN,
               lower = enm$BVN0p025,
               upper = enm$BVN0p975) + 
  scale_y_continuous(limits = c(-200, 400)) +
  ylab('\n\n ') +
  geom_text(x = 3, y = 375, label = 'ENM', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc <- ggplot(abc, aes(x = time, y = BVN))+ 
  geom_boxplot(middle = abc$BVN, 
               lower = abc$BVN0p025,
               upper = abc$BVN0p975)+
  scale_y_continuous(limits = c(-200, 400)) +
  ylab('North-centroid\nBV (m/yr)\n') +
  geom_text(x = 3, y = 375, label = 'Naive ABC', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc_pollen <- ggplot(abc_pollen, aes(x = time, y = BVN))+ 
  geom_boxplot(middle = abc_pollen$BVN, 
               lower = abc_pollen$BVN0p025,
               upper = abc_pollen$BVN0p975)+
  scale_y_continuous(limits = c(-200, 400)) +
  ylab('\n\n ') +
  geom_text(x = 3, y = 375, label = 'ABC-pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc_enm <- ggplot(abc_enm, aes(x = time, y = BVcent))+ 
  geom_boxplot(middle = abc_enm$BVcent, 
               lower = abc_enm$BVcent0p025,
               upper = abc_enm$BVcent0p975)+
  scale_x_discrete(limits = c("21000-20010","20010-19020","19020-18030","18030-17040",
                              "17040-16050","16050-15060","15060-14070","14070-13080",
                              "13080-12090","12090-11100","11100-10110","10110-9120",
                              "9120-8130","8130-7140","7140-6150","6150-5160",
                              "5160-4170","4170-3180","3180-2190","2190-1200","1200-210"))+
  scale_y_continuous(limits = c(-200, 400)) +
  ylab('\n\n ') +
  xlab('Years before present') +
  geom_text(x = 3, y = 375, label = 'ABC-ENM', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

g <- arrangeGrob(plot_pollen, plot_enm, plot_abc, plot_abc_pollen, plot_abc_enm, 
                 ncol = 1, heights = c(1,1,1,1,1.6))
ggsave(file = 'figures/BVNs_all_methods.png', plot = g, height = 10, width = 7, units = 'in')




# PLOT SOUTH CENTROID BVS FOR ALL 5 METHODS
plot_pollen <- ggplot(pollen, aes(x = median_time, y = BVS)) +
  geom_boxplot(middle = pollen$BVS,
               lower = pollen$BVS0p025,
               upper = pollen$BVS0p975) + 
  scale_y_continuous(limits = c(-500,1000)) +
  ylab('\n\n ') +
  geom_text(x = 19, y = -425, label = 'Pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12))

plot_enm <- ggplot(enm, aes(x = median_time, y = BVS)) + 
  geom_boxplot(middle = enm$BVS,
               lower = enm$BVS0p025,
               upper = enm$BVS0p975) + 
  scale_y_continuous(limits = c(-500, 250)) +
  ylab('\n\n ') +
  geom_text(x = 19, y = -450, label = 'ENM', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc <- ggplot(abc, aes(x = time, y = BVS))+ 
  geom_boxplot(middle = abc$BVS, 
               lower = abc$BVS0p025,
               upper = abc$BVS0p975)+
  scale_y_continuous(limits = c(-500, 250)) +
  ylab('South-centroid\nBV (m/yr)\n') +
  geom_text(x = 19, y = -450, label = 'Naive ABC', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc_pollen <- ggplot(abc_pollen, aes(x = time, y = BVS))+ 
  geom_boxplot(middle = abc_pollen$BVS, 
               lower = abc_pollen$BVS0p025,
               upper = abc_pollen$BVS0p975)+
  scale_y_continuous(limits = c(-500, 250)) +
  ylab('\n\n ') +
  geom_text(x = 19, y = -450, label = 'ABC-pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc_enm <- ggplot(abc_enm, aes(x = time, y = BVcent))+ 
  geom_boxplot(middle = abc_enm$BVcent, 
               lower = abc_enm$BVcent0p025,
               upper = abc_enm$BVcent0p975)+
  scale_x_discrete(limits = c("21000-20010","20010-19020","19020-18030","18030-17040",
                              "17040-16050","16050-15060","15060-14070","14070-13080",
                              "13080-12090","12090-11100","11100-10110","10110-9120",
                              "9120-8130","8130-7140","7140-6150","6150-5160",
                              "5160-4170","4170-3180","3180-2190","2190-1200","1200-210"))+
  scale_y_continuous(limits = c(-500, 250)) +
  ylab('\n\n ') +
  xlab('Years before present') +
  geom_text(x = 19, y = -450, label = 'ABC-ENM', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

g <- arrangeGrob(plot_pollen, plot_enm, plot_abc, plot_abc_pollen, plot_abc_enm, 
                 ncol = 1, heights = c(1,1,1,1,1.6))
ggsave(file = 'figures/BVSs_all_methods.png', plot = g, height = 10, width = 7, units = 'in')






#### CALCULATE BVS FOR NEWEST MODEL RUN (AUGUST 2, 2021) ####

# READ IN FRAXINUS DATA
# prediction grid, frax predictions for 50 random iterations
version <- '3.1'
locs_grid <- readRDS(paste0('data/grid_', version, '.RDS'))
preds <- readRDS('output/preds_frax_n50_v3.1_latent.RDS')
preds_over <- readRDS('output/preds_frax_n50_v3.1_overdispersed.RDS')
preds_mat <- readRDS('output/preds_frax_n50_v3.1_matern.RDS')

time <- seq(210, 21000, by = 990)
time <- data.frame(timeFrom = time[-22], timeTo = time[-1])
time$time_mid <- (time$timeFrom + time$timeTo) / 2  # to get midpoint of time bins


# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
proj <- proj4string(stack)

# need to find times in glacier masks closest to median pollen times
stack_times <- rev(seq(0, by = 30, length.out = 701))
time$closest <- NA
for(i in 1:nrow(time)){
  time$closest[i] <- which.min(abs(stack_times - time$time_mid[i]))
}
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

preds_over_list <- list()
for(i in 1:n_iter){
  preds_over_list[[i]] <- cbind(locs_grid, data.frame(preds_over[i,,]))
  preds_over_list[[i]] <- pivot_longer(preds_over_list[[i]], cols = 3:ncol(preds_over_list[[i]]),
                                  names_to = 'time', values_to = 'pred')
  preds_over_list[[i]]$time <- as.integer(substr(preds_over_list[[i]]$time, start = 2,
                                            stop = nchar(preds_over_list[[i]]$time)))
  preds_over_list[[i]] <- split(preds_over_list[[i]], f = preds_over_list[[i]]$time)
  preds_over_list[[i]] <- lapply(preds_over_list[[i]], function(x) { x['time'] <- NULL; x })
}

preds_mat_list <- list()
for(i in 1:n_iter){
  preds_mat_list[[i]] <- cbind(locs_grid, data.frame(preds_mat[i,,]))
  preds_mat_list[[i]] <- pivot_longer(preds_mat_list[[i]], cols = 3:ncol(preds_mat_list[[i]]),
                                       names_to = 'time', values_to = 'pred')
  preds_mat_list[[i]]$time <- as.integer(substr(preds_mat_list[[i]]$time, start = 2,
                                                 stop = nchar(preds_mat_list[[i]]$time)))
  preds_mat_list[[i]] <- split(preds_mat_list[[i]], f = preds_mat_list[[i]]$time)
  preds_mat_list[[i]] <- lapply(preds_mat_list[[i]], function(x) { x['time'] <- NULL; x })
}


# takes 1-2 minutes to run
rasterstack_list <- preds_list
for(i in 1:n_iter){
  for(j in 1:n_times){
    rasterstack_list[[i]][[j]] <- rasterFromXYZ(rasterstack_list[[i]][[j]])
    proj4string(rasterstack_list[[i]][[j]]) <- proj
    rasterstack_list[[i]][[j]] <- raster::resample(rasterstack_list[[i]][[j]], y = stack_sub[[j]])
    rasterstack_list[[i]][[j]] <- mask(rasterstack_list[[i]][[j]], mask = stack_sub[[j]])
  }
}

rasterstack_over_list <- preds_over_list
for(i in 1:n_iter){
  for(j in 1:n_times){
    rasterstack_over_list[[i]][[j]] <- rasterFromXYZ(rasterstack_over_list[[i]][[j]])
    proj4string(rasterstack_over_list[[i]][[j]]) <- proj
    rasterstack_over_list[[i]][[j]] <- raster::resample(rasterstack_over_list[[i]][[j]], y = stack_sub[[j]])
    rasterstack_over_list[[i]][[j]] <- mask(rasterstack_over_list[[i]][[j]], mask = stack_sub[[j]])
  }
}

rasterstack_mat_list <- preds_mat_list
for(i in 1:n_iter){
  for(j in 1:n_times){
    rasterstack_mat_list[[i]][[j]] <- rasterFromXYZ(rasterstack_mat_list[[i]][[j]])
    proj4string(rasterstack_mat_list[[i]][[j]]) <- proj
    rasterstack_mat_list[[i]][[j]] <- raster::resample(rasterstack_mat_list[[i]][[j]], y = stack_sub[[j]])
    rasterstack_mat_list[[i]][[j]] <- mask(rasterstack_mat_list[[i]][[j]], mask = stack_sub[[j]])
  }
}


# convert list of rasters into raster stack
for(i in 1:n_iter){
  rasterstack_list[[i]] <- raster::stack(rasterstack_list[[i]])
}

for(i in 1:n_iter){
  rasterstack_over_list[[i]] <- raster::stack(rasterstack_over_list[[i]])
}

for(i in 1:n_iter){
  rasterstack_mat_list[[i]] <- raster::stack(rasterstack_mat_list[[i]])
}


##############################################
## CALCULATE AND PLOT BIOTIC VELOCITIES
##############################################
# iterate bioticVelocity function through each iteration (takes ~6 minutes)
# first need to reverse the order of rasters in each rasterstack to get the BV function to work
rasterstack_rev <- list()
for(i in 1:n_iter){
  rasterstack_rev[[i]] <- subset(rasterstack_list[[i]], order(21:1))
}
saveRDS(rasterstack_rev, 'output/rasterstacks_for_bv_calc_n50_latent_v3.1.RDS')
rasterstack_rev <- readRDS('output/rasterstacks_for_bv_calc_n50_latent_v3.1.RDS')

rasterstack_over_rev <- list()
for(i in 1:n_iter){
  rasterstack_over_rev[[i]] <- subset(rasterstack_over_list[[i]], order(21:1))
}
saveRDS(rasterstack_over_rev, 'output/rasterstacks_for_bv_calc_n50_overdispersed_v3.1.RDS')
rasterstack_over_rev <- readRDS('output/rasterstacks_for_bv_calc_n50_overdispersed_v3.1.RDS')

rasterstack_mat_rev <- list()
for(i in 1:n_iter){
  rasterstack_mat_rev[[i]] <- subset(rasterstack_mat_list[[i]], order(21:1))
}
saveRDS(rasterstack_mat_rev, 'output/rasterstacks_for_bv_calc_n50_matern_v3.1.RDS')
rasterstack_mat_rev <- readRDS('output/rasterstacks_for_bv_calc_n50_matern_v3.1.RDS')


# might need to force negative values to be 0, otherwise BV function won't work
for(i in 1:n_iter){
  rasterstack_rev[[i]] <- calc(rasterstack_rev[[i]], fun = function(x) ifelse(x < 0, 0, x))
}

for(i in 1:n_iter){
  rasterstack_over_rev[[i]] <- calc(rasterstack_over_rev[[i]], fun = function(x) ifelse(x < 0, 0, x))
}

for(i in 1:n_iter){
  rasterstack_mat_rev[[i]] <- calc(rasterstack_mat_rev[[i]], fun = function(x) ifelse(x < 0, 0, x))
}

# calculate BVs
times <- seq(-20505, by = 990, length.out = 21)

bv_list <- list()
for(i in 1:n_iter){
  bv_list[[i]] <- bioticVelocity(rasterstack_rev[[i]],
                                 # times = rev(time$time_mid) * -1,
                                 times = times,
                                 onlyInSharedCells = TRUE)
}

bv_over_list <- list()
for(i in 1:n_iter){
  bv_over_list[[i]] <- bioticVelocity(rasterstack_over_rev[[i]],
                                 # times = rev(time$time_mid) * -1,
                                 times = times,
                                 onlyInSharedCells = TRUE)
}

bv_mat_list <- list()
for(i in 1:n_iter){
  bv_mat_list[[i]] <- bioticVelocity(rasterstack_mat_rev[[i]],
                                      # times = rev(time$time_mid) * -1,
                                      times = times,
                                      onlyInSharedCells = TRUE)
}

# convert list of dataframes to single dataframe
for(i in 1:n_iter){
  bv_list[[i]]$iter <- i
  bv_list[[i]]$mod_type <- 'latent'
}
bvs <- bind_rows(bv_list)
saveRDS(bvs, 'output/bvs_n50_latent_v3.1.RDS')

for(i in 1:n_iter){
  bv_over_list[[i]]$iter <- i
  bv_over_list[[i]]$mod_type <- 'overdispersed'
}
bvs_over <- bind_rows(bv_over_list)
saveRDS(bvs_over, 'output/bvs_n50_overdispersed_v3.1.RDS')

for(i in 1:n_iter){
  bv_mat_list[[i]]$iter <- i
  bv_mat_list[[i]]$mod_type <- 'matern'
}
bvs_mat <- bind_rows(bv_mat_list)
saveRDS(bvs_mat, 'output/bvs_n50_matern_v3.1.RDS')

# combine BVs from all model types
bvs_all <- rbind(bvs, bvs_over, bvs_mat)

#### PLOT BVS WITH UNCERTAINTY ####
ggplot(bvs_all, aes(x = factor(timeFrom), y = centroidVelocity, color = mod_type)) +
  geom_boxplot() + 
  coord_flip() +
  xlab('Years before present\n') +
  ylab('\nCentroid velocity (m/yr)') +
  labs(color = 'Model type') +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
ggsave('figures/BVs_all_model_types.png', width = 6, height = 6, unit = 'in')

# summarize across 50 iterations to get 97.5 and 2.5 quantiles
bvs <- readRDS('output/bvs_n50_latent_v3.1.RDS')
bvs_summ <- bvs %>% 
  dplyr::group_by(timeFrom) %>% 
  dplyr::summarize(BVcent = quantile(centroidVelocity, 0.5),
            BVcent0p025 = quantile(centroidVelocity, 0.025),
            BVcent0p975 = quantile(centroidVelocity, 0.975),
            BVN = quantile(nCentroidVelocity, 0.5),
            BVN0p025 = quantile(nCentroidVelocity, 0.025),
            BVN0p975 = quantile(nCentroidVelocity, 0.975),
            BVS = quantile(sCentroidVelocity, 0.5),
            BVS0p025 = quantile(sCentroidVelocity, 0.025),
            BVS0p975 = quantile(sCentroidVelocity, 0.975)
  )
time <- seq(210, 21000, by = 990)
time <- data.frame(timeFrom = time[-22], timeTo = time[-1])
time$timeFrame <- paste0(time$timeTo, ' - ', time$timeFrom)

bvs_summ$timeFrame <- rev(time$timeFrame)

saveRDS(pollen_summ,'pollen_bvs_quants_for_plotting.RDS')

pollen <- readRDS('pollen_bvs_quants_for_plotting.RDS')
pollen <- pollen %>% arrange(desc(median_time))  # arrange so time goes from LGM to modern (top to bottom)

plot_pollen <- ggplot(pollen, aes(x = median_time, y = BVcent)) +
  geom_boxplot(middle = pollen$BVcent,
               lower = pollen$BVcent0p025,
               upper = pollen$BVcent0p975) + 
  scale_y_continuous(limits = c(0,1500)) +
  ylab('\n ') +
  geom_text(x = 3, y = 1450, label = 'Pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 12))









#### CALCULATE BVS FOR NEWEST MODEL RUN (AUGUST 23, 2021) ####

# READ IN FRAXINUS DATA
# prediction grid, frax predictions for 50 random iterations
locs_grid <- readRDS('data/grid_3.1.RDS')
preds <- readRDS('output/frax_preds_n50_3.2_latent_overdispersed.RDS')
preds_no_mod <- readRDS('output/frax_preds_n50_3.2_latent_overdispersed_no_modern.RDS')

time <- c(-70, seq(705, by = 990, length.out = 22))
time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
time$time_mid <- (time$timeFrom + time$timeTo) / 2  # to get midpoint of time bins


# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
proj <- proj4string(stack)

# need to find times in glacier masks closest to median pollen times
stack_times <- rev(seq(0, by = 30, length.out = 701))
time$closest <- NA
for(i in 1:nrow(time)){
  time$closest[i] <- which.min(abs(stack_times - time$time_mid[i]))
}
stack_list <- unstack(stack)
stack_sub <- stack_list[time$closest]

for(i in 1:length(stack_sub)){
  stack_sub[[i]][stack_sub[[i]] == 1] <- NA  # convert areas completely covered by ice to NAs
}

stack_sub_no_mod <- stack_sub[-1]


#### DATA PREPARATION ####
# convert prediction array to list of rasterstacks
n_iter <- dim(preds)[1]
n_times <- dim(preds)[3]
n_times_no_mod <- dim(preds_no_mod)[3]

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

preds_no_mod_list <- list()
for(i in 1:n_iter){
  preds_no_mod_list[[i]] <- cbind(locs_grid, data.frame(preds_no_mod[i,,]))
  preds_no_mod_list[[i]] <- pivot_longer(preds_no_mod_list[[i]], cols = 3:ncol(preds_no_mod_list[[i]]),
                                       names_to = 'time', values_to = 'pred')
  preds_no_mod_list[[i]]$time <- as.integer(substr(preds_no_mod_list[[i]]$time, start = 2,
                                                 stop = nchar(preds_no_mod_list[[i]]$time)))
  preds_no_mod_list[[i]] <- split(preds_no_mod_list[[i]], f = preds_no_mod_list[[i]]$time)
  preds_no_mod_list[[i]] <- lapply(preds_no_mod_list[[i]], function(x) { x['time'] <- NULL; x })
}


# takes 1-2 minutes to run
rasterstack_list <- preds_list
for(i in 1:n_iter){
  for(j in 1:n_times){
    rasterstack_list[[i]][[j]] <- rasterFromXYZ(rasterstack_list[[i]][[j]])
    proj4string(rasterstack_list[[i]][[j]]) <- proj
    rasterstack_list[[i]][[j]] <- raster::resample(rasterstack_list[[i]][[j]], y = stack_sub[[j]])
    rasterstack_list[[i]][[j]] <- mask(rasterstack_list[[i]][[j]], mask = stack_sub[[j]])
  }
}

rasterstack_no_mod_list <- preds_no_mod_list
for(i in 1:n_iter){
  for(j in 1:n_times_no_mod){
    rasterstack_no_mod_list[[i]][[j]] <- rasterFromXYZ(rasterstack_no_mod_list[[i]][[j]])
    proj4string(rasterstack_no_mod_list[[i]][[j]]) <- proj
    rasterstack_no_mod_list[[i]][[j]] <- raster::resample(rasterstack_no_mod_list[[i]][[j]], y = stack_sub_no_mod[[j]])
    rasterstack_no_mod_list[[i]][[j]] <- mask(rasterstack_no_mod_list[[i]][[j]], mask = stack_sub_no_mod[[j]])
  }
}

# weight relative abundance by proportion of cell covered by land (for cells with partial glacial coverage)
stack_sub_revalue <- stack_sub
for(i in 1:n_times){
  values(stack_sub_revalue[[i]]) <- 1 - values(stack_sub_revalue[[i]])
}

stack_sub_revalue_no_mod <- stack_sub_revalue[-1]

rasterstack_revalue <- rasterstack_list
for(i in 1:n_iter){
  for(j in 1:n_times){
    rasterstack_revalue[[i]][[j]] <- overlay(rasterstack_list[[i]][[j]], stack_sub_revalue[[j]],
                                             fun = function(x, y){(x * y)})
  }
}

rasterstack_revalue_no_mod <- rasterstack_no_mod_list
for(i in 1:n_iter){
  for(j in 1:n_times_no_mod){
    rasterstack_revalue_no_mod[[i]][[j]] <- overlay(rasterstack_no_mod_list[[i]][[j]], stack_sub_revalue_no_mod[[j]],
                                             fun = function(x, y){(x * y)})
  }
}


# convert list of rasters into raster stack
for(i in 1:n_iter){
  rasterstack_revalue[[i]] <- raster::stack(rasterstack_revalue[[i]])
}

for(i in 1:n_iter){
  rasterstack_revalue_no_mod[[i]] <- raster::stack(rasterstack_revalue_no_mod[[i]])
}


##############################################
## CALCULATE AND PLOT BIOTIC VELOCITIES
##############################################
# iterate bioticVelocity function through each iteration (takes ~6 minutes)
# first need to reverse the order of rasters in each rasterstack to get the BV function to work
rasterstack_rev <- list()
for(i in 1:n_iter){
  rasterstack_rev[[i]] <- subset(rasterstack_revalue[[i]], order(n_times:1))
}
saveRDS(rasterstack_rev, 'output/rasterstacks_for_bv_calc_n50_latent_overdispersed_v3.2.RDS')
rasterstack_rev <- readRDS('output/rasterstacks_for_bv_calc_n50_latent_overdispersed_v3.2.RDS')

rasterstack_no_mod_rev <- list()
for(i in 1:n_iter){
  rasterstack_no_mod_rev[[i]] <- subset(rasterstack_revalue_no_mod[[i]], order(n_times_no_mod:1))
}
saveRDS(rasterstack_no_mod_rev, 'output/rasterstacks_for_bv_calc_n50_latent_overdispersed_no_mod_v3.2.RDS')
rasterstack_no_mod_rev <- readRDS('output/rasterstacks_for_bv_calc_n50_latent_overdispersed_no_mod_v3.2.RDS')


# might need to force negative values to be 0, otherwise BV function won't work
for(i in 1:n_iter){
  rasterstack_rev[[i]] <- calc(rasterstack_rev[[i]], fun = function(x) ifelse(x < 0, 0, x))
}

for(i in 1:n_iter){
  rasterstack_no_mod_rev[[i]] <- calc(rasterstack_no_mod_rev[[i]], fun = function(x) ifelse(x < 0, 0, x))
}

# calculate BVs
bv_list <- list()
for(i in 1:n_iter){
  bv_list[[i]] <- bioticVelocity(rasterstack_rev[[i]],
                                 times = rev(time$time_mid) * -1,
                                 onlyInSharedCells = TRUE)
}

bv_no_mod_list <- list()
for(i in 1:n_iter){
  bv_no_mod_list[[i]] <- bioticVelocity(rasterstack_no_mod_rev[[i]],
                                      times = rev(time$time_mid[-1]) * -1,
                                      onlyInSharedCells = TRUE)
}


# convert list of dataframes to single dataframe
for(i in 1:n_iter){
  bv_list[[i]]$type <- 'modern'
}
bvs <- bind_rows(bv_list)
saveRDS(bvs, 'output/bvs_n50_latent_overdispersed_v3.2.RDS')

for(i in 1:n_iter){
  bv_no_mod_list[[i]]$type <- 'no modern'
}
bvs_no_mod <- bind_rows(bv_no_mod_list)
saveRDS(bvs_no_mod, 'output/bvs_n50_latent_overdispersed_no_modern_v3.2.RDS')


# combine BVs from all model types
bvs_all <- rbind(bvs, bvs_no_mod)

#### PLOT BVS WITH UNCERTAINTY ####
ggplot(bvs_all, aes(x = factor(timeFrom), y = centroidVelocity, color = type)) +
  geom_boxplot() + 
  xlab('\nYears before present') +
  ylab('Centroid velocity (m/yr)\n') +
  labs(color = '') +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12))
ggsave('figures/BVs_modern_vs_no_modern_v3.2.png', 
       width = 9, height = 4, unit = 'in')





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





#### FOR POLLEN MODEL PAPER ####
#### CALCULATING BV FOR ALL TAXA ACROSS ALL (?) ITERATIONS
# here, using 100 (of the total 1000) iterations
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
version <- '3.1'


#### CALCULATE BVS FOR ALL ITERATIONS AND ALL TAXA

# READ IN FRAXINUS DATA
# prediction grid, frax predictions for 50 random iterations
locs_grid <- readRDS(paste0('data/grid_', version, '.RDS'))
preds <- readRDS('output/preds_all_taxa_n100_latent.RDS')

# generate time bins, get midpoint of bins
time <- seq(210, 21000, by = 990)
time <- data.frame(timeFrom = time[-1], timeTo = time[-22])
time$time_mid <- (time$timeFrom + time$timeTo) / 2  # to get midpoint of time bins

# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
proj <- proj4string(stack)

# need to find times in glacier masks closest to median pollen times
stack_times <- rev(seq(0, by = 30, length.out = 701))
time$closest <- NA
for(i in 1:nrow(time)){
  time$closest[i] <- which.min(abs(stack_times - time$time_mid[i]))
}
stack_list <- unstack(stack)
stack_sub <- stack_list[time$closest]

for(i in 1:length(stack_sub)){
  stack_sub[[i]][stack_sub[[i]] == 1] <- NA  # convert areas completely covered by ice to NAs
}


#### DATA PREPARATION ####
# convert prediction array to list of rasterstacks
n_iter <- dim(preds)[1]
n_taxa <- dim(preds)[3]
n_times <- dim(preds)[4]

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
# takes ~23 minutes to run for dimensions [100, 3951, 13, 21]
rasterstack_list <- preds_list
for(i in 1:n_taxa){
  for(j in 1:n_iter){
    for(k in 1:n_times){
      rasterstack_list[[i]][[j]][[k]] <- rasterFromXYZ(rasterstack_list[[i]][[j]][[k]])
      proj4string(rasterstack_list[[i]][[j]][[k]]) <- proj
      rasterstack_list[[i]][[j]][[k]] <- raster::resample(rasterstack_list[[i]][[j]][[k]], y = stack_sub[[k]])
      rasterstack_list[[i]][[j]][[k]] <- raster::mask(rasterstack_list[[i]][[j]][[k]], mask = stack_sub[[k]])
    }
  }
}

# convert list of rasters into rasterstack
for(i in 1:n_taxa){
  for(j in 1:n_iter){
    rasterstack_list[[i]][[j]] <- raster::stack(rasterstack_list[[i]][[j]])
  }
}


##############################################
## CALCULATE AND PLOT BIOTIC VELOCITIES
##############################################
# iterate bioticVelocity function through each iteration
# first need to reverse the order of rasters in each rasterstack to get the BV function to work
rasterstack_rev <- rasterstack_list
for(i in 1:n_taxa){
  for(j in 1:n_iter){
    rasterstack_rev[[i]][[j]] <- subset(rasterstack_list[[i]][[j]], order(21:1))
  }
}
saveRDS(rasterstack_list, 'output/preds_rasterstacks_for_bv_calc_all_taxa_n100_latent_v3.1.RDS')
rasterstack_rev <- readRDS('output/preds_rasterstacks_for_bv_calc_all_taxa_n100_latent_v3.1.RDS')

# need to force negative values to be 0, otherwise BV function won't work
# takes ~3 minutes to run
for(i in 1:n_taxa){
  for(j in 1:n_iter){
    rasterstack_rev[[i]][[j]] <- calc(rasterstack_rev[[i]][[j]], fun = function(x) ifelse(x < 0, 0, x))
  }
}

# calculate BVs
bv_list <- rep(list(list()), n_taxa)
for(i in 1:n_taxa){
  for(j in 1:n_iter){
    bv_list[[i]][[j]] <- bioticVelocity(rasterstack_rev[[i]][[j]],
                                        times = rev(time$time_mid) * -1,
                                        onlyInSharedCells = TRUE)
  }
}

# combine iterations into single dataframe for each taxon
for(i in 1:n_taxa){
  bv_list[[i]] <- bind_rows(bv_list[[i]])
}

# combine taxa dataframes into single dataframe
taxa_vec <- readRDS('data/taxa_3.2.RDS')
taxa_vec <- sort(taxa_vec)

for(i in 1:n_taxa){
  bv_list[[i]]$taxa <- taxa_vec[i]
}
bv_df <- bind_rows(bv_list)
saveRDS(bv_df, 'output/bvs_n100_all_taxa_latent_v3.1.RDS')


#### PLOT BVS WITH UNCERTAINTY ####
time_bv <- data.frame(timeFrom = time$time_mid[-1], 
                      timeTo = time$time_mid[-21])
time_bv$timeFrame <- paste0(time_bv$timeFrom, ' - ', time_bv$timeTo)

ggplot(bv_df, aes(x = factor(timeFrom), y = centroidVelocity)) +
  geom_boxplot() + 
  scale_x_discrete(labels = rev(time_bv$timeFrame)) +
  xlab('\nYears before present') +
  ylab('Centroid velocity (m/yr)\n') +
  facet_wrap(~ taxa) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave('figures/BVs_all_taxa_n100_latent.png', width = 12, height = 6, unit = 'in')

ggplot(bv_df, aes(x = factor(timeFrom), y = sCentroidVelocity)) +
  geom_boxplot() + 
  scale_x_discrete(labels = rev(time_bv$timeFrame)) +
  xlab('\nYears before present') +
  ylab('South centroid velocity (m/yr)\n') +
  facet_wrap(~ taxa) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave('figures/BVs_south_all_taxa_n100_latent.png', width = 12, height = 6, unit = 'in')

ggplot(bv_df, aes(x = factor(timeFrom), y = nCentroidVelocity)) +
  geom_boxplot() + 
  scale_x_discrete(labels = rev(time_bv$timeFrame)) +
  xlab('\nYears before present') +
  ylab('North centroid velocity (m/yr)\n') +
  facet_wrap(~ taxa) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave('figures/BVs_north_all_taxa_n100_latent.png', width = 12, height = 6, unit = 'in')

# DURING WHAT TIME PERIODS DO WE SEE THE HIGHEST BVS FOR EACH TAXON?
summ <- bv_df %>% 
  dplyr::group_by(taxa, timeFrom) %>% 
  dplyr::summarise(median_BV = median(centroidVelocity),
            median_nBV = median(nCentroidVelocity),
            median_sBV = median(sCentroidVelocity))

summ_max <- summ %>%
  arrange(desc(median_BV)) %>% 
  group_by(taxa) %>%
  slice(1:3)  # examine highest 3 BVs for each taxon

# How many taxa have max BVs at each time period? 
table(summ_max$timeFrom)
# 20,505: 4
# 19,515: 4
# 18,525: 4
# 17,535: 3
# 16,545: 5
# 1695: 12




#### PLOTTING N/S VELOCITIES AS SCATTERPLOTS
# can see overall trend (e.g., whether things are expanding northward overall)
bv <- readRDS('output/bvs_n100_all_taxa_latent_v3.1.RDS')

# plot overall N/S movement of all taxa; take mean of BVs across time for each taxon
bv_summ <- bv %>% 
  dplyr::select(taxa, nCentroidVelocity, sCentroidVelocity) %>% 
  group_by(taxa) %>% 
  summarize(across(c(nCentroidVelocity, sCentroidVelocity), mean))

ggplot(bv_summ, aes(x = sCentroidVelocity, y = nCentroidVelocity, color = taxa)) +
  geom_point() + 
  scale_x_continuous(limits = c(0, 310)) +
  scale_y_continuous(limits = c(0, 310)) +
  geom_abline(slope = 1, intercept = 0, lty = 'dashed') +
  theme_classic()

# plot Fraxinus movement for different time periods
frax_summ <- bv %>% 
  filter(taxa == 'Fraxinus') %>% 
  dplyr::select(timeFrom, nCentroidVelocity, sCentroidVelocity) %>% 
  group_by(timeFrom) %>% 
  summarize(across(c(nCentroidVelocity, sCentroidVelocity), mean))
frax_summ$timeFrom <- frax_summ$timeFrom * -1

ggplot(frax_summ, aes(x = sCentroidVelocity, y = nCentroidVelocity, color = timeFrom)) +
  geom_point(size = 4) + 
  scale_x_continuous(limits = c(0, 455)) +
  scale_y_continuous(limits = c(0, 455)) +
  geom_abline(slope = 1, intercept = 0, lty = 'dashed') +
  labs(color = 'Years before\npresent\n') +
  xlab('\nSouth centroid velocity (m/yr)') +
  ylab('North centroid velocity (m/yr)\n') +
  theme_classic() + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

## COULD MAKE ABOVE PLOTS BETTER BY ADDING ERROR BARS
## COULD ADD BV SUMMARIES FROM OTHER 4 METHODS ONTO ONE PLOT TO COMPARE