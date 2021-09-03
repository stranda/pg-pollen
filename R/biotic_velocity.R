setwd('C:/Users/abrow/Documents/pg-pollen')
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
require(gridExtra)
require(ggrepel)

################################################
## CALCULATE AND PLOT BIOTIC VELOCITIES
################################################

# read prediction output (50 iterations randomly sampled from all 1000 model iterations)
locs_grid <- readRDS('data/grid_3.1.RDS')
preds <- readRDS('output/frax_preds_n50_3.2_latent_overdispersed.RDS')

# specify time bins; use median of time bins for subsequent work
time <- c(-70, seq(705, by = 990, length.out = 22))
time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
time$time_mid <- (time$timeFrom + time$timeTo) / 2

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

# convert areas completely covered by ice to NAs
for(i in 1:length(stack_sub)){
  stack_sub[[i]][stack_sub[[i]] == 1] <- NA
}

# convert prediction array to list of list of dataframes
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

# convert list of list of dataframes to list of rasters
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

# weight relative abundance by proportion of cell covered by land (for cells with partial glacial coverage)
# partial glacial coverage is indicated as proportion of cell covered by glacier
# first, convert this value to proportion of cell covered by land
stack_sub_revalue <- stack_sub
for(i in 1:n_times){
  values(stack_sub_revalue[[i]]) <- 1 - values(stack_sub_revalue[[i]])
}

# multiply relative abundance predictions by proportion of cell covered by land
rasterstack_revalue <- rasterstack_list
for(i in 1:n_iter){
  for(j in 1:n_times){
    rasterstack_revalue[[i]][[j]] <- overlay(rasterstack_list[[i]][[j]], stack_sub_revalue[[j]],
                                             fun = function(x, y){(x * y)})
  }
}

# convert list of rasters into rasterstack
for(i in 1:n_iter){
  rasterstack_revalue[[i]] <- raster::stack(rasterstack_revalue[[i]])
}

# need to reverse the order of rasters in each rasterstack to get the BV function to work
rasterstack_rev <- list()
for(i in 1:n_iter){
  rasterstack_rev[[i]] <- subset(rasterstack_revalue[[i]], order(n_times:1))
}

# might need to force negative values to be 0, otherwise BV function won't work
for(i in 1:n_iter){
  rasterstack_rev[[i]] <- calc(rasterstack_rev[[i]], fun = function(x) ifelse(x < 0, 0, x))
}

#### CALCULATE BVs
# iterate bioticVelocity function through each iteration (takes ~6 minutes)
bv_list <- list()
for(i in 1:n_iter){
  bv_list[[i]] <- bioticVelocity(rasterstack_rev[[i]],
                                 times = rev(time$time_mid) * -1,
                                 onlyInSharedCells = TRUE)
}

# convert list of dataframes to single dataframe
bvs <- bind_rows(bv_list)
saveRDS(bvs, 'output/bvs_n50_latent_overdispersed_v3.2.RDS')

#### PLOT BVS WITH UNCERTAINTY
time$time_label <- NA
for(i in 2:nrow(time)){
  time$time_label[i] <- paste0(time$time_mid[i], ' - ', time$time_mid[i-1])
}

ggplot(bvs, aes(x = factor(timeFrom), y = centroidVelocity)) +
  geom_boxplot() + 
  scale_x_discrete(labels = rev(time$time_label[-1])) +
  xlab('\nYears before present') +
  ylab('Centroid velocity (m/yr)\n') +
  labs(color = '') +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12))
ggsave('figures/BVs.png', width = 9, height = 4, unit = 'in')




################################################
## CALCULATE AND PLOT REFUGE LOCATIONS
################################################

# read in 'simulated' dataframe
sim <- raster('data/map-data/study_region_resampled_to_genetic_demographic_simulation_resolution.tif')
proj <- proj4string(sim)

# read in fraxinus predictions for each of the 50 iterations
# subset to only LGM predictions
grid <- readRDS('data/grid_3.0.RDS')
frax <- readRDS('output/frax_preds_n50_3.2_latent_overdispersed.RDS')
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
                                                       sim = sim, threshold = 0.025)
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
saveRDS(refuge_for_allan, 'output/refuge_rasters_n50.RDS')


# PLOT 50 REFUGIA AS ANIMATED GIF
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
    scale_fill_manual(values = my_colors, labels = breaks_vec) +
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
ggsave('figures/frax_refugia_latent_n50.pdf', plots, 
       height = 15, width = 8.5, units = 'in')




#################################################
# PLOT CENTROIDS OF FRAX RANGE OVER TIME ON MAP
#################################################
# require(ggrepel)
# use biotic velocity dataframe
bvs <- readRDS('output/bvs_n50_latent_overdispersed_v3.2.RDS')

# take average BV across iterations
bvs_summ <- bvs %>% 
  dplyr::group_by(timeFrom) %>% 
  dplyr::summarize(x = mean(centroidLong), y = mean(centroidLat))

ggplot(data = bvs_summ) +
  geom_point(aes(x = x, y = y, fill = timeFrom),
             size = 4, alpha = 0.9, pch = 21) +
  geom_label_repel(aes(x = x, y = y, 
                       label = ifelse(timeFrom > -6000, as.character(timeFrom),'')),
                   box.padding   = 0.5,
                   point.padding = 0.5,
                   segment.color = 'red',
                   max.overlaps = 30) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(-700000,1500000)) +
  scale_x_continuous(limits = c(-500000,2500000)) +
  labs(title = 'Fraxinus range centroids over time') +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16)) +
  coord_equal()




##########################################################
# PLOT BIOTIC VELOCITIES FOR ALL 5 METHODS
##########################################################
# setwd('C:/Users/abrow/Documents/green_ash')
# require(gridExtra)
# read bv data
abc_enm <- read.csv('bvs_enm-abc.csv', stringsAsFactors = FALSE)
abc <- read.csv('bvs_naive_abc.csv', stringsAsFactors = FALSE)
abc_pollen <- read.csv('bvs_pollen-abc.csv', stringsAsFactors = FALSE)
load('enm_biotic_velocities.rda')
enm <- velocities
rm(velocities)
pollen <- readRDS('bvs_n50_latent_overdispersed_v3.2.RDS')

# DATA PREPARATION - POLLEN
# summarize across 50 iterations to get 97.5 and 2.5 quantiles
pollen_summ <- pollen %>%
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
pollen_summ$timeFrom <- pollen_summ$timeFrom * -1
saveRDS(pollen_summ,'pollen_bvs_quants_for_plotting.RDS')

plot_pollen <- ggplot(pollen_summ, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVcent)) +
  geom_boxplot(middle = pollen_summ$BVcent,
               lower = pollen_summ$BVcent0p025,
               upper = pollen_summ$BVcent0p975) + 
  scale_y_continuous(limits = c(0,1500)) +
  ylab('\n ') +
  geom_text(x = 3, y = 1450, label = 'Pollen', size = 6) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 12))

# DATA PREPARATION - ENM
# enm dataframe contains BVs for each model type; filter out BVs you don't need
enm_sub <- enm %>% filter(timeSpan == 990, 
                          onlyInSharedCells == TRUE, 
                          onlyInContinuouslyExposedLand == FALSE)

# summarize across 24 models to get 97.5 and 2.5 quantiles
enm_summ <- enm_sub %>% 
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
enm_summ$timeFrom <- enm_summ$timeFrom * -1
saveRDS(enm_summ,'enm_bvs_quants_for_plotting.RDS')

plot_enm <- ggplot(enm_summ, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVcent)) + 
  geom_boxplot(middle = enm_summ$BVcent,
               lower = enm_summ$BVcent0p025,
               upper = enm_summ$BVcent0p975) + 
  scale_y_continuous(limits = c(0, 400)) +
  ylab('\n ') +
  geom_text(x = 3, y = 375, label = 'ENM', size = 6) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# DATA PREPARATION - NAIVE ABC
abc$timeFrom <- as.numeric(gsub('-.*', '', abc$time))

plot_abc <- ggplot(abc, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVcent))+ 
  geom_boxplot(middle = abc$BVcent, 
               lower = abc$BVcent0p025,
               upper = abc$BVcent0p975)+
  scale_y_continuous(limits = c(0, 400)) +
  ylab('Centroid BV (m/yr)\n') +
  geom_text(x = 3, y = 375, label = 'Naive ABC', size = 6) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# DATA PREPARATION - ENM-ABC
abc_enm$timeFrom <- as.numeric(gsub('-.*', '', abc_enm$time))

plot_abc_enm <- ggplot(abc_enm, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVcent))+ 
  geom_boxplot(middle = abc_enm$BVcent, 
               lower = abc_enm$BVcent0p025,
               upper = abc_enm$BVcent0p975)+
  # scale_x_discrete(limits = c("21000-20010","20010-19020","19020-18030","18030-17040",
  #                             "17040-16050","16050-15060","15060-14070","14070-13080",
  #                             "13080-12090","12090-11100","11100-10110","10110-9120",
  #                             "9120-8130","8130-7140","7140-6150","6150-5160",
  #                             "5160-4170","4170-3180","3180-2190","2190-1200","1200-210"))+
  scale_y_continuous(limits = c(0, 400)) +
  ylab('\n ') +
  xlab('Years before present') +
  geom_text(x = 3, y = 375, label = 'ABC-ENM', size = 6) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

# DATA PREPARATION - POLLEN-ABC
abc_pollen$timeFrom <- as.numeric(gsub('-.*', '', abc_pollen$time))

plot_abc_pollen <- ggplot(abc_pollen, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVcent))+ 
  geom_boxplot(middle = abc_pollen$BVcent, 
               lower = abc_pollen$BVcent0p025,
               upper = abc_pollen$BVcent0p975)+
  scale_y_continuous(limits = c(0, 400)) +
  ylab('\n ') +
  geom_text(x = 3, y = 375, label = 'ABC-pollen', size = 6) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

g <- arrangeGrob(plot_pollen, plot_enm, plot_abc, plot_abc_pollen, plot_abc_enm, 
                 ncol = 1, heights = c(1,1,1,1,1))
ggsave(file = 'figures/BVs_all_methods.png', plot = g, height = 10, width = 7, units = 'in')




# PLOT NORTH CENTROID BVS FOR ALL 5 METHODS
plot_pollen <- ggplot(pollen_summ, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVN)) +
  geom_boxplot(middle = pollen_summ$BVN,
               lower = pollen_summ$BVN0p025,
               upper = pollen_summ$BVN0p975) + 
  scale_y_continuous(limits = c(-200,1700)) +
  ylab('\n\n ') +
  geom_text(x = 3, y = 1600, label = 'Pollen', size = 6) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12))

plot_enm <- ggplot(enm_summ, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVN)) + 
  geom_boxplot(middle = enm_summ$BVN,
               lower = enm_summ$BVN0p025,
               upper = enm_summ$BVN0p975) + 
  scale_y_continuous(limits = c(-200, 400)) +
  ylab('\n\n ') +
  geom_text(x = 3, y = 375, label = 'ENM', size = 6) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc <- ggplot(abc, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVN))+ 
  geom_boxplot(middle = abc$BVN, 
               lower = abc$BVN0p025,
               upper = abc$BVN0p975)+
  scale_y_continuous(limits = c(-200, 400)) +
  ylab('North-centroid\nBV (m/yr)\n') +
  geom_text(x = 3, y = 375, label = 'Naive ABC', size = 6) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc_pollen <- ggplot(abc_pollen, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVN))+ 
  geom_boxplot(middle = abc_pollen$BVN, 
               lower = abc_pollen$BVN0p025,
               upper = abc_pollen$BVN0p975)+
  scale_y_continuous(limits = c(-200, 400)) +
  ylab('\n\n ') +
  geom_text(x = 3, y = 375, label = 'ABC-pollen', size = 6) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc_enm <- ggplot(abc_enm, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVcent))+ 
  geom_boxplot(middle = abc_enm$BVcent, 
               lower = abc_enm$BVcent0p025,
               upper = abc_enm$BVcent0p975)+
  # scale_x_discrete(limits = c("21000-20010","20010-19020","19020-18030","18030-17040",
  #                             "17040-16050","16050-15060","15060-14070","14070-13080",
  #                             "13080-12090","12090-11100","11100-10110","10110-9120",
  #                             "9120-8130","8130-7140","7140-6150","6150-5160",
  #                             "5160-4170","4170-3180","3180-2190","2190-1200","1200-210"))+
  scale_y_continuous(limits = c(-200, 400)) +
  ylab('\n\n ') +
  xlab('Years before present') +
  geom_text(x = 3, y = 375, label = 'ABC-ENM', size = 6) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

g <- arrangeGrob(plot_pollen, plot_enm, plot_abc, plot_abc_pollen, plot_abc_enm, 
                 ncol = 1, heights = c(1,1,1,1,1))
ggsave(file = 'figures/BVNs_all_methods.png', plot = g, height = 10, width = 7, units = 'in')




# PLOT SOUTH CENTROID BVS FOR ALL 5 METHODS
plot_pollen <- ggplot(pollen_summ, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVS)) +
  geom_boxplot(middle = pollen_summ$BVS,
               lower = pollen_summ$BVS0p025,
               upper = pollen_summ$BVS0p975) + 
  scale_y_continuous(limits = c(-500,1000)) +
  ylab('\n\n ') +
  geom_text(x = 19, y = -425, label = 'Pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12))

plot_enm <- ggplot(enm_summ, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVS)) + 
  geom_boxplot(middle = enm_summ$BVS,
               lower = enm_summ$BVS0p025,
               upper = enm_summ$BVS0p975) + 
  scale_y_continuous(limits = c(-500, 250)) +
  ylab('\n\n ') +
  geom_text(x = 19, y = -450, label = 'ENM', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc <- ggplot(abc, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVS))+ 
  geom_boxplot(middle = abc$BVS, 
               lower = abc$BVS0p025,
               upper = abc$BVS0p975)+
  scale_y_continuous(limits = c(-500, 250)) +
  ylab('South-centroid\nBV (m/yr)\n') +
  geom_text(x = 19, y = -450, label = 'Naive ABC', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc_pollen <- ggplot(abc_pollen, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVS))+ 
  geom_boxplot(middle = abc_pollen$BVS, 
               lower = abc_pollen$BVS0p025,
               upper = abc_pollen$BVS0p975)+
  scale_y_continuous(limits = c(-500, 250)) +
  ylab('\n\n ') +
  geom_text(x = 19, y = -450, label = 'ABC-pollen', size = 6) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))

plot_abc_enm <- ggplot(abc_enm, aes(x = reorder(factor(timeFrom), -timeFrom), y = BVcent))+ 
  geom_boxplot(middle = abc_enm$BVcent, 
               lower = abc_enm$BVcent0p025,
               upper = abc_enm$BVcent0p975)+
  # scale_x_discrete(limits = c("21000-20010","20010-19020","19020-18030","18030-17040",
  #                             "17040-16050","16050-15060","15060-14070","14070-13080",
  #                             "13080-12090","12090-11100","11100-10110","10110-9120",
  #                             "9120-8130","8130-7140","7140-6150","6150-5160",
  #                             "5160-4170","4170-3180","3180-2190","2190-1200","1200-210"))+
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
                 ncol = 1, heights = c(1,1,1,1,1))
ggsave(file = 'figures/BVSs_all_methods_V2.png', plot = g, height = 10, width = 7, units = 'in')
