
require(tidyr)
require(dplyr)
require(ggplot2)
require(rasterVis)
require(fields)
require(rgdal)
require(raster)
require(enmSdm)

setwd('C:/Users/abrow/Documents/pg-pollen')
version <- 'Dec15'

################################################
# DATA PREPARATION
################################################

#### READ POLLEN DATA AND PREDICTIONS ####
# read raw data
dat <- readRDS(paste0('output/polya-gamma-dat_', version, '.RDS'))
taxa <- data.frame(taxon_num = 1:14, taxon = c(dat$taxa.keep, 'Other'))

# read prediction grid locations
locs_grid <- readRDS(paste0('data/grid_', version, '.RDS'))

# read mean and sds preds (or median and IQR)
means <- readRDS(paste0('output/preds_', version, '_all_taxa_means.RDS'))
sds <- readRDS(paste0('output/preds_', version, '_all_taxa_sds.RDS'))
# medians <- readRDS(paste0('output/preds_', version, '_all_taxa_medians.RDS'))
# iqrs <- readRDS(paste0('output/preds_', version, '_all_taxa_iqr.RDS'))

n_times <- dim(means)[3]
n_taxa <- dim(means)[2]


#### READ MAP DATA ####
# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
names(stack) <- 1:701
sub <- seq(1, 701, by = 33)
stack_sub <- subset(stack, subset = paste0('X', sub))  # only want mask every 990 years
stack_sub[stack_sub > 0.6] <- NA
proj <- proj4string(stack)

na_shp <- readOGR("data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("data/map-data/Great_Lakes.shp", "Great_Lakes")
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


#### PREPARE RAW DATAFRAMES ####
# create list of dataframes, grouped by time period
dat_list <- list()
for(i in 1:n_times){
  dat_list[[i]] <- cbind(dat$locs, data.frame(dat$y[,,i]))
  dat_list[[i]] <- pivot_longer(dat_list[[i]], cols = 3:ncol(dat_list[[i]]),
                                 names_to = 'taxon_num', values_to = 'count')
  dat_list[[i]]$taxon_num <- as.integer(substr(dat_list[[i]]$taxon_num, start = 2, 
                                                stop = nchar(dat_list[[i]]$taxon_num)))
  dat_list[[i]] <- left_join(dat_list[[i]], taxa, by = 'taxon_num')
  dat_list[[i]] <- dat_list[[i]] %>% dplyr::select(-taxon_num)
}

# create list of taxon-specific dataframes within time-specific list items
dat_list <- lapply(dat_list, function(x) split(x, f = x$taxon))
for(i in 1:n_times){
  dat_list[[i]] <- lapply(dat_list[[i]], 
                           function(x) x[!(names(x) %in% c('taxon'))])
}

#### PREPARE MEAN PREDICTION DATAFRAMES ####
# create list of dataframes, grouped by time period
mean_list <- list()
for(i in 1:n_times){
  mean_list[[i]] <- cbind(locs_grid, data.frame(means[,,i]))
  mean_list[[i]] <- pivot_longer(mean_list[[i]], cols = 3:ncol(mean_list[[i]]),
                                 names_to = 'taxon_num', values_to = 'mean')
  mean_list[[i]]$taxon_num <- as.integer(substr(mean_list[[i]]$taxon_num, start = 2, 
                                           stop = nchar(mean_list[[i]]$taxon_num)))
  mean_list[[i]] <- left_join(mean_list[[i]], taxa, by = 'taxon_num')
  mean_list[[i]] <- mean_list[[i]] %>% dplyr::select(-taxon_num)
}

# create list of taxon-specific dataframes within time-specific list items
mean_list <- lapply(mean_list, function(x) split(x, f = x$taxon))
for(i in 1:n_times){
    mean_list[[i]] <- lapply(mean_list[[i]], 
                               function(x) x[!(names(x) %in% c('taxon'))])
  }

#### PREPARE SD PREDICTION DATAFRAMES ####
# create list of dataframes, grouped by time period
sd_list <- list()
for(i in 1:n_times){
  sd_list[[i]] <- cbind(locs_grid, data.frame(sds[,,i]))
  sd_list[[i]] <- pivot_longer(sd_list[[i]], cols = 3:ncol(sd_list[[i]]),
                                 names_to = 'taxon_num', values_to = 'sd')
  sd_list[[i]]$taxon_num <- as.integer(substr(sd_list[[i]]$taxon_num, start = 2, 
                                                stop = nchar(sd_list[[i]]$taxon_num)))
  sd_list[[i]] <- left_join(sd_list[[i]], taxa, by = 'taxon_num')
  sd_list[[i]] <- sd_list[[i]] %>% dplyr::select(-taxon_num)
}

# create list of taxon-specific dataframes within time-specific list items
sd_list <- lapply(sd_list, function(x) split(x, f = x$taxon))
for(i in 1:n_times){
  sd_list[[i]] <- lapply(sd_list[[i]], 
                           function(x) x[!(names(x) %in% c('taxon'))])
}


#### PREPARE MEAN PREDICTION DATA FOR PLOTTING ####
# convert rasterstack masks into list of rasters, reverse order of items to correspond with time
mask_list <- unstack(stack_sub)
mask_names <- NA
for(i in 1:n_times){
  mask_names[i] <- names(mask_list[[i]])
  mask_names[i] <- as.integer(substr(mask_names[i], start = 2, stop = nchar(mask_names[i]))) * 30
}
mask_names <- rev(mask_names)
names(mask_list) <- mask_names
mask_list <- rev(mask_list)

# convert predictions to rasters so you can mask out glaciers/ocean/lakes
raster_list <- mean_list
for(i in 1:n_times){
  for(j in 1:n_taxa){
    raster_list[[i]][[j]] <- rasterFromXYZ(mean_list[[i]][[j]])
    proj4string(raster_list[[i]][[j]]) <- proj
    raster_list[[i]][[j]] <- resample(raster_list[[i]][[j]], y = mask_list[[i]])
    raster_list[[i]][[j]] <- mask(raster_list[[i]][[j]], mask = mask_list[[i]])
  }
}

# convert masked rasters to dataframes for plotting
mask_df_list <- raster_list
for(i in 1:n_times){
  for(j in 1:n_taxa){
    mask_df_list[[i]][[j]] <- as.data.frame(raster_list[[i]][[j]], xy = TRUE)
  }
}

#### PREPARE SD PREDICTION DATA FOR PLOTTING ####
# convert predictions to rasters so you can mask out glaciers/ocean/lakes
raster_list <- sd_list
for(i in 1:n_times){
  for(j in 1:n_taxa){
    raster_list[[i]][[j]] <- rasterFromXYZ(sd_list[[i]][[j]])
    proj4string(raster_list[[i]][[j]]) <- proj
    raster_list[[i]][[j]] <- resample(raster_list[[i]][[j]], y = mask_list[[i]])
    raster_list[[i]][[j]] <- mask(raster_list[[i]][[j]], mask = mask_list[[i]])
  }
}

# convert masked rasters to dataframes for plotting
sd_mask_df_list <- raster_list
for(i in 1:n_times){
  for(j in 1:n_taxa){
    sd_mask_df_list[[i]][[j]] <- as.data.frame(raster_list[[i]][[j]], xy = TRUE)
  }
}

#### COMBINE MEAN AND SD PREDICTIONS, CONVERT BACK TO DATAFRAME FORMAT ####
preds_list <- mask_df_list
for(i in 1:n_times){
  for(j in 1:n_taxa){
    preds_list[[i]][[j]] <- full_join(mask_df_list[[i]][[j]], sd_mask_df_list[[i]][[j]],
                                      by = c('x','y'))
  }
}
# saveRDS(preds_list, 'output/preds_list_mean_sd_all_taxa_Dec15.RDS')

for(i in 1:n_times){
  for(j in 1:n_taxa){
    preds_list[[i]][[j]]$taxon <- names(preds_list[[i]][j])
  }
}
preds_list <- lapply(preds_list, function(x) bind_rows(x))

for(i in 1:n_times){
    preds_list[[i]]$time <- i
}
preds_df <- bind_rows(preds_list)

preds_df <- pivot_longer(preds_df, cols = c('mean','sd'),
                             names_to = 'type', values_to = 'preds')
# saveRDS(preds_df, 'output/preds_long_df_mean_sd_all_taxa_Dec15.RDS')


##########################################################
# PLOT PREDICTION MEANS AND SD BY TAXON X TIME
##########################################################

#### CREATE SPECIES-SPECIFIC QUANTILE BREAKS FOR BETTER VISUALIZATION ####
n_times <- nrow(preds_df %>% dplyr::select(time) %>% distinct)
n_taxa <- nrow(preds_df %>% dplyr::select(taxon) %>% distinct)

preds_list <- split(preds_df, f = preds_df$taxon)
for(i in 1:n_taxa){
  preds_list2 <- preds_list[[i]]
  preds_list2 <- preds_list2[preds_list2$type == 'mean', ]
  quants <- quantile(preds_list2$preds, probs = seq(0, 1, 0.1), na.rm = TRUE)
  quants <- c(0, quants)
  breaks <- cut(preds_list[[i]]$preds, breaks = quants, include.lowest = TRUE, labels = FALSE)
  preds_list[[i]]$breaks <- breaks
}
preds_plot <- bind_rows(preds_list)

# SELECT WHICH TAXA/TIME PERIODS YOU WANT TO PLOT (IF DESIRED)
preds_plot_sub <- preds_plot[preds_plot$taxon %in% c('Fraxinus','Pinus') & 
                               preds_plot$time %in% c(1,2), ]
n_times <- nrow(preds_plot_sub %>% dplyr::select(time) %>% distinct)

pdf('figures/preds_mean_sd_with_time.pdf')
for (i in 1:n_times){
  time_sub <- preds_plot_sub[which(preds_plot_sub$time == i), ]
  
  p <- ggplot() +
    geom_tile(data = time_sub, aes(x = x, y = y, fill = as.factor(breaks))) +
    scale_fill_manual(values = tim.colors(11), name='Preds', drop=FALSE) +
    facet_wrap(taxon ~ type) +
    geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
    geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(limits = xlim) +
    ggtitle(paste0(i-1, ',000 ybp')) +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          # legend.title = element_text(size = 12),
          # legend.text = element_text(size = 10),
          plot.title = element_text(size = 12)) +
    coord_equal()
  print(p)
}
dev.off()
