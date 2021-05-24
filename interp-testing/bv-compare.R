library(sp)
library(enmSdm)
library(ggplot2)
require(raster)

BVs = FALSE

dat_ash = readRDS('output/preds_paleo_mean_frax.RDS')
coords = readRDS('data/grid_3.0.RDS')/1e3

coords_all = coords
sp::coordinates(coords_all) <- ~x + y
proj4string(coords_all) <- CRS("+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# years = rep(1950, 4)
years = seq(from=150, by=500, length.out=42)#c(0, 50, 100, 150)
years_pred = seq(from=150, to=21000, by=30)#

proj <- '+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs'

#########################################################################################################################################
## 
#########################################################################################################################################

occupied = function(dat_rast, presence=0.03){
  n_times = dim(dat_rast)[3]
  n_cells = rep(NA, n_times)
  for (i in 1:n_times){
    print(i)
    vals = dat_rast[[i]][]
    n_cells[i] = sum(vals>=presence, na.rm=TRUE)
  }
  
  return(n_cells)

}

value_sum = function(dat_rast){
  n_times = dim(dat_rast)[3]
  value = rep(NA, n_times)
  for (i in 1:n_times){
    print(i)
    vals = dat_rast[[i]][]
    value[i] = sum(vals, na.rm=TRUE)
  }
  
  return(value)
  
}

#########################################################################################################################################


dat = data.frame(n_cells= n_cells, time=times)


ash_kriged_df = readRDS('data/ash_kriged_df.RDS')
kriged_out = ash_kriged_df[,c('x', 'y', 'timeIndex', 'var1.pred')]
colnames(kriged_out) = c('x', 'y', 'ybp', 'value')
kriged_out[,c('x', 'y')] = kriged_out[,c('x', 'y')]*1e3
years_pred = seq(from=150, to=21000, by=30)#
kriged_out[,'ybp'] = years_pred[kriged_out[,'ybp']]
kriged_out[,c('value')] = kriged_out[,c('value')]/100

saveRDS(kriged_out, 'data/krige_out_df_paleo.RDS')

times = unique(ash_kriged_df$timeIndex)
n_times = length(times)
modern_list <- list()
for(i in 1:n_times){
  dat_sub = ash_kriged_df[which(ash_kriged_df$timeIndex == times[i]),]
  modern_list[[i]] <- as.data.frame(dat_sub[,'var1.pred']/100)
  names(modern_list[[i]]) <- 'pred'
  modern_list[[i]] <- cbind(dat_sub[,c('x', 'y')]*1e3, modern_list[[i]])
  modern_list[[i]] <- rasterFromXYZ(modern_list[[i]], crs = CRS(proj))
}

all_list = modern_list

# MASK OUT SEAS, LAKES, ICE
# rasterstack masks are arranged LGM (item 1) to modern (item 701)
# I think modern = 2000 A.D.
if (!file.exists('data/mask_list.RDS')){
stack <- stack('../data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
stack_list <- unstack(stack)
stack_list <- rev(stack_list)

stack_list_c <- stack_list
for(i in 1:length(stack_list_c)){
  stack_list_c[[i]][stack_list_c[[i]] == 1] <- NA
}
saveRDS(stack_list_c, 'data/mask_list.RDS')

} else {
  stack_list_c = readRDS('data/mask_list.RDS')
}

# Specify pollen times (46 total; 4 modern and 42 paleo)
times <- years_pred#seq(0, 150, by = 30)

# create vector of ice mask times
adam_times <- (seq(0, 700, 1) * 30)
#adam_times = adam_times[6:length(adam_times)]

# # interpolate pollen rasters to temporal resolution of 30yrs
# all_stack <- stack(all_list)
# interp <- interpolateRasters(all_stack, interpFrom = times, interpTo = adam_times)
# interp_list <- unstack(interp)

n_times <- length(times)#length(stack_list_c)
resamp_list <- all_list
for(i in 1:n_times){
  idx_time = which(adam_times == times[i])
  resamp_list[[i]] <- raster::resample(resamp_list[[i]], y = stack_list_c[[idx_time]])
  # mask_list[[i]] <- mask(mask_list[[i]], mask = stack_list_c[[idx_time]])
}
resamp_rast <- stack(resamp_list)

# # GET BIOTIC VELOCITY FROM INTERPOLATED RASTERS
# interp_list <- unstack(interp)
# interp_list <- rev(interp_list)  # oldest needs to be on top

#bv_times <- adam_times * -1
#bv_times <- bv_times[15:697]  # just for now, to get the function to work
#bvs <- bioticVelocity(masked_rast, times = bv_times, onlyInSharedCells = TRUE)
if (BVs){
bvs_krige <- bioticVelocity(resamp_rast, times = times, onlyInSharedCells = TRUE)
saveRDS(bvs_krige, 'data/BVs_krige_paleo.RDS')
}

occ_krige = occupied(resamp_rast)
value_sum_krige = value_sum(resamp_rast)
range_summary = data.frame(n_cells = occ_krige, value_sum = value_sum_krige, time = times, type = rep('krige', length(times)))



# mask out ice/water
n_times <- length(times)#length(stack_list_c)
mask_list <- resamp_list
for(i in 1:n_times){
  idx_time = which(adam_times == times[i])
  # mask_list[[i]] <- raster::resample(mask_list[[i]], y = stack_list_c[[idx_time]])
  mask_list[[i]] <- mask(mask_list[[i]], mask = stack_list_c[[idx_time]])
}
mask_rast <- stack(mask_list)

# # GET BIOTIC VELOCITY FROM INTERPOLATED RASTERS
# interp_list <- unstack(interp)
# interp_list <- rev(interp_list)  # oldest needs to be on top

#bv_times <- adam_times * -1
#bv_times <- bv_times[15:697]  # just for now, to get the function to work
#bvs <- bioticVelocity(masked_rast, times = bv_times, onlyInSharedCells = TRUE)
if (BVs){
bvs_krige_mask <- bioticVelocity(mask_rast, times = times, onlyInSharedCells = TRUE)
saveRDS(bvs_krige_mask, 'data/BVs_krige_mask_paleo.RDS')
}

occ_krige_mask = occupied(mask_rast)
value_sum_krige_mask = value_sum(mask_rast)

range_summary = rbind(range_summary,
                   data.frame(n_cells = occ_krige_mask, 
                              value_sum= value_sum_krige_mask, 
                              time = times, type = rep('krige-mask', length(times))))






####################################################################################################################################

####################################################################################################################################
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
if (!file.exists('data/mask_list.RDS')){
  stack <- stack('../data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
  stack_list <- unstack(stack)
  stack_list <- rev(stack_list)
  
  stack_list_c <- stack_list
  for(i in 1:length(stack_list_c)){
    stack_list_c[[i]][stack_list_c[[i]] == 1] <- NA
  }
  saveRDS(stack_list_c, 'data/mask_list.RDS')
  
} else {
  stack_list_c = readRDS('data/mask_list.RDS')
}

# Specify pollen times (46 total; 4 modern and 42 paleo)
times <- c(seq(0, 150, by = 50), seq(from = 650, by = 500, length.out = n_times))

# create vector of ice mask times
adam_times <- seq(0, 700, 1) * 30

# interpolate pollen rasters to temporal resolution of 30yrs
all_stack <- stack(all_list)
interp <- interpolateRasters(all_stack, interpFrom = times, interpTo = adam_times)
interp_list <- unstack(interp)

n_times <- length(adam_times)#length(stack_list_c)
resamp_list <- interp_list
for(i in 1:n_times){
  # idx_time = which(adam_times == times[i])
  resamp_list[[i]] <- raster::resample(resamp_list[[i]], y = stack_list_c[[i]])
  # mask_list[[i]] <- mask(mask_list[[i]], mask = stack_list_c[[idx_time]])
}
resamp_rast <- stack(resamp_list)

# # GET BIOTIC VELOCITY FROM INTERPOLATED RASTERS
# interp_list <- unstack(interp)
# interp_list <- rev(interp_list)  # oldest needs to be on top

#bv_times <- adam_times * -1
#bv_times <- bv_times[15:697]  # just for now, to get the function to work
#bvs <- bioticVelocity(masked_rast, times = bv_times, onlyInSharedCells = TRUE)
if (BVs){
bvs_interp <- bioticVelocity(resamp_rast, times = adam_times, onlyInSharedCells = TRUE)
saveRDS(bvs_interp, 'data/BVs_interp_paleo.RDS')
}

occ_interp = occupied(resamp_rast)
value_sum_interp = value_sum(resamp_rast)
range_summary = rbind(range_summary,
                   data.frame(n_cells = occ_interp, 
                              value_sum = value_sum_interp,
                              time = adam_times, 
                              type = rep('interp', length(adam_times))))



# mask out ice/water
n_times <- length(stack_list_c)
mask_list <- resamp_list
for(i in 1:n_times){
  # mask_list[[i]] <- raster::resample(mask_list[[i]], y = stack_list_c[[i]])
  mask_list[[i]] <- mask(mask_list[[i]], mask = stack_list_c[[i]])
}

masked_rast <- stack(mask_list)

# THE PLOTS DON'T LOOK LIKE FRAXINUS DISTRIBUTIONS. 
# I DOUBLE CHECKED THE PULL POLLEN CODE TO MAKE SURE MY TAXON LIST WAS CORRECT. LOOKS FINE.

# interp <- stack(interp_list)
# bv_times <- adam_times * -1
# bv_times <- bv_times[15:697]  # just for now, to get the function to work
if (BVs){
bvs_interp_mask <- bioticVelocity(masked_rast, times = adam_times, onlyInSharedCells = TRUE)
saveRDS(bvs_interp_mask, 'data/BVs_interp_mask_paleo.RDS')
}

occ_interp_mask = occupied(masked_rast)
value_sum_interp_mask = value_sum(masked_rast)

range_summary = rbind(range_summary,
                   data.frame(n_cells = occ_interp_mask, 
                              value_sum = value_sum_interp_mask,
                              time = adam_times, 
                              type = rep('interp-mask', length(adam_times))))


saveRDS(range_summary, 'data/range_summary_paleo.RDS')


####################################################################################################################################
## No interp
####################################################################################################################################

# MASK OUT SEAS, LAKES, ICE
# rasterstack masks are arranged LGM (item 1) to modern (item 701)
# I think modern = 2000 A.D.
if (!file.exists('data/mask_list.RDS')){
  stack <- stack('../data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
  stack_list <- unstack(stack)
  stack_list <- rev(stack_list)
  
  stack_list_c <- stack_list
  for(i in 1:length(stack_list_c)){
    stack_list_c[[i]][stack_list_c[[i]] == 1] <- NA
  }
  saveRDS(stack_list_c, 'data/mask_list.RDS')
  
} else {
  stack_list_c = readRDS('data/mask_list.RDS')
}

locs <- readRDS('data/grid_3.0.RDS')
proj <- '+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs'

paleo <- readRDS('output/preds_paleo_mean_frax.RDS')
modern <- readRDS('output/preds_modern_mean_frax.RDS')

# Specify pollen times (46 total; 4 modern and 42 paleo)
# times <- c(seq(0, 150, by = 50), seq(from = 650, by = 500))
times <- c(seq(0, 150, by = 50), seq(from = 650, by = 500, length.out = 42))

# create vector of ice mask times
adam_times <- seq(0, 700, 1) * 30

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

n_times = length(times)

# interpolate pollen rasters to temporal resolution of 30yrs
all_stack <- stack(all_list)

occ_raw = occupied(all_stack)
value_sum_raw = value_sum(all_stack)

range_raw = data.frame(n_cells = occ_raw, value_sum = value_sum_raw, time = times, type = rep('raw', length(times)))

ggplot(data=range_raw) + geom_line(aes(x=time, y=value_sum, color=type))

ggplot(data=range_raw) + geom_line(aes(x=time, y=n_cells, color=type))


bvs_raw <- bioticVelocity(all_stack, times = times, onlyInSharedCells = TRUE)
saveRDS(bvs_raw, 'data/BVs_raw_paleo.RDS')

bvs_raw = readRDS('data/BVs_raw_paleo.RDS')

ggplot(data = bvs_raw, aes(x = timeFrom, y = nCentroidVelocity)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot(data = bvs_raw, aes(x = centroidLong, y = centroidLat, color=timeFrom)) +
  geom_point() +
  labs(x = 'Long', y = 'Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

library(fields)

n_times <- length(times)#length(stack_list_c)
resamp_list <- all_list
dist_mat = rdist(times, adam_times)

for(i in 1:n_times){
  idx_time = which.min(dist_mat[i,])
  print(idx_time)
  resamp_list[[i]] <- raster::resample(resamp_list[[i]], y = stack_list_c[[idx_time]])
  # mask_list[[i]] <- mask(mask_list[[i]], mask = stack_list_c[[idx_time]])
}
resamp_rast <- stack(resamp_list)

bvs_raw_resamp <- bioticVelocity(resamp_rast, times = times, onlyInSharedCells = TRUE)
saveRDS(bvs_raw_resamp, 'data/BVs_raw_resamp_paleo.RDS')

ggplot(data = bvs_raw_resamp, aes(x = timeFrom, y = nCentroidVelocity)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot(data = bvs_raw_resamp, aes(x = centroidLong, y = centroidLat, color=timeFrom)) +
  geom_point() +
  labs(x = 'Long', y = 'Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))


ggplot() +
  geom_point(data = bvs_raw_resamp, aes(x = timeFrom, y = centroidLat)) +
  geom_point(data = bvs_raw, aes(x = timeFrom, y = centroidLat), color='blue') +
  geom_line(data = bvs_raw_resamp, aes(x = timeFrom, y = centroidLat), alpha=0.1) +
  geom_line(data = bvs_raw, aes(x = timeFrom, y = centroidLat), color='blue', alpha=0.1) +
  # labs(x = 'Long', y = 'Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot() +
  geom_point(data = bvs_raw_resamp, aes(x = timeFrom, y = centroidLong)) +
  geom_point(data = bvs_raw, aes(x = timeFrom, y =centroidLong), color='blue') +
  geom_line(data = bvs_raw_resamp, aes(x = timeFrom, y = centroidLong), alpha=0.1) +
  geom_line(data = bvs_raw, aes(x = timeFrom, y = centroidLong), color='blue', alpha=0.1) +
  # labs(x = 'Long', y = 'Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))


ggplot() +
  geom_point(data = bvs_raw_resamp, aes(x = timeFrom, y = centroidVelocity)) +
  geom_point(data = bvs_raw, aes(x = timeFrom, y =centroidVelocity), color='blue') +
  geom_line(data = bvs_raw_resamp, aes(x = timeFrom, y = centroidVelocity), alpha=0.1) +
  geom_line(data = bvs_raw, aes(x = timeFrom, y = centroidVelocity), color='blue', alpha=0.1) +
  # labs(x = 'Long', y = 'Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

# mask out ice/water
n_times <- length(times)
mask_list <- resamp_list
for(i in 1:n_times){
  idx_time = which.min(dist_mat[i,])
  print(idx_time)
  mask_list[[i]] <- mask(mask_list[[i]], mask = stack_list_c[[idx_time]])
}

masked_rast <- stack(mask_list)

bvs_raw_resamp_masked <- bioticVelocity(masked_rast, times = times, onlyInSharedCells = TRUE)
saveRDS(bvs_raw_resamp_masked, 'data/BVs_raw_resamp_masked_paleo.RDS')

ggplot() +
  geom_point(data = bvs_raw_resamp, aes(x = timeFrom, y = centroidVelocity)) +
  geom_point(data = bvs_raw_resamp_masked, aes(x = timeFrom, y = centroidVelocity), color='red') +
  geom_point(data = bvs_raw, aes(x = timeFrom, y =centroidVelocity), color='blue') +
  geom_line(data = bvs_raw_resamp, aes(x = timeFrom, y = centroidVelocity), alpha=0.1) +
  geom_line(data = bvs_raw, aes(x = timeFrom, y = centroidVelocity), color='blue', alpha=0.1) +
  geom_line(data = bvs_raw_resamp_masked, aes(x = timeFrom, y = centroidVelocity), color='red', alpha=0.1) +
  # labs(x = 'Long', y = 'Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))


####################################################################################################################################
## BV plots
####################################################################################################################################

ggplot(data=range_summary) + geom_line(aes(x=time, y=value_sum, color=type))

ggplot(data=range_summary) + geom_line(aes(x=time, y=n_cells, color=type))

bvs_interp = readRDS('data/BVs_interp_paleo.RDS')

ggplot(data = bvs_interp, aes(x = timeFrom*30, y = centroidVelocity)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot(data = bvs_interp, aes(x = timeFrom, y = nCentroidVelocity)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot(data = bvs_interp, aes(x = centroidLong, y = centroidLat, color=timeFrom)) +
  geom_point() +
  labs(x = 'Long', y = 'Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

bvs_interp_mask = readRDS('data/BVs_interp_mask_paleo.RDS')

ggplot(data = bvs_interp_mask, aes(x = timeFrom, y = nCentroidVelocity)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))


ggplot(data = bvs_interp_mask, aes(x = timeFrom, y = centroidLat)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot(data = bvs_interp_mask, aes(x = timeFrom, y = centroidLong)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid Long') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot(data = bvs_interp_mask, aes(x = centroidLong, y = centroidLat, color=timeFrom)) +
  geom_point() +
  labs(x = 'Long', y = 'Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

bvs_krige = readRDS('data/BVs_krige_paleo.RDS')

ggplot(data = bvs_krige, aes(x = timeFrom, y = nCentroidVelocity)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))


ggplot(data = bvs_krige, aes(x = timeFrom, y = centroidLat)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot(data = bvs_krige, aes(x = timeFrom, y = centroidLong)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid Long') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot(data = bvs_krige, aes(x = centroidLong, y = centroidLat, color=timeFrom)) +
  geom_point() +
  labs(x = 'Long', y = 'Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

bvs_krige_mask = readRDS('data/BVs_krige_mask_paleo.RDS')

ggplot(data = bvs_krige_mask, aes(x = timeFrom, y = nCentroidVelocity)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))


ggplot(data = bvs_krige_mask, aes(x = timeFrom, y = centroidLat)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot(data = bvs_krige_mask, aes(x = timeFrom, y = centroidLong)) +
  geom_line() +
  labs(x = 'Time from (years before present)', y = 'Centroid Long') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

ggplot(data = bvs_krige_mask, aes(x = centroidLong, y = centroidLat, color=timeFrom)) +
  geom_point() +
  labs(x = 'Long', y = 'Lat') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        text = element_text(size = 14))

# # minor_breaks <- rev(times$from) * -1
# # breaks <- minor_breaks[seq(1,n_times, by = 1)]
# 
# ggplot(data = bvs, aes(x = timeFrom, y = centroidVelocity)) +
#   geom_line() + 
#   #scale_x_continuous(minor_breaks = minor_breaks, breaks = breaks) +
#   labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
#         text = element_text(size = 14))
# ggsave('figures/bvs_after_30yr_krige_paleo.png')
# 
# 
# ##################################################################################################################################
# 
# dat_centroid = data.frame(long=bvs$centroidLong, lat=bvs$centroidLat, time=bvs$timeFrom)
# 
# ggplot() +
#   geom_point(data=dat_centroid, aes(x=long, y=lat, colour=time))
# 
# 
# 
# ggplot(data = dat, aes(x = times, y = n_cells)) +
#   geom_line() + 
#   #scale_x_continuous(minor_breaks = minor_breaks, breaks = breaks) +
#   labs(x = 'Time (years before present)', y = 'Range size (cells)') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
#         text = element_text(size = 14))
# ggsave('figures/range_size_after_30yr_krige_paleo.png')
# 
# 
# 
# ##############################################################################
# bvs_land <- bioticVelocity(stack, times = seq(0, 700, 1)*30, onlyInSharedCells = TRUE)
# saveRDS(bvs_land, 'data/BVs_land_paleo.RDS')
# 
# 
# # minor_breaks <- rev(times$from) * -1
# # breaks <- minor_breaks[seq(1,n_times, by = 1)]
# 
# ggplot(data = bvs_land, aes(x = timeFrom*30, y = centroidVelocity)) +
#   geom_line() + 
#   #scale_x_continuous(minor_breaks = minor_breaks, breaks = breaks) +
#   labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
#         text = element_text(size = 14))
# ggsave('figures/bvs_land_30yr_interp_paleo.png')
# 
# 
# dat_centroid_land = data.frame(long=bvs_land$centroidLong, lat=bvs_land$centroidLat, time=bvs_land$timeFrom*30)
# 
# ggplot() +
#   geom_point(data=dat_centroid, aes(x=long, y=lat, colour=time))
# ggsave('figures/bvs_land_30yr_centroid_movement_paleo.png')
# 
# ggplot() + 
#   geom_line(data = bvs, aes(x = timeFrom, y = centroidVelocity)) +
#   geom_line(data = bvs_land, aes(x = timeFrom*30, y = centroidVelocity), colour='blue') + 
#   #scale_x_continuous(minor_breaks = minor_breaks, breaks = breaks) +
#   labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
#         text = element_text(size = 14))
# ggsave('figures/bvs_land_30yr_interp_paleo.png')
# 
# dat_centroid_land$time = dat_centroid_land$time*30 
# centroids = rbind(data.frame(dat_centroid_land, type="land"), 
#                   data.frame(dat_centroid, type="frax-krige"))
# 
# ggplot(data=centroids) +
#   geom_point(aes(x=long, y=lat, group=type, colour=type))
#   
# ggsave('figures/bvs_land_frax_30yr_centroid_movement_paleo.png')
# 
# 
# 
# bvs_krige = readRDS('data/BVs_krige_paleo.RDS')
# bvs_interp = readRDS('data/interp_bvs.RDS')
# ggplot() + 
#   geom_line(data = bvs_krige, aes(x = timeFrom, y = centroidVelocity)) +
#   geom_line(data = bvs_interp, aes(x = -timeFrom, y = centroidVelocity), colour='blue') + 
#   #scale_x_continuous(minor_breaks = minor_breaks, breaks = breaks) +
#   labs(x = 'Time from (years before present)', y = 'Centroid velocity (m/yr)') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
#         text = element_text(size = 14))
# ggsave('figures/bvs_30yr_both_paleo.png')
