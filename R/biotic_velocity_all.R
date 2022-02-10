#setwd('C:/Users/abrow/Documents/pg-pollen')
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
require(parallel)

################################################
## CALCULATE AND PLOT BIOTIC VELOCITIES
################################################

# read prediction output (50 iterations randomly sampled from all 1000 model iterations)
locs_grid <- readRDS('data/grid_4.1.RDS')
preds <- readRDS('output/polya-gamma-predictions_4.1_overdispersed.RDS')
taxa = readRDS("data/taxa_4.1.RDS")

# read raster masks (provides masks for spatial domain/resolution of genetic + ENM data)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
proj <- proj4string(stack)


#tx = 1  #taxon number


for (tx in 1:dim(preds[["pi"]])[3])
if (!file.exists(paste0('output/',taxa[tx],'_bvs_n200_v4.1.RDS')))
{
txpreds = preds[["pi"]][,,tx,]

# specify time bins; use median of time bins for subsequent work
time <- c(-70, seq(705, by = 990, length.out = 22))
time <- data.frame(timeFrom = time[-1], timeTo = time[-23])
time$time_mid <- (time$timeFrom + time$timeTo) / 2


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
n_iter <- dim(txpreds)[1]
n_times <- dim(txpreds)[3]

preds_list <- list() 
for(i in 1:n_iter){
  preds_list[[i]] <- cbind(locs_grid, data.frame(txpreds[i,,]))
  preds_list[[i]] <- pivot_longer(preds_list[[i]], cols = 3:ncol(preds_list[[i]]),
                                  names_to = 'time', values_to = 'pred')
  preds_list[[i]]$time <- as.integer(substr(preds_list[[i]]$time, start = 2,
                                            stop = nchar(preds_list[[i]]$time)))
  preds_list[[i]] <- split(preds_list[[i]], f = preds_list[[i]]$time)
  preds_list[[i]] <- lapply(preds_list[[i]], function(x) { x['time'] <- NULL; x })
}

if (FALSE)
    {
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
}

# convert list of list of dataframes to list of rasters
# takes 1-2 minutes to run

rasterstack_list <- mclapply(1:n_iter,mc.cores=12,function(i)
{
    lapply(1:n_times, function(j)
        {
            l <- rasterFromXYZ(preds_list[[i]][[j]])
            proj4string(l) <- proj
            l <- raster::resample(l, y = stack_sub[[j]])
            mask(l, mask = stack_sub[[j]])
  })
})


### weight relative abundance by proportion of cell covered by land (for cells with partial glacial coverage)
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

##parallelized calculate BVs
bv_list = mclapply(rasterstack_rev, mc.cores=12, function(x) {
    enmSdm::bioticVelocity(x, times=(rev(time$time_mid)* -1), onlyInSharedCells=FALSE,
                           quants=c(0.05,0.25,0.5,0.75,0.95)) #false to match climate
    })

# convert list of dataframes to single dataframe
bvs <- bind_rows(bv_list)
saveRDS(bvs, paste0('output/',taxa[tx],'_bvs_n200_v4.1.RDS'))

}


#do some plotting
#first make a mega dataset containing all taxa
abv  = NULL #going to increment rows in a df, slow but not many to do
for (tx in 1:dim(preds[["pi"]])[3])
{
    bv = readRDS(paste0('output/',taxa[tx],'_bvs_n200_v4.1.RDS'))
    bv$taxon=taxa[tx]
    abv = rbind(abv,bv)
}

pdf ("figures/biotic_velocities.pdf",width=12,height=12)

ggplot(abv, aes(x=as.factor(timeFrom), y=centroidVelocity)) + geom_boxplot() + facet_wrap(~taxon) + theme(axis.text.x = element_text(angle=90,hjust=1))


ggplot(abv, aes(x=as.factor(timeFrom), y=nsQuantVelocity_quant0p05)) + geom_boxplot() + facet_wrap(~taxon) + theme(axis.text.x = element_text(angle=90,hjust=1))

ggplot(abv, aes(x=as.factor(timeFrom), y=nsQuantVelocity_quant0p95)) + geom_boxplot() + facet_wrap(~taxon) + theme(axis.text.x = element_text(angle=90,hjust=1))


dev.off()




