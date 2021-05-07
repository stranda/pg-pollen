
require(raster)
require(sp)
require(rgdal)
require(ggplot2)

# setwd('C:/Users/abrow/Documents/pg-pollen')
version <- '3.0'

# READ IN DATA
dat <- readRDS(paste0('data/', 'modern_pollen_locs_', version, '.RDS'))
locs_pollen <- dat$locs

# Get spatial domain and projection you want to use
# (Adam Smith's .tif from: NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
proj <- proj4string(stack)

# build a raster that identifies gric cells that were land at any point since LGM
# 1 if land, 0 if not land
land = stack[[1]]
nrow = dim(land)[1]
ncol = dim(land)[2]
for (i in 1:nlayers(stack)){
  print(paste0("layer :", i))
  raster_sub = stack[[i]]
  land <- overlay(land, raster_sub, fun = function(x, y) {
    x[!is.na(y[])] <- 1
    return(x)
  })
}


# Adam's spatial domain
coords <- data.frame(x = c(stack@extent@xmin,stack@extent@xmax),
                     y = c(stack@extent@ymin,stack@extent@ymax))
sp::coordinates(coords) <- ~x + y
sp::proj4string(coords) <- CRS(proj)
pol_bbox <- bbox(coords)

# construct the grid using Adam's bounding box, 40km resolution, Adam's projection
reconst_grid <- raster::raster(xmn = pol_bbox[1,1],
                               xmx = pol_bbox[1,2],
                               ymn = pol_bbox[2,1],
                               ymx = pol_bbox[2,2],
                               resolution = 60000,
                               crs = proj)
vals = rep(1, ncell(reconst_grid))
reconst_grid = setValues(reconst_grid, vals)

land_resamp <- resample(land,reconst_grid,method='bilinear')

# land_array = as.array(land_resamp)[,,1]

extent(land_resamp) == extent(reconst_grid)

reconst_grid_land <- overlay(reconst_grid, land_resamp, fun = function(x, y) {
  x[is.na(y[])] <- NA
  return(x)
})

# want to work with a data frame, not a raster
reconst_land_df = raster::as.data.frame(reconst_grid_land, xy=TRUE)
reconst_land_df = reconst_land_df[which(!is.na(reconst_land_df$layer)),]

locs_grid = reconst_land_df[,1:2]

saveRDS(locs_grid, paste0('data/grid_', version, '.RDS'))


# visualize grid
na_shp <- readOGR("data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("data/map-data/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj)

ggplot(data = locs_grid, aes(x = x, y = y)) + 
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
  geom_point(alpha=0.3) +
  coord_equal()

