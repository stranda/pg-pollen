library(reshape2)
library(rdist)
library(rgeos)
library(ggplot2)
library(sp)
library(rgdal)
library(raster)
library(fields)

#### READ IN DATA ####
dat = readRDS('data/12taxa_pollen_dat_v1.0.RDS')

# note that locations were scaled to fit the model
# unscaling to think in meters, then will rescale again before prediction
locs_pollen <- data.frame(dat[,c('x', 'y')] )

taxa = colnames(dat)[!(colnames(dat) %in% c('x', 'y'))]
y = dat[,..taxa]

#### CONSTRUCT GRID ####
# function to construct a raster grid
# take in bounding box for grid, resolution in m, and projection
build_grid <- function(veg_box, resolution = 24000, proj = '+init=epsg:3175') {
  raster::raster(xmn = veg_box[1],
                 xmx = veg_box[3],
                 ymn = veg_box[2],
                 ymx = veg_box[4],
                 resolution = resolution,
                 crs = proj)
}

bbox_tran <- function(x, coord_formula = '~ x + y', from, to) {
  
  sp::coordinates(x) <- formula(coord_formula)
  sp::proj4string(x) <- sp::CRS(from)
  bbox <- as.vector(sp::bbox(sp::spTransform(x, CRSobj = sp::CRS(to))))
  return(bbox)
}

# get bounding box from pollen record coordinates
pol_box <- bbox_tran(locs_pollen, '~ x + y',
                     proj_out,
                     proj_out)
xlim = c(pol_box[1]-24000, pol_box[3]+24000)
ylim = c(pol_box[2]-24000, pol_box[4]+24000)

# build the raster grid
# 40 km grid cells
reconst_grid <- build_grid(pol_box,
                           resolution = 60000,
                           proj = proj_out)

# want to work with a data frame not a raster
reconst_grid = as.data.frame(reconst_grid, xy=TRUE)

locs_grid = reconst_grid[,1:2]

saveRDS(locs_grid, 'data/grid.RDS')
