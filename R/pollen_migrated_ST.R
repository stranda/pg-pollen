
library(neotoma)
library(analogue)
library(Bchron)
library(ggplot2)
library(data.table)
library(sp)
library(raster)
library(mapproj)
library(ggplot2)
library(dplyr)
library(tidyr)

# setwd('C:/Users/abrow/Documents/pg-pollen')
version = '3.0'

# READ SHAPEFILES AND ESTABLISH SPATIAL DOMAIN
# USA Contiguous albers equal area
# OLD PROJECTION
# proj_out <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 
# +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

# Get raster masks and spatial domain you want to use
# (Adam Smith's .tif from: NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
proj <- proj4string(stack)

# Adam's spatial domain
coords <- data.frame(x = c(stack@extent@xmin,stack@extent@xmax),
                     y = c(stack@extent@ymin,stack@extent@ymax))
sp::coordinates(coords) <- ~x + y
sp::proj4string(coords) <- CRS(proj)
coords <- sp::spTransform(coords, CRS("+init=epsg:4326"))

# limits
lat_hi  = coords@coords[2,2]
lat_lo  = coords@coords[1,2]
long_west = coords@coords[1,1]
long_east = coords@coords[2,1]

# load pollen data 
compiled.cores <- read.csv('data/pollen_north_america_v6.0.csv', stringsAsFactors = FALSE)
load('data/sites_north_america.rdata')
load('data/pollen.equiv.rda')

sites.cores = sites_north_america

compiled.cores = data.frame(long = sites.cores[match(compiled.cores$dataset_id, sites.cores$datasetid), 'longitude'],
                           lat = sites.cores[match(compiled.cores$dataset_id, sites.cores$datasetid), 'latitude'],
                           compiled.cores)

# plot pollen core locations within bounding box
map <- map_data("world")
ggplot(data = data.frame(map), aes(long, lat)) + 
  geom_polygon(aes(group=group), color = "steelblue", alpha = 0.2) +
  geom_point(data = sites.cores,
             aes(x = longitude, y = latitude), color = 2) +
  xlab("Longitude West") + 
  ylab("Latitude North") +
  coord_map(projection = "albers",
            lat0 = lat_lo,
            lat1 = lat_hi, 
            xlim = c(long_west, long_east),
            ylim = c(lat_lo, lat_hi))


# until there's an answer as to why we get age NAs after running 'compile_downloads'...
# for now just remove rows with age NAs
compiled.cores <- compiled.cores[!is.na(compiled.cores$age), ]

# # translate any dates in radiocarbon years to calendar years
# radio.years <- (compiled.cores$agetype %in% "Radiocarbon years BP") &
#   (compiled.cores$age > 95 ) &
#   (compiled.cores$age < 50193)
# sryears <- sum(radio.years, na.rm = TRUE)
# 
# # BChronCalibrate is in the BChron package:
# calibrated <- BchronCalibrate(compiled.cores$age[radio.years],
#                               ageSds = rep(100, sryears),
#                               calCurves = rep("intcal20", sryears))
# #  we want the weighted means from "calibrated"
# wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
# compiled.cores$age[radio.years] <- sapply(calibrated, wmean.date)
# # saveRDS(compiled.cores, 'data/compiled_cores_P25_all_times.RDS')

# remove any samples with ages greater than 50000 YBP
compiled.cores <- compiled.cores[which(compiled.cores$age < 21500), ]
hist(compiled.cores$age)

# only keep pollen from within domain of interest
compiled.cores <- compiled.cores[compiled.cores$lat < lat_hi & 
                            compiled.cores$lat > lat_lo & 
                            compiled.cores$long > long_west & 
                            compiled.cores$long < long_east ,]

# visualize subsetted data
map <- map_data("world")
ggplot(data = data.frame(map), aes(long, lat)) + 
  geom_polygon(aes(group=group), color = "steelblue", alpha = 0.2) +
  geom_point(data = compiled.cores,
             aes(x = long, y = lat), color = 2, size = 2) +
  xlab("Longitude West") + 
  ylab("Latitude North") +
  coord_map(projection = "albers",
            lat0 = lat_lo,
            lat1 = lat_hi, 
            xlim = c(long_west, -59),
            ylim = c(lat_lo, lat_hi))

# remove pollen cores from isolated locations (islands, mexico)
bermuda <- compiled.cores[compiled.cores$lat < 34 & compiled.cores$long > -78, ]
bahamas_cuba <- with(compiled.cores, compiled.cores[lat < 28 & long > -80, ])
mexico <- with(compiled.cores, compiled.cores[lat < 23, ])
remove.sites <- rbind(bermuda, bahamas_cuba, mexico)
remove.sites <- remove.sites %>% distinct()
compiled.cores <- compiled.cores[!(compiled.cores$siteid) %in% remove.sites$siteid,]

# convert projection to proj
sp::coordinates(compiled.cores) <- ~long+lat
sp::proj4string(compiled.cores) <- sp::CRS('+init=epsg:4326')
albers <- sp::spTransform(compiled.cores, proj)
xy = data.frame(coordinates(albers))
colnames(xy) = c('x', 'y')

# construct data frame with re-projected coordinates
compiled.cores = data.frame(xy, compiled.cores)
compiled.cores = compiled.cores[,which(colnames(compiled.cores)!= 'optional')]
# save this for later when you'll need to calculate relative proportions for model validation
# saveRDS(compiled.cores, 'data/compiled_cores_P25.RDS')

# remove non-tree taxa
taxa.nontree <- c('Other', 'Prairie.Forbs', 'Poaceae', 'Cyperaceae')
tree.cores <- compiled.cores[, which(!(colnames(compiled.cores) %in% taxa.nontree))]

# remove sites with no tree pollen counts
# i.e., rows that contain only zeros/NAs across all taxa
start_col <- which(colnames(tree.cores) == 'Abies')
tree.cores$sum <- apply(tree.cores[,start_col:ncol(tree.cores)], 1, function(x) sum(x, na.rm=TRUE))
tree.cores <- tree.cores[tree.cores$sum > 0, ]
tree.cores <- tree.cores %>% dplyr::select(-sum)

# specify which tree taxa to keep
# use 13 taxa with highest relative proportions
props <- tree.cores[,start_col:ncol(tree.cores)]/rowSums(tree.cores[,13:ncol(tree.cores)], na.rm = TRUE)
props <- pivot_longer(props, cols = colnames(props), names_to = 'taxon', values_to = 'props')
props <- props %>% group_by(taxon) %>% summarise(props = sum(props))
props <- props %>% arrange(desc(props))
taxa.keep <- props[1:13, ]
taxa.keep <- as.character(taxa.keep$taxon)

# ASSIGN TIME CHUNKS
# Modern time = AD 1800 - 2000 (i.e., 150 to -50 YBP)
# Paleo time = 21,000 - 150 YBP

# MODERN TIME
modern <- tree.cores[tree.cores$age >= -50 & tree.cores$age < 150, ]
modern_bins <- seq(-50, 150, by = 50)
modern_cut <- cut(modern$age, include.lowest = TRUE, breaks = modern_bins)
modern$cut <- as.integer(modern_cut)

# assign unique ID to each distinct x, y coordinate
modern_xyid <- modern %>% dplyr::select(x,y) %>% distinct()
modern_xyid$id <- as.character(seq(1, nrow(modern_xyid), by = 1))

# separate dataframe into list of dataframes, one for each time chunk, remove excess columns
modern_time <- split(modern, f = modern$cut)

# PALEO TIME
paleo <- tree.cores[tree.cores$age >= 150, ]
paleo_bins <- seq(150, max(paleo$age), by = 500)
paleo_cut <- cut(paleo$age, include.lowest = TRUE, breaks = paleo_bins)
paleo$cut <- as.integer(paleo_cut)
paleo <- paleo[!is.na(paleo$cut),] # remove ages > highest time bin

# assign unique ID to each distinct x, y coordinate
paleo_xyid <- paleo %>% dplyr::select(x,y) %>% distinct()
paleo_xyid$id <- as.character(seq(1, nrow(paleo_xyid), by = 1))

# separate dataframe into list of dataframes, one for each time chunk, remove excess columns
paleo_time <- split(paleo, f = paleo$cut)


# CREATE 'OTHER' TAXON BY COMBINING TREE COUNTS FROM TREES NOT BEING MODELED INDIVIDUALLY
# (if you want to retain 'Other tree' column, you'll need to use length(n_taxa) + 1 in later code)
# SUM POLLEN COUNTS BY TIME PERIOD/SITE
# MERGE WITH SITE IDENTIFIER; EACH TIME PERIOD SHOULD HAVE SAME # ROWS (SITES)
# CONVERT VALUES TO INTEGERS

# MODERN
n_times <- length(modern_time)
for(i in 1:n_times){
  compiled.meta = modern_time[[i]][,c('x', 'y')]
  compiled.counts = modern_time[[i]][,13:ncol(modern_time[[i]])]
  compiled.other = compiled.counts[, which((!(colnames(compiled.counts) %in% taxa.keep)) & 
                                          (!(colnames(compiled.counts) %in% taxa.nontree)))]
  compiled.counts.sub = data.frame(compiled.counts[, which(colnames(compiled.counts) %in% taxa.keep)],
                                   Other = rowSums(compiled.other, na.rm=TRUE))
  modern_time[[i]] = data.frame(compiled.meta, compiled.counts.sub)
}

for(i in 1:n_times){
  modern_time[[i]] <- modern_time[[i]] %>% group_by(x, y) %>% summarise_all(.funs = sum, na.rm = TRUE)
  modern_time[[i]] <- ungroup(modern_time[[i]])
  modern_time[[i]] <- left_join(modern_xyid, modern_time[[i]], by = c('x','y'))
  modern_time[[i]][,4:ncol(modern_time[[i]])] <- apply(modern_time[[i]][,4:ncol(modern_time[[i]])], c(1,2), as.integer)
}


# PALEO
n_times <- length(paleo_time)
for(i in 1:n_times){
  compiled.meta = paleo_time[[i]][,c('x', 'y')]
  compiled.counts = paleo_time[[i]][,13:ncol(paleo_time[[i]])]
  compiled.other = compiled.counts[, which((!(colnames(compiled.counts) %in% taxa.keep)) & 
                                             (!(colnames(compiled.counts) %in% taxa.nontree)))]
  compiled.counts.sub = data.frame(compiled.counts[, which(colnames(compiled.counts) %in% taxa.keep)],
                                   Other = rowSums(compiled.other, na.rm=TRUE))
  paleo_time[[i]] = data.frame(compiled.meta, compiled.counts.sub)
}

for(i in 1:n_times){
  paleo_time[[i]] <- paleo_time[[i]] %>% group_by(x, y) %>% summarise_all(.funs = sum, na.rm = TRUE)
  paleo_time[[i]] <- ungroup(paleo_time[[i]])
  paleo_time[[i]] <- left_join(modern_xyid, paleo_time[[i]], by = c('x','y'))
  paleo_time[[i]][,4:ncol(paleo_time[[i]])] <- apply(paleo_time[[i]][,4:ncol(paleo_time[[i]])], c(1,2), as.integer)
}


# CONSTRUCT DATA LIST, INCLUDING LOCATIONS, INTEGER ARRAY, AND TAXA.KEEP

# MODERN
modern_locs <- modern_time[[1]][,c('x','y','id')]
n_locs <- nrow(modern_locs)
n_taxa <- length(taxa.keep)
n_times <- length(modern_time)

modern_dat_array <- array(0, dim = c(n_locs, n_taxa, n_times))
for (i in 1:n_times){
  for (j in 1:n_locs){
    site_id = as.numeric(modern_time[[i]]$id[j])
    modern_dat_array[site_id,,i] = as.numeric(unname(modern_time[[i]][j,4:(4+n_taxa-1)]))
  }
}
modern_locs <- modern_locs[,c('x','y')]

saveRDS(modern_dat_array, paste0('data/', 'modern_pollen_dat_', version, '.RDS'))
saveRDS(modern_locs, paste0('data/', 'modern_pollen_locs_', version, '.RDS'))
saveRDS(taxa.keep, paste0('data/', 'pollen_taxa_', version, '.RDS'))


# PALEO
paleo_locs <- paleo_time[[1]][,c('x','y','id')]
n_locs <- nrow(paleo_locs)
n_times <- length(paleo_time)

paleo_dat_array <- array(0, dim = c(n_locs, n_taxa, n_times))
for (i in 1:n_times){
  for (j in 1:n_locs){
    site_id = as.numeric(paleo_time[[i]]$id[j])
    paleo_dat_array[site_id,,i] = as.numeric(unname(paleo_time[[i]][j,4:(4+n_taxa-1)]))
  }
}
paleo_locs <- paleo_locs[,c('x','y')]

saveRDS(paleo_dat_array, paste0('data/', 'paleo_pollen_dat_', version, '.RDS'))
saveRDS(paleo_locs, paste0('data/', 'paleo_pollen_locs_', version, '.RDS'))



#ANDRIA'S EXTRA CODE - TRYING TO FIND SITE UNDER GLACIER

dat_array = readRDS('data/pollen_dat_1.0.RDS')
locs = readRDS('data/pollen_locs_1.0.RDS')
taxa.keep = readRDS('data/pollen_taxa_1.0.RDS')

dat_lgm = dat_array[,,22]
colnames(dat_lgm) = c(taxa.keep,"Other")
dat = data.frame(locs, dat_lgm)
dat = dat[!is.na(dat[,3]),]

test_sub = test[which((test$cut==22)&(test$lat>49)),]


na_shp <- readOGR("data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("data/map-data/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj)

p <- ggplot() +
  geom_point(data = dat, aes(x = x, y = y), col = "red") +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) #+
  scale_y_continuous(limits = ylim) +
  scale_x_continuous(limits = xlim) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        plot.title = element_text(size = 12)) +
  coord_equal()
print(p)
