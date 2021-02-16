
library(neotoma)
library(analogue)
library(Bchron)
library(ggplot2)
library(data.table)
library(sp)
library(dplyr)
library(raster)
library(mapproj)

# setwd('C:/Users/abrow/Documents/pg-pollen')
version = '2.0'

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
load('pollen_north_america_p25.rdata')
load('sites_north_america.rdata')
load('pollen.equiv.rda')

compiled.cores = pollen_north_america
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

# translate any dates in radiocarbon years to calendar years
radio.years <- (compiled.cores$agetype %in% "Radiocarbon years BP") &
  (compiled.cores$age > 95 ) &
  (compiled.cores$age < 50193)
sryears <- sum(radio.years, na.rm = TRUE)

# BChronCalibrate is in the BChron package:
calibrated <- BchronCalibrate(compiled.cores$age[radio.years],
                              ageSds = rep(100, sryears),
                              calCurves = rep("intcal20", sryears))
#  we want the weighted means from "calibrated"
wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
compiled.cores$age[radio.years] <- sapply(calibrated, wmean.date)

# remove any samples with ages greater than 50000 YBP
compiled.cores <- compiled.cores[which(compiled.cores$age < 21000), ]

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
compiled.cores <- compiled.cores[!(row.names(compiled.cores) %in% row.names(remove.sites)),]

# convert projection to proj
sp::coordinates(compiled.cores) <- ~long+lat
sp::proj4string(compiled.cores) <- sp::CRS('+init=epsg:4326')
albers <- sp::spTransform(compiled.cores, proj)
xy = data.frame(coordinates(albers))
colnames(xy) = c('x', 'y')

# construct data frame with re-projected coordinates
compiled.cores = data.frame(xy, compiled.cores)
compiled.cores = compiled.cores[,which(colnames(compiled.cores)!= 'optional')]

# remove non-tree taxa
# taxa.nontree <- c('Other', 'Prairie.Forbs', 'Poaceae')
taxa.nontree <- c('Other', 'Prairie.Forbs', 'Poaceae')
tree.cores <- compiled.cores[, which(!(colnames(compiled.cores) %in% taxa.nontree))]

# remove sites with no tree pollen counts
# i.e., rows that contain only zeros/NAs across all taxa
tree.cores$sum <- apply(tree.cores[,13:ncol(tree.cores)], 1, function(x) sum(x, na.rm=TRUE))
tree.cores <- tree.cores[tree.cores$sum > 0, ]

# ASSIGN TIME CHUNKS
# QUESTIONS
# PLACE ALL NEGATIVE AGES IN SAME BIN AS TIME = 0?
# For now, just bin everything that's <0 together

time_bins = c(min(tree.cores$age), seq(0, 21000, by=990))
N_times = length(time_bins)-1

test.cut <- cut(tree.cores$age, include.lowest = TRUE, breaks = time_bins)
test <- tree.cores
test$cut <- as.integer(test.cut)
test <- test[!is.na(test$cut),] # remove rows with records >21k years old

# assign unique ID to each distinct x, y coordinate
xyid <- test[,c('x','y')] %>% distinct()
xyid$id <- as.character(seq(1, nrow(xyid), by = 1))

# separate dataframe into list of dataframes, one for each time chunk, remove excess columns
time <- split(test, f = test$cut)

# QUESTION: NOW THAT WE AREN'T TESTING THE MODEL ANYMORE, SHOULD WE ADD MORE TAXA?
# FOR NOW, USE TAXA WITH >200,000 TOTAL POLLEN COUNT
# extract pollen data from specified taxa
# and combine the other tree taxa into 'other' column

props <- tree.cores[,13:35]/rowSums(tree.cores[,13:35], na.rm = TRUE)

summ <- data.frame(pollen_sum = apply(tree.cores[,13:35], 2, function(x) sum(x, na.rm=TRUE)))
summ$gtlt <- ifelse(summ$pollen_sum<200000,"less","greater")
taxa.keep <- rownames(summ[summ$pollen_sum >= 200000, ])
taxa.nontree = c('Other', 'Prairie.Forbs', 'Poaceae')

for (i in 1:N_times){
  compiled.meta = time[[i]][,c('x', 'y')]
  compiled.counts = time[[i]][,13:35]
  compiled.other = compiled.counts[, which((!(colnames(compiled.counts) %in% taxa.keep)) & 
                                          (!(colnames(compiled.counts) %in% taxa.nontree)))]
  compiled.counts.sub = data.frame(compiled.counts[, which(colnames(compiled.counts) %in% taxa.keep)],
                                   Other = rowSums(compiled.other, na.rm=TRUE))
  time[[i]]= data.frame(compiled.meta, compiled.counts.sub)
}

# sum pollen counts by time chunk/site
time_sum <- list()
for(i in 1:N_times){
  time_sum[[i]] <- time[[i]] %>% group_by(x, y) %>% summarise_all(.funs = sum, na.rm = TRUE)
  time_sum[[i]] <- ungroup(time_sum[[i]])
}

# add unique x, y coordinate identifier to each dataframe; merge so that all 
# dataframes have the same # rows (nrows in xyid object)
for(i in 1:N_times){
  time_sum[[i]] <- left_join(xyid, time_sum[[i]], by = c('x','y'))
}

# convert numeric pollen counts to integer class
for(i in 1:N_times){
  time_sum[[i]][,4:17] <- apply(time_sum[[i]][,4:17], c(1,2), as.integer)
}

# extract site coordinates
locs <- time_sum[[1]][,c('x','y','id')]
N_locs = nrow(locs)
N_taxa = length(taxa.keep) + 1

dat_array = array(0, dim = c(N_locs, N_taxa, N_times))
for (i in 1:N_times){
  for (j in 1:N_locs){
    site_id = as.numeric(time_sum[[i]]$id[j])
    dat_array[site_id,,i] = as.numeric(unname(time_sum[[i]][j,4:(4+N_taxa-1)]))
  }
}

locs <- locs[,c('x','y')]

saveRDS(dat_array, paste0('data/', 'pollen_dat_', version, '.RDS'))
saveRDS(locs, paste0('data/', 'pollen_locs_', version, '.RDS'))
saveRDS(taxa.keep, paste0('data/', 'pollen_taxa_', version, '.RDS'))

#

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
