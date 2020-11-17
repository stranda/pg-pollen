
library(neotoma)
library(analogue)
library(Bchron)
library(ggplot2)
library(data.table)
library(sp)
library(dplyr)
library(raster)

# setwd('C:/Users/abrow/Documents/pg-pollen')

# READ SHAPEFILES AND ESTABLISH SPATIAL DOMAIN
# USA Contiguous albers equal area
proj_out <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 
+ellps=GRS80 +datum=NAD83 +units=m +no_defs'

# Get raster masks and spatial domain you want to use
# (Adam Smith's .tif from: NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks)
stack <- stack('data/map-data/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')

# Adam's spatial domain
coords <- data.frame(x = c(stack@extent@xmin,stack@extent@xmax),
                     y = c(stack@extent@ymin,stack@extent@ymax))
sp::coordinates(coords) <- ~x + y
sp::proj4string(coords) <- CRS(proj_out)
coords <- sp::spTransform(coords, CRS("+init=epsg:4326"))

# limits
lat_hi  = coords@coords[2,2]
lat_lo  = coords@coords[1,2]
long_west = coords@coords[1,1]
long_east = coords@coords[2,1]

# version
version = 'Oct12'

# DOWNLOAD DATA withing bounded region from NEOTOMA
all.datasets <- get_dataset(loc = c(long_west, lat_lo, long_east, lat_hi),
                            datasettype = "pollen")

# plot pollen core locations within bounding box
map <- map_data("world")
ggplot(data = data.frame(map), aes(long, lat)) + 
  geom_polygon(aes(group=group), color = "steelblue", alpha = 0.2) +
  geom_point(data = get_site(all.datasets),
             aes(x = long, y = lat), color = 2) +
  xlab("Longitude West") + 
  ylab("Latitude North") +
  coord_map(projection = "albers",
            lat0 = lat_lo,
            lat1 = lat_hi, 
            xlim = c(long_west, long_east),
            ylim = c(lat_lo, lat_hi))

# check if existing download contains all datasets
if (file.exists("data/all.downloads.RDS")) {
  all.downloads.old <- readRDS("data/all.downloads.RDS")
  if (all(names(all.downloads.old) == names(all.datasets))){
    all.downloads <- all.downloads.old
  }
} else {
  all.downloads <- get_download(all.datasets, 
                                verbose = FALSE)
  saveRDS(all.downloads, "data/all.downloads.RDS")
}

# are there any age NAs?
names <- names(all.cores)
test <- lapply(names, function(x) anyNA(all.cores$x$chronologies$`palEON-STEPPS`))
test <- data.frame(test)
any(test == TRUE)
test <- lapply(names, function(x) anyNA(all.cores$x$chronologies$`NAPD 1`))
test <- data.frame(test)
any(test == TRUE)

####
# START HERE IF YOU ALREADY HAVE POLLEN DATA DOWNLOADED
####

all.downloads <- readRDS('data/all.downloads.RDS')

# standardize the pollen taxonomy
all.cores <- compile_taxa(all.downloads, "P25")

# convert list to data frame
compiled.cores  <- compile_downloads(all.cores)

# translate any dates in radiocarbon years to calendar years
radio.years <- (compiled.cores$date.type %in% "Radiocarbon years BP") &
  (compiled.cores$age > 71 ) &
  (compiled.cores$age < 46401)
sryears <- sum(radio.years, na.rm = TRUE)

# BChronCalibrate is in the BChron package:
calibrated <- BchronCalibrate(compiled.cores$age[radio.years],
                              ageSds = rep(100, sryears),
                              calCurves = rep("intcal13", sryears))
#  we want the weighted means from "calibrated"
wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))

compiled.cores$age[radio.years] <- sapply(calibrated, wmean.date)

hist(compiled.cores$age)

# visualize subsetted data
map <- map_data("world")
ggplot(data = data.frame(map), aes(long, lat)) + 
  geom_polygon(aes(group=group), color = "steelblue", alpha = 0.2) +
  geom_point(data = compiled.cores,
             aes(x = long, y = lat), color = 2) +
  xlab("Longitude West") + 
  ylab("Latitude North") +
  coord_map(projection = "albers",
            lat0 = lat_lo,
            lat1 = lat_hi, 
            xlim = c(long_west, long_east),
            ylim = c(lat_lo, lat_hi))

# convert projection to proj_out
sp::coordinates(compiled.cores) <- ~long+lat
sp::proj4string(compiled.cores) <- sp::CRS('+init=epsg:4326')
albers <- sp::spTransform(compiled.cores, proj_out)
xy = data.frame(coordinates(albers))
colnames(xy) = c('x', 'y')

# construct data frame with re-projected coordinates
compiled.cores = data.frame(xy, compiled.cores)
compiled.cores = compiled.cores[,which(colnames(compiled.cores)!= 'optional')]

# remove non-tree taxa
taxa.nontree <- c('Other', 'Prairie.Forbs', 'Poaceae')
tree.cores <- compiled.cores[, which(!(colnames(compiled.cores) %in% taxa.nontree))]

# remove sites with no pollen counts (IS THIS RIGHT????)
tree.cores[,13:35]=apply(tree.cores[,13:35], MARGIN = c(1,2), 
                         function(x) if (is.na(x)){x=0}else{x})

# standardize pollen counts so every site has same total # pollen grains across taxa
# (only need to do this if you're analyzing a subset of all taxa in dataset)
stand.cores <- tree.cores
# stand.cores[ ,13:ncol(stand.cores)] <- t(apply(stand.cores[ ,13:ncol(stand.cores)],
#                                          1,
#                                          function(x) x/sum(x, na.rm=TRUE)*500))

# stand.cores[ ,13:ncol(stand.cores)] <- t(apply(stand.cores[ ,13:ncol(stand.cores)],
#                                                1,
#                                                function(x) if (sum(x)==0) {x} else {x/sum(x, na.rm=TRUE)*500}))

# check your work
# all(round(rowSums(stand.cores[, 13:ncol(stand.cores)], na.rm=TRUE)) == 500) 

# create time chunks and add to dataframe
# dat.cut <- cut(stand.cores$age, breaks = c(min(stand.cores$age, na.rm = TRUE), seq(1000, 21000, by = 1000)),
#                labels = FALSE)
# stand.cores$cut <- dat.cut
# stand.cores <- stand.cores[!is.na(stand.cores$cut),]

####
# ASSIGN TIME CHUNKS
####

time_bins = seq(0, 21000, by=990)
N_times = length(time_bins)-1

test.cut <- cut(stand.cores$age, 
                breaks = time_bins)
test <- stand.cores
test$cut <- as.integer(test.cut)
test <- test[!is.na(test$cut),]

# assign unique ID to each distinct x, y coordinate
xyid <- test[,c('x','y')] %>% distinct()
xyid$id <- as.character(seq(1, nrow(xyid), by = 1))

# separate dataframe into list of dataframes, one for each time chunk, remove excess columns
time <- split(test, f = test$cut)

# extract pollen data from specified taxa
taxa.keep = c('Acer', 'Alnus','Betula', 'Cyperaceae', 'Fagus', 'Fraxinus',
              'Ostrya.Carpinus', 'Picea', 'Pinus', 'Quercus', 'Tsuga', 'Ulmus')
taxa.nontree = c('Other', 'Prairie.Forbs', 'Poaceae')

for (i in 1:N_times){
  compiled.meta = time[[i]][,c('x', 'y')]
  compiled.counts = time[[i]][,13:(ncol(time[[i]])-1)]
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

# DO NOT CHANGE NAS TO ZEROS! KEEP THEM AS NA
# change NAs to zeros
# for(i in 1:N_times){
#   time_sum[[i]][,c(4:ncol(time_sum[[i]]))] <- apply(time_sum[[i]][,c(4:ncol(time_sum[[i]]))], c(1,2), 
#                                   function(x) 
#                                     if (is.na(x)){x=0}
#                                   else{x})
# }

# use rounding to get integer counts
counts <- time_sum
for(i in 1:N_times){
  counts[[i]][,c(4:ncol(time_sum[[i]]))] <- round(counts[[i]][,c(4:ncol(time_sum[[i]]))])
}

# extract site coordinates
locs <- counts[[1]][,c('x','y')]
N_locs = nrow(locs)
N_taxa = length(taxa.keep) + 1
locs = cbind(locs, id = seq(1, N_locs))

dat_array = array(0, dim = c(N_locs, N_taxa, N_times))
for (i in 1:N_times){
  for (j in 1:nrow(counts[[i]])){
    site_id = as.numeric(counts[[i]]$id[j])
    dat_array[site_id,,i] = as.numeric(unname(counts[[i]][j,4:(4+N_taxa-1)]))
  }
}

locs <- locs[,c('x','y')]

saveRDS(dat_array, paste0('data/', 'pollen_dat_v', version, '.RDS'))
saveRDS(locs, paste0('data/', 'pollen_locs_v', version, '.RDS'))
