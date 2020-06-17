
# NEXT STEPS: create 3D array dataframe so you can test temporal addition to new model
# but this code aggregates counts by site (across all years! rather than subsetting
# for using 'modern' counts)
# So figure out where to stop using old code and where to start writing new code
# This time, after calibration, aggregate by site and time period
# Start with a few time chunks (0-200, 200-400, 400-600 YBP)

library(neotoma)
library(analogue)
library(Bchron)
library(ggplot2)
library(data.table)
library(sp)

# USA Contiguous albers equal area
proj_out <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 
+ellps=GRS80 +datum=NAD83 +units=m +no_defs'

# limits
lat_hi  = 49
lat_lo  = 44
long_west = -100
long_east = -70

# version
version = 'time_1.0'

# get datasets withing bounded region from NEOTOMA
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

# remove sites with no pollen counts (BELOW CODE ISN'T WORKING)
# tree.cores <- tree.cores[rowSums(tree.cores[,13:ncol(tree.cores)] > 0, na.rm = TRUE), ]



tree.cores[,13:35] <- ifelse(is.na(tree.cores[,13:35]), 0, tree.cores[,13:35])


# standardize pollen counts so every site has same total # pollen grains across taxa
# (only need to do this if you're analyzing a subset of all taxa in dataset)
stand.cores <- tree.cores
stand.cores[ ,13:ncol(stand.cores)] <- t(apply(stand.cores[ ,13:ncol(stand.cores)],
                                         1,
                                         function(x) x/sum(x, na.rm=TRUE)*500))
# check your work
all(round(rowSums(stand.cores[, 13:ncol(stand.cores)], na.rm=TRUE)) == 500) 

# create time chunks and add to dataframe
dat.cut <- cut(stand.cores$age, breaks = c(min(stand.cores$age, na.rm = TRUE), seq(1000, 21000, by = 1000)),
               labels = FALSE)
stand.cores$cut <- dat.cut
stand.cores <- stand.cores[!is.na(stand.cores$cut),]



####
# TO TEST TEMPORAL MODEL, FIRST BREAK INTO A FEW SMALLER TIME CHUNKS
####

test.cut <- cut(stand.cores$age, 
                breaks = c(min(stand.cores$age, na.rm = TRUE), 100, 300, 500))
test <- stand.cores
test$cut <- test.cut
test <- test[!is.na(test$cut),]

# separate dataframe into list of dataframes, one for each time chunk, remove excess columns
time <- split(test, f = test$cut)
# time <- lapply(time, "[", c(1,2,4, 13:(ncol(test)-1)))  # USE FOR NON-TEST DATASET
# (below, for testing partial dataset)
time <- lapply(time, "[", c('x','y','Acer','Betula','Ostrya.Carpinus','Ulmus'))

# remove sites that contain NA or 0 pollen counts for each taxon/time chunk
time_sum <- list()
for(i in 1:3){
    time_sum[[i]] <- time[[i]] %>% group_by(x, y) %>% summarise_all(.funs = sum, na.rm = TRUE)
}

# check to see if any sites have 0 pollen counts; if so, remove those sites from the data
# any(rowSums(time_sum[[1]][,3:6]) == 0)

# separate pollen counts from locations
counts <- list()
locs <- list()
for(i in 1:3){
counts[[i]] <- time_sum[[i]][,3:6]
locs[[i]] <- time_sum[[i]][,c('x', 'y')]
}


# extract pollen data from specified taxa
# taxa.keep = c('Acer', 'Alnus','Betula', 'Cyperaceae', 'Fagus', 'Ostrya.Carpinus', 
#               'Picea', 'Pinus', 'Quercus', 'Tsuga', 'Ulmus')
# taxa.nontree = c('Other', 'Prairie.Forbs', 'Poaceae')
# 
# compiled.other = compiled.counts[, which((!(colnames(compiled.counts) %in% taxa.keep)) & 
#                                            (!(colnames(compiled.counts) %in% taxa.nontree)))]
# 
# compiled.counts.sub = data.frame(compiled.counts[, which(colnames(compiled.counts) %in% taxa.keep)], 
#                                  Other = rowSums(compiled.other, na.rm=TRUE))

# calculate proportion from pollen counts

props <- counts[[i]]/rowSums(counts[[i]])


counts.locs[[i]] <- data.frame(locs[[i]], counts[[i]])

dat = as.data.table(compiled.counts.locs)[, lapply(.SD, sum, na.rm=TRUE), by = list(x, y)]
dat[,3:ncol(dat)] = round(dat[,3:ncol(dat)])
N_taxa = ncol(compiled.counts.sub)

saveRDS(dat, paste0('data/', N_taxa, 'taxa_pollen_dat_v', version, '.RDS'))
