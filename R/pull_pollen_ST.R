
library(neotoma)
library(analogue)
library(Bchron)
library(ggplot2)
library(data.table)
library(sp)
library(dplyr)
library(raster)

setwd('C:/Users/abrow/Documents/pg-pollen')
version = 'Dec15'

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

# FIX CHRONOLOGIES THAT YIELD FUTURE AGES - adjust 'age' and 'age.type' (there are 5 of these)
# Exception for dataset 328: use "COHMAP chron 1"
id <- '328'
all.downloads[[id]]$sample.meta$age <- all.downloads[[id]]$chronologies[['COHMAP chron 1']]$age
all.downloads[[id]]$sample.meta$age.type <- all.downloads[[id]]$chronologies[['COHMAP chron 1']]$age.type

# Exception for dataset 1002: use "NAPD 1"
id <- '1002'
all.downloads[[id]]$sample.meta$age <- all.downloads[[id]]$chronologies[['NAPD 1']]$age
all.downloads[[id]]$sample.meta$age.type <- all.downloads[[id]]$chronologies[['NAPD 1']]$age.type

# Exception for dataset 1984: use "NAPD 1"
id <- '1984'
all.downloads[[id]]$sample.meta$age <- all.downloads[[id]]$chronologies[['NAPD 1']]$age
all.downloads[[id]]$sample.meta$age.type <- all.downloads[[id]]$chronologies[['NAPD 1']]$age.type

# Exception for dataset 13051: use "Neotoma 1"
id <- '13051'
all.downloads[[id]]$sample.meta$age <- all.downloads[[id]]$chronologies[['Neotoma 1']]$age
all.downloads[[id]]$sample.meta$age.type <- all.downloads[[id]]$chronologies[['Neotoma 1']]$age.type

# Exception for dataset 17404: use "Neotoma 1"
id <- '17404'
all.downloads[[id]]$sample.meta$age <- all.downloads[[id]]$chronologies[['Neotoma 1']]$age
all.downloads[[id]]$sample.meta$age.type <- all.downloads[[id]]$chronologies[['Neotoma 1']]$age.type



####
# START HERE IF YOU ALREADY HAVE POLLEN DATA DOWNLOADED
####

all.downloads <- readRDS('data/all.downloads.RDS')

# standardize the pollen taxonomy
all.cores <- compile_taxa(all.downloads, "P25")

# convert list to data frame
compiled.cores  <- compile_downloads(all.cores)

# until there's an answer as to why we get age NAs after running 'compile_downloads'...
# for now just remove rows with age NAs
compiled.cores <- compiled.cores[!is.na(compiled.cores$age), ]

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
             aes(x = long, y = lat), color = 2, size = 2) +
  xlab("Longitude West") + 
  ylab("Latitude North") +
  coord_map(projection = "albers",
            lat0 = lat_lo,
            lat1 = lat_hi, 
            xlim = c(long_west, -59),
            ylim = c(lat_lo, lat_hi))

# # remove pollen cores from islands far off mainland
# remove.sites <- compiled.cores[compiled.cores$lat <= 35 & compiled.cores$long > -70,]
# compiled.cores <- compiled.cores[!(row.names(compiled.cores) %in% row.names(remove.sites)),]

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
taxa.nontree <- c('Other', 'Prairie.Forbs', 'Poaceae')
tree.cores <- compiled.cores[, which(!(colnames(compiled.cores) %in% taxa.nontree))]

# remove sites with no tree pollen counts
# i.e., rows that contain only zeros/NAs across all taxa
tree.cores$sum <- apply(tree.cores[,13:35], 1, function(x) sum(x, na.rm=TRUE))
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

summ <- data.frame(apply(tree.cores[,13:35], 2, function(x) sum(x, na.rm=TRUE)))
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
