library(neotoma)
library(analogue)
library(Bchron)
library(ggplot2)
library(data.table)

# USA Contiguous albers equal area
proj_out <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 
             +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

# limits
lat_hi  = 49
lat_lo  = 44
long_west = -100
long_east = -70

# version
version = '1.0'

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

# subset to include only samples between 0 and 200 YBP
compiled.cores.sub <- subset(compiled.cores, subset=(age>0)&(age<200))

# visualize subsetted data
map <- map_data("world")
ggplot(data = data.frame(map), aes(long, lat)) + 
  geom_polygon(aes(group=group), color = "steelblue", alpha = 0.2) +
  geom_point(data = compiled.cores.sub,
             aes(x = long, y = lat), color = 2) +
  xlab("Longitude West") + 
  ylab("Latitude North") +
  coord_map(projection = "albers",
            lat0 = lat_lo,
            lat1 = lat_hi, 
            xlim = c(long_west, long_east),
            ylim = c(lat_lo, lat_hi))

# convert projection to proj_out
sp::coordinates(compiled.cores.sub) <- ~long+lat
sp::proj4string(compiled.cores.sub) <- sp::CRS('+init=epsg:4326')
albers <- sp::spTransform(compiled.cores.sub, proj_out)
xy = data.frame(coordinates(albers))
colnames(xy) = c('x', 'y')

# construct data frame with re-projected coordinates
compiled.cores.sub = data.frame(xy, compiled.cores.sub)
compiled.cores.sub = compiled.cores.sub[,which(colnames(compiled.cores.sub)!= 'optional')]

# 
compiled.counts = compiled.cores.sub[,15:ncol(compiled.cores.sub)]
compiled.props  = compiled.cores.sub[,15:ncol(compiled.cores.sub)]/rowSums(compiled.cores.sub[,15:ncol(compiled.cores.sub)])
compiled.locs   = compiled.cores.sub[,c('x', 'y')]

# extract pollen data from specified taxa
taxa.keep = c('Acer', 'Alnus','Betula', 'Cyperaceae', 'Fagus', 'Ostrya.Carpinus', 
              'Picea', 'Pinus', 'Quercus', 'Tsuga', 'Ulmus')
taxa.nontree = c('Other', 'Prairie.Forbs', 'Poaceae')

compiled.other = compiled.counts[, which((!(colnames(compiled.counts) %in% taxa.keep)) & 
                                           (!(colnames(compiled.counts) %in% taxa.nontree)))]

compiled.counts.sub = data.frame(compiled.counts[, which(colnames(compiled.counts) %in% taxa.keep)], 
                 Other = rowSums(compiled.other, na.rm=TRUE))

compiled.props.sub = compiled.counts.sub/rowSums(compiled.counts.sub)
compiled.counts.locs = data.frame(compiled.locs, compiled.counts.sub)

dat = as.data.table(compiled.counts.locs)[, lapply(.SD, sum, na.rm=TRUE), by = list(x, y)]
dat[,3:ncol(dat)] = round(dat[,3:ncol(dat)])
N_taxa = ncol(compiled.counts.sub)

saveRDS(dat, paste0('data/', N_taxa, 'taxa_pollen_dat_v', version, '.RDS'))
