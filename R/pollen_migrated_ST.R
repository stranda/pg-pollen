
library(ggplot2)
library(data.table)
library(sp)
library(raster)
library(mapproj)
library(ggplot2)
library(dplyr)
library(tidyr)
version = '4.0'  # v4.0 uses the dataset with corrected lookup table

# Get raster masks and spatial domain you want to use
# (Adam Smith's .tif from: NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks)
# "+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
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
compiled.cores <- read.csv('data/pollen_north_america_ABI_v7.0.csv', stringsAsFactors = FALSE)
load('data/sites_north_america.rdata')
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
# for now just remove rows with age NAs (143 of them as of 17 Sept 2021)
compiled.cores <- compiled.cores[!is.na(compiled.cores$age), ]

# remove any samples with ages greater than 21500 YBP
compiled.cores <- compiled.cores[which(compiled.cores$age < 21500), ]

# only keep pollen from within domain of interest
compiled.cores <- compiled.cores[compiled.cores$lat < lat_hi & 
                            compiled.cores$lat > lat_lo & 
                            compiled.cores$long > long_west & 
                            compiled.cores$long < long_east ,]

# remove pollen cores from isolated locations (islands, mexico)
bermuda <- compiled.cores[compiled.cores$lat < 34 & compiled.cores$long > -78, ]
# bahamas_cuba <- with(compiled.cores, compiled.cores[lat < 28 & long > -80, ])
# mexico <- with(compiled.cores, compiled.cores[lat < 23, ])
# remove.sites <- rbind(bermuda, bahamas_cuba, mexico)
remove.sites <- bermuda %>% distinct()
compiled.cores <- compiled.cores[!(compiled.cores$siteid) %in% remove.sites$siteid,]

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
taxa.nontree <- c('CYPERACEAE')
tree.cores <- compiled.cores[, which(!(colnames(compiled.cores) %in% taxa.nontree))]

# remove sites with no tree pollen counts
# i.e., rows that contain only zeros/NAs across all taxa
start_col <- which(colnames(tree.cores) == 'ageboundyounger') + 1
tree.cores$sum <- apply(tree.cores[,start_col:ncol(tree.cores)], 1, function(x) sum(x, na.rm=TRUE))
tree.cores <- tree.cores[tree.cores$sum > 0, ]
tree.cores <- tree.cores %>% dplyr::select(-sum)

# specify which tree taxa to keep
# use taxa with highest relative proportions
props <- tree.cores[,start_col:ncol(tree.cores)]/rowSums(tree.cores[,start_col:ncol(tree.cores)], na.rm = TRUE)
props <- pivot_longer(props, cols = colnames(props), names_to = 'taxon', values_to = 'props')
props <- props %>% group_by(taxon) %>% summarise(props = sum(props))
props <- props %>% arrange(desc(props))
taxa.keep <- props[1:15, ]
taxa.keep <- as.character(taxa.keep$taxon)

# ASSIGN TIME CHUNKS
# Modern = 1980 AD (to match ENM data) = -30 ybp in pollen data (so add 30 years to age)
# first time chunk: -70 to 705 ybp [represents 210 ybp]
# time 2: 705 - 1695 ybp [1200 ybp] ... etc.
# every 990 years for 22 time chunks (until (20505, 21495), representing 21000 ybp)
tree.cores$age <- tree.cores$age + 30
paleo <- tree.cores[tree.cores$age >= -70, ]
paleo_bins <- c(-70, seq(705, by = 990, length.out = 22))
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
n_times <- length(paleo_time)
for(i in 1:n_times){
  compiled.meta = paleo_time[[i]][,c('x', 'y')]
  compiled.counts = paleo_time[[i]][,start_col:(ncol(paleo_time[[i]])-1)]
  compiled.other = compiled.counts[, which((!(colnames(compiled.counts) %in% taxa.keep)) & 
                                             (!(colnames(compiled.counts) %in% taxa.nontree)))]
  compiled.counts.sub = data.frame(compiled.counts[, which(colnames(compiled.counts) %in% taxa.keep)],
                                   Other = rowSums(compiled.other, na.rm=TRUE))
  paleo_time[[i]] = data.frame(compiled.meta, compiled.counts.sub)
}

for(i in 1:n_times){
  paleo_time[[i]] <- paleo_time[[i]] %>% group_by(x, y) %>% summarise_all(.funs = sum, na.rm = TRUE)
  paleo_time[[i]] <- ungroup(paleo_time[[i]])
  paleo_time[[i]] <- left_join(paleo_xyid, paleo_time[[i]], by = c('x','y'))
  paleo_time[[i]][,4:ncol(paleo_time[[i]])] <- apply(paleo_time[[i]][,4:ncol(paleo_time[[i]])], c(1,2), as.integer)
}
taxa.keep <- colnames(paleo_time[[1]])
taxa.keep <- taxa.keep[!taxa.keep %in% c('x', 'y', 'id')]

# CONSTRUCT DATA LIST, INCLUDING LOCATIONS, INTEGER ARRAY, AND TAXA.KEEP
paleo_locs <- paleo_time[[1]][,c('x','y','id')]
n_locs <- nrow(paleo_locs)
n_times <- length(paleo_time)
n_taxa <- length(taxa.keep)

paleo_dat_array <- array(0, dim = c(n_locs, n_taxa, n_times))
for (i in 1:n_times){
  for (j in 1:n_locs){
    site_id = as.numeric(paleo_time[[i]]$id[j])
    paleo_dat_array[site_id,,i] = as.numeric(unname(paleo_time[[i]][j,4:(4+n_taxa-1)]))
  }
}
paleo_locs <- paleo_locs[,c('x','y')]

saveRDS(paleo_dat_array, paste0('data/', 'pollen_dat_', version, '.RDS'))
saveRDS(paleo_locs, paste0('data/', 'pollen_locs_', version, '.RDS'))
saveRDS(taxa.keep, paste0('data/', 'taxa_', version, '.RDS'))
