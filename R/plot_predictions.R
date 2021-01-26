library(pgR)
library(ggplot2)
library(fields)
library(rgdal)

run='ST_Aug'
setwd('C:/Users/abrow/Documents/pg-pollen')
out = readRDS(paste0('output/polya-gamma-posts_', run,'.RDS'))
dat = readRDS( paste0('output/polya-gamma-dat_', run,'.RDS'))

#### READ MAP DATA ####
# getting data ready
proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
+lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# WGS84
proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
+towgs84=0,0,0"

na_shp <- readOGR("data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("data/map-data/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj_out)


# note that locations were scaled to fit the model
# unscaling to think in meters, then will rescale again before prediction
rescale = dat$rescale
locs_pollen <- dat$locs*rescale 
names(locs_pollen) <- c("x", "y")

taxa.keep = dat$taxa.keep
# taxa.keep =  as.vector(colnames(dat$y))#[!(colnames(dat) %in% c('x', 'y'))])
# y = as.data.frame(dat$y[,taxa.keep])
y = dat$y
X = dat$X
N_cores = nrow(locs_pollen)

#### DISTANCE MATRICES ####
D_pollen <- fields::rdist(locs_pollen/rescale)# N_cores x N_cores
any(D_pollen == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
D_pollen <- ifelse(D_pollen == 0, 0.007, D_pollen)  # remove them
diag(D_pollen) <- 0

locs_grid = readRDS('data/grid_900locs_ESA.RDS')
X_pred <- matrix(rep(1, nrow(locs_grid)), nrow(locs_grid), 1)
locs = locs_pollen/rescale
locs_pred = locs_grid/rescale

bbox_tran <- function(x, coord_formula = '~ x + y', from, to) {
  
  sp::coordinates(x) <- formula(coord_formula)
  sp::proj4string(x) <- sp::CRS(from)
  bbox <- as.vector(sp::bbox(sp::spTransform(x, CRSobj = sp::CRS(to))))
  return(bbox)
}

grid_box <- bbox_tran(locs_grid, '~ x + y',
                      proj_out,
                      proj_out)
xlim = c(grid_box[1], grid_box[3])
ylim = c(grid_box[2], grid_box[4])

preds = readRDS(paste0('output/polya-gamma-preds_', run,'.RDS'))
library(raster)

pred_iter = preds$pi[1,,,]
pm = reshape2::melt(pred_iter)
colnames(pm) = c('cell', 'taxon', 'time', 'value')
pm$x = locs_pred[pm$cell, 'x']
pm$y = locs_pred[pm$cell, 'y']

taxon = 1
t_vals = unique(pm$time)

x <- stack()
for (t in 1:length(t_vals)){
  pm_sub <- pm[which((pm$taxon==taxon)&(pm$time = t_vals[t])),]
  dfr <- rasterFromXYZ(pm_sub[,c('x', 'y', 'value')]) 
  x <- stack(x, dfr)
}


# START HERE AS OF SEPT 2020
# require(devtools)
# install_github("https://github.com/adamlilith/enmSdm")
require(enmSdm)
require(raster)
run='ST_Aug'
dat <- readRDS( paste0('output/polya-gamma-dat_', run,'.RDS'))
taxa <- dat$taxa.keep

# in hpcc, put all 12 raster stacks in a list
# one stack per taxon; only values for 1st iteration for now
# time chunks are stacked within the raster stacks
stack_list <- readRDS('output/preds_raster_stack_Aug.RDS')

# un-project coordinates before running through bv function
# USA Contiguous albers equal area
proj_out <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 
+ellps=GRS80 +datum=NAD83 +units=m +no_defs'

# specify projection (IS THIS RIGHT? IS IT ALBERS?)
for(i in 1:length(stack_list)){
  sp::proj4string(stack_list[[i]]) <- proj_out
}

bv_list <- list()
# run through bv function
for(i in 1:length(stack_list)){
  bv_list[[i]] <- bioticVelocity(stack_list[[i]])
}

# plot results
par(mfrow = c(3, 4))
for(i in 1:length(bv_list)){
  plot(bv_list[[i]]$timeFrom, bv_list[[i]]$centroidVelocity, type = 'l',
       main = taxa[i], xlab = '1,000 YBP', ylab = 'Centroid velocity (m/yr)')
}

acer <- bv_list[[1]]
acer$centroidLat <- acer$centroidLat * 1e3
acer$centroidLong <- acer$centroidLong * 1e3

ggplot(acer) +
  geom_point(aes(x=centroidLong, y=centroidLat, color = timeFrom), size = 4) +
  # geom_text(aes(x=centroidLong, y=centroidLat, color = timeFrom, label = timeFrom),
  #           size = 8) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group), alpha = 0.5) +
  geom_path(data = lake_shp, aes(x = long, y = lat, group = group), alpha = 0.5) +
  scale_y_continuous(limits = ylim) +
  scale_x_continuous(limits = xlim) +
  theme_classic() +
  # theme(axis.ticks = element_blank(),
  #       axis.text = element_blank(),
  #       axis.title = element_blank(),
  #       line = element_blank(),
  #       legend.title = element_text(size = 16),
  #       legend.text = element_text(size = 14),
  #       plot.title = element_blank()) +
  coord_equal()

###########################################################################3333
pi_mean = apply(preds$pi, c(2,3,4), mean, na.rm=TRUE)
preds_melt = reshape2::melt(pi_mean)#, id.vars=c('x', 'y'))
colnames(preds_melt) = c('loc', 'taxon', 'time', 'value')
preds_melt$x = locs_grid[preds_melt$loc,'x']
preds_melt$y = locs_grid[preds_melt$loc,'y']
preds_melt$taxon = taxa.keep[preds_melt$taxon]

props = apply(y, c(1,3), function(x) if (sum(x)==0){rep(0, length(x))} else if (sum(x)!=0){x/sum(x)})
# props = y/rowSums(y)
# dat_melt = reshape2::melt(data.frame(locs_pollen, props), id.vars=c('x', 'y'))
dat_melt = reshape2::melt(props)#, id.vars=c('x', 'y'))
colnames(dat_melt) = c('taxon', 'loc', 'time', 'value')
dat_melt = dat_melt[, c('loc', 'taxon', 'time', 'value')]
dat_melt$x = locs_pollen[dat_melt$loc,'x']
dat_melt$y = locs_pollen[dat_melt$loc,'y']
dat_melt$taxon = taxa.keep[dat_melt$taxon]

all_melt = rbind(data.frame(dat_melt, type="data"),
                 data.frame(preds_melt, type="preds"))

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
all_melt$value_binned = cut(all_melt$value, breaks, include.lowest=TRUE, labels=FALSE)

breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1,
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

#
# ggplot() +
#   geom_point(data=all_melt, aes(x=x, y=y, color=factor(value_binned), fill=factor(value_binned))) +
#   scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) +
#   scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) +
#   # scale_colour_gradientn(colours = tim.colors(10), limits=c(0,1)) +
#   # scale_fill_gradientn(colours = tim.colors(10), limits=c(0,1)) +
#   facet_grid(variable~type) +
#   geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
#   geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
#   scale_y_continuous(limits = ylim) +
#   scale_x_continuous(limits = xlim) +
#   theme_classic() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         line = element_blank(),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_blank()) +
#   coord_equal()
# #ggsave(paste0("../figs/all_binned_", run, ".png"), device="png", type="cairo")

pdf(paste0("figures/preds_binned_tiled_", run, "_by_time.pdf"))
for (tt in 1:N_times){
  sub_melt = all_melt[which(all_melt$time == tt),]
  p <- ggplot() +
    geom_tile(data=subset(all_melt, type=="preds"), aes(x=x, y=y, color=factor(value_binned), fill=factor(value_binned))) +
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) +
    scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) +
    facet_wrap(.~taxon) +
    geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
    geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(limits = xlim) +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          plot.title = element_blank()) +
    coord_equal()
  print(p)
  # ggsave(paste0("figures/preds_binned_tiled_", run, ".png"), device="png", type="cairo")
}
dev.off()

time_bins = seq(250, 8000, by=500)

times_sub = c(1, 4, 8, 12, 16)
N_times_sub = length(times_sub)

all_sub = all_melt[which(all_melt$taxon %in% c('Fagus', 'Ostrya.Carpinus', 'Picea', 'Pinus')),]
all_sub = all_sub[which(all_sub$time %in% times_sub),]
all_sub$time = time_bins[all_sub$time]
all_sub = all_sub[which(all_sub$type=='preds'),]

pdf(paste0("figures/preds_binned_tiled_", run, "_by_time_subset.pdf"))
# for (tt in 1:N_times_sub){
  # sub_melt = all_sub[which(all_sub$time == times_sub[tt]),]
  p <- ggplot() +
    geom_tile(data=all_sub, aes(x=x, y=y, color=factor(value_binned), fill=factor(value_binned))) +
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) +
    scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) +
    facet_grid(time~taxon) +
    geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
    geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(limits = xlim) +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          plot.title = element_blank()) +
    coord_equal()
  print(p)
  # ggsave(paste0("figures/preds_binned_tiled_", run, ".png"), device="png", type="cairo")
# }
dev.off()


pdf(paste0("figures/preds_binned_tiled_", run, "_by_taxon.pdf"))
# for (tt in 1:N_times){
#   sub_melt = all_melt[which(all_melt$time == tt),]
  p <- ggplot() +
    geom_tile(data=subset(all_melt, type=="preds"), aes(x=x, y=y, color=factor(value_binned), fill=factor(value_binned))) +
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) +
    scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) +
    facet_grid(time~taxon) +
    geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
    geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(limits = xlim) +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          plot.title = element_blank()) +
    coord_equal()
  print(p)
  # ggsave(paste0("figures/preds_binned_tiled_", run, ".png"), device="png", type="cairo")
# }
dev.off()
#
# ggplot() +
#   geom_tile(data=all_melt, aes(x=x, y=y, color=factor(value_binned), fill=factor(value_binned))) +
#   scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) +
#   scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) +
#   # scale_colour_gradientn(colours = tim.colors(10), limits=c(0,1)) +
#   # scale_fill_gradientn(colours = tim.colors(10), limits=c(0,1)) +
#   facet_grid(variable~type) +
#   geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
#   geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
#   scale_y_continuous(limits = ylim) +
#   scale_x_continuous(limits = xlim) +
#   theme_classic() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         line = element_blank(),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_blank()) +
#   coord_equal()
# ggsave(paste0("figures/all_binned_tiled_", run, ".png"), device="png", type="cairo")