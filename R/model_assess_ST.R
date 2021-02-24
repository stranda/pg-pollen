library(reshape2)
library(rdist)
library(rgeos)
library(ggplot2)
library(sp)
library(rgdal)
library(raster)
library(fields)

run='Oct12_matern'

# dir.create(file.path('figures'), showWarnings = FALSE)

# get bounding box
bbox_tran <- function(x, coord_formula = '~ x + y', from, to) {
  sp::coordinates(x) <- formula(coord_formula)
  sp::proj4string(x) <- sp::CRS(from)
  bbox <- as.vector(sp::bbox(sp::spTransform(x, CRSobj = sp::CRS(to))))
  return(bbox)
}

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

#### READ IN MODEL DATA AND OUTPUT ####
out = readRDS(paste0('output/polya-gamma-posts_', run,'.RDS'))
dat = readRDS( paste0('output/polya-gamma-dat_', run,'.RDS'))

# note that locations were scaled to fit the model
# unscaling to think in meters, then will rescale again before prediction
rescale = dat$rescale
locs_pollen <- dat$locs*rescale 
names(locs_pollen) <- c("x", "y")
taxa.keep = dat$taxa.keep
N_cores = nrow(locs_pollen)

pol_box <- bbox_tran(locs_pollen, '~ x + y',
                     proj_out,
                     proj_out)
xlim = c(pol_box[1]-24000, pol_box[3]+24000)
ylim = c(pol_box[2]-24000, pol_box[4]+24000)


#### DISTANCE MATRICES ####
D_pollen <- fields::rdist(locs_pollen/rescale)# N_cores x N_cores
any(D_pollen == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
D_pollen <- ifelse(D_pollen == 0, 0.007, D_pollen)  # remove them
diag(D_pollen) <- 0


#### PREDICTIONS ####
N_iter = dim(out$tau2)[1]
J = dim(out$eta)[3] + 1
N_times = dim(out$eta)[4]
burn = 0
N_keep = N_iter-burn#+1

tau2  = out$tau2[burn:N_iter,]
theta = out$theta[burn:N_iter,,]
# omega = out$omega[burn:N_iter,,]
eta   = out$eta[burn:N_iter,,,]
# mu    = out$mu[burn:N_iter,,]
beta  = out$beta[burn:N_iter,,]

## calculate the Matern correlation using parameters theta on the log scale
correlation_function <- function(D, theta) {
  geoR::matern(D, exp(theta[1]), exp(theta[2]))
}

x = seq(0, max(D_pollen), length=1000)
taxa = taxa.keep
cov_df = data.frame(taxon=character(0),
                    distance=numeric(0),
                    covar=numeric(0))
for (j in 1:(J-1)){
  cov_df = rbind(cov_df, 
                 data.frame(taxon = rep(taxa[j], length(x)),
                            distance = x, 
                            # covar = mean(tau)^2*geoR::matern(x, exp(colMeans(out$theta[j,,]))[1], exp(colMeans(out$theta[j,,]))[2])))
                            covar = geoR::matern(x, exp(colMeans(out$theta[,j,]))[1], 
                                                 exp(colMeans(out$theta[,j,]))[2])))

}

ggplot(data=cov_df) + 
  geom_line(aes(x=distance, y=covar, color=taxon), size=2) +
  theme_bw() +
  xlim(c(0,500)) + 
  xlab("Distance (km)") + 
  ylab("Covariance")
ggsave(paste0("figures/covariance_vs_distance_", run, ".pdf"))#, device="pdf", type="cairo")


###############################################################################################################################
## trace
###############################################################################################################################

# tau
tau_melt = reshape2::melt(tau2)
colnames(tau_melt) = c('iter', 'taxon', 'value')
taxa <- readRDS('data/pollen_taxa_Dec15.RDS')
taxa <- data.frame(taxa.id = 1:13, taxa = taxa)
tau_melt <- merge(tau_melt, taxa, by.x = 'taxon', by.y = 'taxa.id')
ggplot() + geom_line(data=tau_melt, aes(x=iter, y=value, color=factor(taxa)))
ggsave(paste0("figures/trace_tau_", run, ".png"), device="png", type="cairo")

# theta
theta_melt = reshape2::melt(theta)
colnames(theta_melt) = c('iter', 'taxon', 'number', 'value')
theta_melt <- merge(theta_melt, taxa, by.x = 'taxon', by.y = 'taxa.id')

ggplot(data=theta_melt) + 
  geom_line(aes(x=iter, y=exp(value), color=factor(taxa))) +
  theme_bw() +
  facet_grid(number~., scales="free_y")
ggsave(paste0("figures/trace_theta_", run, ".png"), device="png", type="cairo")

# mu
mu_melt = reshape2::melt(beta)
colnames(mu_melt) = c('iter', 'taxon', 'value')
mu_melt$taxon = taxa[mu_melt$taxon]
mu_melt <- merge(mu_melt, taxa, by.x = 'taxon', by.y = 'taxa.id')

ggplot(data=mu_melt) + 
  geom_line(aes(x=iter, y=value, color=taxa)) +
  theme_bw()
ggsave(paste0("figures/trace_mu_", run, ".png"), device="png", type="cairo")

###############################################################################################################################
## maps
###############################################################################################################################

# function to convert eta to pi (proportions)
expit <- function(x) {
  1 / (1 + exp(-x))
}

eta_to_pi <- function(eta) {
  ## convert eta to a probability vector pi
  ## can make this more general by first checking if vector vs. matrix and then
  ## calculating the response
  N <- nrow(eta)
  J <- ncol(eta) + 1
  pi <- matrix(0, N, J)
  stick <- rep(1, N)
  for (j in 1:(J - 1)) {
    pi[, j] <- expit(eta[, j]) * stick
    stick <- stick - pi[, j]
  }
  pi[, J] <- stick
  return(pi)
}

pis = array(NA, dim=c(N_cores, J, N_times, N_keep))
for (tt in 1:N_times){
  for (i in 1:N_keep){
    pis[,,tt,i] <- eta_to_pi(eta[i,,,tt])
    # pis <- pis %>% mutate(sum = rowSums(.))  # check to make sure it worked
  }
}

# GET MEAN OF PIS FOR MAPPING
# NEXT, GET SD OF PIS FOR MAPPING UNCERTAINTY
pi_mean = apply(pis, c(1,2,3), mean, na.rm=TRUE)
preds_melt = reshape2::melt(pi_mean)#, id.vars=c('x', 'y'))
colnames(preds_melt) = c('loc', 'taxon', 'time', 'value')
preds_melt$x = locs_pollen[preds_melt$loc,'x']
preds_melt$y = locs_pollen[preds_melt$loc,'y']
preds_melt$taxon = taxa.keep[preds_melt$taxon]

# calculate proportions from raw data
y <- dat$y
props = apply(X = y, MARGIN = c(1,3), FUN = function(x) 
  if (is.na(sum(x))){rep(NA, length(x))}
  else if (sum(x, na.rm = TRUE)==0){rep(0, length(x))}
  else if (sum(x, na.rm = TRUE)!=0){x/sum(x, na.rm = TRUE)})

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

some <- c('Picea', 'Betula', 'Fraxinus', 'Quercus')
all_melt_plot <- all_melt[all_melt$taxon %in% some, ]


pdf(paste0("figures/all_binned_by_time_Oct12_4taxa_4times.pdf"))
for (tt in seq(1, N_times, by = 6)){
  sub_melt = all_melt_plot[which(all_melt_plot$time == tt),]
  p <- ggplot() + 
    geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
    geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
    
    geom_point(data=sub_melt, aes(x=x, y=y, color=factor(value_binned), fill=factor(value_binned),
                                  alpha = 0.5)) + 
    scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
    scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
    facet_grid(taxon~type) + 
    
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(limits = xlim) +
    ggtitle("Time: ", tt) +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) + #,
          # plot.title = element_blank()) +
    coord_equal()
  print(p)
  # ggsave(paste0("../figs/all_binned_", run, "_time", tt, ".png"), device="png", type="cairo")
}
dev.off()

png(paste0("figures/all_binned_tiled_", run, "_by_time.png"))
for (tt in 1:N_times){
  sub_melt = all_melt[which(all_melt$time == tt),]
p <- ggplot() + 
  geom_tile(data=all_melt, aes(x=x, y=y, color=factor(value_binned), fill=factor(value_binned))) + 
  scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
  scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
  # scale_colour_gradientn(colours = tim.colors(10), limits=c(0,1)) + 
  # scale_fill_gradientn(colours = tim.colors(10), limits=c(0,1)) + 
  facet_grid(taxon~type) + 
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
# ggsave(paste0("figures/all_binned_tiled_", run, ".png"), device="png", type="cairo")
}
dev.off()

###############################################################################################################################
## observed versus predicted
###############################################################################################################################

# foo = merge(preds_melt, dat_melt, by=c("x", "y", "variable"))
foo = merge(preds_melt, dat_melt, by=c("x", "y", "loc", "time", "taxon"))

ggplot(data=foo) +
  geom_point(aes(x=value.y, y=value.x)) +
  theme_bw() + 
  coord_equal() +
  xlim(c(0,1)) + 
  ylim(c(0,1)) +
  xlab("Observed proportions") +
  ylab("Predicted proportions") +
  geom_abline(intercept=0, slope=1) + 
  facet_wrap(~taxon)
ggsave(paste0("figures/props_obs_vs_preds_", run, ".png"), device="png", type="cairo")
