library(pgR)
source('pgR/predict-pgSPLM.R')

out = readRDS(paste0('output/polya-gamma-posts_pgR_', run,'.RDS'))
dat = readRDS( paste0('output/polya-gamma-dat_pgR_', run,'.RDS'))

# note that locations were scaled to fit the model
# unscaling to think in meters, then will rescale again before prediction
rescale = dat$rescale
locs_pollen <- dat$locs*rescale 
names(locs_pollen) <- c("x", "y")
taxa.keep =  as.vector(colnames(dat$y))#[!(colnames(dat) %in% c('x', 'y'))])
y = as.data.frame(dat$y[,taxa.keep])
X = dat$X
N_cores = nrow(locs_pollen)

# 
# pol_box <- bbox_tran(locs_pollen, '~ x + y',
#                      proj_out,
#                      proj_out)
# xlim = c(pol_box[1]-24000, pol_box[3]+24000)
# ylim = c(pol_box[2]-24000, pol_box[4]+24000)


#### DISTANCE MATRICES ####
D_pollen <- fields::rdist(locs_pollen/rescale)# N_cores x N_cores
any(D_pollen == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
D_pollen <- ifelse(D_pollen == 0, 0.007, D_pollen)  # remove them
diag(D_pollen) <- 0

locs_grid = readRDS('data/grid.RDS')




X_pred <- matrix(rep(1, nrow(locs_grid)), nrow(locs_grid), 1)

locs = locs_pollen/rescale
locs_pred = locs_grid/rescale

preds = predict_pgSPLM(
  out,
  X,
  X_pred,
  locs = locs,
  locs_pred = locs_pred,
  corr_fun="matern",
  # shared_covariance_params,
  shared_tau = TRUE,
  shared_theta = FALSE,
  diag_adjust = 1.e-8,
  n_cores = 1L,
  progress = TRUE
)

pi_mean = apply(preds$pi, c(2,3), mean, na.rm=TRUE)
colnames(pi_mean) = taxa.keep

preds = data.frame(locs_pred*rescale, pi_mean)
preds_melt = reshape2::melt(preds, id.vars=c('x', 'y'))

props = y/rowSums(y)
dat_melt = reshape2::melt(data.frame(locs_pollen, props), id.vars=c('x', 'y'))

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

ggplot() + 
  geom_tile(data=subset(all_melt, type=="preds"), aes(x=x, y=y, color=factor(value_binned), fill=factor(value_binned))) + 
  scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
  scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
  # scale_colour_gradientn(colours = tim.colors(10), limits=c(0,1)) + 
  # scale_fill_gradientn(colours = tim.colors(10), limits=c(0,1)) + 
  facet_wrap(.~variable) + 
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
ggsave(paste0("figures/preds_binned_tiled_", run, ".png"), device="png", type="cairo")

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