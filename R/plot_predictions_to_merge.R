

#### MAP FRAXINUS ESTIMATES BY TIME - MEAN/SD FROM VERSION 3.1
require(ggplot2)
require(gridExtra)
require(ggforce)
require(dplyr)
require(rgdal)

# read mapping data
proj <- '+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs'
na_shp <- readOGR("data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("data/map-data/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj)

# read prediction output
mean <- readRDS('output/preds_frax_mean_4.0.RDS')
sd <- readRDS('output/preds_frax_sd_4.0.RDS')
grid <- readRDS('data/grid_3.1.RDS')
n_times <- dim(mean)[2]

# create prediction dataframe, combined mean/SD estimates
mean_list <- list()
for(i in 1:n_times){
  mean_list[[i]] <- data.frame(grid, mean[,i])
  names(mean_list[[i]]) <- c('x','y','pred')
  mean_list[[i]]$time <- i
  mean_list[[i]]$type <- 'mean'
}
mean_df <- bind_rows(mean_list)

sd_list <- list()
for(i in 1:n_times){
  sd_list[[i]] <- data.frame(grid, sd[,i])
  names(sd_list[[i]]) <- c('x','y','pred')
  sd_list[[i]]$time <- i
  sd_list[[i]]$type <- 'sd'
}
sd_df <- bind_rows(sd_list)
df <- rbind(mean_df, sd_df)

# specify times bins
bins <- c(-70, seq(705, by = 990, length.out = 22))
bins_df <- data.frame(timeTo = bins[-23], timeFrom = bins[-1])
bins_df$timeframe <- paste0(bins_df$timeFrom, ' - ', bins_df$timeTo)
bins_df$time <- 1:n_times

# make data breaks
breaks_vec <- c(0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.25)
breaks <- cut(df$pred, breaks = breaks_vec, include.lowest = TRUE, labels = FALSE)
df$breaks <- breaks
breaks_vec <- breaks_vec[-1]
df <- left_join(df, bins_df, by = 'time')

# subset df using time bins you want to map
# times_use <- which(bins_df$mid %in% c(400, 4400, 8400, 12400, 16400, 20400))
# df_sub <- df[df$time %in% times_use, ]

xlim <- c(min(df$x), max(df$x))
ylim <- c(min(df$y), max(df$y))

df_sub <- df[df$time < 3, ]
time_labs <- as.character(bins_df$timeframe)
names(time_labs) <- bins_df$time

plot_list <- list()
for(i in 1:11){
  plot_list[[i]] <- ggplot(data = df) +
    geom_tile(aes(x = x, y = y, fill = as.factor(breaks))) +
    scale_fill_viridis_d(labels = breaks_vec) +
    geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
              alpha = 0.2, color = 'white') +
    geom_path(data = lake_shp, aes(x = long, y = lat, group = group), 
              alpha = 0.3, color = 'white') +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(limits = xlim) +
    labs(fill = 'Rel abund \nprediction') +
    facet_wrap_paginate(time ~ type, nrow = 2, ncol = 2, page = i,
                        labeller = labeller(time = time_labs)) +
    # facet_wrap(time ~ type, labeller = labeller(time = time_labs)) +
    theme_classic() +
    theme(
      strip.text = element_text(size = 10),
      strip.background = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      line = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)) +
    # legend.position = 'none') +
    coord_equal()
}

plots <- marrangeGrob(plot_list, nrow = 1, ncol = 1)  # takes several minutes to run
ggsave('figures/frax_preds_mean_sd_4.0.pdf', plots,
       height = 8, width = 7.5, units = 'in')



#### MEAN PREDICTIONS OVER TIME FOR ALL TAXA - LATENT OVERDISPERSED MODEL ####
setwd('C:/Users/abrow/Documents/pg-pollen')
require(ggplot2)
require(ggforce)
require(gridExtra)
require(viridis)
require(dplyr)
require(raster)
require(rgdal)

# read prediction output
mean <- readRDS('output/preds_all_taxa_mean_4.0.RDS')
sd <- readRDS('output/preds_all_taxa_sd_4.0.RDS')
grid <- readRDS('data/grid_3.1.RDS')
taxa <- readRDS('data/taxa_4.0.RDS')
n_times <- dim(mean)[3]
n_taxa <- dim(mean)[2]

# create prediction dataframe
mean_list <- rep(list(list()), n_taxa)
for(i in 1:n_taxa){
  for(j in 1:n_times){
    mean_list[[i]][[j]] <- data.frame(grid, mean[ , i, j])
    names(mean_list[[i]][[j]]) <- c('x','y','pred')
    mean_list[[i]][[j]]$time <- j
    mean_list[[i]][[j]]$taxa <- taxa[i]
    mean_list[[i]][[j]]$type <- 'mean'
  }
  mean_list[[i]] <- bind_rows(mean_list[[i]])
}
mean_df <- bind_rows(mean_list)

sd_list <- rep(list(list()), n_taxa)
for(i in 1:n_taxa){
  for(j in 1:n_times){
    sd_list[[i]][[j]] <- data.frame(grid, sd[ , i, j])
    names(sd_list[[i]][[j]]) <- c('x','y','pred')
    sd_list[[i]][[j]]$time <- j
    sd_list[[i]][[j]]$taxa <- taxa[i]
    sd_list[[i]][[j]]$type <- 'sd'
  }
  sd_list[[i]] <- bind_rows(sd_list[[i]])
}
sd_df <- bind_rows(sd_list)
df <- rbind(mean_df, sd_df)

# specify times bins
bins <- c(-70, seq(705, by = 990, length.out = 22))
bins_df <- data.frame(timeTo = bins[-23], timeFrom = bins[-1])
bins_df$timeframe <- paste0(bins_df$timeFrom, ' - ', bins_df$timeTo)
bins_df$time <- 1:n_times

# north america shapefiles
proj <- '+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs'
na_shp <- readOGR("data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("data/map-data/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj)

xlim <- c(min(df$x), max(df$x))
ylim <- c(min(df$y), max(df$y))

time_labs <- as.character(bins_df$timeframe)
names(time_labs) <- bins_df$time

## PLOT MEAN/SD OF PREDICTIONS - ALL TAXA, ALL TIMES - SEPARATE PDF FOR EACH TAXON
plot_list <- rep(list(list()), n_taxa)

for(i in 1:n_taxa){
  df_sub <- df[df$taxa == taxa[i], ]
  breaks_vec <- quantile(df_sub$pred, probs = seq(0, 1, 0.1))
  breaks <- cut(df_sub$pred, breaks = breaks_vec, include.lowest = TRUE, labels = FALSE)
  df_sub$breaks <- breaks
  breaks_vec <- round(unname(breaks_vec[-1]), digits = 5)

  for(j in 1:11){
    plot_list[[i]][[j]] <- ggplot(data = df_sub) +
      geom_tile(aes(x = x, y = y, fill = factor(breaks))) +
      scale_fill_viridis_d(labels = breaks_vec) +
      # scale_fill_viridis() +
      geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
                alpha = 0.2, color = 'white') +
      geom_path(data = lake_shp, aes(x = long, y = lat, group = group), 
                alpha = 0.3, color = 'white') +
      scale_y_continuous(limits = ylim) +
      scale_x_continuous(limits = xlim) +
      labs(fill = 'Rel abund \nprediction') +
      facet_wrap_paginate(time ~ type, nrow = 2, ncol = 2, page = j,
                          labeller = labeller(time = time_labs)) +
      # facet_wrap(~ time, labeller = labeller(time = time_labs)) +
      ggtitle(label = taxa[i]) + 
      theme_classic() +
      theme(
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
      # legend.position = 'none') +
      coord_equal()
  }
}
# below code returns error: 
# Error: `device` must be NULL, a string or a function.
for(i in 1:n_taxa){
  plot_i <- marrangeGrob(plot_list[[i]], nrow = 1, ncol = 1)
  ggsave(paste0('figures/', taxa[i], '_preds_mean_sd_4.0.pdf', plot_i))
}

# creating and saving PDF for one taxon at a time
taxa
# [1] "PICEA"           "PINUS"           "BETULA"          "ALNUS"           "SALIX"          
# [6] "ULMUS"           "CARYA"           "FRAXINUS"        "ABIES"           "TSUGA"          
# [11] "CUPRESSACEAE"    "OSTRYA.CARPINUS" "ACER"            "FAGUS"           "QUERCUS"        
# [16] "Other"
plot <- marrangeGrob(plot_list[[3]], nrow = 1, ncol = 1)
ggsave('figures/betula_preds_mean_sd_4.0.pdf', plot, 
       width = 6, height = 7, units = 'in')



## PLOTTING *ONLY* MEAN PREDICTIONS - ALL TAXA, ALL TIMES
plot_list <- list()
for(i in 1:n_taxa){
  df_sub <- df[df$taxa == taxa[i], ]
  
  breaks_vec <- quantile(df_sub$pred, probs = seq(0, 1, 0.1))
  breaks <- cut(df_sub$pred, breaks = breaks_vec, include.lowest = TRUE, labels = FALSE)
  df_sub$breaks <- breaks
  breaks_vec <- round(unname(breaks_vec[-1]), digits = 5)
  
  plot_list[[i]] <- ggplot(data = df_sub) +
    geom_tile(aes(x = x, y = y, fill = factor(breaks))) +
    scale_fill_viridis_d(labels = breaks_vec) +
    # scale_fill_viridis() +
    geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
              alpha = 0.2, color = 'white') +
    geom_path(data = lake_shp, aes(x = long, y = lat, group = group), 
              alpha = 0.3, color = 'white') +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(limits = xlim) +
    labs(fill = 'Rel abund \nprediction') +
    facet_wrap_paginate(~ time, page = i, labeller = labeller(time = time_labs)) +
    # facet_wrap(~ time, labeller = labeller(time = time_labs)) +
    ggtitle(label = taxa[i]) + 
    theme_classic() +
    theme(
      strip.text = element_text(size = 10),
      strip.background = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      line = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)) +
    # legend.position = 'none') +
    coord_equal()
}
plots <- marrangeGrob(plot_list, nrow = 1, ncol = 1)  # takes several minutes to run
ggsave('figures/all_taxa_preds_mean_4.0.pdf', plots)
# height = 8, width = 7.5, units = 'in')



## PLOTTING RAW OBSERVATIONS FOR ALL TAXA
require(rgdal)
require(ggplot2)
require(dplyr)
require(viridis)
require(gridExtra)

dat <- readRDS('data/pollen_dat_Other4.0.RDS')
locs <- readRDS('data/pollen_locs_4.0.RDS')
taxa <- readRDS('data/taxa_Other4.0.RDS')
n_taxa <- dim(dat)[2]
n_times <- dim(dat)[3]

dat_list <- rep(list(list()), n_taxa)
for(i in 1:n_taxa){
  for(j in 1:n_times){
    dat_list[[i]][[j]] <- data.frame(locs, dat[, i, j])
    names(dat_list[[i]][[j]]) <- c('x','y','count')
    dat_list[[i]][[j]]$time <- j
    dat_list[[i]][[j]]$taxa <- taxa[i]
  }
  dat_list[[i]] <- bind_rows(dat_list[[i]])
}
dat <- bind_rows(dat_list)

# convert counts to relative proportions
loc_summ <- dat %>% 
  group_by(x, y, time) %>% 
  summarize(across(count, ~ sum(count, na.rm = TRUE)))
names(loc_summ)[4] <- 'total_count'

dat <- left_join(dat, loc_summ, by = c('x','y','time'))
dat$rel_prop <- ifelse(is.na(dat$count), NA, dat$count / dat$total_count)
saveRDS(dat, 'output/obs_rel_prop_all_taxa.RDS')

# specify time labels for plotting
bins <- c(-70, seq(705, by = 990, length.out = 22))
bins_df <- data.frame(timeTo = bins[-23], timeFrom = bins[-1])
bins_df$time <- paste0(bins_df$timeFrom, ' - ', bins_df$timeTo)
bins_df$time_id <- 1:n_times
time_labs <- as.character(bins_df$time)
names(time_labs) <- bins_df$time_id

# north america shapefiles
proj <- '+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs'
na_shp <- readOGR("data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("data/map-data/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj)

xlim <- c(min(dat$x), max(dat$x))
ylim <- c(min(dat$y), max(dat$y))

plot_list <- list()
for(i in 1:n_taxa){
  dat_sub <- dat[dat$taxa == taxa[i], ]
  
  breaks_vec <- quantile(dat_sub$rel_prop, probs = seq(0, 1, 0.1), na.rm = TRUE)
  breaks_rem <- max(which(breaks_vec == 0))
  breaks_new <- breaks_vec[breaks_rem:length(breaks_vec)]
  breaks <- cut(dat_sub$rel_prop, breaks = breaks_new, include.lowest = TRUE, labels = FALSE)
  dat_sub$breaks <- breaks
  breaks_vec <- round(unname(breaks_new[-1]), digits = 5)
  
  plot_list[[i]] <- ggplot(subset(dat_sub, !is.na(breaks))) +
    geom_path(data = cont_shp, aes(x = long, y = lat, group = group), 
              alpha = 0.2) +
    geom_path(data = lake_shp, aes(x = long, y = lat, group = group), 
              alpha = 0.3) +
    geom_point(aes(x = x, y = y, color = factor(breaks)), alpha = 0.7) +
    scale_color_viridis_d(labels = breaks_vec) +
    # scale_color_viridis(trans = 'reverse') +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(limits = xlim) +
    labs(color = 'Relative\nproportion') +
    facet_wrap(~ time, labeller = labeller(time = time_labs)) +
    ggtitle(taxa[i]) +
    theme_classic() +
    theme(
      strip.text = element_text(size = 10),
      strip.background = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      line = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 13, hjust = 0.5)) +
    coord_equal()
}

plots <- marrangeGrob(plot_list, nrow = 1, ncol = 1)  # takes several minutes to run
ggsave('figures/all_taxa_observed_4.0.pdf', plots, 
       width = 12, height = 9, units = 'in')
