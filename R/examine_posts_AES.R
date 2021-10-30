library(tidyverse)
library(pgR)
library(latex2exp)
source("R/plot-trace.R")

version <- "4.0"

out_latent <- readRDS(here::here('output', paste0('polya-gamma-posts_', version, '_latent_overdispersed.RDS')))
plot_trace_latent(out_latent, base_size = 7, file = "figures/trace-plots_latent.pdf")


