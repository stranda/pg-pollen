#' Trace plots from the output of `pg_stlm()` or `pg_stlm_mra()`
#'
#' @param out The output from `pg_stlm()` or `pg_stlm_mra()`
#' @param base_size The base size for the plot
#' @param file If `file = NULL`, the ggplot object is returned. If `file` is not NULL, an image is saved to the file path specified by `file`
#' @param num_plot is the number of locations to plot
#' @param width If a file path is specified, `width` determines the width of the saved image (in inches)
#' @param height If a file path is specified, `height` determines the height of the saved image (in inches)
#'
#' @return Either a ggplot object of the model trace plots (if `file = NULL`) or a saved image file with no return (`file` is not NULL)
#' @export
#'
#' @import patchwork
#' @import tidyverse
#' @import latex2exp
#'

plot_trace_latent <- function(out, base_size = 12, num_plot = 10, file = NULL, width = 16, height = 9) {
    
    # if (!inherits(out, c("pg_stlm", "pg_stlm_mra")))
    #     stop('out must be of class "pg_stlm" or "pg_stlm_mra"')
    # if (!is_positive_numeric(width, 1))
    #     stop("width must be a positive number")
    # if (!is_positive_numeric(height, 1))
    #     stop("height must be a positive number")
    # if (!is_positive_numeric(base_size, 1))
    #     stop("base_size must be a positive number")
    # if (!is.null(file) & !is.character(file))
    #     stop("file must be a character string")
    
    library(latex2exp)
    
    p_eta <- NULL
    
    if (class(out) %in% c("pg_stlm", "pg_stlm_overdispersed", "pg_stlm_latent_overdispersed")) {
        eta_post <- out$eta
        dims <- dim(eta_post)
        dimnames(eta_post) <- list(iteration = 1:dims[1], location = 1:dims[2], species = 1:dims[3], time = 1:dims[4])
        dat_eta <- as.data.frame.table(eta_post, responseName = "eta") %>%
            mutate(iteration = as.numeric(iteration),
                   time      = as.numeric(time))
                                 
        plot_idx <- sample(1:dims[2], num_plot)
        p_eta <- dat_eta %>%
            filter(location %in% plot_idx) %>%
            ggplot(aes(x = .data$iteration, y = .data$eta, group = .data$location, color = .data$location)) +
            geom_line(alpha = 0.5) +
            facet_wrap( ~ time, ncol = 3) +
            scale_color_viridis_d(end = 0.8) +
            ggtitle(TeX("Trace plots for $\\eta$")) +
            theme_bw(base_size = base_size) +
            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
            ylab(TeX("$\\eta$"))
        
    }
    if (class(out) == "pg_stlm_mra") {
        stop('models of class "pg_stlm_mra" are currently not supported')

    }

    p_pi <- NULL
    
    if (class(out) %in% c("pg_stlm", "pg_stlm_overdispersed", "pg_stlm_latent_overdispersed")) {
        pi_post <- out$pi
        dims <- dim(pi_post)
        dimnames(pi_post) <- list(iteration = 1:dims[1], location = 1:dims[2], species = 1:dims[3], time = 1:dims[4])
        dat_pi <- as.data.frame.table(pi_post, responseName = "pi") %>%
            mutate(iteration = as.numeric(iteration),
                   time      = as.numeric(time))
        
        plot_idx <- sample(1:dims[2], num_plot)
        
        p_pi <- dat_pi %>%
            filter(location %in% plot_idx) %>%
            ggplot(aes(x = .data$iteration, y = .data$pi, group = .data$location, color = .data$location)) +
            geom_line(alpha = 0.5) +
            facet_wrap( ~ time, ncol = 3) +
            scale_color_viridis_d(end = 0.8) +
            ggtitle(TeX("Trace plots for $\\pi$")) +
            theme_bw(base_size = base_size) +
            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
            ylab(TeX("$\\pi$"))
        
    }
    
    if (class(out) == "pg_stlm_mra") {
        stop('models of class "pg_stlm_mra" are currently not supported')
        
    }

    
    betas <- out$beta
    dimnames(betas) <- list(
        iteration = 1:dim(out$beta)[1],
        covariate = 1:dim(out$beta)[2],
        species = 1:dim(out$beta)[3]
    )
    p_beta <- as.data.frame.table(betas, responseName = "beta") %>%
    # p_beta  <- data.frame(
    #     beta      = c(out$beta),
    #     iteration = rep(1:dim(out$beta)[1], times = dim(out$beta)[2] * dim(out$beta)[3]),
    #     covariate = factor(rep(1:dim(out$beta)[2], each = dim(out$beta)[1])),
    #     species   = factor(rep(1:dim(out$beta)[3], each = dim(out$beta)[1] * dim(out$beta)[2]))
    # ) %>%
        ggplot(aes(x = .data$iteration, y = .data$beta, group = .data$species, color = .data$species)) +
        geom_line(alpha = 0.75) +
        scale_color_viridis_d(end = 0.8) +
        ggtitle(TeX("Trace plots for $\\beta$")) +
        theme_bw(base_size = base_size) +
        facet_wrap(~ .data$covariate) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ylab(TeX("$\\beta$"))
    
 
    
    if (is.null(file)) {
        if (class(out) %in% c("pg_stlm", "pg_stlm_overdispersed", "pg_stlm_latent_overdispersed")) {
            return(p_eta  + p_pi)
        } else {
            return(p_eta + p_pi)
        }
    } else {
        if (class(out) %in% c("pg_stlm", "pg_stlm_overdispersed", "pg_stlm_latent_overdispersed")) {
            ggsave(filename = file,
                   plot = p_eta / p_pi,
                   device = "png",
                   width = width,
                   height = height,
                   units = "in")
        } else {
            ggsave(filename = file,
                   plot = p_eta  + p_pi,
                   device = "png",
                   width = width,
                   height = height,
                   units = "in")
        }
    }
}
