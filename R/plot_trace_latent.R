#' Trace plots from the output of `pg_stlm_overdispersed' 'pg_stlm_latent_overdispersed'
#'
#' @param out The output from `pg_stlm()` or `pg_stlm_mra()`
#' @param base_size The base size for the plot
#' @param file If `file = NULL`, the ggplot object is returned. If `file` is not NULL, an image is saved to the file path specified by `file`
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

plot_trace_latent <- function(out, base_size = 12, file = NULL, width = 16, height = 9) {
    
    p_tau2 <- NULL
    

        p_tau2 <- data.frame(
            tau2      = c(out$tau2),
            iteration = rep(1:nrow(out$tau2)),
            species = factor(rep(1:ncol(out$tau2), each = nrow(out$tau2)))
        ) %>%
            ggplot(aes(x = .data$iteration, y = .data$tau2, group = .data$species, color = .data$species)) +
            geom_line(alpha = 0.75) +
            scale_color_viridis_d(end = 0.8) +
            ggtitle("Trace plots for tau2") +
            theme_bw(base_size = base_size) +
            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
            ylab("tau^2")
    
    p_sigma2 <- NULL
    
        p_sigma2 <- data.frame(
            sigma2    = c(out$sigma2),
            iteration = 1:dim(out$sigma2)[1],
            species   = rep(1:dim(out$sigma2)[2], each = dim(out$sigma2)[1])) %>%
            ggplot(aes(x = .data$iteration, y = .data$sigma2)) +
            geom_line(alpha = 0.75) +
            scale_color_viridis_d(end = 0.8) +
            ggtitle("Trace plots for sigma2") +
            theme_bw(base_size = base_size) +
            facet_wrap(~ species) +
            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
            ylab("sigma^2")
    
    
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
        ylab("beta$")
    
        thetas <- exp(out$theta)
        dimnames(thetas) <- list(
            iteration = 1:dim(thetas)[1],
            species = 1:dim(thetas)[2],
            parameter = c("range", "smoothness")
        )
        
        p_theta <- as.data.frame.table(thetas, responseName = "theta") %>%
            mutate(iteration = as.numeric(iteration)) %>%
            ggplot(aes(x = .data$iteration, y = .data$theta, group = .data$species, color = .data$species)) +
            geom_line(alpha = 0.75) +
            scale_color_viridis_d(end = 0.8) +
            ggtitle("Trace plots for thetas") +
            theme_bw(base_size = base_size) +
            facet_wrap(~ parameter, nrow = 2, scales = "free") +
            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
            ylab("Thetas")
    
    
    p_rho <- NULL
    
        p_rho <- data.frame(rho = out$rho) %>%
            mutate(iteration = 1:n()) %>%
            ggplot(aes(x = .data$iteration, y = .data$rho)) +
            geom_line(alpha = 0.75) +
            ggtitle(TeX("Trace plots for $\\rho$")) +
            theme_bw(base_size = base_size) +
            ylab(TeX("$\\rho$"))
    
    if (is.null(file)) {
                return((p_tau2  + p_theta) / (p_beta + p_rho))
    } else {
            ggsave(filename = file,
                   plot = (p_tau2  + p_theta) / (p_beta + p_rho),
                   device = cairo_pdf,
                   width = width,
                   height = height,
                   units = "in")
    }
}
