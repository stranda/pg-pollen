plot_trace <-
function(out, base_size = 12, file = NULL, width = 16, height = 9) {
    
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
    
    p_tau2 <- NULL
    
    if (class(out) == "pg_stlm") {
        p_tau2 <- data.frame(
            tau2      = c(out$tau2),
            iteration = rep(1:nrow(out$tau2)),
            species = factor(rep(1:ncol(out$tau2), each = nrow(out$tau2)))
        ) %>%
            ggplot(aes(x = .data$iteration, y = .data$tau2, group = .data$species, color = .data$species)) +
            geom_line(alpha = 0.75) +
            scale_color_viridis_d(end = 0.8) +
            ggtitle(TeX("Trace plots for $\\tau2$")) +
            theme_bw(base_size = base_size) +
            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
            ylab(TeX("$\\tau^2$"))
    }
    if (class(out) == "pg_stlm_mra") {
        tau2s <- out$tau2
        dimnames(tau2s) <- list(
            iteration  = 1:dim(tau2s)[1],
            resolution = 1:dim(tau2s)[2],
            species    = 1:dim(tau2s)[3]
        )
        dat_tau2 <- as.data.frame.table(tau2s, responseName = "tau2")
        
        
        p_tau2 <- dat_tau2 %>%
            ggplot(aes(x = .data$iteration, y = .data$tau2, group = .data$resolution, color = .data$resolution)) +
            geom_line(alpha = 0.75) +
            scale_color_viridis_d(end = 0.8) +
            ggtitle(TeX("Trace plots for $\\tau2$")) +
            theme_bw(base_size = base_size) +
            facet_wrap(~ species) +
            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
            ylab(TeX("$\\tau^2$"))
    }
    
    p_sigma2 <- NULL
    if (class(out) == "pg_stlm_mra") {
        p_sigma2 <- data.frame(
            sigma2    = c(out$sigma2),
            iteration = 1:dim(out$sigma2)[1],
            species   = rep(1:dim(out$sigma2)[2], each = dim(out$sigma2)[1])) %>%
            ggplot(aes(x = .data$iteration, y = .data$sigma2)) +
            geom_line(alpha = 0.75) +
            scale_color_viridis_d(end = 0.8) +
            ggtitle(TeX("Trace plots for $\\sigma2$")) +
            theme_bw(base_size = base_size) +
            facet_wrap(~ species) +
            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
            ylab(TeX("$\\sigma^2$"))
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
    
    if (class(out) == "pg_stlm") {
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
            ggtitle(TeX("Trace plots for $\\theta$")) +
            theme_bw(base_size = base_size) +
            facet_wrap(~ parameter, nrow = 2, scales = "free") +
            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
            ylab("")
    }
    
    p_rho <- NULL
    if (class(out) == "pg_stlm") {
        p_rho <- data.frame(rho = out$rho) %>%
            mutate(iteration = 1:n()) %>%
            ggplot(aes(x = .data$iteration, y = .data$rho)) +
            geom_line(alpha = 0.75) +
            ggtitle(TeX("Trace plots for $\\rho$")) +
            theme_bw(base_size = base_size) +
            ylab(TeX("$\\rho$"))
    } else {
        p_rho <- data.frame(
            rho       = c(out$rho), 
            iteration = 1:dim(out$rho)[1],
            species   = rep(1:dim(out$rho)[2], each = dim(out$rho)[1])) %>%
            ggplot(aes(x = .data$iteration, y = .data$rho)) +
            geom_line(alpha = 0.75) +
            ggtitle(TeX("Trace plots for $\\rho$")) +
            facet_wrap(~ species) +
            theme_bw(base_size = base_size) +
            ylab(TeX("$\\rho$"))
    }    
    
    if (is.null(file)) {
        if (class(out) == "pg_stlm") {
            return((p_tau2  + p_theta) / (p_beta + p_rho))
        } else {
            return((p_tau2 + p_sigma2) / (p_rho + p_beta))
        }
    } else {
        if (class(out) == "pg_stlm") {
            ggsave(filename = file,
                   plot = (p_tau2  + p_theta) / (p_beta + p_rho),
                   device = "png",
                   width = width,
                   height = height,
                   units = "in")
        } else {
            ggsave(filename = file,
                   plot = (p_tau2 + p_sigma2) / (p_rho + p_beta),
                   device = "png",
                   width = width,
                   height = height,
                   units = "in")
        }
    }
}
