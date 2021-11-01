\name{plot_trace}
\alias{plot_trace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_trace(out, base_size = 12, file = NULL, width = 16, height = 9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{
%%     ~~Describe \code{out} here~~
}
  \item{base_size}{
%%     ~~Describe \code{base_size} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
  \item{height}{
%%     ~~Describe \code{height} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (out, base_size = 12, file = NULL, width = 16, height = 9) 
{
    p_tau2 <- NULL
    if (class(out) == "pg_stlm") {
        p_tau2 <- data.frame(tau2 = c(out$tau2), iteration = rep(1:nrow(out$tau2)), 
            species = factor(rep(1:ncol(out$tau2), each = nrow(out$tau2)))) \%>\% 
            ggplot(aes(x = .data$iteration, y = .data$tau2, group = .data$species, 
                color = .data$species)) + geom_line(alpha = 0.75) + 
            scale_color_viridis_d(end = 0.8) + ggtitle(TeX("Trace plots for $\\tau2$")) + 
            theme_bw(base_size = base_size) + theme(axis.title.y = element_text(angle = 0, 
            vjust = 0.5)) + ylab(TeX("$\\tau^2$"))
    }
    if (class(out) == "pg_stlm_mra") {
        tau2s <- out$tau2
        dimnames(tau2s) <- list(iteration = 1:dim(tau2s)[1], 
            resolution = 1:dim(tau2s)[2], species = 1:dim(tau2s)[3])
        dat_tau2 <- as.data.frame.table(tau2s, responseName = "tau2")
        p_tau2 <- dat_tau2 \%>\% ggplot(aes(x = .data$iteration, 
            y = .data$tau2, group = .data$resolution, color = .data$resolution)) + 
            geom_line(alpha = 0.75) + scale_color_viridis_d(end = 0.8) + 
            ggtitle(TeX("Trace plots for $\\tau2$")) + theme_bw(base_size = base_size) + 
            facet_wrap(~species) + theme(axis.title.y = element_text(angle = 0, 
            vjust = 0.5)) + ylab(TeX("$\\tau^2$"))
    }
    p_sigma2 <- NULL
    if (class(out) == "pg_stlm_mra") {
        p_sigma2 <- data.frame(sigma2 = c(out$sigma2), iteration = 1:dim(out$sigma2)[1], 
            species = rep(1:dim(out$sigma2)[2], each = dim(out$sigma2)[1])) \%>\% 
            ggplot(aes(x = .data$iteration, y = .data$sigma2)) + 
            geom_line(alpha = 0.75) + scale_color_viridis_d(end = 0.8) + 
            ggtitle(TeX("Trace plots for $\\sigma2$")) + theme_bw(base_size = base_size) + 
            facet_wrap(~species) + theme(axis.title.y = element_text(angle = 0, 
            vjust = 0.5)) + ylab(TeX("$\\sigma^2$"))
    }
    betas <- out$beta
    dimnames(betas) <- list(iteration = 1:dim(out$beta)[1], covariate = 1:dim(out$beta)[2], 
        species = 1:dim(out$beta)[3])
    p_beta <- as.data.frame.table(betas, responseName = "beta") \%>\% 
        ggplot(aes(x = .data$iteration, y = .data$beta, group = .data$species, 
            color = .data$species)) + geom_line(alpha = 0.75) + 
        scale_color_viridis_d(end = 0.8) + ggtitle(TeX("Trace plots for $\\beta$")) + 
        theme_bw(base_size = base_size) + facet_wrap(~.data$covariate) + 
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
        ylab(TeX("$\\beta$"))
    if (class(out) == "pg_stlm") {
        thetas <- exp(out$theta)
        dimnames(thetas) <- list(iteration = 1:dim(thetas)[1], 
            species = 1:dim(thetas)[2], parameter = c("range", 
                "smoothness"))
        p_theta <- as.data.frame.table(thetas, responseName = "theta") \%>\% 
            mutate(iteration = as.numeric(iteration)) \%>\% ggplot(aes(x = .data$iteration, 
            y = .data$theta, group = .data$species, color = .data$species)) + 
            geom_line(alpha = 0.75) + scale_color_viridis_d(end = 0.8) + 
            ggtitle(TeX("Trace plots for $\\theta$")) + theme_bw(base_size = base_size) + 
            facet_wrap(~parameter, nrow = 2, scales = "free") + 
            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
            ylab("")
    }
    p_rho <- NULL
    if (class(out) == "pg_stlm") {
        p_rho <- data.frame(rho = out$rho) \%>\% mutate(iteration = 1:n()) \%>\% 
            ggplot(aes(x = .data$iteration, y = .data$rho)) + 
            geom_line(alpha = 0.75) + ggtitle(TeX("Trace plots for $\\rho$")) + 
            theme_bw(base_size = base_size) + ylab(TeX("$\\rho$"))
    }
    else {
        p_rho <- data.frame(rho = c(out$rho), iteration = 1:dim(out$rho)[1], 
            species = rep(1:dim(out$rho)[2], each = dim(out$rho)[1])) \%>\% 
            ggplot(aes(x = .data$iteration, y = .data$rho)) + 
            geom_line(alpha = 0.75) + ggtitle(TeX("Trace plots for $\\rho$")) + 
            facet_wrap(~species) + theme_bw(base_size = base_size) + 
            ylab(TeX("$\\rho$"))
    }
    if (is.null(file)) {
        if (class(out) == "pg_stlm") {
            return((p_tau2 + p_theta)/(p_beta + p_rho))
        }
        else {
            return((p_tau2 + p_sigma2)/(p_rho + p_beta))
        }
    }
    else {
        if (class(out) == "pg_stlm") {
            ggsave(filename = file, plot = (p_tau2 + p_theta)/(p_beta + 
                p_rho), device = "png", width = width, height = height, 
                units = "in")
        }
        else {
            ggsave(filename = file, plot = (p_tau2 + p_sigma2)/(p_rho + 
                p_beta), device = "png", width = width, height = height, 
                units = "in")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
