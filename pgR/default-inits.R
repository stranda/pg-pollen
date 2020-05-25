#' Initialized default intial conditions
#'
#'  A function for setting up the default pgLM priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @noRd
default_inits_pgLM <- function(Y, X, priors) {
    
    inits <- list(
        beta = t(mvnfast::rmvn(ncol(Y) - 1, priors$mu_beta, priors$Sigma_beta))
    )
    return(inits)
}


#' Initialized default intial conditions
#'
#'  A function for setting up the default pgLM priors
#'
#' @param Y is a \eqn{n \times d}{n x d} matrix of multinomial count data.
#' @param X is a \eqn{n \times p}{n x p} matrix of variables.
#' @param priors
#' @param shared_theta
#' @param shared_tau
#' @noRd
default_inits_pgSPLM <- function(Y, X, priors, corr_fun = "exponential", shared_theta, shared_tau) {

    J <- ncol(Y)
    theta_mean <- NULL
    theta_var  <- NULL
    if (corr_fun == "matern") {
        theta_mean <- c(priors$mean_range, priors$mean_nu)
        theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)
    } else if (corr_fun == "exponential") {
        theta_mean <- priors$mean_range
        theta_var  <- priors$sd_range^2
    }
    
    inits <- list(
        beta  = t(mvnfast::rmvn(J-1, priors$mu_beta, priors$Sigma_beta)),
        tau2  = if (shared_tau) {
            min(1 / rgamma(1, priors$alpha_tau, priors$beta_tau), 10)
        } else {
            pmin(1 / rgamma(J-1, priors$alpha_tau, priors$beta_tau), 10)
        },
        theta = if (shared_theta) {
            if (corr_fun == "matern") {
                as.vector(pmin(pmax(mvnfast::rmvn(1, theta_mean, theta_var), -1), 0.1))
            } else if (corr_fun == "exponential") {
                pmin(pmax(rnorm(1, theta_mean, sqrt(theta_var)), -1), 0.1)
            }
        } else {
            if (corr_fun == "matern") {
                pmin(pmax(mvnfast::rmvn(J-1, theta_mean, theta_var), -1), 0.1)
            } else if (corr_fun == "exponential") {
                pmin(pmax(rnorm(J-1, theta_mean, sqrt(theta_var)), -1), 0.1)
            }
        }
    )
    return(inits)
}


