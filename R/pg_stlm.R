#' Bayesian Polya-gamma regression
#' 
#' this function runs the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param Y is a \eqn{N \times J \times T}{N x J x T} array of compositional count data.
#' @param X is a \eqn{N \times p}{n_sites x p} matrix of climate variables.
#' @param locs is a \eqn{n_sites \times 2}{n_sites x 2} matrix of observation locations.
#' @param params is a list of parameter settings. The list
#' \code{params} must contain the following values:
#' * \code{n_adapt}: A positive integer number of adaptive MCMC iterations.
#' * \code{n_mcmc}: A positive integer number of total MCMC iterations
#' post adaptation.
#' * \code{n_thin}: A positive integer number of MCMC iterations per saved
#' sample.
#' * \code{n_message}: A positive integer number of frequency of iterations
#'  to output a progress message. For example, \code{n_message = 50}
#'  outputs progress messages every 50 iterations.
#' @param priors is a list of prior settings. 
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param shared_covariance_params is a logicial input that determines whether to fit the spatial process with component specifice parameters. If TRUE, each component has conditionally independent Gaussian process parameters theta and tau2. If FALSE, all components share the same Gaussian process parameters theta and tau2. 
#' @param inits is the list of intial values if the user wishes to specify initial values. If these values are not specified, then the intital values will be randomly sampled from the prior.
#' @param config is the list of configuration values if the user wishes to specify initial values. If these values are not specified, then default a configuration will be used.
#' @param n_chain is the MCMC chain id. The default is 1.
#' @param progress is a logicial input that determines whether to print a progress bar.
#' @param verbose is a logicial input that determines whether to print more detailed messages.
#' @importFrom stats rmultinom
#' @importFrom hms as_hms
#' @export

## polya-gamma spatial linear regression model
pgSTLM <- function(
    Y, 
    X,
    locs, 
    params,
    priors,
    corr_fun                 = "exponential",
    n_cores                  = 1L,
    shared_covariance_params = TRUE,
    inits                    = NULL,
    config                   = NULL,
    n_chain                  = 1,
    progress                 = FALSE,
    verbose                  = FALSE
) {
    
    start <- Sys.time()
    
    ##
    ## Run error checks
    ## 
    
    pgR:::check_input_pg_stlm(Y, X, locs)
    pgR:::check_params(params)
    pgR:::check_corr_fun(corr_fun)
    # check_inits_pgLM(params, inits)
    # check_config(params, config)
    
    ## add in faster parallel cholesky as needed
    
    ## 
    ## setup config
    ##
    
    ## do we sample the functional relationship parameters? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_beta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_beta']])) {
            sample_beta <- config[['sample_beta']]
        }
    }
    
    ## do we sample the climate autocorrelation parameter? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_rho <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_rho']])) {
            sample_rho <- config[['sample_rho']]
        }
    }
    
    ## do we sample the climate variance parameter? This is primarily 
    ## used to troubleshoot model fitting using simulated data
    sample_tau2 <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_tau2']])) {
            sample_tau2 <- config[['sample_tau2']]
        }
    }
    
    ## do we sample the climate spatial covariance parameters? This is
    ## primarily used to troubleshoot model fitting using simulated data
    sample_theta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_theta']])) {
            sample_theta <- config[['sample_theta']]
        }
    }
    
    ## do we sample the latent intensity parameter eta
    sample_eta <- TRUE
    if (!is.null(config)) {
        if (!is.null(config[['sample_eta']])) {
            sample_eta <- config[['sample_eta']]
        }
    }
    
    
    ##
    ## setup constats
    ##
    
    N      <- nrow(Y)
    J      <- ncol(Y)
    n_time <- dim(Y)[3]
    p      <- ncol(X)
    D      <- fields::rdist(locs)
    
    ## Add in a counter for the number of regularized Cholesky factors.
    ## This is useful in correcting for numerical errors resulting in 
    ## covariance matrices that are not full rank
    num_chol_failures <- 0

    
    ## We assume a partially missing observation is the same as 
    ## fully missing. The index allows for fast accesing of missing
    ## observations
    missing_idx <- matrix(FALSE, N, n_time)
    for (i in 1:N) {
        for (tt in 1:n_time) {
            missing_idx[i, tt] <- any(is.na(Y[i, , tt]))
        }
    }
    
    message("There are ", ifelse(any(missing_idx), sum(missing_idx), "no"), " observations with missing count vectors")
    
    ## Calculate Mi and kappa
    Mi    <- array(0, dim = c(N, J - 1, n_time))
    kappa <- array(0, dim = c(N, J - 1, n_time))
    for (i in 1:N){
        for (tt in 1:n_time) {
            if (missing_idx[i, tt]) {
                Mi[i, , tt]    <- 0
                kappa[i, , tt] <- 0
            } else {
                Mi[i, , tt] <- sum(Y[i, , tt]) - c(0, cumsum(Y[i, , tt][1:(J - 2)]))
                kappa[i, , tt] <- Y[i, 1:(J - 1), tt] - Mi[i, , tt] / 2
            }
        }
    }
    
    ## create an index for nonzero values
    nonzero_idx <- Mi != 0
    
    ##
    ## initial values
    ##
    
    ## default priors
    
    mu_beta        <- rep(0, p)
    Sigma_beta     <- 10 * diag(p)
    
    ## check if priors for mu_beta are specified
    if (!is.null(priors[['mu_beta']])) {
        if (all(!is.na(priors[['mu_beta']]))) {
            mu_beta <- priors[['mu_beta']]
        }
    }
    
    ## check if priors for Sigma_beta are specified
    if (!is.null(priors[['Sigma_beta']])) {
        if (all(!is.na(priors[['Sigma_beta']]))) {
            Sigma_beta <- priors[['Sigma_beta']]
        }
    }
    Sigma_beta_chol <- tryCatch(
        chol(Sigma_beta),
        error = function(e) {
            if (verbose)
                message("The Cholesky decomposition of the prior covariance Sigma_beta was ill-conditioned and mildy regularized.")
            chol(Sigma_beta + 1e-8 * diag(N))                    
        }
    )
    Sigma_beta_inv  <- chol2inv(Sigma_beta_chol)
    
    ##
    ## initialize beta
    ##
    
    beta <- t(mvnfast::rmvn(J-1, mu_beta, Sigma_beta_chol, isChol = TRUE))
    ## check if initial value for beta is given
    if (!is.null(inits[['beta']])) {
        if (all(!is.na(inits[['beta']]))) {
            beta <- inits[['beta']]
        }
    }
    Xbeta <- X %*% beta
    
    ##
    ## initialize spatial Gaussian process -- share parameters across the different components
    ##    can generalize to each component getting its own covariance
    ##
    ## assume the GP parameters don't change through time
    
    theta_mean <- NULL
    theta_var  <- NULL
    if (corr_fun == "matern") {
        theta_mean <- c(priors$mean_range, priors$mean_nu)
        theta_var  <- diag(c(priors$sd_range, priors$sd_nu)^2)
    } else if (corr_fun == "exponential") {
        theta_mean <- priors$mean_range
        theta_var  <- priors$sd_range^2
    }
    
    theta <- NULL
    if (shared_covariance_params) {
        if (corr_fun == "matern") {
            theta <- as.vector(pmin(pmax(mvnfast::rmvn(1, theta_mean, theta_var), -2), 0.1))            
        } else if (corr_fun == "exponential") {
            theta <- pmin(pmax(rnorm(1, theta_mean, sqrt(theta_var)), -2), 0.1)
        }
    } else {
        if (corr_fun == "matern") {
            theta <- pmin(pmax(mvnfast::rmvn(J-1, theta_mean, theta_var), -2), 0.1)
        } else if (corr_fun == "exponential") {
            theta <- pmin(pmax(rnorm(J-1, theta_mean, sqrt(theta_var)), -2), 0.1)
        }
        
    }
    
    if (!is.null(inits[['theta']])) {
        if (all(!is.na(inits[['theta']]))) {
            theta <- inits[['theta']]
        }
    }
    
    ## check dimensions of theta
    if (shared_covariance_params) {
        if (corr_fun == "matern")
            if (!is_numeric_vector(theta, 2)) 
                stop('If shared_covariance_params is TRUE, theta must be a numeric vector of length 2 when corr_fun is "matern"')
        if (corr_fun == "exponential")
            if (!is_numeric_vector(theta, 1)) 
                stop('If shared_covariance_params is TRUE, theta must be a numeric of length 1 when corr_fun is "exponential"')
    } else {
        if (corr_fun == "matern")
            if (!is_numeric_matrix(theta, J-1, 2))
                stop('If shared_covariance_params is FALSE, theta must be a J-1 by 2 numeric matrix when corr_fun is "matern"')
        if (corr_fun == "exponential")
            if (!is_numeric_vector(theta, J-1)) 
                stop('If shared_covariance_params is FALSE, theta must be a numeric vector of length J-1 when corr_fun is "exponential"')
    }
    
    tau2 <- NULL
    if (shared_covariance_params) {
        tau2 <- min(1 / rgamma(1, priors$alpha_tau, priors$beta_tau), 10)
    } else {
        tau2 <- pmin(1 / rgamma(J-1, priors$alpha_tau, priors$beta_tau), 10)
    }
    if (!is.null(inits[['tau2']])) {
        if (all(!is.na(inits[['tau2']]))) {
            ## if tau2 passes error checks
            tau2 <- inits[['tau2']]
        }
    }
    
    ## check dimensions of tau2
    if (shared_covariance_params) {
        if (!is_positive_numeric(tau2, 1))
            stop("If shared_covariance_params is FALSE, tau2 must be a numeric scalar")
    } else {
        if (!is_positive_numeric(tau2, J-1))
            stop("If shared_covariance_params is TRUE, tau2 must be a J-1 positive numeric vector ")
    }
    
    
    Sigma <- NULL
    if (shared_covariance_params) {
        Sigma <- tau2 * correlation_function(D, theta, corr_fun = corr_fun)
    } else {
        Sigma <- array(0, dim = c(J - 1, N, N))
        for (j in 1:(J-1)) {
            if (corr_fun == "matern") {
                Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j, ], corr_fun = corr_fun)
            } else if (corr_fun == "exponential") {
                Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j], corr_fun = corr_fun)
            }
        }
    }
    
    Sigma_chol <- NULL
    if (shared_covariance_params) {
        Sigma_chol <- tryCatch(
            chol(Sigma),
            error = function(e) {
                if (verbose)
                    message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                num_chol_failures <- num_chol_failures + 1
                chol(Sigma + 1e-8 * diag(N))                    
            }
        )
    } else {
        Sigma_chol <- array(0, dim = c(J-1, N, N))
        for (j in 1:(J-1)) {
            ## add a warning for the Cholesky function
            Sigma_chol[j, , ] <- tryCatch(
                chol(Sigma[j, , ]),
                error = function(e) {
                    if (verbose)
                        message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                    num_chol_failures <- num_chol_failures + 1
                    chol(Sigma[j, , ] + 1e-8 * diag(N))                    
                }
            )
        }
    }
    
    Sigma_inv  <- NULL
    if (shared_covariance_params) {
        Sigma_inv <- chol2inv(Sigma_chol)   
    } else {
        Sigma_inv <- array(0, dim = c(J-1, N, N))
        for (j in 1:(J-1)) {
            Sigma_inv[j, , ] <- chol2inv(Sigma_chol[j, , ])
        }
    }
    
    ## temporal autocorrelation
    rho <- runif(1, 0, 1)
    if (!is.null(inits[['rho']])) {
        if (all(!is.na(inits[['rho']]))) {
            ## if rho passes error checks
            rho <- inits[['rho']]
        }
    }
    
    eta  <- array(0, dim = c(N, J-1, n_time))
    for (tt in 1:n_time) {
        if (tt == 1) {
            for (j in 1:(J-1)) {
                if (shared_covariance_params) {
                    eta[, j, 1] <- Xbeta[, j] + t(mvnfast::rmvn(1, rep(0, N), Sigma_chol, isChol = TRUE))
                } else{
                    eta[, j, 1] <- Xbeta[, j] + t(mvnfast::rmvn(1, rep(0, N), Sigma_chol[j, , ], isChol = TRUE))
                }
            }
        } else {
            for (j in 1:(J-1)) {
                if (shared_covariance_params) {
                    eta[, j, tt] <- Xbeta[, j] + mvnfast::rmvn(1, rho * eta[, j, tt - 1], Sigma_chol, isChol = TRUE)
                } else {
                    eta[, j, tt] <- Xbeta[, j] + t(mvnfast::rmvn(1, rho * eta[, j, tt - 1], Sigma_chol[j, , ], isChol = TRUE))
                }
            }
        }
    }    
    if (!is.null(inits[['eta']])) {
        if (all(!is.na(inits[['eta']]))) {
            eta <- inits[['eta']]
        }
    }
    
    ##
    ## sampler config options -- to be added later
    ## 
    #
    
    ##
    ## initialize omega
    ##
    
    omega <- array(0, dim = c(N, J-1, n_time))
    omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
    
    ##
    ## setup save variables
    ##
    
    n_save     <- params$n_mcmc / params$n_thin
    beta_save  <- array(0, dim = c(n_save, p, J-1))
    tau2_save  <- NULL
    if (shared_covariance_params) {
        tau2_save <- rep(0, n_save)
    } else {
        tau2_save <- matrix(0, n_save, J-1)
    }
    theta_save <- NULL
    if (shared_covariance_params) {
        if (corr_fun == "matern") {
            theta_save <- matrix(0, n_save, 2)
        } else if (corr_fun == "exponential") {
            theta_save <- rep(0, n_save)
        }
    } else {
        if (corr_fun == "matern") {
            theta_save <- array(0, dim = c(n_save, J-1, 2))
        } else if (corr_fun == "exponential") {
            theta_save <- matrix(0, n_save, J-1)
        }
    }
    eta_save   <- array(0, dim = c(n_save, N, J-1, n_time))
    pi_save    <- array(0, dim = c(n_save, N, J, n_time))
    rho_save   <- rep(0, n_save)
    
    ## 
    ## initialize tuning 
    ##
    
    ##
    ## tuning variables for adaptive MCMC
    ##
    
    theta_batch           <- NULL
    theta_accept          <- NULL
    theta_accept_batch    <- NULL
    lambda_theta          <- NULL
    Sigma_theta_tune      <- NULL
    Sigma_theta_tune_chol <- NULL
    theta_tune            <- NULL
    
    if (shared_covariance_params) {
        
        theta_accept       <- 0
        theta_accept_batch <- 0
        
        if (corr_fun == "matern") {
            theta_batch <- matrix(0, 50, 2) 
            lambda_theta          <- 0.05
            Sigma_theta_tune      <- 1.8 * diag(2) - .8
            Sigma_theta_tune_chol <- tryCatch(
                chol(Sigma_theta_tune),
                error = function(e) {
                    if (verbose)
                        message("The Cholesky decomposition of the Metroplois-Hastings adaptive tuning matrix for Matern parameters theta was ill-conditioned and mildy regularized.")
                    chol(Sigma_theta_tune + 1e-8 * diag(2))                    
                }
            )
        } else if (corr_fun == "exponential") {
            theta_tune <- mean(D) / 2
        }
        
    } else {
        
        theta_accept       <- rep(0, J-1)
        theta_accept_batch <- rep(0, J-1)
        
        if (corr_fun == "matern") {
            theta_batch <- array(0, dim = c(50, 2, J-1))      
            lambda_theta     <- rep(0.05, J-1)
            Sigma_theta_tune <- array(0, dim = c(2, 2, J-1))
            for (j in 1:(J-1)) {
                Sigma_theta_tune[, , j] <- 1.8 * diag(2) - .8
            }
            Sigma_theta_tune_chol <- array(0, dim = c(2, 2, J-1))
            for (j in 1:(J-1)) {
                Sigma_theta_tune_chol[, , j] <- tryCatch(
                    chol(Sigma_theta_tune[, , j]),
                    error = function(e) {
                        if (verbose)
                            message("The Cholesky decomposition of the Metroplois-Hastings adaptive tuning matrix for Matern parameters theta was ill-conditioned and mildy regularized.")
                        chol(Sigma_theta_tune[, , j] + 1e-8 * diag(2))                    
                    }
                )
            }
        } else if (corr_fun == "exponential") {
            theta_tune <- rep(mean(D) / 2, J-1)
        }
    }
    
    ## tuning for rho
    rho_accept       <- 0
    rho_accept_batch <- 0
    rho_tune         <- 0.025
    
    ##
    ## Starting MCMC chain
    ##
    
    message("Starting MCMC for chain ", n_chain, ", running for ", params$n_adapt, " adaptive iterations and ", params$n_mcmc, " fitting iterations \n")
    if (progress) {
        progressBar <- utils::txtProgressBar(style = 3)
    }
    percentage_points <- round((1:100 / 100) * (params$n_adapt + params$n_mcmc))
    
    
    for (k in 1:(params$n_adapt + params$n_mcmc)) {
        if (k == params$n_adapt + 1) {
            message("Starting MCMC fitting for chain ", n_chain, ", running for ", params$n_mcmc, " iterations \n")
        }
        if (k %% params$n_message == 0) {
            if (k <= params$n_adapt) {
                message("MCMC adaptation iteration ", k, " for chain ", n_chain)
            } else {
                message("MCMC fitting iteration ", k - params$n_adapt, " for chain ", n_chain)
            }
        }
        
        ##
        ## sample Omega
        ##
        
        if (verbose)
            message("sample omega")
        
        omega[nonzero_idx] <- pgdraw(Mi[nonzero_idx], eta[nonzero_idx], cores = n_cores)
        
        ##
        ## sample beta
        ##
        
        ## can parallelize this update -- each group of parameters is 
        ## conditionally independent given omega and kappa(y)
        
        if (sample_beta){
            if (verbose)
                message("sample beta")
            
            if (shared_covariance_params) {
                for (j in 1:(J-1)) {
                    tXSigma_inv <- t(X) %*% Sigma_inv
                    A <- n_time * tXSigma_inv %*% X + Sigma_beta_inv
                    ## guarantee a symmetric matrix
                    A <- (A + t(A)) / 2
                    b <- rowSums(tXSigma_inv %*% eta[, j, ]) + Sigma_beta_inv %*% mu_beta
                    beta[, j]   <- rmvn_arma(A, b)
                }
            } else {
                for (j in 1:(J-1)) {
                    tXSigma_inv <- t(X) %*% Sigma_inv[j, , ]
                    A <- n_time * tXSigma_inv %*% X + Sigma_beta_inv
                    ## guarantee a symmetric matrix
                    A <- (A + t(A)) / 2
                    b <- rowSums(rho * tXSigma_inv %*% eta[, j, ]) + Sigma_beta_inv %*% mu_beta
                    beta[, j]   <- rmvn_arma(A, b)
                }
            }        
            Xbeta <- X %*% beta
        }
        
        ##
        ## sample spatial correlation parameters theta
        ##
        
        if(sample_theta) {      
            if (verbose)
                message("sample theta")
            
            if (shared_covariance_params) {
                ## update a common theta for all processes
                theta_star <- NULL
                if (corr_fun == "matern") {
                    theta_star <- as.vector(
                        mvnfast::rmvn( 
                            n      = 1,
                            mu     = theta,
                            sigma  = lambda_theta * Sigma_theta_tune_chol,
                            isChol = TRUE
                        )
                    )
                } else if (corr_fun == "exponential") {
                    theta_star <- rnorm(1, theta, theta_tune)
                }
                Sigma_star       <- tau2 * correlation_function(D, theta_star, corr_fun = corr_fun)
                ## add in faster parallel cholesky as needed
                Sigma_chol_star <- tryCatch(
                    chol(Sigma_star),
                    error = function(e) {
                        if (verbose)
                            message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                        num_chol_failures <- num_chol_failures + 1
                        chol(Sigma_star + 1e-8 * diag(N))                    
                    }
                )
                Sigma_inv_star  <- chol2inv(Sigma_chol_star)
                ## parallelize this
                mh1 <- sum(
                    sapply(1:(J-1), function(j) {
                        mvnfast::dmvn(eta[, j, 1], Xbeta[, j], Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) 
                    })
                ) +
                    sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho * eta[, j, tt - 1], Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) 
                            })
                        }) 
                        
                    ) +
                    ## prior
                    mvnfast::dmvn(theta_star, theta_mean, theta_var, log = TRUE)
                ## parallelize this        
                mh2 <- sum(
                    sapply(1:(J-1), function(j) {
                        mvnfast::dmvn(eta[, j, 1], Xbeta[, j], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores) 
                    })
                ) +
                    sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho * eta[, j, tt - 1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores)
                            })
                        })
                    ) +
                    ## prior
                    mvnfast::dmvn(theta, theta_mean, theta_var, log = TRUE)
                
                mh <- exp(mh1 - mh2)
                if (mh > runif(1, 0, 1)) {
                    theta      <- theta_star
                    Sigma      <- Sigma_star
                    Sigma_chol <- Sigma_chol_star
                    Sigma_inv  <- Sigma_inv_star 
                    if (k <= params$n_adapt) {
                        theta_accept_batch <- theta_accept_batch + 1 / 50
                    } else {
                        theta_accept <- theta_accept + 1 / params$n_mcmc
                    }
                }
                ## adapt the tuning
                if (k <= params$n_adapt) {
                    if (corr_fun == "matern") {
                        save_idx <- k %% 50
                        if ((k %% 50) == 0) {
                            save_idx <- 50
                        } 
                        theta_batch[save_idx, ] <- theta 
                    }
                    if (k %% 50 == 0) {
                        if (corr_fun == "matern") {
                            out_tuning <- update_tuning_mv(
                                k,
                                theta_accept_batch,
                                lambda_theta,
                                theta_batch,
                                Sigma_theta_tune,
                                Sigma_theta_tune_chol
                            )
                            theta_batch           <- out_tuning$batch_samples
                            Sigma_theta_tune      <- out_tuning$Sigma_tune
                            Sigma_theta_tune_chol <- out_tuning$Sigma_tune_chol
                            lambda_theta          <- out_tuning$lambda
                            theta_accept_batch    <- out_tuning$accept
                        } else if (corr_fun == "exponential") {
                            out_tuning         <- update_tuning(k, theta_accept_batch, theta_tune)
                            theta_tune         <- out_tuning$tune
                            theta_accept_batch <- out_tuning$accept
                        }
                    }
                }
            } else {
                ## 
                ## theta varies for each component
                ##
                for (j in 1:(J-1)) {
                    theta_star <- NULL
                    if (corr_fun == "matern") {
                        theta_star <- as.vector(
                            mvnfast::rmvn( 
                                n      = 1,
                                mu     = theta[j, ],
                                sigma  = lambda_theta[j] * Sigma_theta_tune_chol[, , j],
                                isChol = TRUE
                            )
                        )
                    } else if (corr_fun == "exponential") {
                        theta_star <- rnorm(1, theta[j], theta_tune)
                    }
                    
                    Sigma_star      <- tau2[j] * correlation_function(D, theta_star, corr_fun = corr_fun)
                    ## add in faster parallel cholesky as needed
                    Sigma_chol_star <- tryCatch(
                        chol(Sigma_star),
                        error = function(e) {
                            if (verbose)
                                message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                            num_chol_failures <- num_chol_failures + 1
                            chol(Sigma_star + 1e-8 * diag(N))                    
                        }
                    )
                    Sigma_inv_star  <- chol2inv(Sigma_chol_star)
                    
                    ## parallelize this
                    mh1 <- mvnfast::dmvn(eta[, j, 1], Xbeta[, j], Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) +
                      sum(
                        sapply(2:n_time, function(tt) {
                            mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho * eta[, j, tt - 1], Sigma_chol_star, isChol = TRUE, log = TRUE, ncores = n_cores) 
                        })) +
                      ## prior
                      mvnfast::dmvn(theta_star, theta_mean, theta_var, log = TRUE)

                    ## parallelize this        
                    mh2 <- mvnfast::dmvn(eta[, j, 1], Xbeta[, j], Sigma_chol[j, , ], isChol = TRUE, log = TRUE, ncores = n_cores) +
                        sum(sapply(2:n_time, function(tt) {
                                    mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho * eta[, j, tt - 1], Sigma_chol[j, , ], isChol = TRUE, log = TRUE, ncores = n_cores)
                            })) +
                        ## prior
                        mvnfast::dmvn(theta[j, , drop = FALSE], theta_mean, theta_var, log = TRUE)

                    mh <- exp(mh1 - mh2)
                    if (mh > runif(1, 0, 1)) {
                        if (corr_fun == "matern") {
                            theta[j, ] <- theta_star
                        } else if (corr_fun == "exponential") {
                            theta[j] <- theta_star
                        }
                        Sigma[j, , ]      <- Sigma_star
                        Sigma_chol[j, , ] <- Sigma_chol_star
                        Sigma_inv[j, , ]  <- Sigma_inv_star 
                        if (k <= params$n_adapt) {
                            theta_accept_batch[j] <- theta_accept_batch[j] + 1 / 50
                        } else {
                            theta_accept[j] <- theta_accept[j] + 1 / params$n_mcmc
                        }
                    }
                }
                ## adapt the tuning
                if (k <= params$n_adapt) {
                    if (corr_fun == "matern") {
                        save_idx <- k %% 50
                        if ((k %% 50) == 0) {
                            save_idx <- 50
                        } 
                        theta_batch[save_idx, , ] <- theta 
                    }
                    if (k %% 50 == 0) {
                        if (corr_fun == "matern") {
                            out_tuning <- update_tuning_mv_mat(
                                k,
                                theta_accept_batch,
                                lambda_theta,
                                theta_batch,
                                Sigma_theta_tune,
                                Sigma_theta_tune_chol
                            )
                            theta_batch           <- out_tuning$batch_samples
                            Sigma_theta_tune      <- out_tuning$Sigma_tune
                            Sigma_theta_tune_chol <- out_tuning$Sigma_tune_chol
                            lambda_theta          <- out_tuning$lambda
                            theta_accept_batch    <- out_tuning$accept
                        } else if (corr_fun == "exponential") {
                            out_tuning         <- update_tuning_vec(k, theta_accept_batch, theta_tune)
                            theta_tune         <- out_tuning$tune
                            theta_accept_batch <- out_tuning$accept
                        }
                    }   
                }        
            }
        }        
        
        ##
        ## sample spatial process variance tau2
        ##
        
        if (sample_tau2) {
            if (verbose)
                message("sample tau2")
            
            ## can we make this more efficient?
            devs <- array(0, dim = c(N, J-1, n_time))
            devs[, , 1] <- eta[, , 1] - Xbeta
            for (tt in 2:n_time) {
                devs[, , tt] <- eta[, , tt] -  rho * eta[, , tt - 1] - (1 - rho) * Xbeta
            }
            
            if (shared_covariance_params) {
                
                ## double check this math later -- seems right for now
                SS <- sum(
                    sapply(1:n_time, function(tt) {
                        devs[, , tt] * (tau2 * Sigma_inv %*% devs[, , tt])
                    })
                )
                tau2       <- 1 / rgamma(1, N * (J-1) * n_time / 2 + priors$alpha_tau, SS / 2 + priors$beta_tau) 
                Sigma      <- tau2 * correlation_function(D, theta, corr_fun = corr_fun) 
                ## add in faster parallel cholesky as needed
                ## see https://github.com/RfastOfficial/Rfast/blob/master/src/cholesky.cpp
                Sigma_chol <- tryCatch(
                    chol(Sigma),
                    error = function(e) {
                        if (verbose)
                            message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                        num_chol_failures <- num_chol_failures + 1
                        chol(Sigma + 1e-8 * diag(N))                    
                    }
                )
                Sigma_inv  <- chol2inv(Sigma_chol) 
            } else {
                for (j in 1:(J-1)) {
                    SS <- sum(
                        sapply(1:n_time, function(tt) {
                            devs[, j, tt] * (tau2[j] * Sigma_inv[j, , ] %*% devs[, j, tt])
                        })
                    )
                    tau2[j]    <- 1 / rgamma(1, N / 2 * n_time + priors$alpha_tau, SS / 2 + priors$beta_tau) 
                    if (corr_fun == "matern") {
                        Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j, ], corr_fun = corr_fun) 
                    } else {
                        Sigma[j, , ] <- tau2[j] * correlation_function(D, theta[j], corr_fun = corr_fun) 
                    }
                    ## add in faster parallel cholesky as needed
                    ## see https://github.com/RfastOfficial/Rfast/blob/master/src/cholesky.cpp
                    Sigma_chol[j, , ] <- tryCatch(
                        chol(Sigma[j, , ]),
                        error = function(e) {
                            if (verbose)
                                message("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized. If this warning is rare, this should be safe to ignore.")
                            num_chol_failures <- num_chol_failures + 1
                            chol(Sigma[j, , ] + 1e-8 * diag(N))                    
                        }
                    )
                    Sigma_inv[j, , ]  <- chol2inv(Sigma_chol[j, , ])
                }
            }
        }
        
        ##
        ## sample eta
        ##
        
        if (sample_eta) {
            if (verbose)
                message("sample eta")
            
            ## need to add in the autocorrelation process (how does this affect the kappas?)
            for (tt in 1:n_time) {
                for (j in 1:(J-1)) {
                    A <- NULL
                    b <- NULL
                    if (tt == 1) {
                        if (shared_covariance_params) {
                            A <- (1 + rho^2) * Sigma_inv + diag(omega[, j, tt])
                            b <- Sigma_inv %*% ((1 - rho + rho^2) * Xbeta[, j] + rho * eta[, j, tt + 1]) + kappa[, j, tt]
                        } else {
                            A <- (1 + rho^2) * Sigma_inv[j, , ] + diag(omega[, j, tt])
                            b <- Sigma_inv[j,,] %*% ((1 - rho + rho^2) * Xbeta[, j] + rho * eta[, j, tt + 1]) + kappa[, j, tt]
                        }
                    } else if (tt == n_time) {
                        if (shared_covariance_params) {
                            A <- Sigma_inv + diag(omega[, j, tt])
                            b     <- Sigma_inv %*% ((1 - rho) * Xbeta[, j] + rho * eta[, j, tt - 1]) + kappa[, j, tt]
                        } else {
                            A <- Sigma_inv[j, , ] + diag(omega[, j, tt])
                            b <- Sigma_inv[j,,] %*% ((1 - rho) * Xbeta[, j] + rho * eta[, j, tt - 1]) + kappa[, j, tt]
                        }
                    } else {
                        if (shared_covariance_params) {
                            A <- (1 + rho^2) * Sigma_inv + diag(omega[, j, tt])
                            b <- Sigma_inv %*% ((1 - rho)^2 * Xbeta[, j] + rho * (eta[, j, tt - 1] +  eta[, j, tt + 1])) + kappa[, j, tt]
                        } else {
                            A <- (1 + rho^2) * Sigma_inv[j, , ] + diag(omega[, j, tt])
                            b <- Sigma_inv[j,,] %*% ((1 - rho)^2 * Xbeta[, j] + rho * (eta[, j, tt - 1] +  eta[, j, tt + 1])) + kappa[, j, tt]
                        }
                        ## is this (1 - rho) * Xbeta[, j] or (1 + rho) * Xbeta[, j]???
                    }
                    ## guarantee that A is symmetric
                    A            <- (A + t(A)) / 2
                    eta[, j, tt] <- rmvn_arma(A, b)
                }
            }
        }
        
        ##
        ## sample rho
        ##
        
        if (sample_rho) {
            if (verbose)
                message("sample rho")
            
            rho_star <- rnorm(1, rho, rho_tune)
            if (rho_star < 1 & rho_star > -1) {
                mh1 <- NULL
                mh2 <- NULL
                if (shared_covariance_params){
                    mh1 <- sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho_star * eta[, j, tt - 1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores) 
                            })
                        }) 
                        
                    ) 
                    ## parallelize this        
                    mh2 <- sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho * eta[, j, tt - 1], Sigma_chol, isChol = TRUE, log = TRUE, ncores = n_cores) 
                            })
                        }) 
                    ) 
                } else {
                    mh1 <- sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho_star * eta[, j, tt - 1], Sigma_chol[j,,], isChol = TRUE, log = TRUE, ncores = n_cores) 
                            })
                        }) 
                        
                    ) 
                    ## parallelize this        
                    mh2 <- sum(
                        sapply(2:n_time, function(tt) {
                            sapply(1:(J-1), function(j) {
                                mvnfast::dmvn(eta[, j, tt], Xbeta[, j] + rho * eta[, j, tt - 1], Sigma_chol[j,,], isChol = TRUE, log = TRUE, ncores = n_cores) 
                            })
                        }) 
                    )
                }
                
                mh <- exp(mh1 - mh2)
                if (length(mh) > 1)
                    stop("error in mh for rho")
                if (mh > runif(1, 0.0, 1.0)) {
                    rho   <- rho_star
                    if (k <= params$n_adapt) {
                        rho_accept_batch <- rho_accept_batch + 1.0 / 50.0
                    } else {
                        rho_accept <- rho_accept + 1.0 / params$n_mcmc
                    }
                }
                
                ## update tuning
                if (k <= params$n_adapt) {
                    if (k %% 50 == 0){
                        out_tuning <- update_tuning(
                            k,
                            rho_accept_batch, 
                            rho_tune
                        )
                        rho_tune         <- out_tuning$tune
                        rho_accept_batch <- out_tuning$accept
                    }
                }
            }
        }
        
        ##
        ## save variables
        ##
        
        if (k >= params$n_adapt) {
            if (k %% params$n_thin == 0) {
                save_idx                <- (k - params$n_adapt) / params$n_thin
                beta_save[save_idx, , ] <- beta
                if (shared_covariance_params) {
                    if (corr_fun == "matern") {
                        theta_save[save_idx, ]  <- theta
                    } else if (corr_fun == "exponential") {
                        theta_save[save_idx]  <- theta
                    }
                    tau2_save[save_idx]     <- tau2
                } else {
                    if (corr_fun == "matern") {
                        theta_save[save_idx, , ]  <- theta
                    } else if (corr_fun == "exponential") {
                        theta_save[save_idx, ]  <- theta
                    }
                    tau2_save[save_idx, ]     <- tau2
                }
                eta_save[save_idx, , , ] <- eta
                for (tt in 1:n_time) {
                    pi_save[save_idx, , , tt]  <- eta_to_pi(eta[, , tt])
                }
                rho_save[save_idx]       <- rho
            }
          
          out <- list(
            beta  = beta_save,
            theta = theta_save,
            tau2  = tau2_save,
            eta   = eta_save,
            pi    = pi_save,
            rho   = rho_save
          )
          # class(out) <- "pg_stlm"
          saveRDS(out, paste0('output/pg-inter-output.RDS'))
        }
        
        
        
        
        ##
        ## End of MCMC loop
        ##
        
        if (k %in% percentage_points && progress) {
            utils::setTxtProgressBar(progressBar, k / (params$n_adapt + params$n_mcmc))
        }
    }
    
    ## print out acceptance rates -- no tuning in this model
    
    if (num_chol_failures > 0)
        warning("The Cholesky decomposition of the Matern correlation function was ill-conditioned and mildy regularized ", num_chol_failures, " times. If this warning is rare, this should be safe to ignore.")
    
    ## eventually create a model class and include this as a variable in the class
    message("Acceptance rate for theta is ", mean(theta_accept))
    message("Acceptance rate for rho is ", mean(rho_accept))
    
    
    ##
    ## return the MCMC output -- think about a better way to make this a class
    ## 
    
    if (progress) {
        close(progressBar)
    }
    
    stop    <- Sys.time()
    runtime <- stop - start
    
    message("MCMC took ", hms::as_hms(runtime))
    
    
    out <- list(
        beta  = beta_save,
        theta = theta_save,
        tau2  = tau2_save,
        eta   = eta_save,
        pi    = pi_save,
        rho   = rho_save
    )
    class(out) <- "pg_stlm"
    
    return(out)
}


