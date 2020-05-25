#' Bayesian Polya-gamma regression prediction
#' 
#' this function generates predictions from the Bayesian multinomial regression using Polya-gamma data augmentation
#' @param out is a list of MCMC outputs from pgSPLM
#' @param X is a \eqn{n \times p}{n x p} matrix of covariates at the observed locations.
#' @param X_pred is a \eqn{n_{pred} \times p}{n_{pred} x p} matrix of covariates at the locations where predictions are to be made. 
#' @param locs is a \eqn{n \times 2}{n x 2} matrix of locations where observations were taken.
#' @param locs_pred is a \eqn{n_pred \times 2}{n_pred x 2} matrix of locations where predictions are to be made.
#' @param corr_fun is a character that denotes the correlation function form. Current options include "matern" and "exponential".
#' @param shared_covariance_params
#' @param n_cores is the number of cores for parallel computation using openMP.
#' @param progress is a logicial input that determines whether to print a progress bar.
#' 
#' @export 

predict_pgSPLM <- function(
    out,
    X,
    X_pred,
    locs,
    locs_pred,
    corr_fun,
    # shared_covariance_params,
    shared_tau = FALSE,
    shared_theta = FALSE,
    diag_adjust = 0,
    n_cores = 1L,
    progress = TRUE
) {
    check_corr_fun(corr_fun)
    
    beta      <- out$beta
    theta     <- out$theta
    tau2      <- out$tau2
    eta       <- out$eta
    n_samples <- nrow(beta)    
    n_pred    <- nrow(X_pred)
    J         <- dim(beta)[3] + 1
    
    if (n_pred > 10000) {
        stop("Number of prediction points must be less than 10000")
    }
    
    D_obs      <- fields::rdist(locs)
    D_pred     <- fields::rdist(locs_pred)
    D_pred_obs <- fields::rdist(locs_pred, locs)
    
    eta_pred <- array(0, dim = c(n_samples, n_pred, J-1))
    
    if (progress) {
        message("Beginning Kriging estimates")
        progressBar <- txtProgressBar(style = 3)
    }
    percentage_points <- round((1:100 / 100) * n_samples)   
    
    ## parallelize this later
    for (k in 1:n_samples) {
        if ((shared_theta)&(shared_tau)) {
            if (corr_fun == "matern") {
                Sigma           <- tau2[k] * correlation_function(D_obs, theta[k, ], corr_fun = corr_fun)
                Sigma_unobs     <- tau2[k] * correlation_function(D_pred, theta[k, ], corr_fun = corr_fun)
                Sigma_unobs_obs <- tau2[k] * correlation_function(D_pred_obs, theta[k, ], corr_fun = corr_fun)
            } else if (corr_fun == "exponential") {
                Sigma           <- tau2[k] * correlation_function(D_obs, theta[k], corr_fun = corr_fun)
                Sigma_unobs     <- tau2[k] * correlation_function(D_pred, theta[k], corr_fun = corr_fun)
                Sigma_unobs_obs <- tau2[k] * correlation_function(D_pred_obs, theta[k], corr_fun = corr_fun)
            }           
            Sigma_inv       <- chol2inv(chol(Sigma))        
            for (j in 1:(J - 1)) {
                pred_mean <- Sigma_unobs_obs %*% (Sigma_inv %*% (eta[k, , j] - X %*% beta[k, , j])) + X_pred %*% beta[k, , j]
                pred_var  <- Sigma_unobs - (Sigma_unobs_obs %*% Sigma_inv) %*% t(Sigma_unobs_obs) + diag(diag_adjust, nrow=n_pred, ncol=n_pred)  
                eta_pred[k, , j] <- mvnfast::rmvn(1, pred_mean, pred_var)
            } 
        } else if (!shared_theta & shared_tau){
            for (j in 1:(J - 1)) {
                if (corr_fun == "matern") {
                    Sigma           <- tau2[k] * correlation_function(D_obs, theta[k, j, ], corr_fun = corr_fun)
                    Sigma_unobs     <- tau2[k] * correlation_function(D_pred, theta[k, j, ], corr_fun = corr_fun)
                    Sigma_unobs_obs <- tau2[k] * correlation_function(D_pred_obs, theta[k, j, ], corr_fun = corr_fun)
                } else if (corr_fun == "exponential") {
                    Sigma           <- tau2[k] * correlation_function(D_obs, theta[k, j], corr_fun = corr_fun)
                    Sigma_unobs     <- tau2[k] * correlation_function(D_pred, theta[k, j], corr_fun = corr_fun)
                    Sigma_unobs_obs <- tau2[k] * correlation_function(D_pred_obs, theta[k, j], corr_fun = corr_fun)
                }
                
                Sigma_inv       <- chol2inv(chol(Sigma))        
                pred_mean <- Sigma_unobs_obs %*% (Sigma_inv %*% (eta[k, , j] - X %*% beta[k, , j])) + X_pred %*% beta[k, , j]
                eta_pred[k, , j] <- pred_mean
                
                # pred_var  <- Sigma_unobs - (Sigma_unobs_obs %*% Sigma_inv) %*% t(Sigma_unobs_obs)
                # eta_pred[k, , j] <- mvnfast::rmvn(1, pred_mean, pred_var)
            } 
        } else if (!shared_theta & !shared_tau){
          for (j in 1:(J - 1)) {
            if (corr_fun == "matern") {
              Sigma           <- tau2[k, j] * correlation_function(D_obs, theta[k, j, ], corr_fun = corr_fun)
              Sigma_unobs     <- tau2[k, j] * correlation_function(D_pred, theta[k, j, ], corr_fun = corr_fun)
              Sigma_unobs_obs <- tau2[k, j] * correlation_function(D_pred_obs, theta[k, j, ], corr_fun = corr_fun)
            } else if (corr_fun == "exponential") {
              Sigma           <- tau2[k, j] * correlation_function(D_obs, theta[k, j], corr_fun = corr_fun)
              Sigma_unobs     <- tau2[k, j] * correlation_function(D_pred, theta[k, j], corr_fun = corr_fun)
              Sigma_unobs_obs <- tau2[k, j] * correlation_function(D_pred_obs, theta[k, j], corr_fun = corr_fun)
            }
            
            Sigma_inv       <- chol2inv(chol(Sigma))        
            pred_mean <- Sigma_unobs_obs %*% (Sigma_inv %*% (eta[k, , j] - X %*% beta[k, , j])) + X_pred %*% beta[k, , j]
            pred_var  <- Sigma_unobs - (Sigma_unobs_obs %*% Sigma_inv) %*% t(Sigma_unobs_obs)  + diag(diag_adjust, nrow=n_pred, ncol=n_pred)  
            eta_pred[k, , j] <- mvnfast::rmvn(1, pred_mean, pred_var)
          } 
        }
        if (k %in% percentage_points && progress) {
            setTxtProgressBar(progressBar, k / n_samples)
        }
    }
    
    if (progress) {
        close(progressBar)
    }
    
    ## convert from eta to pi
    pi_pred <- sapply(1:n_samples, function(i) eta_to_pi(eta_pred[i, , ]), simplify = "array")
    ## permute to be in order of MCMC samples (rows), 
    ##    observations (columns), components (slices)
    pi_pred <- aperm(pi_pred, c(3, 1, 2))
    

    return(
        list(
            eta = eta_pred, 
            pi  = pi_pred
        )
    )
}

