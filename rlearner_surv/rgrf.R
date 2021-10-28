#' @include utils.R
#'
#' @title R-learner, implemented via xgboost (boosting)
#'
#' @description  R-learner, as proposed by Nie and Wager (2017), implemented via xgboost (boosting)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param k_folds number of folds used for cross fitting and cross validation
#' @param p_hat pre-computed estimates on E[W|X] corresponding to the input x. rboost will compute it internally if not provided.
#' @param m_hat pre-computed estimates on E[Y|X] corresponding to the input x. rboost will compute it internally if not provided.
#' @param ntrees_max the maximum number of trees to grow for xgboost
#' @param num_search_rounds the number of random sampling of hyperparameter combinations for cross validating on xgboost trees
#' @param print_every_n the number of iterations (in each iteration, a tree is grown) by which the code prints out information
#' @param early_stopping_rounds the number of rounds the test error stops decreasing by which the cross validation in finding the optimal number of trees stops
#' @param nthread the number of threads to use. The default is NULL, which uses all available threads
#' @param verbose boolean; whether to print statistic
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rboost_fit = rboost(x, w, y)
#' rboost_est = predict(rboost_fit, x)
#' }
#'
#' @export
rgrf = function(x, w, y, D,
                times = NULL,
                k_folds = NULL,
                p_hat = NULL,
                m_hat = NULL,
                c_hat = NULL,      # censoring weight 
                num.trees = 2000,
                nthread = NULL,
                verbose = FALSE){


  input = sanitize_input(x,w,y,D)
  x = input$x
  w = as.numeric(input$w)
  y = input$y
  D = input$D
  nobs = nrow(x)
  pobs = ncol(x)

  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(y)/4)))
  }

  tempdat <- data.frame(y, D, w, x)
  
  if (is.null(p_hat)){
    w_fit <- regression_forest(x, w, tune.parameters = "all") 
    p_hat <- w_fit$predictions
  }else{
    w_fit = NULL
  }
  
  if (is.null(m_hat)){
    y_fit <- survival_forest(as.matrix(tempdat[, 3:dim(tempdat)[2]]),
                             tempdat$y, 
                             tempdat$D)
    tempdat1 <- tempdat; tempdat1$w <- 1
    tempdat0 <- tempdat; tempdat0$w <- 0
    surf1 <- predict(y_fit, as.matrix(tempdat1[,3:dim(tempdat1)[2]]))$predictions[, which.min(abs(y_fit$failure.times-times))]
    surf0 <- predict(y_fit, as.matrix(tempdat0[,3:dim(tempdat0)[2]]))$predictions[, which.min(abs(y_fit$failure.times-times))]
    m_hat  <- p_hat * surf1 + (1 - p_hat) * surf0 
  }else {
    y_fit = NULL
  }

  if (is.null(c_hat)){
    c_fit <- survival_forest(cbind(w,x),
                             y, 
                             1 - D,
                             prediction.type = "Nelson-Aalen")
    C.hat <- c_fit$predictions
    cent <- y
    cent[D==0] <- times
    C.index <- rep(NA, length(cent))
    for (h in 1:length(cent)){
      C.index[h] <- which.min(abs(sort(unique(y[D==0])) - cent[h]))
    }
    c_hat <- C.hat[cbind(1:length(y), C.index)]
  }else{
    c_fit <- NULL
  }

  # use binary data
  tempdata <- data.frame(y, D, w, m_hat, p_hat, c_hat, x)
  binary_data <- tempdata[tempdata$D==1|tempdata$y > times,]          # remove subjects who got censored before the time of interest t50
  binary_data$D[binary_data$D==1 & binary_data$y > times] <- 0     # recode the event status for subjects who had events after t50
  binary_data <- binary_data[complete.cases(binary_data),]
  
  y_tilde <- (1 - binary_data$D) - binary_data$m_hat
  w_tilde <-  binary_data$w - binary_data$p_hat
  pseudo_outcome <- y_tilde/w_tilde
  weights <- w_tilde^2/binary_data$c_hat  
  
  tau_dat <- data.frame(pseudo_outcome, binary_data[,7:dim(binary_data)[2]])
  tau_fit <- regression_forest(tau_dat[, 2:dim(tau_dat)[2]], 
                               tau_dat$pseudo_outcome, 
                               sample.weights = weights, 
                               tune.parameters = "all")

  ret <- list(tau_fit = tau_fit,
              pseudo_outcome = pseudo_outcome,
              weights = weights,
              w_fit = w_fit,
              y_fit = y_fit,
              c_fit = c_fit,
              p_hat = p_hat,
              m_hat = m_hat,
              c_hat = c_hat)
  class(ret) <- "rgrf"
  ret
}

#' predict for rgrf
#'
#' get estimated tau(x) using the trained rgrf model
#'
#' @param object a rgrf object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param tau_only if set to TRUE, onlly return prediction on tau. Otherwise, return a list including prediction on tau, propensity score, and baseline main effect.
#' @param ... additional arguments (currently not used)
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rboost_fit = rboost(x, w, y)
#' rboost_est = predict(rboost_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.rgrf <- function(object,
                         newx = NULL,
                         tau_only = T,
                          ...) {
  if (!is.null(newx)){
    newx = sanitize_x(newx)
  }
  if (tau_only) {
    return(predict(object$tau_fit, newx)$predictions)
  } else {
    tau <- predict(object$tau_fit, newx)$predictions
    e = predict(object$w_fit, newx)$predictions
    m = predict(object$y_fit, newx)$predictions
    mu1 = m + (1-e) * tau
    mu0 = m - e * tau
    return(list(tau=tau, e=e, m=m, mu1 = mu1, mu0 = mu0))
  }
}
