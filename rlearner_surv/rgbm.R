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
#' @param n.trees_max the maximum number of trees to grow for xgboost
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
rgbm = function(x, w, y, D,
                times = NULL,
                k_folds = NULL,
                p_hat = NULL,
                m_hat = NULL,
                c_hat = NULL,      # censoring weight 
                n.trees = NULL,
                shrinkage = NULL,
                cv.folds = NULL,
                nthread = NULL,
                cf = TRUE, 
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
  foldid <- sample(rep(seq(k_folds), length = length(w)))
  
  if (is.null(p_hat)){
    if (cf) {
      p_hat <- rep(NA, length(w))
      for (k in 1:k_folds){
        testdat <- tempdat[foldid==k, ]
        traindat <- tempdat[!foldid==k, ]
        w_dat <- traindat[, 3:dim(traindat)[2]]
        w_fit <- gbm(w ~ ., 
                     data = w_dat, 
                     distribution = "bernoulli",
                     shrinkage = shrinkage, 
                     n.trees = n.trees,
                     cv.folds = cv.folds) 
        best_iter <- gbm.perf(w_fit, plot.it = F, method = 'cv')
        log_odds <- predict(w_fit, testdat[, 4:dim(testdat)[2]], n.trees = best_iter)
        p_hat[foldid==k] <- 1/(1 + exp(-log_odds))
      }
    } else {
      w_dat <- tempdat[, 3:dim(tempdat)[2]]
      w_fit <- gbm(w ~ ., 
                   data = w_dat, 
                   distribution = "bernoulli",
                   shrinkage = shrinkage, 
                   n.trees = n.trees,
                   cv.folds = cv.folds, 
                   n.cores = detectCores()) 
      best_iter <- gbm.perf(w_fit, plot.it = F, method = 'cv')
      log_odds <- predict(w_fit, w_dat[, 2:dim(w_dat)[2]], n.trees = best_iter)
      p_hat <- 1/(1 + exp(-log_odds))
    }
  }else{
    w_fit = NULL
  }
  
  if (is.null(m_hat)){
    if (cf){
      survt1 <- survt0 <- rep(NA, length(w))
      for (k in 1:k_folds){
        traindat <- tempdat[!foldid==k, ]
        testdat  <- tempdat[foldid==k, ]
        testdat1 <- testdat; testdat1$w <- 1
        testdat0 <- testdat; testdat0$w <- 0
        
        y_fit <- gbm(Surv(y, D)~.,
                     data = traindat,
                     distribution = "coxph",
                     shrinkage = shrinkage, 
                     n.trees = n.trees,
                     cv.folds = cv.folds, 
                     n.cores = detectCores()) 
        best_iter <- gbm.perf(y_fit, plot.it = F, method = 'cv')
        pred.train <- predict(y_fit, traindat, n.trees = best_iter)
        basehaz.cum <- basehaz.gbm(traindat$y, traindat$D, pred.train, t.eval = times, cumulative = TRUE)
        
        # survival probs
        pred.test1 <- predict(y_fit, testdat1, n.trees = best_iter)
        pred.test0 <- predict(y_fit, testdat0, n.trees = best_iter)
        survt1[foldid==k] <- exp(-exp(pred.test1)*basehaz.cum)
        survt0[foldid==k] <- exp(-exp(pred.test0)*basehaz.cum)
      }
      m_hat  <- p_hat * survt1 + (1 - p_hat) * survt0
    } else {
      testdat  <- traindat <- tempdat
      testdat1 <- testdat; testdat1$w <- 1
      testdat0 <- testdat; testdat0$w <- 0
      
      y_fit <- gbm(Surv(y, D)~.,
                   data = traindat,
                   distribution = "coxph",
                   shrinkage = shrinkage, 
                   n.trees = n.trees,
                   cv.folds = cv.folds, 
                   n.cores = detectCores()) 
      best_iter <- gbm.perf(y_fit, plot.it = F, method = 'cv')
      pred.train <- predict(y_fit, traindat, n.trees = best_iter)
      basehaz.cum <- basehaz.gbm(traindat$y, traindat$D, pred.train, t.eval = times, cumulative = TRUE)
      
      # survival probs
      pred.test1 <- predict(y_fit, testdat1, n.trees = best_iter)
      pred.test0 <- predict(y_fit, testdat0, n.trees = best_iter)
      survt1 <- exp(-exp(pred.test1)*basehaz.cum)
      survt0 <- exp(-exp(pred.test0)*basehaz.cum)
      
      m_hat  <- p_hat * survt1 + (1 - p_hat) * survt0
    }
  }else {
    y_fit = NULL
  }

  if (is.null(c_hat)){
    if (cf) {
      c_hat <- rep(NA, length(w))
      for (k in 1:k_folds){
        testdat <- tempdat[foldid==k, ]
        traindat <- tempdat[!foldid==k, ]
        c_fit <- gbm(Surv(y, 1-D)~.,
                     data = traindat, 
                     distribution = "coxph",
                     shrinkage = shrinkage, 
                     n.trees = n.trees,
                     cv.folds = cv.folds, 
                     n.cores = detectCores()) 
        best_iter <- gbm.perf(c_fit, plot.it = F, method = 'cv')
        pred.train <- predict(c_fit, traindat, n.trees = best_iter)
        cent <- testdat$y
        cent[testdat$D==0] <- times
        basehaz.cum <- rep(NA, length(cent))
        for (z in 1:length(cent)){
          basehaz.cum[z] <- basehaz.gbm(traindat$y, 1-traindat$D, pred.train, t.eval = cent[z], cumulative = TRUE)
        }
        pred.test <- predict(c_fit, testdat, n.trees = best_iter)
        c_hat[foldid==k] <- exp(-exp(pred.test)*basehaz.cum)
      }
    } else {
      c_fit <- gbm(Surv(y, 1 - D)~.,
                   data = tempdat, 
                   distribution = "coxph",
                   shrinkage = shrinkage, 
                   n.trees = n.trees,
                   cv.folds = cv.folds, 
                   n.cores = detectCores()) 
      best_iter <- gbm.perf(c_fit, plot.it = F, method = 'cv')
      pred.train <- predict(c_fit, tempdat, n.trees = best_iter)
      cent <- tempdat$y
      cent[tempdat$D==0] <- times
      basehaz.cum <- rep(NA, length(cent))
      for (z in 1:length(cent)){
        basehaz.cum[z] <- basehaz.gbm(tempdat$y, 1 - tempdat$D, pred.train, t.eval = cent[z], cumulative = TRUE)
      }
      pred.test <- predict(c_fit, tempdat, n.trees = best_iter)
      c_hat <- exp(-exp(pred.test) * basehaz.cum)
    }
  }else{
    c_fit = NULL
  }

  # use binary data
  tempdat <- data.frame(y, D, w, m_hat, p_hat, c_hat, x)
  binary_data <- tempdat[tempdat$D==1|tempdat$y > times,]          # remove subjects who got censored before the time of interest t50
  binary_data$D[binary_data$D==1 & binary_data$y > times] <- 0     # recode the event status for subjects who had events after t50
  binary_data <- binary_data[complete.cases(binary_data),]
  
  y_tilde = (1 - binary_data$D) - binary_data$m_hat
  w_tilde =  binary_data$w - binary_data$p_hat
  pseudo_outcome = y_tilde/w_tilde
  weights = w_tilde^2/binary_data$c_hat  
  
  tau_dat <- data.frame(pseudo_outcome, binary_data[,7:dim(binary_data)[2]])
  tau_fit <- gbm(pseudo_outcome ~ ., 
                 data = tau_dat, 
                 distribution = "gaussian",
                 weights = weights,
                 shrinkage = shrinkage, 
                 n.trees = n.trees,
                 cv.folds = cv.folds, 
                 n.cores = detectCores()) 

  ret = list(tau_fit = tau_fit,
             pseudo_outcome = pseudo_outcome,
             weights = weights,
             w_fit = w_fit,
             y_fit = y_fit,
             c_fit = c_fit,
             p_hat = p_hat,
             m_hat = m_hat,
             c_hat = c_hat)
  class(ret) = "rgbm"
  ret
}

#' predict for rgbm
#'
#' get estimated tau(x) using the trained rboost model
#'
#' @param object a rboost object
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
predict.rgbm <- function(object,
                         newx = NULL,
                         tau_only = T,
                          ...) {
  if (!is.null(newx)){
    newx = sanitize_x(newx)
  }
  if (tau_only) {
    best_iter_tau <- gbm.perf(object$tau_fit, plot.it = F, method = 'cv')
    return(predict(object$tau_fit, data.frame(newx), n.trees = best_iter_tau))
  } else {
    best_iter_tau <- gbm.perf(object$tau_fit, plot.it = F, method = 'cv')
    tau <- predict(object$tau_fit, data.frame(newx), n.trees = best_iter_tau)
    e = predict(object$w_fit, newx = data.frame(newx))
    m = predict(object$y_fit, newx = data.frame(newx))
    mu1 = m + (1-e) * tau
    mu0 = m - e * tau
    return(list(tau=tau, e=e, m=m, mu1 = mu1, mu0 = mu0))
  }
}
