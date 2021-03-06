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
                failure.times = NULL, 
                num.trees = 2000,
                sample.weights = NULL,
                clusters = NULL,
                equalize.cluster.weights = FALSE,
                sample.fraction = 0.5,
                mtry = min(ceiling(sqrt(ncol(x)) + 20), ncol(x)),
                min.node.size = 5,
                honesty = TRUE,
                honesty.fraction = 0.5,
                honesty.prune.leaves = TRUE,
                alpha = 0.05,
                imbalance.penalty = 0,
                stabilize.splits = TRUE,
                ci.group.size = 2,
                tune.parameters = "none",
                compute.oob.predictions = TRUE,
                num.threads = NULL,
                seed = runif(1, 0, .Machine$integer.max),
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
  
  if (is.null(p_hat)){
    w_fit <- regression_forest(x, w, num.trees = max(50, num.trees / 4),
                               sample.weights = sample.weights, clusters = clusters,
                               equalize.cluster.weights = equalize.cluster.weights,
                               sample.fraction = sample.fraction, mtry = mtry,
                               min.node.size = 5, honesty = TRUE,
                               honesty.fraction = 0.5, honesty.prune.leaves = TRUE,
                               alpha = alpha, imbalance.penalty = imbalance.penalty,
                               ci.group.size = 1, compute.oob.predictions = TRUE,
                               num.threads = num.threads, seed = seed) 
    p_hat <- predict(w_fit)$predictions
  }else if (length(p_hat) == 1) {
    w_fit = NULL
    p_hat <- rep(p_hat, nrow(x))
  }else if (length(p_hat) != nrow(x)){
    stop("p_hat has incorrect length.")
  }
  
  args.nuisance <- list(failure.times = failure.times,
                        num.trees = max(50, num.trees / 4),
                        sample.weights = sample.weights,
                        clusters = clusters,
                        equalize.cluster.weights = equalize.cluster.weights,
                        sample.fraction = sample.fraction,
                        mtry = mtry,
                        min.node.size = 15,
                        honesty = TRUE,
                        honesty.fraction = 0.5,
                        honesty.prune.leaves = TRUE,
                        alpha = alpha,
                        prediction.type = "Nelson-Aalen",
                        compute.oob.predictions = FALSE,
                        num.threads = num.threads,
                        seed = seed)
  
  if (is.null(m_hat)){
    y_fit <- do.call(survival_forest, c(list(X = cbind(x, w), Y = y, D = D), args.nuisance))
    y_fit[["X.orig"]][, ncol(x) + 1] <- rep(1, nrow(x))
    S1.hat <- predict(y_fit)$predictions
    y_fit[["X.orig"]][, ncol(x) + 1] <- rep(0, nrow(x))
    S0.hat <- predict(y_fit)$predictions
    y_fit[["X.orig"]][, ncol(x) + 1] <- w
    
    times.index <- findInterval(times, y_fit$failure.times)
    surf1 <- S1.hat[, times.index]
    surf0 <- S0.hat[, times.index]
    m_hat  <- p_hat * surf1 + (1 - p_hat) * surf0 
  }else {
    y_fit = NULL
  }
  
  if (is.null(failure.times)) {
    Y.grid <- sort(unique(y))
  } else {
    Y.grid <- failure.times
  }

  args.nuisance$compute.oob.predictions <- TRUE
  if (is.null(c_hat)){
    c_fit <- do.call(survival_forest, c(list(X = cbind(x, w), Y = y, D = 1 - D), args.nuisance))
    C.hat <- predict(c_fit, failure.times = Y.grid)$predictions
    cent <- y
    cent[D==0] <- times
    cen.times.index <- findInterval(cent, Y.grid)
    c_hat <- C.hat[cbind(1:length(y), cen.times.index)]
    c_hat[c_hat==0] <- min(c_hat[c_hat!=0])
  }else{
    c_fit <- NULL
  }
  
  if (any(c_hat <= 0.05)) {
    warning(paste("Estimated censoring probabilites go as low as:", min(c_hat),
                  "- an identifying assumption is that there exists a fixed positive constant M",
                  "such that the probability of observing an event time past the maximum follow-up time Y.max",
                  "is at least M. Formally, we assume: P(Y >= Y.max | X) > M.",
                  "This warning appears when M is less than 0.05, at which point causal survival forest",
                  "can not be expected to deliver reliable estimates."))
  } else if (any(c_hat < 0.2 & c_hat > 0.05)) {
    warning(paste("Estimated censoring probabilites are lower than 0.2.",
                  "An identifying assumption is that there exists a fixed positive constant M",
                  "such that the probability of observing an event time past the maximum follow-up time Y.max",
                  "is at least M. Formally, we assume: P(Y >= Y.max | X) > M."))
  }

  # create binary data
  tempdata <- data.frame(y, D, w, m_hat, p_hat, c_hat, x)
  binary_data <- tempdata[tempdata$D==1|tempdata$y > times,]       # remove subjects who got censored before the time of interest t50
  binary_data$D[binary_data$D==1 & binary_data$y > times] <- 0     # recode the event status for subjects who had events after t50
  binary_data <- binary_data[complete.cases(binary_data),]
  
  y_tilde <- (1 - binary_data$D) - binary_data$m_hat
  w_tilde <-  binary_data$w - binary_data$p_hat
  pseudo_outcome <- y_tilde/w_tilde
  weights <- w_tilde^2/binary_data$c_hat  
  
  tau_dat <- data.frame(pseudo_outcome, binary_data[,7:dim(binary_data)[2]])
  tau_fit <- regression_forest(tau_dat[, 2:dim(tau_dat)[2]], 
                               tau_dat$pseudo_outcome, 
                               num.trees = num.trees,
                               sample.weights = weights,
                               clusters = clusters,
                               sample.fraction = sample.fraction,
                               mtry = mtry,
                               min.node.size = min.node.size,
                               honesty = honesty,
                               honesty.fraction = honesty.fraction,
                               honesty.prune.leaves = honesty.prune.leaves,
                               alpha = alpha,
                               imbalance.penalty = imbalance.penalty,
                               ci.group.size = ci.group.size,
                               compute.oob.predictions = compute.oob.predictions,
                               num.threads = num.threads,
                               seed = seed)

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
