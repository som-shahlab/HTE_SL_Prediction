#' @include utils.R
#'
#' @title S-learner, implemented via glmnet (lasso)
#'
#' @description  S-learner, as proposed by Imai and Ratkovic (2013), implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds number of folds for cross validation
#' @param foldid user-supplied foldid. Must have length equal to length(w). If provided, it overrides the k_folds option.
#' @param lambda user-supplied lambda sequence for cross validation
#' @param lambda_choice how to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty_factor user-supplied penalty factor, must be of length the same as number of features in x
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' slasso_fit = slasso(x, w, y)
#' slasso_est = predict(slasso_fit, x)
#' }
#' @export
slasso_surv = function(x, w, y, D, times,
                       alpha = 1,
                       k_folds = NULL,
                       foldid = NULL,
                       lambda = NULL,
                       lambda_choice = c("lambda.min", "lambda.1se"),
                       penalty_factor = NULL){

  input = sanitize_input(x,w,y,D)
  x = input$x
  w = input$w
  y = input$y
  D = input$D
  
  x_int = data.frame(as.numeric(w - 0.5) * x)
  names(x_int) = paste0(names(data.frame(x)), "_int") 
  x_all = cbind(as.matrix(x_int), x)
  standardization = caret::preProcess(x_all, method=c("center", "scale")) # get the standardization params
  x_scl = predict(standardization, x_all)							 # standardize the input
  x_scl = x_scl[,!is.na(colSums(x_scl)), drop = FALSE]
  x_scl_tilde = cbind(w - 0.5, x_scl)
  
  x_int1 = data.frame(0.5 * x)
  names(x_int1) = paste0(names(data.frame(x)), "_int")
  x_all1 = cbind(as.matrix(x_int1), x)
  standardization = caret::preProcess(x_all1, method=c("center", "scale")) # get the standardization params
  x_scl_1 = predict(standardization, x_all1)							 # standardize the input
  x_scl_1 = x_scl_1[,!is.na(colSums(x_scl_1)), drop = FALSE]
  x_scl_tilde_1 = cbind(0.5, x_scl_1)

  x_int0 = data.frame(-0.5 * x)
  names(x_int0) = paste0(names(data.frame(x)), "_int")
  x_all0 = cbind(as.matrix(x_int0), x)
  standardization = caret::preProcess(x_all0, method=c("center", "scale")) # get the standardization params
  x_scl_0 = predict(standardization, x_all0)							 # standardize the input
  x_scl_0 = x_scl_0[,!is.na(colSums(x_scl_0)), drop = FALSE]
  x_scl_tilde_0 = cbind(-0.5, x_scl_0)
  
  lambda_choice = match.arg(lambda_choice)
  
  nobs = nrow(x_scl)
  pobs = ncol(x_scl)
  
  
  if (is.null(foldid) || length(foldid) != length(w)) {
    
    if (!is.null(foldid) && length(foldid) != length(w)) {
      warning("supplied foldid does not have the same length ")
    }
    
    if (is.null(k_folds)) {
      k_folds = floor(max(3, min(10,length(w)/4)))
    }
    
    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(k_folds), length = length(w)))
    
  }
  
  if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
    if (!is.null(penalty_factor) && length(penalty_factor) != pobs + 1) {
      warning("penalty_factor supplied is not 1 plus 2 times the number of columns in x. Using all ones instead.")
    }
    penalty_factor = c(0, rep(1, pobs))
  }

  s_fit <- glmnet::cv.glmnet(x_scl_tilde, 
                             Surv(y, D),
                             family = "cox",
                             foldid = foldid, 
                             lambda = lambda,
                             penalty.factor = penalty_factor,
                             standardize = FALSE, 
                             alpha = alpha)

  s_beta <- as.vector(t(coef(s_fit, s = lambda_choice)))
  
  link1 <- exp(x_scl_tilde_1 %*% s_beta)
  link0 <- exp(x_scl_tilde_0 %*% s_beta)
  
  d <- data.frame(table(y[D == 1]))[,2]  # number of events at each unique time                               
  h0 <- rep(NA, length(sort(unique(y[D == 1]))))
  for(l in 1:length(sort(unique(y[D == 1])))){
    h0[l] <- d[l] / sum(exp(x_scl_tilde[y >= sort(unique(y[D == 1]))[l], ] %*% s_beta))
  }    
  S0 <- exp(-cumsum(h0))
  S0_t <- S0[sort(unique(y[D == 1]))>=times][1]
  
  surv1 <- S0_t^exp(link1)
  surv0 <- S0_t^exp(link0)
 
  tau_hat <- as.numeric(surv1 - surv0)  

  ret = list(s_fit = s_fit,
             x_org = x_scl_tilde, 
             y_org = y, 
             D_org = D,
             s_beta = s_beta,
             tau_hat = tau_hat,
             lambda_choice = lambda_choice)

  class(ret) <- "slasso_surv"
  ret
}

#' predict for slasso
#'
#' get estimated tau(x) using the trained slasso model
#'
#' @param object a slasso object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
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
#' slasso_fit = slasso(x, w, y)
#' slasso_est = predict(slasso_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.slasso_surv <- function(object,
                                newx = NULL,
                                times,
                                ...) {
  if (!is.null(newx)) {
    newx = sanitize_x(newx)
    x_int1 <- cbind(newx, newx)
    x_int0 <- cbind(0 * newx, newx)
    x_pred1 = predict(object$standardization, x_int1)
    x_pred1 = x_pred1[,!is.na(colSums(x_pred1)), drop = FALSE]
    x_pred0 = predict(object$standardization, x_int0)
    x_pred0 = x_pred0[,!is.na(colSums(x_pred0)), drop = FALSE]
    newx_scl_pred1 = cbind(1, x_pred1)
    newx_scl_pred0 = cbind(0, x_pred0)

    link1 <- exp(newx_scl_pred1 %*% object$s_beta)
    link0 <- exp(newx_scl_pred0 %*% object$s_beta)
    
    d <- data.frame(table(object$y_org[object$D_org == 1]))[,2]                                
    h0 <- rep(NA, length(sort(unique(object$y_org[object$D_org == 1]))))
    for(l in 1:length(sort(unique(object$y_org[object$D_org == 1])))){
      h0[l] <- d[l] / sum(exp(object$x_org[object$y_org >= sort(unique(object$y_org[object$D_org == 1]))[l], ] %*% object$s_beta))
    }    
    S0 <- exp(-cumsum(h0))
    S0_t <- S0[sort(unique(object$y_org[object$D_org == 1]))>=times][1]
    
    surv1 <- S0_t^exp(link1)
    surv0 <- S0_t^exp(link0)
    
    tau_hat <- as.numeric(surv1 - surv0)
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}
