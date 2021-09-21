# Simulations based on the four scenarios in the causal survival forest paper 
# Ref: Estimating heterogeneous treatment effects with right-censored data via causal survival forests (Cui 2020)
library(grf)
n <- 4000
p <- 3
nsim <- 200
nmethods <- 8

## Scenario 1: AFT for T and Cox PH for C
mse <- cr <- matrix(NA, nsim, nmethods)
for (i in 1: nsim){
  X <- matrix(runif(n * p), n, p)
  beta <- dbeta(X[, 1], 2, 4)
  pi <- (1 + beta) / 4
  W <- rbinom(n, 1, pi)
  epsilon <- rbinom(n, 0, 1)
  tau <- 1.5
  
  # Failure time 
  failure.time.log <- -1.85 - 0.8 * (X[, 1] < 0.5) + 0.7 * X[, 2]^(0.5) + 0.2 * X[, 3] + (0.7 - 0.4 * (X[, 1] < 0.5) - 0.4 * X[, 2]^(0.5)) * W + epsilon
  failure.time <- exp(failure.time.log)
  failure.time.r <- pmin(failure.time, tau)
  
  # Truth 
  failure.time.r1 <- pmin(exp(-1.85 - 0.8 * (X[, 1] < 0.5) + 0.7 * X[, 2]^(0.5) + 0.2 * X[, 3] + (0.7 - 0.4 * (X[, 1] < 0.5) - 0.4 * X[, 2]^(0.5)) + epsilon), tau)
  failure.time.r0 <- pmin(exp(-1.85 - 0.8 * (X[, 1] < 0.5) + 0.7 * X[, 2]^(0.5) + 0.2 * X[, 3] + epsilon), tau)
  delta0 <- failure.time.r1 - failure.time.r0
  
  # Censoring time 
  log.risk <- -1.75 - 0.5 * X[, 2]^(0.5) + 0.2 * X[, 3] + (1.15 + 0.5 * (X[, 1] < 0.5) - 0.3 * X[, 2]^(0.5)) * W
  v <- runif(n)
  censor.time <- (-1 * log(v) * exp(-1*log.risk))^(0.5)  # H0 inverse function = t^(1/2)

  # Observed data 
  Y <- pmin(failure.time.r, censor.time)
  D <- as.integer(failure.time.r <= censor.time)
  sim_dat1 <- data.frame(X, Y, W, D, delta0)
  
  index <- sample(1:n, 2000, replace = F)
  sim_dat1_train <- sim_dat1[index,]
  sim_dat1_test <- sim_dat1[-index,]
  
  cs.forest <- causal_survival_forest(sim_dat1_train[,1:3], sim_dat1_train$Y, sim_dat1_train$W, sim_dat1_train$D)
  out <- predict(cs.forest, sim_dat1_test[,1:3])$predictions
  
  # MSE 
  mse[i,1] <- sum((out - sim_dat1_test$delta0)^2)/(n/2)
  
  # classification error 
  cr[i,1] <- 1-sum((sign(out)==sign(sim_dat1_test$delta0)))/(n/2)
}

# MSE 
mse.ave <- colMeans(mse)

# Excess MSE 


# classification error 
cr.ave <- colMeans(cr)



## Scenario 2: Cox PH for T and AFT for C
mse <- cr <- matrix(NA, nsim, nmethods)
for (i in 1:nsim){
  X <- matrix(runif(n * p), n, p)
  beta <- dbeta(X[, 2], 2, 4)
  pi <- (1 + beta) / 4
  W <- rbinom(n, 1, pi)
  epsilon <- rbinom(n, 0, 1)
  tau <- 2
  
  # Failure time 
  log.risk <- X[, 1] + (-0.5 + X[, 2]) * W
  v <- runif(n)
  failure.time <- (-1 * log(v) / exp(log.risk))^2  # H0 inverse function = t^2
  failure.time.r <- pmin(failure.time, tau)
  
  # Truth 
  failure.time.r1 <- pmin((-1 * log(v) / exp(X[, 1] + (-0.5 + X[, 2]) * 1))^2, tau)
  failure.time.r0 <- pmin((-1 * log(v) / exp(X[, 1] + (-0.5 + X[, 2]) * 0))^2, tau)
  delta0 <- failure.time.r1 - failure.time.r0
  
  # Censoring time 
  censor.time <- runif(n, 0, 3)
  
  # Observed data 
  Y <- pmin(failure.time.r, censor.time)
  D <- as.integer(failure.time.r <= censor.time)
  sim_dat2 <- data.frame(X, Y, W, D, delta0)
  
  index <- sample(1:n, 2000, replace = F)
  sim_dat2_train <- sim_dat2[index,]
  sim_dat2_test <- sim_dat2[-index,]
  
  cs.forest <- causal_survival_forest(sim_dat2_train[,1:3], sim_dat2_train$Y, sim_dat2_train$W, sim_dat2_train$D)
  out <- predict(cs.forest, sim_dat2_test[,1:3])$predictions
  
  # MSE 
  mse[i,1] <- sum((out - sim_dat2_test$delta0)^2)/(n/2) # 0.06
  
  # classification error 
  cr[i,1] <- 1-sum((sign(out)==sign(sim_dat2_test$delta0)))/(n/2)  # 0.137
}


## Scenario 3: Poisson dist. for both T and C 
mse <- cr <- matrix(NA, nsim, nmethods)
for (i in 1:nsim){
  X <- matrix(runif(n * p), n, p)
  beta <- dbeta(X[, 1], 2, 4)
  pi <- (1 + beta) / 4
  W <- rbinom(n, 1, pi)
  epsilon <- rbinom(n, 0, 1)
  tau <- 15
  
  # Failure time
  mean.T <- X[, 2]^2 + X[, 3] + 6 + 2 * (X[, 1]^(1/2) - 0.3) * W
  failure.time <- rpois(n, mean.T)
  failure.time.r <- pmin(failure.time, tau)
  
  # Truth 
  failure.time.r1 <- pmin(rpois(n, X[, 2]^2 + X[, 3] + 6 + 2 * (X[, 1]^(1/2) - 0.3) * 1), tau)
  failure.time.r0 <- pmin(rpois(n, X[, 2]^2 + X[, 3] + 6 + 2 * (X[, 1]^(1/2) - 0.3) * 0), tau)
  delta0 <- failure.time.r1 - failure.time.r0
  
  # Censoring time 
  mean.C <- 12 + log(1 + exp(X[, 3]))
  censor.time <- rpois(n, mean.C)
  
  # Observed data 
  Y <- pmin(failure.time.r, censor.time)
  D <- as.integer(failure.time.r <= censor.time)
  sim_dat3 <- data.frame(X, Y, W, D, delta0)
  
  index <- sample(1:n, 2000, replace = F)
  sim_dat3_train <- sim_dat3[index,]
  sim_dat3_test <- sim_dat3[-index,]
  
  cs.forest <- causal_survival_forest(sim_dat3_train[,1:3], sim_dat3_train$Y, sim_dat3_train$W, sim_dat3_train$D)
  out <- predict(cs.forest, sim_dat3_test[,1:3])$predictions
  
  # MSE 
  mse[i,1] <- sum((out - sim_dat3_test$delta0)^2)/(n/2)  # 15.19
  
  # classification error 
  cr[i,1] <- 1-sum((sign(out[sim_dat3_test$delta0!=0])==sign(sim_dat3_test$delta0[sim_dat3_test$delta0!=0])))/(length(out[sim_dat3_test$delta0!=0])) 
}



## Scenario 4: Poisson dist. for both T and C 
mse <- cr <- matrix(NA, nsim, nmethods)
for (i in 1:nsim){
  X <- matrix(runif(n * p), n, p)
  pi <- ((1 + exp(-1 * X[, 1])) * (1 + exp(-1 * X[, 2])))^(-1)
  W <- rbinom(n, 1, pi)
  epsilon <- rbinom(n, 0, 1)
  tau <- 3
  
  # Failure time
  mean.T <- X[, 2] + X[, 3] + pmax(0, X[, 1] - 0.3) * W
  failure.time <- rpois(n, mean.T)
  failure.time.r <- pmin(failure.time, tau)
  
  # Truth 
  failure.time.r1 <- pmin(rpois(n, X[, 2] + X[, 3] + pmax(0, X[, 1] - 0.3) * 1), tau)
  failure.time.r0 <- pmin(rpois(n, X[, 2] + X[, 3] + pmax(0, X[, 1] - 0.3) * 0), tau)
  delta0 <- failure.time.r1 - failure.time.r0
  
  # Censoring time 
  mean.C <- 1 + log(1 + exp(X[, 3]))
  censor.time <- rpois(n, mean.C)
  
  # Observed data 
  Y <- pmin(failure.time.r, censor.time)
  D <- as.integer(failure.time.r <= censor.time)
  sim_dat4 <- data.frame(X, Y, W, D, delta0)
  
  index <- sample(1:n, 2000, replace = F)
  sim_dat4_train <- sim_dat4[index,]
  sim_dat4_test <- sim_dat4[-index,]
  
  cs.forest <- causal_survival_forest(sim_dat4_train[,1:3], sim_dat4_train$Y, sim_dat4_train$W, sim_dat4_train$D)
  out <- predict(cs.forest, sim_dat4_test[,1:3])$prediction
  
  # MSE 
  mse[i,1] <- sum((out - sim_dat4_test$delta0)^2)/(n/2)
  
  # classification error 
  # how to deal with delta0 = 0
  cr[i,1] <- 1-sum((sign(out)==sign(sim_dat4_test$delta0)))/(n/2)
  
  cr[i,1] <- 1-sum((sign(out[sim_dat4_test$delta0!=0])==sign(sim_dat4_test$delta0[sim_dat4_test$delta0!=0])))/(length(out[sim_dat4_test$delta0!=0]))
}

head(cbind(out, sim_dat4_test$delta0),30)

# min.node.sizes <- c(15, ceiling(p^(1/2)), ceiling(log(p)))
# mtrys <- c(ceiling(p/3), ceiling(2*p/3), ceiling(p^(1/2)))
# cs.forest <- causal_survival_forest(X, Y, W, D, num.trees=500)














