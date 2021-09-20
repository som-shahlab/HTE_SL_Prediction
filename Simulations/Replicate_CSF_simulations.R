# Simulations based on the four scenarios in the causal survival forest paper 
# Ref: Estimating heterogeneous treatment effects with right-censored data via causal survival forests (Cui 2020)

nsim <- 500

## Scenario 1: AFT for T and Cox PH for C
set.seed(3)
n <- 1600
p <- 3
X <- matrix(runif(n * p), n, p)
#W <- rbinom(n, 1, 0.5)
beta <- dbeta(X[, 1], 2, 4)
pi <- (1 + beta) / 4
W <- rbinom(n, 1, pi)
epsilon <- rbinom(n, 0, 1)
tau <- 1.5

# Failure time 
failure.time.log <- -1.85 - 0.8 * (X[, 1] < 0.5) + 0.7 * X[, 2]^(0.5) + 0.2 * X[, 3] + (0.7 - 0.4 * (X[, 1] < 0.5) - 0.4 * X[, 2]^(0.5)) * W + epsilon
failure.time <- exp(failure.time.log)
failure.time.r <- pmin(failure.time + W, tau)

# Censoring time 
log.risk <- -1.75 - 0.5 * X[, 2]^(0.5) + 0.2 * X[, 3] + (1.15 + 0.5 * (X[, 1] < 0.5) - 0.3 * X[, 2]^(0.5)) * W
v <- runif(n)
censor.time <- (-1 * log(v) / exp(log.risk))^(0.5)  # H0 inverse function = t^(1/2)

# Observed data 
Y <- pmin(failure.time.r, censor.time)
D <- as.integer(failure.time.r <= censor.time)
sim_dat1 <- data.frame(X, Y, W, D)

## Scenario 2: Cox PH for T and AFT for C
set.seed(3)
n <- 1600
p <- 3
X <- matrix(runif(n * p), n, p)
beta <- dbeta(X[, 2], 2, 4)
pi <- (1 + beta) / 4
W <- rbinom(n, 1, pi)
epsilon <- rbinom(n, 0, 1)
tau <- 2

# Failure time 
log.risk <- X[, 1] + (-0.4 + X[, 2] * W)
v <- runif(n)
failure.time <- (-1 * log(v) / exp(log.risk))^2  # H0 inverse function = t^2
failure.time.r <- pmin(failure.time + W, tau)

# Censoring time 
censor.time.log <- X[, 1] - X[, 3] * W + epsilon
censor.time <- exp(censor.time.log)

# Observed data 
Y <- pmin(failure.time.r, censor.time)
D <- as.integer(failure.time.r <= censor.time)
sim_dat2 <- data.frame(X, Y, W, D)


## Scenario 3: Poisson dist. for both T and C 
set.seed(3)
n <- 1600
p <- 3
X <- matrix(runif(n * p), n, p)
beta <- dbeta(X[, 1], 2, 4)
pi <- (1 + beta) / 4
W <- rbinom(n, 1, pi)
epsilon <- rbinom(n, 0, 1)
tau <- 15

# Failure time
mean.T <- X[, 2]^2 + X[, 3] + 6 + 2 * (X[, 1]^(1/2) - 0.3) * W
failure.time <- rpois(n, mean.T)
failure.time.r <- pmin(failure.time + W, tau)

# Censoring time 
mean.C <- 12 + log(1 + exp(X[, 3]))
censor.time <- rpois(n, mean.C)

# Observed data 
Y <- pmin(failure.time.r, censor.time)
D <- as.integer(failure.time.r <= censor.time)
sim_dat3 <- data.frame(X, Y, W, D)

## Scenario 4: Poisson dist. for both T and C 
set.seed(3)
n <- 1600
p <- 3
X <- matrix(runif(n * p), n, p)
pi <- ((1 + exp(-1 * X[, 1])) * (1 + exp(-1 * X[, 2])))^(-1)
W <- rbinom(n, 1, pi)
epsilon <- rbinom(n, 0, 1)
tau <- 5

# Failure time
mean.T <- pmax(0, X[, 2] + X[, 3]) + pmax(0, X[, 1] - 0.3) * W
failure.time <- rpois(n, mean.T)
failure.time.r <- pmin(failure.time + W, tau)

# Censoring time 
mean.C <- 1 + log(1 + exp(X[, 3]))
censor.time <- rpois(n, mean.C)

# Observed data 
Y <- pmin(failure.time.r, censor.time)
D <- as.integer(failure.time.r <= censor.time)
sim_dat4 <- data.frame(X, Y, W, D)


library(grf)
min.node.sizes <- c(15, ceiling(p^(1/2)), ceiling(log(p)))
mtrys <- c(ceiling(p/3), ceiling(2*p/3), ceiling(p^(1/2)))
cs.forest <- causal_survival_forest(X, Y, W, D, num.trees=500)














