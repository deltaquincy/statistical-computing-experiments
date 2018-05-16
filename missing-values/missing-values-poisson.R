# -----------------------------------------------------------
# Statistical Computing Experiments
# -----------------------------------------------------------
# EM Algorithm for Poisson Mixture Distribution
# Author: Zilong Liang
# Date: 2018-04-04
# -----------------------------------------------------------

# -----------------------------
# Main function of EM Algorithm
# -----------------------------

poismm <- function(x, k, tol = 1e-6, iter.max = 200) {
  # Check the samples
  n <- nrow(x)
  
  # Initialize the parameters
  tau <- runif(k - 1, 0, 1 / k); tau <- c(tau, 1 - sum(tau))
  lambda <- runif(k)
  
  # EM iterations
  iter <- 0  # Iteration index
  repeat {
    # E-Step
    f <- matrix(nrow = n, ncol = k)
    for (j in 1:k) {  # TODO: vectorization?
      f[, j] <- dpois(x, lambda[j])
    }
    e.weighted <- f %*% tau
    e <- matrix(0, nrow = n, ncol = k)
    for (j in 1:k) {  # TODO: vectorization?
      e[, j] <- f[, j] * tau[j] 
    }
    e <- e / rowSums(e.weighted)

    # M-Step
    sum.e <- colSums(e)
    tau.new <- sum.e / n
    lambda.new <- c(t(e) %*% x / sum.e)
    
    # Judge convergence
    iter <- iter + 1
    if (iter > iter.max) { break }
    err.tau <- norm(rbind(tau.new - tau), "I") / norm(rbind(tau), "I")
    err.lambda <- norm(rbind(lambda.new - lambda), "I") / 
                  norm(rbind(lambda), "I")
    err.max <- max(c(err.tau, err.lambda))
    if (err.max < tol) { break }
    
    # Iterate the parameters
    tau <- tau.new
    lambda <- lambda.new
  }
  
  return (list(tau, lambda))
}


# ----------
# Experiment
# ----------

# Read sample data
x <- read.csv("data.csv")
x <- as.matrix(x[, 2])

# Estimate parameters
estimates <- poismm(x, 2)
tau <- estimates[[1]]
lambda <- estimates[[2]]

# Prepare plotting
pfunc <- function(xx, tau, lambda) {
  # PDF of Gaussian Mixture Distribution
  k = length(tau)
  yy <- 0
  for (j in 1:k) {
    yy <- yy + tau[j] * dpois(xx, lambda[j])
  }
  return (yy)
}

# Plotting
xx <- seq(0, 11)
yy = pfunc(xx, tau, lambda)
hist(x, freq = FALSE,
     breaks = c(-0.5:11.5),
     xlab = "k", 
     ylab =  "Density or Possibility", 
     main = "EM Algorithm on Poisson Mixture Distribution",
     family = "serif")
points(xx, yy, pch = 18, cex = 1.5, col = "#b28fce")
