# -----------------------------------------------------------
# Statistical Computing Experiments
# -----------------------------------------------------------
# EM Algorithm for Multivariate Gaussian Mixture Distribution
# Author: Zilong Liang
# Date: 2018-03-23
# -----------------------------------------------------------

library("mvtnorm")

# -----------------------------
# Main function of EM Algorithm
# -----------------------------

gmm <- function(x, k, tol = 1e-6, iter.max = 200) {
  # Check the samples
  d <- ncol(x)
  n <- nrow(x)
  
  # Initialize the parameters
  tau <- runif(k - 1, 0, 1 / k); tau <- c(tau, 1 - sum(tau))
  mu <- matrix(runif(d * k, 0, 5), nrow = k, ncol = d)
  sigma <- array(rep(diag(d), k), dim = c(d, d, k))
  
  # EM iterations
  iter <- 0  # Iteration index
  repeat {
    # E-Step
    f <- apply(as.matrix(1:k), 1, function(j) {
      return(dmvnorm(x, mu[j, ], sigma[, , j]))
    })
    e.weighted <- f %*% tau
    e <- apply(as.matrix(1:k), 1, function(j) {
      return(f[, j] * tau[j])
    })
    e <- e / rowSums(e.weighted)

    # M-Step
    sum.e <- colSums(e)
    tau.new <- sum.e / n
    mu.new <- t(e) %*% x / sum.e
    sigma.new <- array(rep(0, d * d* k), dim = c(d, d, k))
    for (j in 1:k) {
      for (i in 1:n) {
        sigma.new[, , j] <- sigma.new[, , j] + e[i, j] *
                            ((x[i, ] - mu.new[j, ]) %o% (x[i, ] - mu.new[j, ]))
      }
      sigma.new[, , j] <- sigma.new[, , j] / sum.e[j]
    }
    
    # Judge convergence and iterate the parameters
    iter <- iter + 1
    if (iter > iter.max) { break }
    err.tau <- norm(rbind(tau.new - tau), "I") / rbind(tau)
    err.mu <- norm(mu.new - mu, "I") / norm(mu, "I")
    err.sigma <- norm(colSums(sigma.new - sigma), "I") / norm(colSums(sigma), "I")
    err.max <- max(c(err.tau, err.mu, err.sigma))
    
    tau <- tau.new
    mu <- mu.new
    sigma <- sigma.new
    
    if (err.max < tol) { break }
  }
  
  return (list(tau, mu, sigma, iter))
}


# ----------
# Experiment
# ----------

# Read sample data
x <- read.csv("data.csv")
x <- as.matrix(x[, 2:3])

# Estimate parameters
estimates <- gmm(x, 3)
tau <- estimates[[1]]
mu <- estimates[[2]]
sigma <- estimates[[3]]

# Prepare plotting
pfunc <- function(pp, tau, mu, sigma) {
  # PDF of Gaussian Mixture Distribution
  k = length(tau)
  zz <- 0
  for (j in 1:k) {
    zz <- zz + tau[j] * dmvnorm(pp, mu[j, ], sigma[, , j])
  }
  return (zz)
}
pnum <- 300
xx <- seq(-6, 10, length.out = pnum)
yy <- seq(-6, 5, length.out = pnum)
pp <- cbind(rep(xx, pnum), sort(rep(yy, pnum)))  # Plotting points
zz <- matrix(pfunc(pp, tau, mu, sigma), pnum, pnum)

# Plotting
contour(xx, yy, zz, xlab = "x", ylab = "y", family = "serif", col = "#2fa9df")
points(x, pch = 20, cex = 0.4, col = "#b28fce")
points(estimates[[2]], pch = 20, cex = 1, col ="#4e4f97")
title("Experiment of EM Algorithm", family = "serif")
