# -----------------------------------------------------------
# Statistical Computing Experiments
# -----------------------------------------------------------
# EM Algorithm for Bivariate Normal Data with Missing Values
# Author: Zilong Liang
# Date: 2018-04-04
# -----------------------------------------------------------

# Samples
x0 <- cbind(c(1, 1, -1, -1), c(1, -1, 1, -1))
x1 <- cbind(c(NA, NA, NA, NA), c(2, 2, -2, -2))
x2 <- cbind(c(2, 2, -2, -2), c(NA, NA, NA, NA))
n <- 12
m1 <- 5
m2 <- 9

# Parameters initialization
mu <- c(1, 1)
sigma <- matrix(c(1, 0, 0, 1), ncol = 2, nrow = 2)
tol <- 1e-10
iter.max <- 500

iter <- 0
repeat {
  # E-Step
  rho <- sigma[1, 2] / sqrt(sigma[1, 1] * sigma[2, 2])
  x1[, 1] <- mu[1] + (sigma[1, 2]/sigma[2, 2]) * (x1[, 2] - mu[2])
  x2[, 2] <- mu[2] + (sigma[1, 2]/sigma[1, 1]) * (x2[, 1] - mu[1])
  x1square <- x1[, 1] ^ 2 + sigma[1, 1] * (1 - rho^2)
  x2square <- x2[, 2] ^ 2 + sigma[2, 2] * (1 - rho^2)
  x <- rbind(x0, x1, x2)
  xsquare <- x ^ 2
  xsquare[m1:(m2-1), 1] <- x1square
  xsquare[m2:n, 2] <- x2square

  # M-Step
  mu.new <- colSums(x) / n
  s11 <- sum(xsquare[, 1]) / n - mu.new[1] ^ 2
  s12 <- sum(x[, 1] * x[, 2]) / n - mu.new[1] * mu.new[2]
  s22 <- sum(xsquare[, 2]) / n - mu.new[2] ^ 2
  sigma.new <- matrix(c(s11, s12, s12, s22), ncol = 2, nrow = 2)

  # Convergence Judgement
  iter <- iter + 1
  if (iter > iter.max) { break }
  err.sigma <- norm(sigma.new - sigma) / norm(sigma)
  if (err.sigma < tol) { break }

  # Parameters iteration
  mu <- mu.new
  sigma <- sigma.new
}

print(mu)
print(sigma)
