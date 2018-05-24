# -----------------------------------------------------------
# Statistical Computing Experiments
# -----------------------------------------------------------
# MCMC Simulation
# Author: Zilong Liang
# Date: 2018-05-23
# -----------------------------------------------------------

# Clear
rm(list = ls())

# ---------------
# Helper function
# ---------------

# Continuous random variable generator
# Acceptance-rejection technique
generator.continuous.accrej <- function(y.gen, Pf, Qf, c, n) {
  # Initialize
  i <- 0
  x <- rep(0, n)
  
  # Produce
  while (i < n) {
    y <- y.gen(1)
    u <- runif(1)
    if (u < Pf(y) / (c * Qf(y))) {
      i <- i + 1
      x[i] <- y
    }
  }
  
  return (x)
}

# ---------------
# MCMC simulators
# ---------------

# Gibbs sampler
sampler.gibbs <- function(r.list, range.list, k, n) {
  # Initialize
  nn <- n + floor(0.3 * n)  # Generate more and pick the end
  x <- matrix(0, ncol = k, nrow = nn)
  for (j in 1:k) {
    x[1, j] <- runif(1, min = range.list[[j]][1], max = range.list[[j]][2])
  }
  
  # Sample
  for (i in 2:nn) {
    x[i, ] <- x[i-1, ]
    for (j in 1:k) {
      x[i, j] <- r.list[[j]](x[i, ])
    }
  }
  return (x[(nn-n+1):nn, ])
}

# Metropolis sampler
sampler.metropolis <- function(x0, f, n, walk.gen, 
                               proposal = function(x1, x2) {return (1)}) {
  # Check inputs and initialize
  k <- length(x0)
  nn <- n + floor(0.5 * n)  # Generate more and pick the end
  x <- matrix(0, ncol = k, nrow = nn)
  x[1, ] <- x0
  
  # Sample
  for (i in 2:nn) {
    x.new <- x[i-1, ] + walk.gen()
    alpha <- (f(x.new) * proposal(x.new, x[i-1, ])) /
             (f(x[i-1, ]) * proposal(x[i-1, ], x.new))
    if (runif(1) <= alpha) {
      x[i, ] <- x.new
    } else {
      x[i, ] <- x[i-1, ]
    }
  }
  
  return (x[(nn-n+1):nn, ])
}

# ------------------------------
# Experiments in the assignments
# ------------------------------

# Generation size
n <- 3000

# Assignment 1
print("Processing: assignment 1")
B <- 10
range.list <- list(c(0, B), c(0, B))
y.gen1 <- function(n) {
  return (runif(n, min = 0, max = B))
}
Qf1 <- function(x) {
  return (dunif(x, min = 0, max = B))
}
r1 <- function(x.array) {
  y <- x.array[2]
  Pf1 <- function(x) {
    return ((-1/y * (1 - exp(-B*y))) * exp(-x * y))
  }
  c1 <- (-1/y * (1 - exp(-B*y))) * B
  return (generator.continuous.accrej(y.gen1, Pf1, Qf1, c1, 1))
}
y.gen2 <- function(n) {
  return (runif(n, min = 0, max = B))
}
Qf2 <- function(x) {
  return (dunif(x, min = 0, max = B))
}
r2 <- function(x.array) {
  x <- x.array[1]
  Pf2 <- function(y) {
    return ((-1/x * (1 - exp(-B*x))) * exp(-x * y))
  }
  c2 <- (-1/x * (1 - exp(-B*x))) * B
  return (generator.continuous.accrej(y.gen2, Pf2, Qf2, c2, 1))
}
r.list <- list(r1, r2)
k <- 2
x <- sampler.gibbs(r.list, range.list, k, n)
Ex <- apply(x, 2, mean)[1]
Exy <- mean(x[, 1] * x[, 2])
print(sprintf("E(X):  %.6f", Ex))
print(sprintf("E(XY): %.6f", Exy))
pdf("pic01.pdf", width = 7, height = 5)
  plot(x, pch = 20, cex = 0.4, col = "#6497b1",
       xlab = "x", ylab = "y")
  title("Assignment 1")
dev.off()
print("The plot has been saved as 'pic01.pdf'")

# Experiment 2
print("Processing: assignment 2")
x <- matrix(rexp(3 * n), ncol = 3, nrow = n)
x <- apply(x, 1, function(xx) {
  return (xx[1] + 2*xx[2] + 3*xx[3])
})
x1 <- x[x > 15]
E1 <- mean(x1)
x2 <- x[x < 1]
E2 <- mean(x2)
print(sprintf("E(>15): %.6f", E1))
print(sprintf("E(<1):  %.6f", E2))


# Experiment 3
print("Processing: assignment 3")
k <- 3
x0 <- rexp(k)
walk.gen <- function() {return (rnorm(k))}
f <- function(x) {
  return (ifelse(all(x > 0),
                 exp(-(sum(x)+x[1]*x[2]+x[1]*x[3]+x[2]*x[3])),
                 0))
}
x <- sampler.metropolis(x0, f, n, walk.gen)
Exyz <- mean(apply(x, 1, prod))
print(sprintf("E(XYZ): %.6f", Exyz))
pdf("pic03.pdf", width = 13, height = 5)
  par(mfrow = c(1, 3))
  hist(x[, 1], xlab = "x", ylab = "frequency", 
       freq = FALSE, main = "Assignment 3 (x)",
       col = "#005b96")
  hist(x[, 2], xlab = "y", ylab = "frequency", 
       freq = FALSE, main = "Assignment 3 (y)",
       col = "#6497b1")
  hist(x[, 3], xlab = "z", ylab = "frequency", 
       freq = FALSE, main = "Assignment 3 (z)",
       col = "#011f4b")
dev.off()
print("The plot has been saved as 'pic03.pdf'")

# Experiment 4
print("Processing: assignment 4")
k <- 3
x0 <- c(2, 0.5, 3)  # Fixed initialization
walk.gen <- function() {
  rx <- ifelse(runif(1) <= 1/2, 1, -1)
  ry <- runif(1, min = -0.2, max = 0.2)
  rz <- ifelse(runif(1) <= 1/2, 1, -1)

  return (c(rx, ry, rz))
}
f <- function(x) {
  i <- x[1]
  y <- x[2]
  n <- x[3]
  if (i < 0 || i > n || n < 0 || y < 0) {
    return (0)
  } else {
    return (choose(n, i) * y^(i+1) * (1-y)^(n-i+2) * dpois(n, 4))
  }
}
x <- sampler.metropolis(x0, f, n, walk.gen)
E <- apply(x, 2, mean)
print(sprintf("E(X): %.6f", E[1]))
print(sprintf("E(Y): %.6f", E[2]))
print(sprintf("E(Z): %.6f", E[3]))
pdf("pic04.pdf", width = 13, height = 5)
  par(mfrow = c(1, 3))
  barplot(table(x[, 1])/length(x[, 1]),
          xlab = "x", ylab = "frequency", 
          main = "Assignment 4 (x)",
          col = "#83adb5")
  hist(x[, 2], xlab = "y", ylab = "frequency", 
       freq = FALSE, main = "Assignment 4 (y)",
       col = "#c7bbc9")
  barplot(table(x[, 3])/length(x[, 3]),
          xlab = "z", ylab = "frequency", 
          main = "Assignment 4 (z)",
          col = "#bfb5b2")
dev.off()
print("The plot has been saved as 'pic04.pdf'")

# Assignment 5
print("Processing: assignment 5")
k <- 2
x0 <- rnorm(k)
walk.gen <- function() {return (rnorm(k))}
mu1 <- c(1, 4)
mu2 <- c(-2, -1)
sigma1 <- matrix(c(1, 0.3, 0.3, 2), c(2, 2))
sigma2 <- matrix(c(6, 0.9, 0.9, 1), c(2, 2))
sigma1.inv <- solve(sigma1)
sigma2.inv <- solve(sigma2)
sigma1.det <- sqrt(det(sigma1.inv))
sigma2.det <- sqrt(det(sigma2.inv))
f <- function(theta) {
  b1 <- (theta - mu1) %*% sigma1.inv %*% (theta - mu1)
  b2 <- (theta - mu2) %*% sigma2.inv %*% (theta - mu2)
  r <- sigma1.det * exp(-1/2 * b1)
  r <- r + sigma2.det * exp(-1/2 * b2)
  r <- r / (4 * pi)
  return (r)
}
x <- sampler.metropolis(x0, f, n, walk.gen)
pdf("pic05.pdf", width = 7, height = 7)
  # Plotting points
  pnum <- 300
  xx <- seq(-8, 5, length.out = pnum)
  yy <- seq(-5, 9, length.out = pnum)
  pp <- cbind(rep(xx, pnum), sort(rep(yy, pnum)))
  zz <- matrix(apply(pp, 1, f), pnum, pnum)
  # Plotting
  contour(xx, yy, zz,
          xlab = "x", ylab = "y", col = "#7670ab")
  points(x, pch = 20, cex = 0.4, col = "#54b2a9")
  points(c(mu1[1], mu2[1]), c(mu1[2], mu2[2]),
         pch = 20, cex = 2, col = "#7670ab")
  title("Assignment 5")
dev.off()
print("The plot has been saved as 'pic05.pdf'")