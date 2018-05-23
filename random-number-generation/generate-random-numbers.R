# -----------------------------------------------------------
# Statistical Computing Experiments
# -----------------------------------------------------------
# Generate Random Numbers
# Author: Zilong Liang
# Date: 2018-05-22
# -----------------------------------------------------------

# Clear
rm(list = ls())

# --------------------------
# Random variable generators
# --------------------------

# Discrete random variable generator
# Inverse transform method
generator.discrete.inv <- function(range, P, n) {
  # Check inputs
  k <- length(P) - 1
  P.cum <- cumsum(P)
  
  # Initialize
  r <- runif(n)
  x <- rep(0, n)
  
  # Generate
  x[r < P.cum[1]] <- range[1]
  for (i in 1:k) {
    x[r >= P.cum[i] & r < P.cum[i+1]] <- range[i+1]
  }
  
  return (x)
}

# Discrete random variable generator
# Acceptance-rejection technique
generator.discrete.accrej <- function(y.gen, P, Q, n) {
  # Check inputs
  c <- max(P / Q)
  
  # Initialize
  i <- 0
  x <- rep(0, n)
  
  # Produce
  while (i < n) {
    y <- y.gen(1)
    u <- runif(1)
    if (u < P[y+1] / (c * Q[y+1])) {
      i <- i + 1
      x[i] <- y
    }
  }
  
  return (x)
}

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

# ------------------------------
# Experiments in the assignments
# ------------------------------

# Hyper-parameters
n <- 5000

# Assignment 1
print("Processing: assignment 1")
lambda <- 2
k <- 15
range <- 0:k
P <- dpois(0:k, lambda = lambda)
P <- P / sum(P)
y.gen <- function(n) {
  return (sample(range, n))
}
Q <- rep(1/(k+1), k+1)
x1 <- generator.discrete.inv(range, P, n)
x2 <- generator.discrete.accrej(y.gen, P, Q, n)
pdf("pic01.pdf", width = 7, height = 5)
  data1 <- prop.table(table(x1))
  data2 <- prop.table(table(x2))
  min.num <- min(length(data1), length(data2))
  data <- rbind(data1[1:min.num], data2[1:min.num])
  barplot(data, beside = TRUE, 
          xlab = "n", ylab = "frequency", 
          legend = c("Inverse transform", "Acceptance-rejection"),
          col = c("#37629c", "#8da1e7"))
  lines(seq(2, 30, 3), P[1:10], lty = 5, lwd = 2)
  title("Assignment 1")
dev.off()
print("The plot has been saved as 'pic01.pdf'")


# Assignment 2
print("Processing: assignment 2")
range <- 5:14
P <- rep(c(0.11, 0.09), 5)
x <- generator.discrete.inv(range, P, n)
data <- prop.table(table(x))
pdf("pic02.pdf", width = 7, height = 5)
  barplot(data,
          xlab = "n", ylab = "frequency",
          col = "#37629c")
  lines(seq(0.7, 11.5, length.out = 10), P, lty = 5, lwd = 2)
  title("Assignment 2")
dev.off()
print("The plot has been saved as 'pic02.pdf'")

# Assignment 3
print("Processing: assignment 3")
y.gen <- rnorm
Qf <- dnorm
Pf <- function(x) {
  return (ifelse(x < 0, exp(2 * x), exp(-2 * x)))
}
c <- sqrt(2 * pi)
x <- generator.continuous.accrej(y.gen, Pf, Qf, c, n)
pdf("pic03.pdf", width = 7, height = 5)
  hist(x, freq = FALSE,
       xlab = "x", ylab = "density",
       col = "#b4eaea",
       main = "Assignment 3")
  min.x <- min(x)
  max.x <- max(x)
  xx <- seq(min.x, max.x, length.out = 300)
  yy <- Pf(xx)
  lines(xx, yy, lty = 5, lwd = 2, col = "#37629c")
dev.off()
print("The plot has been saved as 'pic03.pdf'")

# Assignment 4
print("Processing: assignment 4")
y.gen <- runif
Qf <- dunif
Pf <- function(x) {
  return (ifelse(x < 0 | x > 1, 0, 30 * (x^2 - 2*x^3 + x^4)))
}
c <- 15 / 8
x <- generator.continuous.accrej(y.gen, Pf, Qf, c, n)
pdf("pic04.pdf", width = 7, height = 5)
  hist(x, freq = FALSE,
       xlab = "x", ylab = "density",
       col = "#8da1e7",
       main = "Assignment 4")
  min.x <- min(x)
  max.x <- max(x)
  xx <- seq(min.x, max.x, length.out = 300)
  yy <- Pf(xx)
  lines(xx, yy, lty = 5, lwd = 2, col = "#1f49a1")
dev.off()
print("The plot has been saved as 'pic04.pdf'")

# Assignment 5
print("Processing: assignment 4")
y.gen <- function(n) {
  return (rexp(n, rate = 1/3))
}
Qf <- function(x) {
  return (dexp(x, rate = 1/3))
}
Pf <- function(x) {
  return (ifelse(x >= 0, 1/2 * x^2 * exp(-x), 0))
}
c <- 27 / (2 * exp(2))
x <- generator.continuous.accrej(y.gen, Pf, Qf, c, n)
pdf("pic05.pdf", width = 7, height = 5)
  hist(x, freq = FALSE,
       xlab = "x", ylab = "density",
       col = "#68b6c1",
       main = "Assignment 5")
  min.x <- min(x)
  max.x <- max(x)
  xx <- seq(min.x, max.x, length.out = 300)
  yy <- Pf(xx)
  lines(xx, yy, lty = 5, lwd = 2, col = "#166b46")
dev.off()
print("The plot has been saved as 'pic05.pdf'")
