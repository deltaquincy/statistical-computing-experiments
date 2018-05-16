# -----------------------------------------------------------
# Statistical Computing Experiments
# -----------------------------------------------------------
# Exercise on Hidden Markov Model
# Author: Zilong Liang
# Date: 2018-04-27
# -----------------------------------------------------------

# The implementation of HMM slightly refers to the 'HMM'
# library in CRAN, on data structures and tricks to avoid
# underflow problems.

forward <- function(A, B, pi, obs) {
  # Helper function: Forward algorithm
  
  # Inputs checking
  states <- colnames(A)
  n <- length(states)
  t <- length(obs)
  
  # Initialization
  a <- matrix(0, nrow = t, ncol = n)
  colnames(a) <- states
  
  # Forward iterations
  a[1, ] <- log(pi * B[, obs[1]])
  for (k in 2:t) {
    for (state in states) {
      logsum.a <- -Inf
      for (state.pre in states) {
        temp <- a[k-1, state.pre] + log(A[state.pre, state])
        if (temp > -Inf) {
          logsum.a <- temp + log(1 + exp(logsum.a - temp))
        }
      }
      a[k, state] <- log(B[state, obs[k]]) + logsum.a
    }
  }
  
  return (a)
}

backward <- function(A, B, pi, obs) {
  # Helper function: Backward algorithm
  
  # Inputs checking
  states <- colnames(A)
  n <- ncol(A)
  t <- length(obs)
  
  # Initialization
  b <- matrix(0, nrow = t, ncol = n)
  colnames(b) <- states
  
  # Backward iterations
  b[t, ] = rep(log(1), n)
  for (k in (t-1):1) {
    for (state in states) {
      logsum.b <- -Inf
      for (state.next in states) {
        temp <- b[k+1, state.next] + 
                log(A[state, state.next]*B[state.next, obs[k+1]])
        if (temp > -Inf) {
          logsum.b <- temp + log(1 + exp(logsum.b - temp))
        }
      }
      b[k, state] <- logsum.b
    }
  }
  
  return (b)
}

evaluation <- function(A, B, pi, obs) {
  # Evaluation task (using Forward algorithm)
  
  # Inputs checking
  states <- colnames(A)
  n <- length(states)
  t <- length(obs)
  
  # Evaluation
  a <- forward(A, B, pi, obs)
  return (sum(exp(a[t, ])))
}

viterbi <- function(A, B, pi, obs) {
  # Decoding task (using Viterbi algorithm)
  
  # Inputs checking
  states <- colnames(A)
  n <- ncol(A)
  t <- length(obs)
  
  # Initialization
  delta <- pi * B[, obs[1]]
  psi <- matrix(0, nrow = t, ncol = n)
  
  # Dynamic programming iterations
  for (k in 2:t) {
    for (i in 1:n) {
      temp <- delta * A[, i]
      delta[i] <- max(temp) * B[i, obs[k]]
      psi[k, i] <- which.max(temp)
    }
  }
  
  # Path backtracking
  path <- rep(0, t)
  path[t] <- which.max(delta)
  for (k in t:2) {
    path[k-1] <- psi[k, path[k]]
  }
  
  return (states[path])
}

baumwelch <- function(obs, states, iter.max = 1000, tol = 1e-6) {
  # Learning parameters (using Baum-Welch algorithm)
  
  # Inputs checking
  states.obs <- levels(obs)
  t <- length(obs)
  n <- length(states)
  m <- length(states.obs)
  
  # Initialization
  A <- matrix(runif(n*n), ncol = n, nrow = n)
  A <- A / apply(A, 1, sum)
  colnames(A) <- states
  rownames(A) <- states
  
  B <- matrix(runif(n*m), ncol = m, nrow = n)
  B <- B / apply(B, 1, sum)
  colnames(B) <- states.obs
  rownames(B) <- states
  
  pi <- runif(n)
  pi <- pi / sum(pi)
  
  # EM iterations
  for (iter in 1:iter.max) {
    
    # ------
    # E-step
    # ------
    
    A.new <- matrix(0, ncol = n, nrow = n)
    colnames(A.new) <- states
    rownames(A.new) <- states
    
    B.new <- matrix(0, ncol = m, nrow = n)
    colnames(B.new) <- states.obs
    rownames(B.new) <- states
    
    a <- forward(A, B, pi, obs)
    b <- backward(A, B, pi, obs)
    
    prob.obs <- a[t, 1]
    for (i in 2:n) {
      prob.obs <- a[t, i] + log(1 + exp(prob.obs - a[t, i]))
    }
    
    # ------
    # M-step
    # ------
    
    # M-step on A (transition matrix)
    for (x in states) {
      for (y in states) {
        temp <- a[1, x] + log(A[x, y]) + log(B[y, obs[1+1]]) + b[1+1, y]
        for (i in 2:(t-1)) {
          temp.2 <- a[i, x] + log(A[x, y]) + log(B[y, obs[i+1]]) + b[i+1, y]
          if (temp.2 > -Inf) {
            temp <- temp.2 + log(1 + exp(temp - temp.2))
          }
        }
        temp <- exp(temp - prob.obs)
        A.new[x, y] <- temp
      }
    }
    A.new <- A.new / apply(A.new, 1, sum)
    
    
    # M-step on B (emmision matrix)
    for (x in states) {
      for (s in states.obs) {
        temp <- -Inf
        for (i in 1:t) {
          if (s == obs[i]) {
            temp.2 <- a[i, x] + b[i, x]
            if (temp.2 > -Inf) {
              temp <- temp.2 + log(1 + exp(temp - temp.2))
            }
          }
        }
        temp <- exp(temp - prob.obs)
        B.new[x, s] <- temp
      }
    }
    B.new <- B.new / apply(B.new, 1, sum)
    
    # M-step on pi
    pi <- exp(a[1, ] + b[1, ])
    pi <- pi / sum(pi)
    
    # Convergence checking
    err <- max(norm(A - A.new), norm(B - B.new))
    A <- A.new
    B <- B.new
    if (err < tol) {
      break
    }
    
    # XXX: Print current step
    print(iter)
    print(list(A, B, pi))
  }
  
  return (list(A, B, pi))
}


# Experiments

# --- Problem 01--03 ---
# A <- matrix(c(0.6, 0.4, 0.4, 0.6), c(2, 2))
# colnames(A) <- c("F", "B")
# rownames(A) <- c("F", "B")
# B <- matrix(c(0.5, 0.8, 0.5, 0.2), c(2, 2))
# colnames(B) <- c("H", "T")
# rownames(B) <- c("F", "B")
# obs <- c("H", "T", "T", "T")
# pi <- c(0.5, 0.5)
# prob <- evaluation(A, B, pi, obs)
# path <- viterbi(A, B, pi, obs)

# --- Problem 04 ---
obs <- read.csv("data.csv")[, 2]
states <- c("A", "B")
para <- baumwelch(obs, states)
A <- para[[1]]
B <- para[[2]]
pi <- para[[3]]
path <- viterbi(A, B, pi, obs)
