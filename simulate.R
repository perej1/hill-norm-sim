#' Hill estimator for readily ordered data
#'
#' @param data_inc Double data vector in increasing order.
#' @param n Integer, sample size.
#' @param k Integer, threshold for tail observations.
#'
#' @returns Double.
estimate_gamma <- function(data_inc, n, k) {
  mean((log(data_inc[(n - k):n]) - log(data_inc[n - k]))[-1])
}


#' Simulate dense discretization from fractional Brownian motion path
#'
#' @param h Hurst index
#'
#' @returns Tibble with two columns, time points and values from the path
sim_fbm_path <- function(h) {
  if (h <= 0 | h >= 1) rlang::abort("Invalid Hurst index.")
  m <- 2^20 # Discretization level, number of increments
  
  # Autocovariance of unit-step fractional Gaussian noise, lags 0..n
  k <- 0:m
  gam <- 0.5 * (abs(k + 1)^(2 * h) - 2 * abs(k)^(2 * h) + abs(k - 1)^(2 * h))
  
  # Circulant embedding: eigenvalues via FFT (non-negative for fGn)
  lam <- Re(fft(c(gam, gam[m:2])))
  lam[lam < 0 & lam > -1e-10] <- 0
  if (any(lam < 0)) rlang::abort("Circulant embedding failed (negative eigenvalue).")
  
  z <- complex(real = rnorm(2 * m), imaginary = rnorm(2 * m))
  fgn <- Re(fft(sqrt(lam / (2 * m)) * z))[1:m] # unit-step fGn
  
  tibble::tibble(
    time = (0:m) / m,
    value = c(0L, cumsum(fgn * (1 / m)^h)) # rescale to step 1/m
  )
}


#' Simulate dense discretization from Poisson process
#'
#' @param lambda Rate
#'
#' @returns Tibble with two columns, time points and values from the path
sim_poisson_path <- function(lambda) {
  if (lambda <= 0) rlang::abort("Invalid rate.")
  m <- 2^20 # Discretization level, number of increments
  n <- m - 1
  incr <- rpois(m, lambda / m)
  tibble::tibble(
    time = (0:m) / m,
    value = c(0L, cumsum(incr))
  )
}


#' Sample points from a realization of a process
#'
#' @param path Tibble representing the realized path
#' @param n Number of sampled points
#'
#' @returns Tibble with four columns, first two columns represent the
#'   deterministic sampling and two last columns represent the random sampling
sample_path <- function(path, n) {
  n <- n - 1
  m <- nrow(path) - 1 # Discretization level, number of increments
  if(m %% n != 0) {
    rlang::abort("Number of sampled points is not multiple of the total discretization level")
  }
  
  # Create time grids
  grid <- (0:n) * m / n + 1
  rand <- sort(sample.int(m + 1, n + 1, replace = TRUE))
  
  tibble::tibble(
    t_grid = path$time[grid],
    value_grid = path$value[grid],
    t_rand = path$time[rand],
    value_rand = path$value[rand]
  )
}


n <- 2^2 + 1

# Test fbm
fbm_path <- sim_fbm_path(0.5)
fbm_sample <- sample_path(fbm_path, n)
plot(fbm_path$time, fbm_path$value, type = "l")
points(fbm_sample$t_rand, fbm_sample$value_rand, cex = 2,  pch = 16, col  = "blue")
points(fbm_sample$t_grid, fbm_sample$value_grid, cex = 2,  pch = 16, col = "red")

# Test Poisson
poisson_path <- sim_poisson_path(5)
poisson_sample <- sample_path(poisson_path, n)
plot(poisson_path$time, poisson_path$value, type = "l")
points(poisson_sample$t_rand, poisson_sample$value_rand, cex = 2,  pch = 16, col  = "blue")
points(poisson_sample$t_grid, poisson_sample$value_grid, cex = 2,  pch = 16, col = "red")
