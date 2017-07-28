# Rcpp samesource test code
#
# Rcpp tests:
#  convergence of Gibbs sampler under H_p (same source)


rm(list = ls())

library(testthat)
library(bayessource)


context('samesource.cpp (Rcpp): Gibbs sampler convergence')
# skip('Skipped: not batch runnable')

p <- 2
# seed <- round(runif(1, 1, 100))
seed <- 1      # extreme sampling error!
set.seed(seed)

# Generate n random integers in {min, ..., max}
randi <- function(n, min, max) { round(runif(n, min, max)) }
# Generate random invertible symmetric matrix
randinvmat <- function(p, M = 3, alpha = 1) {
   mat <- matrix(randi(p^2, -M, M), nrow = p, ncol = p)
   t(mat) %*% mat + diag(p) * alpha
}

# Inverted Wishart RNG ----------------------------------------------------

#' Generate random sample from Inverted Wishart.
#'
#' Using Press parametrization.
#'
#' Uses \pkg{MCMCpack}::riwish.
#'
#' @param v dof (\eqn{> 2p})
#' @param S the scale matrix (pxp)
#' @return a single random variate from IW(v, S)
riwish_Press <- function(v, S){
   p <- nrow(S)
   stopifnot(v > 2*p)
   v.Anderson <- v - p - 1
   return(MCMCpack::riwish(v.Anderson, S))
}

# Bottom-up construction ---------------------------------------------------------------

# Upper hierarchical level
B.exact <- randinvmat(p)
U.exact <- randinvmat(p)
mu.exact <- randi(p, 0, 10)

# Press parametrization
# nw.exact <- 2*(p + 1)     # has no mean!
nw.exact <- 2*(p + 1) + 1

# Middle hierarchical level
theta.exact.1 <- mvtnorm::rmvnorm(1, mu.exact, B.exact)
W.exact.1 <- riwish_Press(nw.exact, U.exact)
# theta.exact.1 <- c(0,0)
# W.exact.1 <- diag(p) + matrix(c(0, 0.2, 0.2, -0.1), nrow = p, ncol = p)

# Hp: same source
theta.exact.2 <- theta.exact.1
W.exact.2 <- W.exact.1

# Lower hierarchical level: data
n <- 200
df.1 <- rmvnorm(n, theta.exact.1, W.exact.1)
df.2 <- rmvnorm(n, theta.exact.2, W.exact.2)
df <- rbind(df.1, df.2)
X <- as.matrix(df)

# Run the Gibbs sampler ---------------------------------------------------------------

burn.in <- 1000
n.iter <- 20000
n.iter.full <- burn.in + n.iter

# Oracle
W.inv <- solve(W.exact.1)
U <- U.exact
B.inv <- solve(B.exact)
mu <- mu.exact
nw <- nw.exact

# Fixed
# W.inv <- 0.01*diag(p)
# U <- 0.01*diag(p)
# B.inv <- 0.01*diag(p)
# mu <- c(0,0)
# nw <- 2*(p + 1) + 1

results <- marginalLikelihood(X, n.iter.full, B.inv, W.inv, U, nw, mu, burn.in, output.mcmc = TRUE, verbose = FALSE, Gibbs_only = TRUE)
postproc <- mcmc_postproc(results$mcmc, compute.ML = TRUE, cumulative = TRUE)

# Subsample the posterior samples (to plot)
theta.samples.mtx <- as.matrix(postproc$theta.samples)
idx.rows <- seq(1, nrow(theta.samples.mtx), length.out = 10000)

# Data sample mean: posterior mean should be close (if n is large enough)
X.mean <- colMeans(X)

# graphics.off()
plot(theta.samples.mtx[idx.rows,1], theta.samples.mtx[idx.rows,2], pch = 16, cex = 0.3, col = '#7F7F7F') #,
     # xlim = mu.exact[1] + c(-1,1) * 6, ylim = mu.exact[2] + c(-1,1) * 6)
points(theta.exact.1[1], theta.exact.1[2], col = 'red', pch = 16, cex = 2)
# points(postproc$theta.samples.ML.cum[1,1], postproc$theta.samples.ML.cum[1,2], col = 'black', pch = 16, cex = 2)
lines(postproc$theta.samples.ML.cum[,1], postproc$theta.samples.ML.cum[,2], lwd = 2)
points(postproc$theta.samples.ML[1], postproc$theta.samples.ML[2], col = 'blue', pch = 16, cex = 2)

# Check accuracy for posterior mean ---------------------------------------


# coda::effectiveSize(postproc$theta.samples)


theta.se <- coda::batchSE(postproc$theta.samples)

# Check if the real mean is contained into the confidence interval
is.lower <- (X.mean[1] < (postproc$theta.samples.ML + 2*theta.se))
is.upper <- (X.mean[2] > (postproc$theta.samples.ML - 2*theta.se))
# expect_true(all(is.lower), 'theta.ML is outside (greater than) the credibility interval.')
# expect_true(all(is.upper), 'theta.ML is outside (smaller than) the credibility interval.')


# Testing full output -----------------------------------------------------

verbose <- FALSE
verbose <- TRUE
n_cores <- 1
set.seed(1)
mlik.target <- -1457.40472756     # digits = 12
results.full <- marginalLikelihood(X, n.iter.full, B.inv, W.inv, U, nw, mu, burn.in, output.mcmc = TRUE, verbose = verbose, Gibbs_only = FALSE, n_cores = n_cores)
results.full$value
expect_equal(results.full$value, mlik.target)

mlik.target <- -1457.41174677
results.full <- marginalLikelihood(X, n.iter.full, B.inv, W.inv, U, nw, mu, burn.in, output.mcmc = TRUE, verbose = FALSE, Gibbs_only = FALSE, n_cores = n_cores)
results.full$value
expect_equal(results.full$value, mlik.target)


# Testing full output w/ OpenMP -------------------------------------------
library(microbenchmark)

verbose <- FALSE
n.times <- 10
f <- function(n.cores) { marginalLikelihood(X, n.iter.full, B.inv, W.inv, U, nw, mu, burn.in, output.mcmc = TRUE, verbose = verbose, Gibbs_only = FALSE, n_cores = n_cores) }
tmp <- microbenchmark(f(1), f(2), times = n.times)
plot(tmp)
