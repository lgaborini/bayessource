# Rcpp samesource test code
#
# Rcpp tests:
#  convergence of Gibbs sampler under H_p (same source)


rm(list = ls())

library(testthat)
library(bayessource)


context('samesource.cpp (Rcpp): Gibbs sampler convergence')

# Data ----------------------------------------------------------------------

p <- 2
seed <- round(runif(1, 1, 100))

# Generate n random integers in {min, ..., max}
randi <- function(n, min, max) { round(runif(n, min, max)) }
# Generate random invertible symmetric matrix
randinvmat <- function(p, M = 3, alpha = 1) {
   mat <- matrix(randi(p^2, -M, M), nrow = p, ncol = p)
   t(mat) %*% mat + diag(p) * alpha
}

theta.1 <- c(0,0)
theta.2 <- c(1,1)
Sigma.1 <- randinvmat(p)
Sigma.2 <- randinvmat(p)


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
theta.exact.1 <- c(0,0)
W.exact.1 <- diag(p) + matrix(c(0, 0.2, 0.2, -0.1), nrow = p, ncol = p)

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

n.iter <- 100000
burn.in <- 1000

# Oracle
W.inv <- solve(W.exact.1)
U <- U.exact
B.inv <- solve(B.exact)
mu <- mu.exact
nw <- nw.exact

# Fixed
W.inv <- diag(p)
U <- diag(p)
B.inv <- diag(p)
mu <- c(0,0)
nw <- 3

results <- marginalLikelihood(X, n.iter, B.inv, W.inv, U, nw, mu, burn.in, output.mcmc = TRUE, verbose = FALSE)
postproc <- mcmc_postproc(results$mcmc, compute.ML = TRUE, cumulative = TRUE)


theta.samples.mtx <- as.matrix(postproc$theta.samples)
plot(theta.samples.mtx[,1], theta.samples.mtx[,2])
points(theta.exact.1[1], theta.exact.1[2], col = 'red')
points(postproc$theta.samples.ML[1], postproc$theta.samples.ML[2], col = 'blue')



# Check accuracy for posterior mean ---------------------------------------


# effectiveSize(postproc$theta.samples)

theta.se <- effectiveSize(postproc$theta.samples)

# Check if the real mean is contained into the confidence interval
is.lower <- (theta.exact.1 < (postproc$theta.samples.ML + 2*theta.se))
is.upper <- (theta.exact.1 > (postproc$theta.samples.ML - 2*theta.se))
expect_true(all(is.lower), 'theta.ML is outside (greater than) the credibility interval.')
expect_true(all(is.upper), 'theta.ML is outside (smaller than) the credibility interval.')
