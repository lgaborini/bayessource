# Rcpp samesource test code
#
# Rcpp tests:
#  mathematical functions (inverses, etc.)
#  logsumexp

library(testthat)
# library(mvtnorm)

context('statistical_functions.cpp: mathematical functions, logSumExp')

# Refresh cpp file without running R code chunks (avoid build-test loop)
# Rcpp::sourceCpp('samesource.cpp', embeddedR = FALSE)

# Data ----------------------------------------------------------------------

p <- round(runif(1, 3, 9))
seed <- round(runif(1, 1, 100))


# Inverses and Cholesky ---------------------------------------------------
# p <- 8

# X is sympd
X <- matrix(rnorm(p^2), ncol = p, nrow = p)
X <- t(X) %*% X + diag(p)

X.inv <- solve(X)
X.chol <- chol(X)
X.inv.chol <- chol(X.inv)

# X.chol is triangular

test_that('Rcpp::inv_triangular is correct.',
	expect_equal(solve(X.chol), bayessource:::inv_triangular(X.chol)))

test_that('Rcpp::chol2inv is correct.', {
	expect_equal(base::chol2inv(X.chol), bayessource:::chol2inv(X.chol))
	expect_equal(bayessource:::chol2inv(X.chol), X.inv)
})

test_that('inv_pd is correct.', {
   expect_equal(base::solve(X), inv_pd(X))
})

X.inv.chol.arma <- bayessource:::inv_Cholesky_from_Cholesky(X.chol)
test_that('Rcpp::inv_Cholesky_from_Cholesky is correct.',
	expect_equal(X.inv.chol, X.inv.chol.arma)
)

ldetX <- log(det(X))
ldetX.chol <- bayessource:::ldet_from_Cholesky(X.chol)

test_that('Rcpp::ldet_from_Cholesky is correct.',
	expect_equal(ldetX, ldetX.chol)
)

test_that('Rcpp::inv_sympd_tol is correct.', {
	expect_equal(X.inv, bayessource:::inv_sympd_tol(X))
})




# Rcpp: logsumexp ---------------------------------------------------------

# Create some data
n <- 1000
x <- rnorm(n)
# densities
probs <- dnorm(x)
lprobs <- log(probs)

# Direct mean: must obtain this
mean.probs <- mean(probs)

# C++
mean.lse.arma <- exp(bayessource:::logSumExp(lprobs) - log(length(lprobs)))
mean.lsem.arma <- exp(bayessource:::logSumExpMean(lprobs))

test_that('Mean with Rcpp::logSumExp is correct', expect_equal(mean.lse.arma, mean.probs))
test_that('Rcpp::logSumExpMean is correct', expect_equal(mean.lsem.arma, mean.probs))
