# Rcpp samesource test code
#
# Rcpp tests:
#	multivariate normal (dmvnorm)

library(testthat)

context('statistical_functions.cpp: multivariate normal')

# Data ----------------------------------------------------------------------

p <- round(runif(1, 3, 9))
seed <- round(runif(1, 1, 100))



# Multivariate normal tests ---------------------------------------------------------

# A symmetric positive definite matrix
Sigma <- matrix(rnorm(p^2), ncol = p, nrow = p)
Sigma <- t(Sigma) %*% Sigma + diag(p)
Sigma.chol <- chol(Sigma)
mu <- seq(p)
x <- rev(seq(p))

dens <- mvtnorm::dmvnorm(x, mu, Sigma, log = FALSE)
ldens <- log(dens)

test_that('dmvnorm is correct: pass full', {
	dens.arma <- bayessource:::dmvnorm(t(x), mu, Sigma, logd = FALSE, is_chol = FALSE)[1]
	ldens.arma <- bayessource:::dmvnorm(t(x), mu, Sigma, logd = TRUE, is_chol = FALSE)[1]
	expect_equal(dens, dens.arma)
	expect_equal(ldens, ldens.arma)
})

test_that('dmvnorm is correct: pass Cholesky', {
	dens.arma <- bayessource:::dmvnorm(t(x), mu, Sigma.chol, logd = FALSE, is_chol = TRUE)[1]
	ldens.arma <- bayessource:::dmvnorm(t(x), mu, Sigma.chol, logd = TRUE, is_chol = TRUE)[1]
	expect_equal(dens, dens.arma)
	expect_equal(ldens, ldens.arma)
})

# Normal RNG
test_that('rmvnorm is correct: full vs Cholesky', {
	set.seed(seed)
	r.full <- bayessource:::rmvnorm(1, mu, Sigma, is_chol = FALSE)
	set.seed(seed)
	r.chol <- bayessource:::rmvnorm(1, mu, Sigma.chol, is_chol = TRUE)
	expect_equal(r.full, r.chol)
})
