# Rcpp samesource test code
#
# Rcpp tests:
#    Wishart/Inverted Wishart equivalence

library(testthat)

context('statistical_functions.cpp: Wishart/Inverted Wishart densities')


# Verify 1D case: Wishart is a Chi^2 ----------------------------------------------------------------------
# p = 1, Sigma = 1:
# X ~ W(n, Sigma) = Chisq(n)

X <- seq(from = 1e-5, to = 8, length.out = 100)
df <- 3
Sigma <- 1

dens.wishart <- sapply(X, function(x) bayessource::dwishart(matrix(x), df, Sigma, log = FALSE, is.chol = FALSE))
dens.wishart.chol <- sapply(X, function(x) bayessource::dwishart(chol(matrix(x)), df, chol(Sigma), log = FALSE, is.chol = TRUE))
dens.chisq <- dchisq(X, df)

test_that('1D Wishart is a Chi^2, pass full', expect_equal(dens.wishart, dens.chisq))
test_that('1D Wishart is a Chi^2, pass Cholesky', expect_equal(dens.wishart.chol, dens.chisq))

# Data ----------------------------------------------------------------------

p <- round(runif(1, 3, 9))
seed <- round(runif(1, 1, 100))

# A symmetric positive definite matrix
Sigma <- matrix(rnorm(p^2), ncol = p, nrow = p)
Sigma <- t(Sigma) %*% Sigma + diag(p)
Sigma.chol <- chol(Sigma)
Sigma.inv <- solve(Sigma)
Sigma.inv.chol <- solve(Sigma.chol)

# Generic p-dimensional case: ----------------------------------------------------------------------
# X ~ Wishart(df, Sigma)    (according to Anderson/Press parametrization: rWishart, dwishart)


# Generate some data
df <- p + round(runif(1, 1, 10))       # Wishart: Press/Anderson, Inverted Wishart: Anderson
df.Anderson <- df
df.Press <- df + p + 1        # Inverted Wishart

set.seed(seed)
X <- bayessource:::rwish(df, Sigma, FALSE, FALSE)
X.inv <- solve(X)
X.chol <- chol(X)
X.inv.chol <- chol(X.inv)




# Inverse Wishart densities ---------------------------------------------------------

dens.iwishart <- bayessource:::diwishart_inverse_R(X.inv, df.Press, Sigma,            log  = FALSE, is.chol = FALSE)
dens.iwishart.arma <- bayessource:::diwishart_inverse(X.inv, df.Press, Sigma,  logd = FALSE, is_chol = FALSE)
ldens.iwishart.arma <- bayessource:::diwishart_inverse(X.inv, df.Press, Sigma, logd = TRUE,  is_chol = FALSE)

test_that('R diwishart_inverse (Press) = Rcpp diwishart_inverse (Press), pass full', {
   expect_equal(dens.iwishart, dens.iwishart.arma)
   # expect_equal(log(dens.iwishart), ldens.iwishart.arma)
})

# Cholesky passes
dens.iwishart.chol <- bayessource:::diwishart_inverse_R(X.inv.chol, df.Press, Sigma.chol,            log  = FALSE, is.chol = TRUE)
dens.iwishart.arma.chol <- bayessource:::diwishart_inverse(X.inv.chol, df.Press, Sigma.chol,  logd = FALSE, is_chol = TRUE)
ldens.iwishart.arma.chol <- bayessource:::diwishart_inverse(X.inv.chol, df.Press, Sigma.chol, logd = TRUE,  is_chol = TRUE)

test_that('R diwishart_inverse (Press) = Rcpp diwishart_inverse (Press), pass Cholesky', {
	expect_equal(dens.iwishart.chol, dens.iwishart.arma.chol)
	# expect_equal(log(dens.iwishart.chol), ldens.iwishart.arma.chol, tolerance = -log(1.5^(-8)), scale = 1)     # larger tolerance
})

# Full / Cholesky passes
test_that('Rcpp diwishart_inverse (Press), pass full = Rcpp diwishart_inverse (Press), pass Cholesky', {
   expect_equal(dens.iwishart.arma, dens.iwishart.arma.chol)
})

# R / R through inverse / Rcpp through inverse equivalence
# Probably duplicated...
dens.iwishart <- bayessource:::diwishart(X, df.Press, Sigma,                               log  = FALSE, is.chol = FALSE)
dens.iwishart.inverse <- bayessource:::diwishart_inverse_R(X.inv, df.Press, Sigma,           log  = FALSE, is.chol = FALSE)
dens.iwishart.inverse.arma <- bayessource:::diwishart_inverse(X.inv, df.Press, Sigma, logd = FALSE, is_chol = FALSE)

test_that('R diwishart (Press) = R diwishart_inverse (Press) = Rcpp diwishart_inverse (Press), pass full', {
	expect_equal(dens.iwishart, dens.iwishart.inverse)
	expect_equal(dens.iwishart, dens.iwishart.inverse.arma)
	expect_equal(dens.iwishart.inverse, dens.iwishart.inverse.arma)
})

dens.iwishart.chol <- bayessource:::diwishart(X.chol, df.Press, Sigma.chol,                               log  = FALSE, is.chol = TRUE)
dens.iwishart.inverse.chol <- bayessource:::diwishart_inverse_R(X.inv.chol, df.Press, Sigma.chol,           log  = FALSE, is.chol = TRUE)
dens.iwishart.inverse.arma.chol <- bayessource:::diwishart_inverse(X.inv.chol, df.Press, Sigma.chol, logd = FALSE, is_chol = TRUE)

test_that('R diwishart (Press) = R diwishart_inverse (Press) = Rcpp diwishart_inverse (Press), pass Choleksy', {
	expect_equal(dens.iwishart.chol, dens.iwishart.inverse.chol)
	expect_equal(dens.iwishart.chol, dens.iwishart.inverse.arma.chol)
	expect_equal(dens.iwishart.inverse.chol, dens.iwishart.inverse.arma.chol)
})

test_that('R diwishart (Press) = R diwishart_inverse (Press) = Rcpp diwishart_inverse (Press), pass Choleksy = pass full', {
	expect_equal(dens.iwishart, dens.iwishart.chol)
	expect_equal(dens.iwishart.inverse, dens.iwishart.inverse.arma.chol)
	expect_equal(dens.iwishart.inverse.arma, dens.iwishart.inverse.chol)
})


# Inverse Wishart density, special case ---------------------------------------------------------
# From https://github.com/stan-dev/math/blob/develop/test/unit/math/mix/mat/prob/inv_wishart_test.cpp
# Also check out https://github.com/stan-dev/math/blob/develop/test/unit/math/mix/mat/prob/wishart_test.cpp

# TODO: to add in statistical_functions.R

# p <- 3
# X <- matrix(c(12.147233, -11.9036079, 1.0910458, -11.903608, 16.7585782, 0.8530256, 1.091046, 0.8530256, 2.5786609), nrow = p, ncol = p, byrow = TRUE)
# Sigma <- matrix(c(7.785215, 3.0597878, 1.1071663, 3.059788, 10.3515035, -0.1232598, 1.107166, -0.1232598, 7.7623386), nrow = p, ncol = p, byrow = TRUE)
# df <- 4
# log_p = log(2.008407e-08)


# Wishart - Inverse Wishart equivalence ----------------------------------------------------------
#
# X ~ Wishart(df, Sigma)    (according to Anderson/Press parametrization: rWishart, dwishart)
# then X^(-1) ~ IWishart(df + p + 1, Sigma^(-1))        (Press: diwishart)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! DO NOT FORGET THE JACOBIAN !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# X ~ Wishart(df, Sigma)
# then X^(-1) ~ IWishart(df + p + 1, Sigma^(-1))
det.jacobian <- det(X.inv)^(-(p + 1))
dens.wishart <- bayessource:::dwishart(X, df, Sigma,                                                      log = FALSE, is.chol = FALSE)
dens.iwishart.Press.corrected <- bayessource:::diwishart(X.inv, df + p + 1, Sigma.inv,                    log = FALSE, is.chol = FALSE) / det.jacobian
test_that('dwishart = diwishart (Press), pass full', expect_equal(dens.wishart, dens.iwishart.Press.corrected))

# Cholesky
dens.wishart.chol <- bayessource:::dwishart(X.chol, df, Sigma.chol,                                       log = FALSE, is.chol = TRUE)
dens.iwishart.Press.corrected.chol <- bayessource:::diwishart(X.inv.chol, df + p + 1, Sigma.inv.chol,     log = FALSE, is.chol = TRUE) / det.jacobian
test_that('dwishart = diwishart (Press), pass Cholesky', expect_equal(dens.wishart.chol, dens.iwishart.Press.corrected.chol))

# Cholesky-full equivalence
test_that('diwishart (Press), pass Cholesky = diwishart (Press), pass full', expect_equal(dens.iwishart.Press.corrected.chol, dens.iwishart.Press.corrected))
test_that('dwishart, pass Cholesky = dwishart, pass full', expect_equal(dens.wishart, dens.wishart.chol))

