#
# Wishart density test code
#
# Checks that Wishart / IWishart densities are correctly computed
#
# R tests

library(testthat)

rm(list = ls())

# source('statistical_functions.R')
# source('samesource_cpp.R')

context('R Wishart/Inverted Wishart')

# Wishart test ------------------------------------------------------------


# Verify 1D case: Wishart is a Chi^2 -----------
# p = 1, Sigma = 1:
# X ~ W(n, Sigma) = Chisq(n)

X <- seq(from = 1e-5, to = 8, length.out = 100)
df <- 3
Sigma <- 1

dens.wishart <- sapply(X, function(x) dwishart(matrix(x), df, Sigma, log = FALSE, is.chol = FALSE))
dens.wishart.chol <- sapply(X, function(x) dwishart(chol(matrix(x)), df, chol(Sigma), log = FALSE, is.chol = TRUE))
dens.chisq <- dchisq(X, df)

test_that('1D Wishart is a Chi^2, pass full', expect_equal(dens.wishart, dens.chisq))
test_that('1D Wishart is a Chi^2, pass Cholesky', expect_equal(dens.wishart.chol, dens.chisq))


# Verify 2D case: ---------------------------------------------------------
# p <- 2
# Sigma <- diag(2)
# df <- 2
# X <- diag(2)
# dwishart(X, df, Sigma, log = FALSE, is.chol = FALSE)
#

# Verify pxp case ----------------------------------------------------------
# X ~ Wishart(df, Sigma)    (according to Anderson/Press parametrization: rWishart, dwishart)
#
# then X^(-1) ~ IWishart(df + p + 1, Sigma^(-1))        (Press: diwishart)
# then X^(-1) ~ IWishart(df, Sigma^(-1))				(Anderson: diwishart_Anderson)
p <- round(runif(1, 3, 6))
df <- p + 2

# Two parametrizations
df.Anderson <- df
df.Press <- df.Anderson + p + 1

Sigma <- matrix(rnorm(p^2), ncol = p, nrow = p)
Sigma <- t(Sigma) %*% Sigma + diag(p)
X <- rWishart(1, df, Sigma)[,,1]

X.inv <- solve(X)
Sigma.inv <- solve(Sigma)
Sigma.chol <- chol(Sigma)
Sigma.inv.chol <- solve(Sigma.chol)
X.chol <- chol(X)
X.inv.chol <- chol(X.inv)

# Wishart / Inverse Wishart equivalence -----------------------------------
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!! DO NOT FORGET THE JACOBIAN !!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# X ~ Wishart(df, Sigma)
# then X^(-1) ~ IWishart(df + p + 1, Sigma^(-1))
det.jacobian <- det(X.inv)^(-(p + 1))
dens.wishart <- dwishart(X, df, Sigma,                                                      log = FALSE, is.chol = FALSE)
dens.iwishart.Anderson.corrected <- diwishart_Anderson(X.inv, df, Sigma.inv,                log = FALSE, is.chol = FALSE) / det.jacobian
dens.iwishart.Press.corrected <- diwishart(X.inv, df + p + 1, Sigma.inv,                    log = FALSE, is.chol = FALSE) / det.jacobian
test_that('dwishart = diwishart (Press), pass full', expect_equal(dens.wishart, dens.iwishart.Press.corrected))
test_that('dwishart = diwishart (Anderson), pass full', expect_equal(dens.wishart, dens.iwishart.Anderson.corrected))

# Cholesky
dens.wishart.chol <- dwishart(X.chol, df, Sigma.chol,                                       log = FALSE, is.chol = TRUE)
dens.iwishart.Anderson.corrected.chol <- diwishart_Anderson(X.inv.chol, df, Sigma.inv.chol, log = FALSE, is.chol = TRUE) / det.jacobian
dens.iwishart.Press.corrected.chol <- diwishart(X.inv.chol, df + p + 1, Sigma.inv.chol,     log = FALSE, is.chol = TRUE) / det.jacobian
test_that('dwishart = diwishart (Press), pass Cholesky', expect_equal(dens.wishart.chol, dens.iwishart.Press.corrected.chol))
test_that('dwishart = diwishart (Anderson), pass Cholesky', expect_equal(dens.wishart.chol, dens.iwishart.Anderson.corrected.chol))

# Cholesky-full equivalence
test_that('diwishart (Press), pass Cholesky = diwishart (Press), pass full', expect_equal(dens.iwishart.Press.corrected.chol, dens.iwishart.Press.corrected))
test_that('diwishart (Anderson), pass Cholesky = diwishart (Anderson), pass full', expect_equal(dens.iwishart.Anderson.corrected.chol, dens.iwishart.Anderson.corrected))
test_that('dwishart, pass Cholesky = dwishart, pass full', expect_equal(dens.wishart, dens.wishart.chol))

# # Log-densities
# ldens.wishart <- log(dens.wishart)
# ldens.iwishart.Anderson.corrected <- log(dens.iwishart.Anderson) - log(det.jacobian)
# ldens.iwishart.Press.corrected <- log(dens.iwishart) - log(det.jacobian)
# test_that('ldwishart (Anderson) = ldiwishart (Anderson), pass full', expect_equal(ldens.wishart, ldens.iwishart.Anderson.corrected))
# test_that('ldwishart (Anderson) = ldiwishart (Press), pass full', expect_equal(ldens.wishart, ldens.iwishart.Press.corrected))


# # The opposite relation
# then X^(-1) ~ IWishart(df + p + 1, Sigma^(-1))        (Press: diwishart)
# then X^(-1) ~ IWishart(df, Sigma^(-1))				(Anderson: diwishart_Anderson)
dens.iwishart.Press <- diwishart(X, df.Press, Sigma,                                        log = FALSE, is.chol = FALSE)
dens.wishart.corrected <- dwishart(X.inv, df, Sigma.inv,                                    log = FALSE, is.chol = FALSE) / det.jacobian
dens.iwishart.Press.chol <- diwishart(X.chol, df.Press, Sigma.chol,                         log = FALSE, is.chol = TRUE)
dens.wishart.corrected.chol <- dwishart(X.inv.chol, df, Sigma.inv.chol,                     log = FALSE, is.chol = TRUE) / det.jacobian
test_that('diwishart (Press) = dwishart, pass full', expect_equal(dens.iwishart.Press, dens.wishart.corrected))
test_that('diwishart (Press) = dwishart, pass Cholesky', expect_equal(dens.iwishart.Press.chol, dens.wishart.corrected.chol))


# Given the inverse
dens.iwishart.Press <- diwishart(X, df.Press, Sigma,                          log = FALSE, is.chol = FALSE)
dens.iwishart.Press.inverse <- diwishart_inverse(X.inv, df.Press, Sigma,      log = FALSE, is.chol = FALSE)
test_that('diwishart (Press) = diwishart_inverse (Press), pass full', expect_equal(dens.iwishart.Press, dens.iwishart.Press.inverse))
dens.iwishart.Press.chol <- diwishart(X.chol, df.Press, Sigma.chol,                     log = FALSE, is.chol = TRUE)
dens.iwishart.Press.inverse.chol <- diwishart_inverse(X.inv.chol, df.Press, Sigma.chol, log = FALSE, is.chol = TRUE)
test_that('diwishart (Press) = diwishart_inverse (Press), pass Cholesky', expect_equal(dens.iwishart.Press.chol, dens.iwishart.Press.inverse.chol))

# test_that('diwishart_inverse (Press), pass full = diwishart_inverse (Press), pass Cholesky', expect_equal(dens.iwishart.Press, dens.iwishart.Press.chol))
test_that('diwishart_inverse (Press), pass full = diwishart_inverse (Press), pass Cholesky', expect_equal(dens.iwishart.Press.inverse, dens.iwishart.Press.inverse.chol))

# Equivalence between Inverted Wishart parametrizations ---------------------------
# Inverted Wishart
# df.Anderson = df.Press - p - 1 = df
# df.Press = df.Anderson + p + 1


# Avoid comparison of logarithms:
# for large p, dens is very small, tolerance needs to be adjusted according to p
# not easy

dens.iwishart.Anderson <- diwishart_Anderson(X.inv, df.Anderson, Sigma,  log = FALSE, is.chol = FALSE)
dens.iwishart.Press <- diwishart_Anderson(X.inv, df.Press, Sigma,        log = FALSE, is.chol = FALSE)
ldens.iwishart.Anderson <- diwishart_Anderson(X.inv, df.Anderson, Sigma, log = TRUE, is.chol = FALSE)
ldens.iwishart.Press <- diwishart_Anderson(X.inv, df.Press, Sigma,       log = TRUE, is.chol = FALSE)
test_that('diwishart (Press) = diwishart (Anderson), pass full', {
	expect_equal(dens.iwishart.Anderson, dens.iwishart.Press)
	# expect_equal(ldens.iwishart.Anderson, ldens.iwishart.Press, tolerance = -log(1.5^(-12)), scale = 1)
})

dens.iwishart.Anderson.chol <- diwishart_Anderson(X.inv.chol, df.Anderson, Sigma.chol,  log = FALSE, is.chol = TRUE)
dens.iwishart.Press.chol <- diwishart_Anderson(X.inv.chol, df.Press, Sigma.chol,        log = FALSE, is.chol = TRUE)
ldens.iwishart.Anderson.chol <- diwishart_Anderson(X.inv.chol, df.Anderson, Sigma.chol, log = TRUE, is.chol = TRUE)
ldens.iwishart.Press.chol <- diwishart_Anderson(X.inv.chol, df.Press, Sigma.chol,       log = TRUE, is.chol = TRUE)
test_that('diwishart (Press) = diwishart (Anderson), pass Cholesky', {
	expect_equal(dens.iwishart.Anderson.chol, dens.iwishart.Press.chol)
	# expect_equal(ldens.iwishart.Anderson.chol, ldens.iwishart.Press.chol, tolerance = -log(1.5^(-12)), scale = 1)
})

# # Rcpp: in its own file ----------------------------------------------------------

# dens.iwishart <- diwishart_inverse(solve(X), df + p + 1, Sigma, log = FALSE)
# dens.iwishart
# dens.iwishart.arma <- diwishart_inverse_arma(solve(X), df + p + 1, Sigma, logd = FALSE)
# dens.iwishart.arma
# test_that('diwishart (Press) = diwishart (Rcpp)', expect_equal(dens.iwishart, dens.iwishart.arma))
