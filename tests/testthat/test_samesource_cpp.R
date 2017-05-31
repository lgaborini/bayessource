# Rcpp samesource test code
#
# Rcpp tests:
#    Wishart, Normal, etc.
#    statistical functions
#    mathematical functions

library(testthat)
# library(mvtnorm)

context('samesource.cpp (Rcpp)')
# source('diwishart_inverse.R')

# Refresh cpp file without running R code chunks (avoid build-test loop)
# Rcpp::sourceCpp('samesource.cpp', embeddedR = FALSE)

# Data ----------------------------------------------------------------------

p <- round(runif(1, 3, 9))
seed <- round(runif(1, 1, 100))

# A symmetric positive definite matrix
Sigma <- matrix(rnorm(p^2), ncol = p, nrow = p)
Sigma <- t(Sigma) %*% Sigma + diag(p)
Sigma.chol <- chol(Sigma)

# Wishart RNG tests ---------------------------------------------------------

df <- p + round(runif(1, 1, 10))			# Wishart: Press/Anderson, Inverted Wishart: Anderson
df.Anderson <- df
df.Press <- df + p + 1			# Inverted Wishart

set.seed(seed)
X <- rwish(df, Sigma, FALSE, FALSE)
X.inv <- solve(X)
X.chol <- chol(X)
X.inv.chol <- chol(X.inv)


test_that('rwish is correct: pass Cholesky, return full', {
	set.seed(seed)
	X.test <- bayessource:::rwish(df, Sigma.chol, TRUE, FALSE)
	expect_equal(X, X.test)
	expect_equal(chol(X), chol(X.test))
})

test_that('rwish is correct: pass Cholesky, return Cholesky', {
	set.seed(seed)
	X.chol.test <- bayessource:::rwish(df, Sigma.chol, TRUE, TRUE)
	expect_equal(X, t(X.chol.test) %*% X.chol.test)
	expect_equal(X.chol, X.chol.test)
})

test_that('rwish is correct: pass full, return Cholesky', {
	set.seed(seed)
	X.chol.test <- bayessource:::rwish(df, Sigma, FALSE, TRUE)
	expect_equal(X, t(X.chol.test) %*% X.chol.test)
	expect_equal(X.chol, X.chol.test)
})


# Wishart RNG: consistency tests ------------------------------------------
# rwish should generate from the same distribution as rWishart does
# Press/Anderson parametrization
#
# Cannot fix the RNG seed, as RNGs are called in different orders

# Exact mean
wishart.mean.exact <- df.Anderson * Sigma
# Exact variances
wishart.variance.exact <- df.Anderson * (Sigma^2 + tcrossprod(diag(Sigma)))

n.samples <- 10000
# set.seed(seed)
rWishart.samples <- rWishart(n.samples, df.Anderson, Sigma)
list.rwish.samples <- replicate(n.samples, rwish(df.Anderson, Sigma), simplify = FALSE)

rWishart.mean <- apply(rWishart.samples, c(1,2), mean)
rwish.mean <- apply(simplify2array(list.rwish.samples), c(1,2), mean)

# wishart.mean.exact
# rWishart.mean

# stopifnot(all.equal(rWishart.mean, wishart.mean.exact, tolerance = 1e-2))
# stopifnot(all.equal(rwish.mean, wishart.mean.exact, tolerance = 1e-2))

# MC error: sigma/N
sd.mc.error <- sqrt(wishart.variance.exact / n.samples)

n.comp <- p^2
alpha.0 <- 0.01
alpha <- alpha.0 / n.comp       # Bonferroni correction
# alpha <- alpha.0
z <- qnorm(alpha/2)

# Confidence interval
# (mean_sample - mean_true)/sd.mc.error ~ N(0, 1)
#
# H0: MC converges to the mean
# H1: MC converges to some other value
is.upper <- (rwish.mean + z * sd.mc.error) <  wishart.mean.exact
is.lower <- (rwish.mean - z * sd.mc.error) >  wishart.mean.exact
test.result <- all(is.upper & is.lower)

test_that(paste('rwish converges to the mean, confidence', 1 - alpha),
          expect_true(test.result))



# Inverse Wishart densities ---------------------------------------------------------

# from "statistical_functions.R"
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
dens.iwishart <- diwishart(X, df.Press, Sigma,                               log  = FALSE, is.chol = FALSE)
dens.iwishart.inverse <- bayessource:::diwishart_inverse_R(X.inv, df.Press, Sigma,           log  = FALSE, is.chol = FALSE)
dens.iwishart.inverse.arma <- bayessource:::diwishart_inverse(X.inv, df.Press, Sigma, logd = FALSE, is_chol = FALSE)

test_that('R diwishart (Press) = R diwishart_inverse (Press) = Rcpp diwishart_inverse (Press), pass full', {
	expect_equal(dens.iwishart, dens.iwishart.inverse)
	expect_equal(dens.iwishart, dens.iwishart.inverse.arma)
	expect_equal(dens.iwishart.inverse, dens.iwishart.inverse.arma)
})

dens.iwishart.chol <- diwishart(X.chol, df.Press, Sigma.chol,                               log  = FALSE, is.chol = TRUE)
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




# Multivariate normal tests ---------------------------------------------------------

# p <- 8
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


# Inverses and Cholesky ---------------------------------------------------
# p <- 8
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

X.inv.chol.arma <- bayessource:::inv_Cholesky_from_Cholesky(X.chol)
test_that('Rcpp::inv_Cholesky_from_Cholesky is correct.',
	expect_equal(X.inv.chol, X.inv.chol.arma)
)

ldetX <- log(det(X))
ldetX.chol <- bayessource:::ldet_from_Cholesky(X.chol)

test_that('Rcpp::ldet_from_Cholesky is correct.',
	expect_equal(ldetX, ldetX.chol)
)



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
