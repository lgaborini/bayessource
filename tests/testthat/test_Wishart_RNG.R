# Rcpp samesource test code
#
# Rcpp tests:
#    Wishart RNG, consistency with standard library

library(testthat)

context('statistical_functions.cpp: Wishart/Inverted Wishart RNG')

# Refresh cpp file without running R code chunks (avoid build-test loop)
# Rcpp::sourceCpp('samesource.cpp', embeddedR = FALSE)

# Data ----------------------------------------------------------------------

p <- round(runif(1, 3, 9))
seed <- round(runif(1, 1, 100))

# A symmetric positive definite matrix
Sigma <- matrix(rnorm(p^2), ncol = p, nrow = p)
Sigma <- t(Sigma) %*% Sigma + diag(p)
Sigma.chol <- chol(Sigma)
Sigma.inv <- solve(Sigma)
Sigma.inv.chol <- solve(Sigma.chol)

# Wishart RNG tests ---------------------------------------------------------
# X ~ Wishart(df, Sigma)    (according to Anderson/Press parametrization: rWishart, dwishart)

df <- p + round(runif(1, 2, 10))			# Wishart: Press/Anderson, Inverted Wishart: Anderson
df.Anderson <- df
df.Press <- df + p + 1			# Inverted Wishart

set.seed(seed)
X <- bayessource:::rwish(df, Sigma, FALSE, FALSE)
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




# Inverted Wishart RNG: consistency tests ------------------------------------------
# Press parametrization
#
# Cannot fix the RNG seed, as RNGs are called in different orders
#
# Generate from W ~ IW(v, Sigma)
#   v > 2*(p + 1)
#   (v_Anderson = v - p - 1)
#
# Then:
#
#   E[W] = Sigma / (v - 2(p + 1))
#
# Variances from:
#   https://en.wikipedia.org/wiki/Inverse-Wishart_distribution#Moments

# Exact mean
iwishart.mean.exact <- Sigma / (df.Press - 2*(p + 1))
# Exact variances (only the diagonal)
iwishart.variances.exact <- 2*diag(Sigma)^2 / ((df.Anderson - p - 1)^2 * (df.Anderson - p - 3))

n.samples <- 10000

riwish.samples <- replicate(n.samples, riwish_Press(df.Press, Sigma), simplify = FALSE)
riwish.mean <- apply(simplify2array(riwish.samples), c(1,2), mean)

# MC error: sigma/N
sd.mc.error <- sqrt(iwishart.variances.exact / n.samples)

n.comp <- p

alpha.0 <- 0.01
alpha <- alpha.0 / n.comp       # Bonferroni correction
# alpha <- alpha.0
z <- qnorm(alpha/2)

# Confidence interval
# (mean_sample - mean_true)/sd.mc.error ~ N(0, 1)
#
# H0: MC converges to the mean
# H1: MC converges to some other value
is.upper <- (riwish.mean + z * sd.mc.error) <  iwishart.mean.exact
is.lower <- (riwish.mean - z * sd.mc.error) >  iwishart.mean.exact
test.result <- all(is.upper & is.lower)

test_that(paste('riwish converges to the mean, confidence', 1 - alpha),
          expect_true(test.result))


