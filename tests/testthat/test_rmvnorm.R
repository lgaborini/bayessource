# Rcpp samesource test code
#
# Rcpp tests:
#  multivariate normal RNG


rm(list = ls())

library(testthat)
library(bayessource)


context('samesource.cpp (Rcpp): multivariate normal RNG')


p <- 3

# Generate test data ------------------------------------------------------

seed.gen <- 1337
set.seed(seed.gen)

mu <- seq(p)

# A symmetric positive definite matrix
Sigma <- matrix(rnorm(p^2), ncol = p, nrow = p)
Sigma <- t(Sigma) %*% Sigma + diag(p)
Sigma.chol <- chol(Sigma)



# Test --------------------------------------------------------------------


# seed <- round(runif(1, 1, 100))
seed <- 1

# Test case
# set.seed(1)
# out.verified <- bayessource::rmvnorm(1, mu, Sigma, is_chol = FALSE)
# dput(out.verified)

# With seed=1 set just before calling bayessource::rmvnorm
out.verified <- structure(c(-1.50053441582829, -0.66596029349688, 4.70205092462118), .Dim = c(1L, 3L))




# Full pass
set.seed(seed)
out.sample.full <- bayessource::rmvnorm(1, mu, Sigma, is_chol = FALSE)

# Full replicates
set.seed(seed)
out.sample.full.2 <- bayessource::rmvnorm(1, mu, Sigma, is_chol = FALSE)
set.seed(seed)
out.sample.full.3 <- bayessource::rmvnorm(1, mu, Sigma, is_chol = FALSE)

# Cholesky pass
set.seed(seed)
out.sample.chol <- bayessource::rmvnorm(1, mu, Sigma.chol, is_chol = TRUE)



test_that('bayessource::rmvnorm is correct: pass full', {
          expect_equal(out.sample.full, out.verified)
          expect_equal(out.sample.full.2, out.verified)
          expect_equal(out.sample.full.3, out.verified)
          })
test_that('bayessource::rmvnorm is correct: pass Cholesky', expect_equal(out.sample.chol, out.verified))
