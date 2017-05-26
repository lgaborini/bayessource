#
# R logSumExp test code
#
# Checks that means are correcly computed in logspace

library(testthat)

context('R/Rcpp logSumExp functions')

# source('statistical_functions.R')

# Create some data without overflow -----------------------------------------------------------
n <- 10
x <- abs(rnorm(n)*5)
lx <- log(x)

# Verify matrixStats::logSumExp
sum.x <- log(sum(exp(lx)))
sum.logSumExp <- matrixStats::logSumExp(lx)
test_that('matrixStats::logSumExp is correct (no overflow)', expect_equal(sum.x, sum.logSumExp))

# The mean
# exp(logSumExp(lx) - log(length(lx))) == mean(x)

# Verify bigger case ----------------------------------------------------
n <- 1000
x <- rnorm(n)

# Densities
probs <- dnorm(x)
lprobs <- log(probs)

# Direct mean: must obtain this
mean.probs <- mean(probs)
mean.lprobs <- sum(exp(lprobs))/length(lprobs)       # may overflow

mean.lse <- exp(matrixStats::logSumExp(lprobs) - log(length(lprobs)))
mean.lsem <- exp(logSumExpMean(lprobs))

test_that('matrixStats::logSumExp is correct', {
   expect_equal(mean.lprobs, mean.lse)
   expect_equal(mean.probs, mean.lprobs)
})
test_that('R logSumExpMean is correct', expect_equal(mean.lprobs, mean.lsem))


# C++: in its own file
# mean.lse.arma <- exp(.logSumExp_arma(lprobs) - log(length(lprobs)))
# mean.lsem.arma <- exp(.logSumExpMean_arma(lprobs))

# test_that('Rcpp::logSumExp is correct', expect_equal(mean.lse.arma, mean.probs))
# test_that('Rcpp::logSumExpMean is correct', expect_equal(mean.lsem.arma, mean.probs))
