#
# Wishart sampling test code
#
# Test RNG Wishart samplers in statistical_functions.R
# Rcpp RNG tests in test_samesource_cpp.R


context('Wishart RNG libraries')

# rm(list = ls())
# source('statistical_functions.R')
library(MCMCpack)


# Make data ---------------------------------------------------------------


p <- 4
Sigma <- matrix(rnorm(p^2), nrow = p, ncol = p)
Sigma <- t(Sigma) %*% Sigma + diag(p) * 0.01


# Test mean equality --------------------------------------------------------------


df.Anderson <- p + 3
df.Press <- df.Anderson + p + 1



riwish.Anderson.mean.exact <- Sigma / (df.Anderson - p - 1)
riwish.Press.mean.exact <- Sigma / (df.Press - 2*(p + 1))

test_that('Exact mean of MCMCpack::riwish = Exact mean of riwish_Press',
          expect_equal(riwish.Anderson.mean.exact, riwish.Press.mean.exact))

riwish.mean.exact <- riwish.Anderson.mean.exact

# Test one RNG ------------------------------------------------------------

seed <- round(runif(1, 1, 1000))

set.seed(seed)
IW.Anderson <- MCMCpack::riwish(df.Anderson, Sigma)
set.seed(seed)
IW.Press <- riwish_Press(df.Press, Sigma)

test_that('One run of MCMCpack::riwish = One run of riwish_Press',
   expect_equal(IW.Anderson, IW.Press))

# Test convergence to the mean --------------------------------------------

test_that('riwish_Press converges to the mean', {
   skip('Not run automatically')
})

if (FALSE) {

   set.seed(seed)
   n.samples <- 10000
   list.Anderson <- replicate(n.samples, MCMCpack::riwish(df.Anderson, Sigma), simplify = FALSE)
   # set.seed(seed)
   list.Press <- replicate(n.samples, riwish_Press(df.Press, Sigma), simplify = FALSE)

   riwish.Anderson.mean <- apply(simplify2array(list.Anderson), c(1,2), mean)
   riwish.Press.mean <- apply(simplify2array(list.Press), c(1,2), mean)
   riwish.Anderson.cummean <- apply(simplify2array(list.Anderson), c(1,2), cummean)
   riwish.Press.cummean <- apply(simplify2array(list.Press), c(1,2), cummean)

   norm(riwish.Anderson.mean - riwish.mean.exact, type = 'f')
   norm(riwish.Press.mean - riwish.mean.exact, type = 'f')

   fcn.norm <- function(x){ norm(x - riwish.mean.exact, type = 'f')}
   norm.Anderson.cummean <- apply(riwish.Anderson.cummean, 1, fcn.norm)
   norm.Press.cummean <- apply(riwish.Press.cummean, 1, fcn.norm)

   plot(seq(n.samples), norm.Anderson.cummean, type = 'l', log = 'x')
   lines(seq(n.samples), norm.Press.cummean, type = 'l', col = 'blue')
   abline(h = 0, col = 'red')
   title('Convergence in norm')
}

