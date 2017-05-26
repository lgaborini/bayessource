#
# Statistical distributions R functions
#
# Wishart density
# multivariate normal density (faster)
# utilities


library(matrixStats)
loadNamespace('MCMCpack')

# Accurately compute log(mean(x)), where x = log(sum(exp(y)))
# https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/
logSumExpMean <- function(x) {
   ret <- matrixStats::logSumExp(x) - log(length(x))
   return(ret)
}

logMultivariateGamma <- function(theta, p) {
   p * (p - 1)/4 * log(pi) + sum(log(gamma(theta + (1 - (1:p))/2)))
}

# Normalization constant, Press parametrization
iwishart.logconstant <- function(df, p) {
   -(df - p - 1) * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((df - p - (1:p))/2)))
}


# Press parametrization
diwishart <- function(X, df, Sigma, log = FALSE, is.chol = FALSE) {
   p <- ncol(X)
   
   # The constant
   lc0 <- -(df - p - 1) * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((df - p - (1:p))/2)))
   
   if (!is.chol) {
      X.chol <- chol(X)
      Sigma.chol <- chol(Sigma)
   } else {
      X.chol <- X
      Sigma.chol <- Sigma
   }
   # ldetSigma <- log(det(Sigma))
   # ldetX <- log(det(X))
   ldetSigma <- 2 * sum(log(diag(Sigma.chol)))
   ldetX <- 2 * sum(log(diag(X.chol)))
   
   # Inverse of Cholesky factor 
   # NOT inverse from Cholesky: chol2inv(X.chol) = solve(X)!
   # NOT Cholesky of inverse
   X.chol.inv <- backsolve(r = X.chol, x = diag(p))
   
   # The trace: all equivalent expressions
   # ltrace <- sum(diag(solve(X) %*% Sigma))
   # ltrace <- sum(solve(X) * Sigma)
   # ltrace <- sum(diag(chol2inv(X.chol) %*% crossprod(Sigma.chol)))
   # ltrace <- sum((Sigma.chol %*% solve(X.chol))^2)
   ltrace <- sum((Sigma.chol %*% X.chol.inv)^2)
   
   ldens <- lc0 + (df - p - 1)/2 * ldetSigma - (df / 2) * ldetX - ltrace / 2
   
   if (log == FALSE) { return(exp(ldens)) }
   return(ldens)
}



# Inverse Wishart density(df, Sigma), by supplying (X^(-1), df, Sigma) rather than (X, df, Sigma)
# Press parametrization
diwishart_inverse <- function(X.inv, df, Sigma, log = FALSE, is.chol = FALSE) {
   p <- ncol(X.inv)
   
   # The constant
   lc0 <- -(df - p - 1) * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((df - p - (1:p))/2)))
   
   if (!is.chol) {
      X.inv.chol <- chol(X.inv)
      Sigma.chol <- chol(Sigma)
   } else {
      X.inv.chol <- X.inv
      Sigma.chol <- Sigma
   }
   
   # Cholesky factor of X: inverse-transpose of Cholesky factor of X.inv
   # X.chol <- t(backsolve(r = X.inv.chol, x = diag(p)))
   
   # ldetSigma <- log(det(Sigma))
   # ldetX <- log(det(X))
   ldetSigma <- 2 * sum(log(diag(Sigma.chol)))
   ldetX.inv <- 2 * sum(log(diag(X.inv.chol)))
   ldetX <- -ldetX.inv
   
   # The trace:
   # ltrace <- sum(diag(X.inv %*% Sigma))
   # ltrace <- sum(X.inv * Sigma)
   ltrace <- sum((Sigma.chol %*% t(X.inv.chol))^2)          # how to avoid transpose?
   
   ldens <- lc0 + (df - p - 1)/2 * ldetSigma - (df / 2) * ldetX - ltrace / 2
   
   if (log == FALSE) { return(exp(ldens)) }
   return(ldens)
}


# Anderson parametrization
diwishart_Anderson <- function(X, df, Sigma, log = FALSE, is.chol = FALSE) {
   p <- ncol(X)
   
   # The constant
   lc0 <- -df * p/2 * log(2) - logMultivariateGamma(df/2, p)
   
   if (!is.chol) {
      X.chol <- chol(X)
      Sigma.chol <- chol(Sigma)
   } else {
      X.chol <- X
      Sigma.chol <- Sigma
   }
   # ldetSigma <- log(det(Sigma))
   # ldetX <- log(det(X))
   ldetSigma <- 2 * sum(log(diag(Sigma.chol)))
   ldetX <- 2 * sum(log(diag(X.chol)))
   
   X.chol.inv <- backsolve(r = X.chol, x = diag(p))      # inverse of Cholesky factor (NOT chol2inv(X.chol) = solve(X)!)
   
   # The log(trace):
   # ltrace <- sum(diag(Sigma %*% solve(X)))
   # ltrace <- sum(solve(X) * Sigma)
   # ltrace <- sum(diag(chol2inv(X.chol) %*% crossprod(Sigma.chol)))
   # ltrace <- sum((Sigma.chol %*% solve(X.chol))^2)
   ltrace <- sum((Sigma.chol %*% X.chol.inv)^2)
   
   ldens <- lc0 + df/2 * ldetSigma - ((df + p + 1) / 2) * ldetX - ltrace / 2
   
   if (log == FALSE) { return(exp(ldens)) }
   return(ldens)
}

# Press/Anderson/... parametrization
dwishart <- function(X, df, Sigma, log = FALSE, is.chol = FALSE) {
   p <- ncol(X)

   # The constant
   lc0 <- -df * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((df + 1 - (1:p))/2)))

   if (!is.chol) {
      X.chol <- chol(X)
      Sigma.chol <- chol(Sigma)
   } else {
      X.chol <- X
      Sigma.chol <- Sigma
   }
   ldetSigma <- 2 * sum(log(diag(Sigma.chol)))
   ldetX <- 2 * sum(log(diag(X.chol)))
   Sigma.chol.inv <- backsolve(r = Sigma.chol, x = diag(p))      # inverse of Cholesky factor (NOT chol2inv(Sigma.chol) = solve(X)!)
   
   # The log(trace):
   # ltrace <- sum(diag(chol2inv(Sigma.chol) %*% crossprod(X.chol)))
   # ltrace <- sum((Sigma.chol %*% solve(X.chol))^2)
   ltrace <- sum((X.chol %*% Sigma.chol.inv)^2)

   ldens <- lc0 + (df - p - 1)/2 * ldetX - (df / 2) * ldetSigma - ltrace / 2

   if (log == FALSE) { return(exp(ldens)) }
   return(ldens)
}


# Multivariate normal, without symmetry check
dmvnorm_fast <- function(x, mean = rep(0, p), sigma = diag(p), log = FALSE, is.chol = FALSE) 
{
   if (is.vector(x)) 
      x <- matrix(x, ncol = length(x))
   p <- ncol(x)
   if (!missing(mean)) {
      if (!is.null(dim(mean))) 
         dim(mean) <- NULL
      if (length(mean) != p) 
         stop("mean and sigma have non-conforming size")
   }
   if (!missing(sigma)) {
      if (p != ncol(sigma)) 
         stop("x and sigma have non-conforming size")
   }
   if (is.chol == FALSE) {
      dec <- tryCatch(chol(sigma), error = function(e) e)  
   } else {
      dec <- sigma
   }
   if (inherits(dec, "error")) {
      x.is.mu <- colSums(t(x) != mean) == 0
      logretval <- rep.int(-Inf, nrow(x))
      logretval[x.is.mu] <- Inf
   }
   else {
      tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
      rss <- colSums(tmp^2)
      logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
   }
   names(logretval) <- rownames(x)
   if (log) 
      return(logretval)
   else
      return(exp(logretval))
}

# Generate random Inverted Wishart
#
# Using Press parametrization
# v > 2p
#
# X ~ IW(v, S)
# E[X] = S / (v - 2(p + 1))
riwish_Press <- function(v, S){
   p <- nrow(S)
   stopifnot(v > 2*p)
   v.Anderson <- v - p - 1
   return(MCMCpack::riwish(v.Anderson, S))
}
