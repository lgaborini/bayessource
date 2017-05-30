#
# Statistical distributions R functions
#
# Wishart density
# multivariate normal density (faster)
# utilities


# library(matrixStats)
# loadNamespace('MCMCpack')

#' Accurately compute the log-mean of exponentials
#'
#' Accurately compute log(mean(x)), where x = log(sum(exp(y))).
#' \url{https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/}
#'
#' Uses \pkg{matrixStats}::logSumExp
#'
#' @param x x
#' @return log(mean(x))
#' @keywords internal
logSumExpMean <- function(x) {
   ret <- matrixStats::logSumExp(x) - log(length(x))
   return(ret)
}

#' Computes the log(gamma(theta)) where theta is a vector from R^p (multivariate gamma)
#'
#' @param theta the variable
#' @param p dimension of theta
#'
#' @keywords internal
logMultivariateGamma <- function(theta, p) {
   p * (p - 1)/4 * log(pi) + sum(log(gamma(theta + (1 - (1:p))/2)))
}

#' Normalization constant, Press parametrization
#'
#' @param df dof
#' @param p size
#'
#' @keywords internal
iwishart.logconstant <- function(df, p) {
   -(df - p - 1) * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((df - p - (1:p))/2)))
}


#' Inverted Wishart density, parametrization according to Press.
#'
#' @param X observation
#' @param df degrees of freedom
#' @param Sigma scale matrix
#' @param log if TRUE, return the log-density
#' @param is.chol if TRUE, Sigma and X are the Cholesky factors of Sigma and X
#' @return the density in X
#' @export
#' @seealso \code{\link{diwishart_Anderson}, \code{\link{diwishart_inverse}}}
#'
#' @template InverseWishart_Press
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



#' Inverted Wishart density from the inverse
#'
#' Computes the density of an Inverted Wishart (df, Sigma) in X, by supplying (X^(-1), df, Sigma) rather than (X, df, Sigma).
#' Avoids a matrix inversion.
#' 
#' Computes the pdf p_X(x) by knowing x^(-1)
#'
#' @param X.inv inverse of X (the observation)
#' @param df degrees of freedom
#' @param Sigma scale matrix
#' @param log if TRUE, return the log-density
#' @param is.chol if TRUE, Sigma and X.inv are the Cholesky factors of Sigma and X.inv
#' @return the density in X
#' @export
#' @seealso \code{\link{diwishart}}
#'
#' @template InverseWishart_Press
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


#' Inverted Wishart density, parametrization according to Anderson
#'
#' @param X the observation
#' @param df degrees of freedom
#' @param Sigma scale matrix
#' @param log if TRUE, return the log-density
#' @param is.chol if TRUE, Sigma and X are the Cholesky factors of Sigma and X
#' @return the density in X
#' @export
#' @seealso \code{\link{diwishart}}
#'
#' @template InverseWishart_Anderson
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

#' Wishart density.
#'
#' Parametrization according to Press/Anderson.
#'
#' @param X the observation
#' @param df degrees of freedom
#' @param Sigma scale matrix
#' @param log if TRUE, return the log-density
#' @param is.chol if TRUE, Sigma and X are the Cholesky factors of Sigma and X
#' @return the density in X
#' @export
#' @seealso \code{\link{riwish_Press}}
#'
#' @template Wishart_eqn
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


#' Multivariate normal density. Faster, assumes symmetry.
#'
#' @param x the observation
#' @param mean mean vector
#' @param sigma covariance matrix
#' @param log if TRUE, return the log-density
#' @param is.chol if TRUE, sigma is the Cholesky factor of sigma
#' @return the density in X
#' @export
#'
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

#' Generate random sample from Inverted Wishart.
#'
#' Using Press parametrization.
#'
#' Uses \pkg{MCMCpack}::riwish.
#'
#' @param v dof (\eqn{> 2p})
#' @param S the scale matrix (pxp)
#' @return a single random variate from IW(v, S)
#' @export
#'
#' @template InverseWishart_Press
riwish_Press <- function(v, S){
   p <- nrow(S)
   stopifnot(v > 2*p)
   v.Anderson <- v - p - 1
   return(MCMCpack::riwish(v.Anderson, S))
}
