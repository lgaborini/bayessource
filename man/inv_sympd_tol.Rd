% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{inv_sympd_tol}
\alias{inv_sympd_tol}
\title{Compute the inverse of a symmetric positive definite matrix}
\usage{
inv_sympd_tol(U_sympd)
}
\description{
Compute the inverse of a symmetric positive definite matrix.
Does not output warnings on default symmetry tolerance: basically, symmetry is forced.
}
\references{
https://github.com/RcppCore/RcppArmadillo/issues/257
}
\keyword{internal}
