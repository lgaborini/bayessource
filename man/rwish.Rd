% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rwish}
\alias{rwish}
\title{Generate random sample from Wishart (faster).}
\usage{
rwish(v, S, is_chol = FALSE, return_chol = FALSE)
}
\arguments{
\item{v}{dof}

\item{S}{the scale matrix (pxp)}

\item{is_chol}{if TRUE, S is the upper Cholesky factor of S}

\item{return_chol}{if TRUE, the upper Cholesky factor is returned}
}
\value{
a single random variate from W(v, S)
}
\description{
Same code as \code{\link[stats]{rWishart}} function in package \pkg{stats}.
}
\section{Wishart parametrization}{


Uses \insertCite{Press2012Applied,Anderson2003introduction}{bayessource} parametrization.

\deqn{X ~ W(v, S)}
with \eqn{S = pxp} matrix, \eqn{v >= p} (the degrees of freedom).

Then:
\deqn{E[X] = v * S}
}

\references{
\insertAllCited{}
}
