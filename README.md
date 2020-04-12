# bayessource

This R package evaluates through the Bayes Factor whether two sets of samples come from the same Multivariate Gaussian distribution or not.

Currently implements a Gibbs sampler for the Multivariate Normal - Inverse Wishart model.

## Installation

The packages is not released on CRAN yet.   
It must be installed using `devtools` or `remotes` from this repository:

```
# install.packages('remotes')
remotes::install_github('lgaborini/bayessource')
```

<!-- badges: start -->
[![R build status](https://github.com/lgaborini/bayessource/workflows/R-CMD-check/badge.svg)](https://github.com/lgaborini/bayessource/actions) [![Codecov test coverage](https://codecov.io/gh/lgaborini/bayessource/branch/master/graph/badge.svg)](https://codecov.io/gh/lgaborini/bayessource?branch=master) [![DOI](https://zenodo.org/badge/163959964.svg)](https://zenodo.org/badge/latestdoi/163959964)

<!-- badges: end -->


## Documentation

Documentation is available on [GitHub pages](https://lgaborini.github.io/bayessource).

## References

Bozza, Taroni, Marquis, Schmittbuhl, “Probabilistic Evaluation of Handwriting Evidence: Likelihood Ratio for Authorship.” Journal of the Royal Statistical Society: Series C (Applied Statistics) 57, no. 3 (June 1, 2008): 329–41. https://doi.org/10.1111/j.1467-9876.2007.00616.x.
