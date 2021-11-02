# bayessource

This R package evaluates through the Bayes Factor whether two sets of samples come from the same Multivariate Gaussian distribution or not.

It currently implements a fast Gibbs sampler for the Multivariate Normal - Inverse Wishart model.

## Installation

The package is not on CRAN yet.   
It must be installed using `devtools` or `remotes` from this repository:

```r
# install.packages('remotes')
remotes::install_github('lgaborini/bayessource')
```

<!-- badges: start -->
[![R build status](https://github.com/lgaborini/bayessource/workflows/R-CMD-check/badge.svg)](https://github.com/lgaborini/bayessource/actions) [![Codecov test coverage](https://codecov.io/gh/lgaborini/bayessource/branch/master/graph/badge.svg)](https://codecov.io/gh/lgaborini/bayessource?branch=master) [![DOI](https://zenodo.org/badge/163959964.svg)](https://zenodo.org/badge/latestdoi/163959964)

<!-- badges: end -->


## Documentation

Documentation is available on [GitHub pages](https://lgaborini.github.io/bayessource), or in the [docs/index.html]([docs/index.html]) file of the repository.

Also see the vignettes:

- [Introduction](https://lgaborini.github.io/bayessource/articles/introduction.html)
- [Prior elicitation](https://lgaborini.github.io/bayessource/articles/prior_elicitation.html)

### Main functions

- `make_priors_and_init()`: obtain hyperpriors and initialization from a background dataset
- `marginalLikelihood()`: fast computation of the marginal likelihood
- `samesource_C()`: fast computation of the Bayes Factor (same source vs. different sources)
- `mcmc_postproc()`: collect and tidy posterior samples from this package

## Extending

### Writing Rd documentation

The documentation uses some roxygen2 Rd templates to enter parametrization/model details.
These are stored in the directory `man-roxygen`.
When updating Rd templates, one must pay attenton that:

- LaTeX is supported only through the Rd `\eqn{latex}{ascii}` and `\deqn{latex}{ascii}` tags.
- it is best to write plain Rd or roxygen2 tags rather than Markdown tags
- sections must start with the `\@section title:` and end up after the Details.
  Do not forget the `:` at the end.


## References

Bozza, Taroni, Marquis, Schmittbuhl, “Probabilistic Evaluation of Handwriting Evidence: Likelihood Ratio for Authorship.” Journal of the Royal Statistical Society: Series C (Applied Statistics) 57, no. 3 (June 1, 2008): 329–41. https://doi.org/10.1111/j.1467-9876.2007.00616.x.
