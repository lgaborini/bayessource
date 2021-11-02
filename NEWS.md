# bayessource 0.3.2

* Improved documentation.
* Add examples.
* Pkgdown site: adding function families.
* Improved vignettes.

# bayessource 0.3.1

* Fig bug in `make_priors_and_init` when a background item is seen only once.    
  Now warns or fails gracefully.

# bayessource 0.3

* Add prior and init elicitation functions.
* Add vignette for prior elicitation.
* Using Roxygen with Markdown. Add Rdpack dependency to generate references with Bibtex.
* Add runtime checks for matrix invertibility.
* CRITICAL: fix prior matrices estimation.
* Removed dplyr Suggests.
* Fix incompatibility with tibbles.
* Updated to conform to tolerance in new RcppArmadillo `inv_sympd`: introduced tolerant function `inv_sympd_chol`.
* `samesource_C` gains parameter `marginals`: if `TRUE`, return also all log-marginals in the BF calculation.

# bayessource 0.2.6

* Removed dependencies magrittr, matrixStats.

# bayessource 0.2.5

* Removed OpenMP support.
* Modified vignette.
* Added unit test to check sign of LR in known conditions.

# bayessource 0.2.4

* Added information if Wishart is degenerate.
* Added 1D Wishart test.
* Added `bayessource::rmvnorm` RNG test.

# bayessource 0.2.3

* Fixed OpenMP support (now RcppArmadillo does not complain anymore).
* Fixed dependencies.
* Added `Rcpp::checkUserInterrupt()` to MCMC loops.

# bayessource 0.2.2

* Added a `NEWS.md` file to track changes to the package.



