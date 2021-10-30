#ifndef _CONFIG_H
#define _CONFIG_H

// Global configuration file
//
// Store constants, package-wide stuff
//-----------------------------------------------

// Enable C++11 (also look for flags in Makefiles)
// [[Rcpp::plugins(cpp11)]]



// Use Cholesky factorization when possible.
// Matrices are propagated through their Cholesky factors.
#ifndef USE_CHOLESKY
#define USE_CHOLESKY false
// #define USE_CHOLESKY true
#endif

//' Check whether Cholesky factorization is used.
//' @keywords internal
//' @family C++ functions
// [[Rcpp::export(rng = false)]]
bool isCholeskyOn();



#endif
