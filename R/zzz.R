#
# Package on-attach function

# Currently prints out Cholesky flag
.onLoad <- function(...) {
	s <- ifelse(bayessource:::isCholeskyOn(),
			'Using Cholesky decomposition (experimental).',
			'Not using Cholesky decomposition (safe).')
	packageStartupMessage(paste0('Package bayessource loaded. ', s))
}

