#
# Package on-attach function

.onLoad <- function(...) {

   # isCholeskyOn()

	# s <- ifelse(bayessource:::isCholeskyOn(),
	# 		'Using Cholesky decomposition (experimental).',
	# 		'Not using Cholesky decomposition (safe).')
	# packageStartupMessage(paste0('Package bayessource loaded. ', s))
}

# Detach the DLL when unloading the library
.onUnload <- function(libpath) {
   library.dynam.unload('bayessource', libpath)
}
