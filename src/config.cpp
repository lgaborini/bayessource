#include "config.h"

// Macro to return Cholesky status. Must be callable (hidden) from R.
bool isCholeskyOn(){
   return (USE_CHOLESKY);
}
