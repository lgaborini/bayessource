# knitr markup functions
requireNamespace('xtable')

bmatrix <- function(x, digits = NULL, display = TRUE, pre = '', post = '', ...) {

   default_args = list(
      include.colnames = FALSE,
      only.contents = TRUE,
      include.rownames = FALSE,
      hline.after = NULL,
      comment = FALSE,
      print.results = FALSE
   )
   passed_args = list(...)
   calling_args = c(list(x = xtable::xtable(x, digits = digits)),
                    c(passed_args,
                      default_args[setdiff(names(default_args), names(passed_args))]))
   if (display) {
      pre.tag <- '$$'
      post.tag <- '$$'
   } else {
      pre.tag <- ''
      post.tag <- ''
   }

   return(cat(
      pre.tag,
      pre,
      "\\begin{bmatrix}\n",
      do.call(xtable::print.xtable, calling_args),
      "\\end{bmatrix}",
      post,
      post.tag
   ))
}

# #' Example:
# #+results='asis', echo=FALSE
# # bmatrix(diag(2), pre='A =')
