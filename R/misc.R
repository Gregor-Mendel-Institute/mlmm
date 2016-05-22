##' @import emma

.onAttach <- function(libname, pkgname) {
  if(! requireNamespace("utils", quietly=TRUE))
    stop("Pkg utils needed for this function to work. Please install it.",
         call.=FALSE)
  msg <- paste0("package '", pkgname,
                "' (version ", utils::packageVersion(pkgname), ")",
                " is loaded",
                "\ndev at https://github.com/Gregor-Mendel-Institute/mlmm")
  packageStartupMessage(msg)
}
