library(Rcpp)
library(BH)
Sys.setenv("PKG_CXXFLAGS"="-std=c++17")
sourceCpp("../src/PermFDRAdjust.cpp")

#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

# This function controls FDR using the permutation method described in our manuscript. Like the BH method above, it corrects the rejection threshold rather than the p-values themselves.
# It returns the new threshold for P-value rejection.
# It uses the Rcpp and BH packages to leverage fast C++ code.
permFDP.adjust.threshold = function(pVals, threshold, myDesign, intOnly, nPerms, nc, nt) {
  pVals = pVals[order(pVals)]
  intMatrix = as.matrix(intOnly)
  return(permFDRAdjustCpp(pVals, threshold, myDesign, intMatrix, nPerms, nc, nt))
}
