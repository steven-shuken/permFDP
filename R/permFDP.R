library(Rcpp)
library(BH)
Sys.setenv("PKG_CXXFLAGS"="-std=c++17")
sourceCpp("src/permFDP.cpp")

#' Permutation-Based FDP Method for Rejection Threshold Correction
#'
#' This function controls FDR using the permutation method described in our manuscript. Like the BH method above, it corrects the rejection threshold rather than the p-values themselves.
#' It returns the new threshold for P-value rejection.
#' It uses the Rcpp and BH packages to leverage fast C++ code.
#' @param pVals Vector of p-values. The length of this vector must be the same as the number of rows in the intOnly data frame
#' @param threshold Original threshold that will be adjusted.
#' @param myDesign Vector of 1s and 2s specifying which columns in intOnly belong to the control group (1) and which belong to the test group (2).
#' @param intOnly A data frame. Each column is a sample, each row is an analyte.
#' @param nPerms The number of permutations to perform. At least 100 is recommended.
#' @keywords p-values FDP FDR permutation
#' @export
#' @examples
#' controlVals = matrix(rnorm(300), ncol = 3, nrow = 100)
#' testVals = matrix(rnorm(300, mean = 3), ncol = 3, nrow = 100)
#' intOnly = data.frame(cbind(controlVals, testVals))
#' myDesign = c(1,1,1,2,2,2)
#' pVals = c()
#' for (row in 1:nrow(intOnly)) {
#' pVals = c(pVals, t.test(intOnly[row, 1:3], intOnly[row, 4:6])$p.value)
#' }
#' threshold = 0.05
#' corrThreshold = permFDP::permFDP.adjust.threshold(pVals, threshold, myDesign, intOnly, 100)
#' corrThreshold

permFDP.adjust.threshold = function(pVals, threshold, myDesign, intOnly, nPerms) {
  pVals = pVals[order(pVals)]
  intMatrix = as.matrix(intOnly)
  nc = length(which(myDesign == 1))
  nt = length(which(myDesign == 2))
  return(permFDRAdjustCpp(pVals, threshold, myDesign, intMatrix, nPerms, nc, nt))
}
