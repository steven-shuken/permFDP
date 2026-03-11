#' @useDynLib permFDP, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

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
#' controlVals <- matrix(rnorm(300), ncol = 3, nrow = 100)
#' testVals <- matrix(rnorm(300, mean = 3), ncol = 3, nrow = 100)
#' intOnly <- data.frame(cbind(controlVals, testVals))
#' myDesign <- c(1, 1, 1, 2, 2, 2)
#' pVals <- c()
#' for (row in 1:nrow(intOnly)) {
#'   pVals <- c(pVals, t.test(intOnly[row, 1:3], intOnly[row, 4:6])$p.value)
#' }
#' threshold <- 0.05
#' corrThreshold <- permFDP::permFDP.adjust.threshold(pVals, threshold, myDesign, intOnly, 100)
#' corrThreshold
permFDP.adjust.threshold <- function(
  pVals,
  threshold,
  myDesign,
  intOnly,
  nPerms
) {
  # Validate pVals
  if (missing(pVals) || is.null(pVals)) {
    stop("pVals is required and cannot be NULL")
  }
  if (!is.numeric(pVals)) {
    stop("pVals must be a numeric vector")
  }
  if (length(pVals) == 0) {
    stop("pVals cannot be empty")
  }
  if (any(is.na(pVals))) {
    stop("pVals contains NA values")
  }
  if (any(pVals < 0 | pVals > 1)) {
    stop("All p-values must be between 0 and 1")
  }

  # Validate threshold
  if (missing(threshold) || is.null(threshold)) {
    stop("threshold is required and cannot be NULL")
  }
  if (length(threshold) == 1 && is.na(threshold)) {
    stop("threshold cannot be NA")
  }
  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("threshold must be a single numeric value")
  }
  if (threshold <= 0 || threshold >= 1) {
    stop("threshold must be between 0 and 1 (exclusive)")
  }

  # Validate myDesign
  if (missing(myDesign) || is.null(myDesign)) {
    stop("myDesign is required and cannot be NULL")
  }
  if (!is.numeric(myDesign)) {
    stop("myDesign must be a numeric vector")
  }
  if (length(myDesign) == 0) {
    stop("myDesign cannot be empty")
  }
  if (any(is.na(myDesign))) {
    stop("myDesign contains NA values")
  }
  if (!all(myDesign %in% c(1, 2))) {
    stop("myDesign must contain only values 1 (control) and 2 (test)")
  }

  # Validate intOnly
  if (missing(intOnly) || is.null(intOnly)) {
    stop("intOnly is required and cannot be NULL")
  }
  if (!is.data.frame(intOnly) && !is.matrix(intOnly)) {
    stop("intOnly must be a data frame or matrix")
  }
  if (nrow(intOnly) == 0 || ncol(intOnly) == 0) {
    stop("intOnly cannot be empty")
  }

  # Validate nPerms
  if (missing(nPerms) || is.null(nPerms)) {
    stop("nPerms is required and cannot be NULL")
  }
  if (length(nPerms) == 1 && is.na(nPerms)) {
    stop("nPerms cannot be NA")
  }
  if (!is.numeric(nPerms) || length(nPerms) != 1) {
    stop("nPerms must be a single numeric value")
  }
  if (nPerms != as.integer(nPerms)) {
    stop("nPerms must be an integer")
  }
  if (nPerms < 1) {
    stop("nPerms must be at least 1 (100 or more recommended)")
  }
  if (nPerms < 100) {
    warning(
      "nPerms < 100 may produce unreliable results. At least 100 permutations are recommended."
    )
  }

  # Validate consistency between inputs
  if (length(pVals) != nrow(intOnly)) {
    stop(sprintf(
      "Length of pVals (%d) must match number of rows in intOnly (%d)",
      length(pVals),
      nrow(intOnly)
    ))
  }

  if (length(myDesign) != ncol(intOnly)) {
    stop(sprintf(
      "Length of myDesign (%d) must match number of columns in intOnly (%d)",
      length(myDesign),
      ncol(intOnly)
    ))
  }

  nc <- length(which(myDesign == 1))
  nt <- length(which(myDesign == 2))

  if (nc < 2 || nt < 2) {
    stop(sprintf(
      "Both control (n=%d) and test (n=%d) groups must have at least 2 samples",
      nc,
      nt
    ))
  }

  # All validations passed - proceed with computation
  pVals <- pVals[order(pVals)]
  intMatrix <- as.matrix(intOnly)
  return(permFDRAdjustCpp(
    pVals,
    threshold,
    myDesign,
    intMatrix,
    nPerms,
    nc,
    nt
  ))
}
