test_that("function returns numeric threshold", {
  # Create valid test data
  set.seed(123)
  controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
  testVals <- matrix(rnorm(30, mean = 1), ncol = 3, nrow = 10)
  intOnly <- data.frame(cbind(controlVals, testVals))
  myDesign <- c(1, 1, 1, 2, 2, 2)
  pVals <- runif(10, 0, 1)
  threshold <- 0.05
  nPerms <- 100

  result <- permFDP.adjust.threshold(
    pVals,
    threshold,
    myDesign,
    intOnly,
    nPerms
  )

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(is.numeric(result))
  expect_false(is.na(result))
})

test_that("function accepts boundary values correctly", {
  set.seed(456)
  controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
  testVals <- matrix(rnorm(30, mean = 1), ncol = 3, nrow = 10)
  intOnly <- data.frame(cbind(controlVals, testVals))
  myDesign <- c(1, 1, 1, 2, 2, 2)
  pVals <- runif(10, 0, 1)
  nPerms <- 100

  # Test with very small threshold
  expect_silent(
    permFDP.adjust.threshold(pVals, 0.001, myDesign, intOnly, nPerms)
  )

  # Test with large threshold
  expect_silent(
    permFDP.adjust.threshold(pVals, 0.999, myDesign, intOnly, nPerms)
  )

  # Test with p-values at boundaries
  pVals_boundary <- c(0, 0.5, 1, runif(7, 0, 1))
  expect_silent(
    permFDP.adjust.threshold(pVals_boundary, 0.05, myDesign, intOnly, nPerms)
  )
})

test_that("function handles different sample sizes", {
  set.seed(789)
  nPerms <- 100
  threshold <- 0.05

  # Small dataset (2 samples per group)
  controlVals <- matrix(rnorm(20), ncol = 2, nrow = 10)
  testVals <- matrix(rnorm(20, mean = 1), ncol = 2, nrow = 10)
  intOnly <- data.frame(cbind(controlVals, testVals))
  myDesign <- c(1, 1, 2, 2)
  pVals <- runif(10, 0, 1)

  expect_silent(
    permFDP.adjust.threshold(pVals, threshold, myDesign, intOnly, nPerms)
  )

  # Medium dataset (5 samples per group)
  controlVals <- matrix(rnorm(50), ncol = 5, nrow = 10)
  testVals <- matrix(rnorm(50, mean = 1), ncol = 5, nrow = 10)
  intOnly <- data.frame(cbind(controlVals, testVals))
  myDesign <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
  pVals <- runif(10, 0, 1)

  expect_silent(
    permFDP.adjust.threshold(pVals, threshold, myDesign, intOnly, nPerms)
  )

  # Unbalanced design (3 control, 5 test)
  controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
  testVals <- matrix(rnorm(50, mean = 1), ncol = 5, nrow = 10)
  intOnly <- data.frame(cbind(controlVals, testVals))
  myDesign <- c(1, 1, 1, 2, 2, 2, 2, 2)
  pVals <- runif(10, 0, 1)

  expect_silent(
    permFDP.adjust.threshold(pVals, threshold, myDesign, intOnly, nPerms)
  )
})

test_that("function handles different numbers of analytes", {
  set.seed(321)
  myDesign <- c(1, 1, 1, 2, 2, 2)
  threshold <- 0.05
  nPerms <- 100

  # Very few analytes
  controlVals <- matrix(rnorm(6), ncol = 3, nrow = 2)
  testVals <- matrix(rnorm(6, mean = 1), ncol = 3, nrow = 2)
  intOnly <- data.frame(cbind(controlVals, testVals))
  pVals <- runif(2, 0, 1)

  expect_silent(
    permFDP.adjust.threshold(pVals, threshold, myDesign, intOnly, nPerms)
  )

  # Many analytes
  controlVals <- matrix(rnorm(1500), ncol = 3, nrow = 500)
  testVals <- matrix(rnorm(1500, mean = 1), ncol = 3, nrow = 500)
  intOnly <- data.frame(cbind(controlVals, testVals))
  pVals <- runif(500, 0, 1)

  expect_silent(
    permFDP.adjust.threshold(pVals, threshold, myDesign, intOnly, nPerms)
  )
})

test_that("function is reproducible with same input", {
  # Create test data
  set.seed(999)
  controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
  testVals <- matrix(rnorm(30, mean = 1), ncol = 3, nrow = 10)
  intOnly <- data.frame(cbind(controlVals, testVals))
  myDesign <- c(1, 1, 1, 2, 2, 2)
  pVals <- runif(10, 0, 1)
  threshold <- 0.05
  nPerms <- 100

  # Run twice with same seed - results should be identical
  # Note: The function uses random permutations internally, so we can't guarantee
  # exact reproducibility without controlling the C++ RNG, but the results should
  # be in the same general range
  result1 <- permFDP.adjust.threshold(
    pVals,
    threshold,
    myDesign,
    intOnly,
    nPerms
  )
  result2 <- permFDP.adjust.threshold(
    pVals,
    threshold,
    myDesign,
    intOnly,
    nPerms
  )

  # expect_equal allows for small numerical differences (default tolerance ~1.5e-8)
  expect_equal(result1, result2)
})

test_that("function handles matrix input for intOnly", {
  set.seed(654)
  controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
  testVals <- matrix(rnorm(30, mean = 1), ncol = 3, nrow = 10)
  intMatrix <- cbind(controlVals, testVals)
  myDesign <- c(1, 1, 1, 2, 2, 2)
  pVals <- runif(10, 0, 1)
  threshold <- 0.05
  nPerms <- 100

  # Should work with matrix input
  result_matrix <- permFDP.adjust.threshold(
    pVals,
    threshold,
    myDesign,
    intMatrix,
    nPerms
  )

  # Should also work with data.frame input
  result_df <- permFDP.adjust.threshold(
    pVals,
    threshold,
    myDesign,
    as.data.frame(intMatrix),
    nPerms
  )

  expect_type(result_matrix, "double")
  expect_type(result_df, "double")
})
