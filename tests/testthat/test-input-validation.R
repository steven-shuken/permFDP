test_that("pVals validation works correctly", {
  # Helper function to create valid test data
  create_valid_data <- function() {
    controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
    testVals <- matrix(rnorm(30, mean = 1), ncol = 3, nrow = 10)
    intOnly <- data.frame(cbind(controlVals, testVals))
    myDesign <- c(1, 1, 1, 2, 2, 2)
    pVals <- runif(10, 0, 1)
    threshold <- 0.05
    nPerms <- 100
    list(
      pVals = pVals,
      threshold = threshold,
      myDesign = myDesign,
      intOnly = intOnly,
      nPerms = nPerms
    )
  }

  # Test: pVals is missing
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "pVals is required and cannot be NULL"
  )

  # Test: pVals is NULL
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = NULL,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "pVals is required and cannot be NULL"
  )

  # Test: pVals is not numeric
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = c("a", "b", "c"),
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly[1:3, ],
      nPerms = data$nPerms
    ),
    "pVals must be a numeric vector"
  )

  # Test: pVals is empty
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = numeric(0),
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly[integer(0), ],
      nPerms = data$nPerms
    ),
    "pVals cannot be empty"
  )

  # Test: pVals contains NA values
  data <- create_valid_data()
  pVals_with_na <- data$pVals
  pVals_with_na[1] <- NA
  expect_error(
    permFDP.adjust.threshold(
      pVals = pVals_with_na,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "pVals contains NA values"
  )

  # Test: pVals has values < 0
  data <- create_valid_data()
  pVals_negative <- data$pVals
  pVals_negative[1] <- -0.1
  expect_error(
    permFDP.adjust.threshold(
      pVals = pVals_negative,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "All p-values must be between 0 and 1"
  )

  # Test: pVals has values > 1
  data <- create_valid_data()
  pVals_large <- data$pVals
  pVals_large[1] <- 1.5
  expect_error(
    permFDP.adjust.threshold(
      pVals = pVals_large,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "All p-values must be between 0 and 1"
  )
})

test_that("threshold validation works correctly", {
  create_valid_data <- function() {
    controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
    testVals <- matrix(rnorm(30, mean = 1), ncol = 3, nrow = 10)
    intOnly <- data.frame(cbind(controlVals, testVals))
    myDesign <- c(1, 1, 1, 2, 2, 2)
    pVals <- runif(10, 0, 1)
    threshold <- 0.05
    nPerms <- 100
    list(
      pVals = pVals,
      threshold = threshold,
      myDesign = myDesign,
      intOnly = intOnly,
      nPerms = nPerms
    )
  }

  # Test: threshold is missing
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "threshold is required and cannot be NULL"
  )

  # Test: threshold is NULL
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = NULL,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "threshold is required and cannot be NULL"
  )

  # Test: threshold is not numeric
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = "0.05",
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "threshold must be a single numeric value"
  )

  # Test: threshold is a vector
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = c(0.05, 0.10),
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "threshold must be a single numeric value"
  )

  # Test: threshold is NA
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = NA,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "threshold cannot be NA"
  )

  # Test: threshold is a vector with NA
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = c(NA, 0.05),
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "threshold must be a single numeric value"
  )

  # Test: threshold is a vector of NAs
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = c(NA, NA),
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "threshold must be a single numeric value"
  )

  # Test: threshold is <= 0
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = 0,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "threshold must be between 0 and 1 \\(exclusive\\)"
  )

  # Test: threshold is >= 1
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = 1,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "threshold must be between 0 and 1 \\(exclusive\\)"
  )

  # Test: threshold is negative
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = -0.05,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "threshold must be between 0 and 1 \\(exclusive\\)"
  )
})

test_that("myDesign validation works correctly", {
  create_valid_data <- function() {
    controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
    testVals <- matrix(rnorm(30, mean = 1), ncol = 3, nrow = 10)
    intOnly <- data.frame(cbind(controlVals, testVals))
    myDesign <- c(1, 1, 1, 2, 2, 2)
    pVals <- runif(10, 0, 1)
    threshold <- 0.05
    nPerms <- 100
    list(
      pVals = pVals,
      threshold = threshold,
      myDesign = myDesign,
      intOnly = intOnly,
      nPerms = nPerms
    )
  }

  # Test: myDesign is missing
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "myDesign is required and cannot be NULL"
  )

  # Test: myDesign is NULL
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = NULL,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "myDesign is required and cannot be NULL"
  )

  # Test: myDesign is not numeric
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = c("1", "1", "1", "2", "2", "2"),
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "myDesign must be a numeric vector"
  )

  # Test: myDesign is empty
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = numeric(0),
      intOnly = data$intOnly[, integer(0)],
      nPerms = data$nPerms
    ),
    "myDesign cannot be empty"
  )

  # Test: myDesign contains NA values
  data <- create_valid_data()
  myDesign_with_na <- data$myDesign
  myDesign_with_na[1] <- NA
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = myDesign_with_na,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "myDesign contains NA values"
  )

  # Test: myDesign contains values other than 1 and 2
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = c(1, 1, 1, 3, 3, 3),
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "myDesign must contain only values 1 \\(control\\) and 2 \\(test\\)"
  )

  # Test: myDesign contains values other than 1 and 2 (0)
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = c(0, 1, 1, 2, 2, 2),
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "myDesign must contain only values 1 \\(control\\) and 2 \\(test\\)"
  )
})

test_that("intOnly validation works correctly", {
  create_valid_data <- function() {
    controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
    testVals <- matrix(rnorm(30, mean = 1), ncol = 3, nrow = 10)
    intOnly <- data.frame(cbind(controlVals, testVals))
    myDesign <- c(1, 1, 1, 2, 2, 2)
    pVals <- runif(10, 0, 1)
    threshold <- 0.05
    nPerms <- 100
    list(
      pVals = pVals,
      threshold = threshold,
      myDesign = myDesign,
      intOnly = intOnly,
      nPerms = nPerms
    )
  }

  # Test: intOnly is missing
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      nPerms = data$nPerms
    ),
    "intOnly is required and cannot be NULL"
  )

  # Test: intOnly is NULL
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = NULL,
      nPerms = data$nPerms
    ),
    "intOnly is required and cannot be NULL"
  )

  # Test: intOnly is not a data frame or matrix
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = c(1, 2, 3),
      nPerms = data$nPerms
    ),
    "intOnly must be a data frame or matrix"
  )

  # Test: intOnly is empty (no rows)
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly[0, ],
      nPerms = data$nPerms
    ),
    "intOnly cannot be empty"
  )

  # Test: intOnly works with matrix input
  data <- create_valid_data()
  expect_silent(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = as.matrix(data$intOnly),
      nPerms = data$nPerms
    )
  )
})

test_that("nPerms validation works correctly", {
  create_valid_data <- function() {
    controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
    testVals <- matrix(rnorm(30, mean = 1), ncol = 3, nrow = 10)
    intOnly <- data.frame(cbind(controlVals, testVals))
    myDesign <- c(1, 1, 1, 2, 2, 2)
    pVals <- runif(10, 0, 1)
    threshold <- 0.05
    nPerms <- 100
    list(
      pVals = pVals,
      threshold = threshold,
      myDesign = myDesign,
      intOnly = intOnly,
      nPerms = nPerms
    )
  }

  # Test: nPerms is missing
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly
    ),
    "nPerms is required and cannot be NULL"
  )

  # Test: nPerms is NULL
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = NULL
    ),
    "nPerms is required and cannot be NULL"
  )

  # Test: nPerms is not numeric
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = "100"
    ),
    "nPerms must be a single numeric value"
  )

  # Test: nPerms is a vector
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = c(10, 20)
    ),
    "nPerms must be a single numeric value"
  )

  # Test: nPerms is NA
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = NA
    ),
    "nPerms cannot be NA"
  )

  # Test: nPerms is a vector with NA
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = c(NA, 10)
    ),
    "nPerms must be a single numeric value"
  )

  # Test: nPerms is a vector of NAs
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = c(NA, NA)
    ),
    "nPerms must be a single numeric value"
  )

  # Test: nPerms is not an integer
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = 10.5
    ),
    "nPerms must be an integer"
  )

  # Test: nPerms is < 1
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = 0
    ),
    "nPerms must be at least 1"
  )

  # Test: nPerms < 100 gives warning
  data <- create_valid_data()
  expect_warning(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = 50
    ),
    "nPerms < 100 may produce unreliable results"
  )

  # Test: nPerms >= 100 does not give warning
  data <- create_valid_data()
  expect_silent(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = 100
    )
  )
})

test_that("cross-validation between inputs works correctly", {
  create_valid_data <- function() {
    controlVals <- matrix(rnorm(30), ncol = 3, nrow = 10)
    testVals <- matrix(rnorm(30, mean = 1), ncol = 3, nrow = 10)
    intOnly <- data.frame(cbind(controlVals, testVals))
    myDesign <- c(1, 1, 1, 2, 2, 2)
    pVals <- runif(10, 0, 1)
    threshold <- 0.05
    nPerms <- 100
    list(
      pVals = pVals,
      threshold = threshold,
      myDesign = myDesign,
      intOnly = intOnly,
      nPerms = nPerms
    )
  }

  # Test: pVals length doesn't match intOnly rows
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals[1:5],
      threshold = data$threshold,
      myDesign = data$myDesign,
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "Length of pVals \\(5\\) must match number of rows in intOnly \\(10\\)"
  )

  # Test: myDesign length doesn't match intOnly columns
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals,
      threshold = data$threshold,
      myDesign = c(1, 1, 2, 2),
      intOnly = data$intOnly,
      nPerms = data$nPerms
    ),
    "Length of myDesign \\(4\\) must match number of columns in intOnly \\(6\\)"
  )

  # Test: Control group has < 2 samples
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals[1:5],
      threshold = data$threshold,
      myDesign = c(1, 2, 2, 2, 2),
      intOnly = data$intOnly[1:5, 1:5],
      nPerms = data$nPerms
    ),
    "Both control \\(n=1\\) and test \\(n=4\\) groups must have at least 2 samples"
  )

  # Test: Test group has < 2 samples
  data <- create_valid_data()
  expect_error(
    permFDP.adjust.threshold(
      pVals = data$pVals[1:5],
      threshold = data$threshold,
      myDesign = c(1, 1, 1, 1, 2),
      intOnly = data$intOnly[1:5, 1:5],
      nPerms = data$nPerms
    ),
    "Both control \\(n=4\\) and test \\(n=1\\) groups must have at least 2 samples"
  )

  # Test: Both groups have exactly 2 samples (valid edge case)
  data <- create_valid_data()
  expect_silent(
    permFDP.adjust.threshold(
      pVals = data$pVals[1:5],
      threshold = data$threshold,
      myDesign = c(1, 1, 2, 2),
      intOnly = data$intOnly[1:5, 1:4],
      nPerms = data$nPerms
    )
  )
})
