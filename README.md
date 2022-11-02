# permFDP
An R package for correcting p-values for multiple hypothesis testing in comparative quantitative omics experiments using permutation-based FDP estimation.

## Installation

Once you have the `devtools` package installed, you can install `permFDP` by using the `install_github()` command:

```
install_github('steven-shuken/permFDP')
```

Expect to see two warnings.

## Usage

The only function currently intended for use is `permFDP.adjust.threshold()`. Use `?permFDP::permFDP.adjust.threshold()` for instructions. Briefly, you need:
1. A vector of p-values
2. The uncorrected rejection threshold
3. A data frame of quantities (for generating the permutated hit counts). Note: the row order does not have to match the p-value order, i.e., the first p-value does not have to correspond to the first analyte/row, etc.
4. A vector of 1s and 2s specifying which columns in the data frame are in which group
5. The number of permutations to perform (at least 100 is recommended).

The function will return the corrected rejection threshold to control the FDR according to the rejection threshold you supplied. E.g., if you supply an uncorrected rejection threshold of 0.05, the estimated FDR will now be 5% using the new corrected rejection threshold.
