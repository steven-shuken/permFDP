// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <random>
#include <boost/math/statistics/t_test.hpp>

/*
 * Performs a T test on the measurements according to the design (1s and 2s)
 * and returns a P value.
 */
// [[Rcpp::export]]
double designTTest(std::vector<double> ints, std::vector<int> design) {
  if (ints.size() != design.size()) Rcpp::stop("ERROR: DESIGN VECTOR AND MEASUREMENT VECTOR NOT EQUAL IN LENGTH!");

  std::vector<double> cIntensities, tIntensities;

  for (int i = 0; i < design.size(); i++) {
    if (design[i] == 1) {
      cIntensities.push_back(ints[i]);
    } else if (design[i] == 2) {
      tIntensities.push_back(ints[i]);
    } else {
      Rcpp::stop("DESIGN SYMBOL IS NOT 1 OR 2!");
    }
  }

  auto [t, p] = boost::math::statistics::two_sample_t_test(cIntensities, tIntensities);
  return p;
}

/*
 * This function randomizes an experimental design while keeping the test and control samples as balanced as possible.
 */
std::vector<int> randBalDesign(int nc, int nt) {

  // Calculate the maximally balanced number of control samples to put into the randomized control group
  int cIntoC = round(nc * nc / (nc + nt));

  // Put control and test samples into each group accordingly
  int tIntoC = nc - cIntoC;
  int cIntoT = nc - cIntoC;
  int tIntoT = nt - cIntoT;

  std::vector<int> controlGroup;
  std::vector<int> testGroup;
  for (int cc = 0; cc < cIntoC; cc++) {
    controlGroup.push_back(1);
  }
  for (int ct = 0; ct < cIntoT; ct++) {
    testGroup.push_back(1);
  }
  for (int tc = 0; tc < tIntoC; tc++) {
    controlGroup.push_back(2);
  }
  for (int tt = 0; tt < tIntoT; tt++) {
    testGroup.push_back(2);
  }

  // obtain a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

  // Randomly shuffle each group
  std::shuffle(controlGroup.begin(), controlGroup.end(), std::default_random_engine(seed));
  std::shuffle(testGroup.begin(), testGroup.end(), std::default_random_engine(seed));

  // Return new design
  std::vector<int> output;
  for (int i = 0; i < controlGroup.size(); i++) {
    output.push_back(controlGroup[i]);
  }
  for (int i = 0; i < testGroup.size(); i++) {
    output.push_back(testGroup[i]);
  }

  return output;
}

/*
 * Returns the number of p values that are below the threshold.
 * The p values passed in must be sorted in ascending order.
 */
int countHits(std::vector<double> sortedPVals, double threshold) {
  for (int i=0; i < sortedPVals.size(); i++) {
    if (sortedPVals[i] > threshold) {
      return i;
    }
  }
  return sortedPVals.size();
}

/*
 * Returns the index of the highest fdp that is at or below the threshold.
 */
int    getHighestPositionBelowThresh(std::vector<double> fdp, double threshold) {
  std::vector<int> positions;
  for (int i=0; i < fdp.size(); i++) {
    if (fdp[i] <= threshold) {
      positions.push_back(i);
    }
  }

  if (positions.size() < 1)
    return -1;

  return *max_element(positions.begin(), positions.end());
}

/*
 * This function controls the FDR using the permutation method described in our manuscript.
 * Like the BH method, it corrects the rejection threshold rather than the p-values themselves.
 * It returns the new threshold for p-value rejection.
 * intOnly is a vector of samples. Each sample is a vector of proteins.
 */
// [[Rcpp::export]]
double permFDRAdjust(std::vector<double> expPs, double threshold, std::vector<int> design,
                     std::vector<std::vector<double>> intOnly, int nPerms, int nc, int nt) {
  // Instantiate the lists of random p-values to estimate false hits
  std::vector<double> firstSample = intOnly[0];
  int nProts = firstSample.size();
  std::vector<std::vector<double>> randPVals(nPerms);

  // Get experimentally observed p-values
  if (expPs.size() < 1) {
    for (int i_prot = 0; i_prot < nProts; i_prot++) {
      std::vector<double> measurements;
      int nSamples = intOnly.size();
      for (int i_sample = 0; i_sample < nSamples; i_sample++) {
        measurements.push_back(intOnly[i_sample][i_prot]);
      }
      expPs[i_prot] = designTTest(measurements, design);
    }
  }
  std::sort(expPs.begin(), expPs.end()); // sort ascending

  // For each permutation, perform t-tests and save the p-values
  for (int i_perm = 0; i_perm < nPerms; i_perm++) {

    // Get randomized balanced design
    std::vector<int> balDesign = randBalDesign(nc, nt);

    // Instantiate vector of randomized T test P values
    std::vector<double> pVals_i(intOnly[0].size());

    // Pull out the protein intensities and do the t tests
    for (int i_prot = 0; i_prot < intOnly[0].size(); i_prot++) {
      std::vector<double> cIntensities;
      std::vector<double> tIntensities;

      for (int i_sample = 0; i_sample < balDesign.size(); i_sample++) {
        if (balDesign[i_sample] == 1) {
          cIntensities.push_back(intOnly[i_sample][i_prot]);
        } else if (balDesign[i_sample] == 2) {
          tIntensities.push_back(intOnly[i_sample][i_prot]);
        } else {
          Rcpp::stop("ERROR: DESIGN SYMBOL IS NOT 1 OR 2!");
        }
      }

      auto [t, p] = boost::math::statistics::two_sample_t_test(cIntensities, tIntensities);

      pVals_i[i_prot] = p;
    }
    randPVals[i_perm] = pVals_i;

    // Sort the p values in this permutation in ascending order
    std::sort(randPVals[i_perm].begin(), randPVals[i_perm].end());
  }

  // Make rank vector: 1 to m, and the estimated FDP vector.
  std::vector<double> rank(expPs.size());
  std::vector<double> fdp(expPs.size());
  for (int i = 0; i < expPs.size(); i++) {
    rank[i] = 1.0 + i;
  }

  // At every p value (rank), get the # hits for each permutation, then store the average.
  for (int i = 0; i< rank.size(); i++) {
    double thresh_i = expPs[i];
    std::vector<int> hitCounts(nPerms);

    // Calculate the average # hits at this threshold
    for (int i_perm = 0; i_perm < randPVals.size(); i_perm++) {
      hitCounts[i_perm] = countHits(randPVals[i_perm], thresh_i);
    }

    double sum = std::accumulate(hitCounts.begin(), hitCounts.end(), 0.0);
    double v = sum / hitCounts.size();
    fdp[i] = v / rank[i];
  }

  // Return the highest p value for which est FDP <= threshold
  int bestIndex = getHighestPositionBelowThresh(fdp, threshold);

  // index = -1: nothing is below the threshold.
  // return a new threshold below all the P values.
  if (bestIndex < 0)
    return expPs[0] / 2;

  // index is at the end: all p-values are below the threshold.
  // Return something just above the highest or halfway between the highest and 1.
  if (bestIndex == fdp.size() - 1) {
    double worstP = expPs[expPs.size() - 1];
    if (worstP + 0.05 <= 1)
      return worstP + 0.05;

    return (worstP + 1) / 2;
  }

  // Otherwise: return a threshold between the p-value at the highest index and that at the next one.
  return (expPs[bestIndex] + expPs[bestIndex + 1]) / 2;
}

// [[Rcpp::export]]
double permFDRAdjustCpp(Rcpp::NumericVector expPs, double threshold, Rcpp::NumericVector design,
                        Rcpp::NumericMatrix intMatrix, int nPerms, int nc, int nt) {
  // Convert NumericMatrix to vector of vectors
  std::vector<std::vector<double>> intVecVec(intMatrix.ncol());
  for (int i_col = 0; i_col < intVecVec.size(); i_col++) {
    std::vector<double> column(intMatrix.nrow());
    for (int i_row = 0; i_row < column.size(); i_row++) {
      column[i_row] = intMatrix(i_row, i_col);
    }
    intVecVec[i_col] = column;
  }

  // Convert everything else
  std::vector<double> expPVec(expPs.begin(), expPs.end());
  std::vector<int> designVec(design.begin(), design.end());

  return permFDRAdjust(expPVec, threshold, designVec, intVecVec, nPerms, nc, nt);
}
