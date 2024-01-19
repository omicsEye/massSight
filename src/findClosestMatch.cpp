#include <Rcpp.h>
using namespace Rcpp;

// Helper function to calculate rms
double calculateRMS(NumericVector query, NumericVector hit, NumericVector stds) {
  double rt_diff = query[0] - hit[0];
  double mz_diff = query[1] - hit[1];
  double intensity_diff = query[2] - hit[2];

  // Calculate the rms score based on the provided stds
  double rms_score = sqrt(pow(rt_diff/stds[0], 2) + pow(mz_diff/stds[1], 2) + pow(intensity_diff/stds[2], 2));
  return rms_score;
}

// [[Rcpp::export]]
int findClosestMatch(NumericVector query, DataFrame ref, NumericVector stds) {
  NumericVector refRT = ref["RT"];
  NumericVector refMZ = ref["MZ"];
  NumericVector refIntensity = ref["Intensity"];
  IntegerVector refCompoundID = ref["Compound_ID"];

  double mass_plus, mass_minus;
  int best_index = -1;
  double min_score = R_PosInf; // Initialize to positive infinity

  for (int i = 0; i < refRT.size(); ++i) {
    // Apply cutoffs for RT
    if ((refRT[i] <= query[0] + stds[0]) && (refRT[i] >= query[0] - stds[0])) {
      // Apply cutoffs for MZ
      mass_plus = query[1] * stds[1] / 10000;
      mass_minus = query[1] * stds[1] / (10000 + stds[1]);
      if ((refMZ[i] < query[1] + mass_plus) && (refMZ[i] > query[1] - mass_minus)) {
        // Apply cutoffs for Intensity if there is a third standard deviation provided
        if (stds.size() > 2 && stds[2] > 0) {
          double cutoffs_2 = pow(10.0, stds[2]);
          if ((refIntensity[i] < query[2] * cutoffs_2) && (refIntensity[i] > query[2] / cutoffs_2)) {
            // Calculate rms score
            NumericVector hit = NumericVector::create(refRT[i], refMZ[i], refIntensity[i]);
            double score = calculateRMS(query, hit, stds);

            // Check if this is the best score so far
            if (score < min_score) {
              min_score = score;
              best_index = refCompoundID[i];
            }
          }
        }
      }
    }
  }

  // Return the index of the best match, or -1 if no match is found
  return best_index + 1;
}
