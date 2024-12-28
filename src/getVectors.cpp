#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector getVectors(DataFrame df, double rt_sim, double mz_sim) {  // mz_sim is now in ppm
  NumericVector RT = df["RT"];
  NumericVector MZ = df["MZ"];
  CharacterVector Compound_ID = df["Compound_ID"];

  int n = RT.size();
  std::vector<std::string> rt_metabolites;
  std::vector<std::string> mz_metabolites;
  std::vector<std::string> metabolites;

  // Sort by RT and MZ using indices
  IntegerVector rt_order = match(clone(RT).sort(), RT);
  IntegerVector mz_order = match(clone(MZ).sort(), MZ);

  // RT logic remains the same
  for(int i = 0; i < n - 1; ++i) {
    int idx = rt_order[i] - 1;
    double rt = RT[idx];
    double diff_rt = RT[rt_order[i+1] - 1] - rt;

    if (i == 0 && diff_rt > rt_sim) {
      rt_metabolites.push_back(Rcpp::as<std::string>(Compound_ID[idx]));
    } else if (i > 0) {
      double diff_up_rt = rt - RT[rt_order[i-1] - 1];
      if (diff_rt >= rt_sim && diff_up_rt >= rt_sim) {
        rt_metabolites.push_back(Rcpp::as<std::string>(Compound_ID[idx]));
      }
    }
  }

  // MZ logic modified to use ppm
  for(int i = 0; i < n - 1; ++i) {
    int idx = mz_order[i] - 1;
    double mz = MZ[idx];
    // Calculate ppm difference with next mass
    double diff_mz_ppm = (MZ[mz_order[i+1] - 1] - mz) / mz * 1e6;

    if (i == 0 && diff_mz_ppm > mz_sim) {
      mz_metabolites.push_back(Rcpp::as<std::string>(Compound_ID[idx]));
    } else if (i > 0) {
      // Calculate ppm difference with previous mass
      double diff_up_mz_ppm = (mz - MZ[mz_order[i-1] - 1]) / MZ[mz_order[i-1] - 1] * 1e6;
      if (diff_mz_ppm > mz_sim && diff_up_mz_ppm > mz_sim) {
        mz_metabolites.push_back(Rcpp::as<std::string>(Compound_ID[idx]));
      }
    }
  }

  // Rest of the code remains the same
  std::set<std::string> unique_metabolites(rt_metabolites.begin(), rt_metabolites.end());
  unique_metabolites.insert(mz_metabolites.begin(), mz_metabolites.end());

  for (auto& id : unique_metabolites) {
    metabolites.push_back(id);
  }

  CharacterVector result(metabolites.size());
  for (size_t i = 0; i < metabolites.size(); ++i) {
    result[i] = metabolites[i];
  }

  if (result.size() == 0) {
    stop("Not enough features isolated. Consider lowering RT and MZ thresholds.");
  }

  return result;
}