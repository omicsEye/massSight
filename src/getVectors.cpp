#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector getVectors(DataFrame df, double rt_sim, double mz_sim) {
  NumericVector RT = df["RT"];
  NumericVector MZ = df["MZ"];
  CharacterVector Compound_ID = df["Compound_ID"]; // Use CharacterVector for string IDs

  int n = RT.size();
  std::vector<std::string> rt_metabolites;
  std::vector<std::string> mz_metabolites;
  std::vector<std::string> metabolites;

  // Sort by RT and MZ using indices
  IntegerVector rt_order = match(clone(RT).sort(), RT);
  IntegerVector mz_order = match(clone(MZ).sort(), MZ);

  // RT logic
  for(int i = 0; i < n - 1; ++i) {
    int idx = rt_order[i] - 1; // Adjusting for 0-based indexing in C++
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

  // MZ logic
  for(int i = 0; i < n - 1; ++i) {
    int idx = mz_order[i] - 1; // Adjusting for 0-based indexing in C++
    double mz = MZ[idx];
    double diff_mz = MZ[mz_order[i+1] - 1] - mz;

    if (i == 0 && diff_mz > mz_sim * mz / 1e6) {
      mz_metabolites.push_back(Rcpp::as<std::string>(Compound_ID[idx]));
    } else if (i > 0) {
      double diff_up_mz = mz - MZ[mz_order[i-1] - 1];
      if (diff_mz > mz_sim * mz / 1e6 && diff_up_mz > mz_sim * mz / 1e6) {
        mz_metabolites.push_back(Rcpp::as<std::string>(Compound_ID[idx]));
      }
    }
  }

  // Merge RT and MZ metabolites
  std::set<std::string> unique_metabolites(rt_metabolites.begin(), rt_metabolites.end());
  unique_metabolites.insert(mz_metabolites.begin(), mz_metabolites.end());

  for (auto& id : unique_metabolites) {
    metabolites.push_back(id);
  }

  // Convert std::vector<std::string> to CharacterVector before returning
  CharacterVector result(metabolites.size());
  for (size_t i = 0; i < metabolites.size(); ++i) {
    result[i] = metabolites[i];
  }

  if (result.size() == 0) {
    stop("Not enough features isolated. Consider lowering RT and MZ thresholds.");
  }

  return result;
}
