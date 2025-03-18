# massSight 0.3.0

* Changed `auto_combine()` to `mass_combine()`
* Fixed bug in `auto_combine()` that prevented combining `MSObjects` with no labeled metabolites
* Fixed bug in `auto_combine()` that prevented combining `MSObjects` with no intensity column
* Added linear smoothing for MZ drift
* Added linear smoothing for RT drift
* Added `pref` argument to `mass_combine()` to allow for pool reference based alignment
* Enhanced `get_unique_matches()` with additional metadata columns (MZ, RT, Metabolite)
* Added Shiny app interface via `run_massSight_app()`
* Improved Bayesian optimization in mass_combine
* Added comprehensive test suite for mass_combine functionality
* Code reformatting across multiple files for better readability and maintainability
* Enhanced `ml_match()` with advanced machine learning-based metabolite matching features
* Added improved feature engineering for more accurate metabolite matching
* Added semi-supervised learning capabilities to `ml_match()` for better performance with limited labeled data
* Added configurable options for intensity-based feature usage in matching algorithms

# massSight 0.2.2

* Working `auto_scale()` function for the scaling of metabolite features
* Fixed logging feature for `auto_combine()`

# massSight 0.2.1

* Fixed bug in `final_plots()` that prevented the display of all matched metabolites
* Removed pre-isolation matching from `auto_combine()`

# massSight 0.2.0

* Added more robust plotting features
* `auto_combine()` output can now be used as input for future `auto_combine()` calls
* Ported `get_vectors()` to C++

# massSight 0.1.0

* Initial GitHub release
