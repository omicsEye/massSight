# massSight (development version)

* Changed `auto_combine()` to `mass_combine()`
* Fixed bug in `auto_combine()` that prevented combining `MSObjects` with no labeled metabolites
* Fixed bug in `auto_combine()` that prevented combining `MSObjects` with no intensity column
* Added linear smoothing for MZ drift
* Added linear smoothing for RT drift
* Added `pref` argument to `mass_combine()` to allow for pool reference based alignment

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
