from mass_sight.workbench import (
    extract_workbench_metadata,
    resolve_use_intensity_mode,
)


def test_extract_workbench_metadata_normalizes_core_fields_and_scale():
    mwtab = {
        "METABOLOMICS WORKBENCH": {"ANALYSIS_ID": "AN000123", "STUDY_ID": "ST000999"},
        "STUDY": {"STUDY_TITLE": "Example study", "INSTITUTE": "Example Institute"},
        "SUBJECT": {"SUBJECT_SPECIES": "Homo sapiens", "SUBJECT_TYPE": "Human"},
        "COLLECTION": {"SAMPLE_TYPE": "Blood plasma"},
        "CHROMATOGRAPHY": {"CHROMATOGRAPHY_TYPE": "Reversed phase LC"},
        "ANALYSIS": {"ANALYSIS_TYPE": "MS Untargeted"},
        "MS": {
            "INSTRUMENT_NAME": "Thermo Q Exactive",
            "INSTRUMENT_TYPE": "Orbitrap",
            "MS_TYPE": "ESI",
            "ION_MODE": "NEGATIVE",
            "MS_COMMENTS": "Data were LOESS normalized to QC samples and missing values were imputed.",
            "MS_RESULTS_FILE": "ST000999_AN000123_Results.txt UNITS:Normalized peak area Has m/z:Yes Has RT:Yes RT units:Minutes",
        },
        "MS_METABOLITE_DATA": {"Units": "Normalized peak area", "Metabolites": [{"a": 1}], "Data": [{"x": 1}]},
    }

    meta = extract_workbench_metadata(mwtab)
    assert meta["analysis_id"] == "AN000123"
    assert meta["study_id"] == "ST000999"
    assert meta["ion_mode_norm"] == "negative"
    assert meta["chromatography_norm"] == "reversed_phase"
    assert meta["ms_results_has_mz_semantics"] == "yes"
    assert meta["ms_results_has_rt"] == "yes"
    assert meta["ms_results_rt_units"] == "minutes"
    assert meta["value_scale"] == "normalized"
    assert "normalized" in str(meta["value_scale_hints"])
    assert "imputed" in str(meta["value_scale_hints"])


def test_extract_workbench_metadata_does_not_treat_normalized_collision_energy_as_data_normalization():
    mwtab = {
        "METABOLOMICS WORKBENCH": {"ANALYSIS_ID": "AN000124"},
        "CHROMATOGRAPHY": {"CHROMATOGRAPHY_TYPE": "HILIC"},
        "MS": {
            "ION_MODE": "POSITIVE",
            "INSTRUMENT_TYPE": "Orbitrap",
            "MS_COMMENTS": "The normalized collision energy was set to 20, 30 and 40 for fragmentation.",
            "MS_RESULTS_FILE": "ST_AN_Results.txt UNITS:Peak area Has m/z:Yes Has RT:Yes RT units:Seconds",
        },
    }
    meta = extract_workbench_metadata(mwtab)
    assert meta["ion_mode_norm"] == "positive"
    assert meta["chromatography_norm"] == "hilic"
    assert meta["value_scale"] == "raw_like"
    assert "normalized" not in str(meta["value_scale_hints"])


def test_resolve_use_intensity_mode_auto_is_conservative_with_scaled_data():
    rows = [
        {"value_scale": "raw_like", "instrument_type": "Orbitrap"},
        {"value_scale": "scaled", "instrument_type": "Orbitrap"},
    ]
    use_int, standardize, reason = resolve_use_intensity_mode(rows, mode="auto")
    assert use_int is False
    assert standardize == "none"
    assert reason == "auto_disabled_scaled_or_log_values"
