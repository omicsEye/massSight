# core.py

import numpy as np
import pandas as pd
import polars as pl
from scipy.interpolate import interp1d
from sklearn.cluster import DBSCAN
from statsmodels.gam.api import GLMGam


def mad_based_outlier(points: int, thresh: float = 5.0) -> np.ndarray:
    """
    Identifies outliers based on the Median Absolute Deviation (MAD).
    """
    if len(points) == 0:
        return np.array([], dtype=bool)
    points = points[~np.isnan(points)]
    median = np.median(points)
    diff = np.sqrt((points - median) ** 2)
    med_abs_deviation = np.median(diff)
    if med_abs_deviation == 0:
        modified_z_score = np.zeros_like(diff)
    else:
        modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score > thresh


def get_vectors(df: pd.DataFrame, rt_sim: float = 0.01, mz_sim: float = 2) -> pd.Series:
    """
    This is a conceptual translation of the getVectors logic.
    A more direct translation might require a specific graph-based implementation.
    This version finds high-density regions by sorting and checking neighbors.
    """
    df = df.sort_values(by=["MZ", "RT"]).reset_index(drop=True)

    # Calculate differences with the next feature
    df["mz_diff"] = df["MZ"].diff().abs()
    df["rt_diff"] = df["RT"].diff().abs()

    # Identify features that are "close" to their next neighbor
    # and part of a potential cluster
    is_close = (df["mz_diff"] < mz_sim) & (df["rt_diff"] < rt_sim)

    # A simple heuristic: keep points that are part of a close pair
    # This captures the start of each "cluster"
    # To get all members, we'd need to propagate the inclusion
    close_indices = df.index[is_close].tolist()

    # Include the point after the close pair as well
    for i in df.index[is_close]:
        if i + 1 < len(df):
            close_indices.append(i + 1)

    return df.loc[list(set(close_indices)), "Compound_ID"]


def iso_dbscan(df: pd.DataFrame, eps: float = 0.1, min_samples: int = 2) -> pd.DataFrame:
    """
    Performs DBSCAN clustering to isolate features.
    """
    features = df[["RT", "MZ"]].values
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(features)
    # Return features that are not noise
    return df[db.labels_ != -1]


def align_isolated_compounds(ms1_iso: pd.DataFrame, ms2_iso: pd.DataFrame, rt_delta: float = 0.5, mz_delta: float = 15) -> pd.DataFrame:
    """
    Aligns isolated compounds from two datasets based on RT and MZ windows.
    """
    ms1_iso = ms1_iso.copy()
    ms2_iso = ms2_iso.copy()

    ms1_iso["rt_lower"] = ms1_iso["RT"] - rt_delta
    ms1_iso["rt_upper"] = ms1_iso["RT"] + rt_delta
    ms1_iso["mz_lower"] = ms1_iso["MZ"] - (ms1_iso["MZ"] * mz_delta / 1e6)
    ms1_iso["mz_upper"] = ms1_iso["MZ"] + (ms1_iso["MZ"] * mz_delta / 1e6)

    # Perform a cross join and then filter
    ms1_iso["key"] = 1
    ms2_iso["key"] = 1

    merged = pd.merge(ms1_iso, ms2_iso, on="key", suffixes=("_1", "_2")).drop(
        "key", axis=1
    )

    aligned = merged[
        (merged["RT_2"] >= merged["rt_lower"])
        & (merged["RT_2"] <= merged["rt_upper"])
        & (merged["MZ_2"] >= merged["mz_lower"])
        & (merged["MZ_2"] <= merged["mz_upper"])
    ]

    return aligned.drop(["rt_lower", "rt_upper", "mz_lower", "mz_upper"], axis=1)


def scale_smooth(query_values: np.ndarray, smooth_x: np.ndarray, smooth_y: np.ndarray) -> np.ndarray:
    """
    Apply smoothing correction.
    """
    interp_func = interp1d(smooth_x, smooth_y, kind="linear", fill_value="extrapolate")
    corrections = interp_func(query_values)
    return query_values - corrections


def smooth_drift(iso_matched_df: pd.DataFrame, df2_raw: pd.DataFrame, smooth_method: str = "gam") -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Models and corrects for systematic drift in RT and MZ.
    """
    results = iso_matched_df.sort_values(by="RT_1").copy()
    results["delta_RT"] = results["RT_2"] - results["RT_1"]

    # Smooth RT
    if smooth_method == "gam":
        gam_rt = GLMGam.from_formula("delta_RT ~ bs(RT_1, df=10)", data=results).fit()
        smooth_x_rt = results["RT_1"]
        smooth_y_rt = gam_rt.predict(pd.DataFrame({"RT_1": smooth_x_rt}))
    else:
        raise NotImplementedError(
            f"Smoothing method '{smooth_method}' not implemented."
        )

    scaled_rts = scale_smooth(df2_raw["RT"], smooth_x_rt + smooth_y_rt, smooth_y_rt)

    # Smooth MZ
    results["mass_error_ppm"] = (
        (results["MZ_2"] - results["MZ_1"]) / results["MZ_1"] * 1e6
    )
    if smooth_method == "gam":
        gam_mz = GLMGam.from_formula(
            "mass_error_ppm ~ bs(MZ_1, df=10)", data=results
        ).fit()
        smooth_x_mz = results["MZ_1"]
        smooth_y_mz = gam_mz.predict(pd.DataFrame({"MZ_1": smooth_x_mz}))
    else:
        raise NotImplementedError(
            f"Smoothing method '{smooth_method}' not implemented."
        )

    # Apply MZ corrections
    interp_func_mz = interp1d(
        smooth_x_mz, smooth_y_mz, kind="linear", fill_value="extrapolate"
    )
    mz_corrections = interp_func_mz(df2_raw["MZ"])
    scaled_mzs = df2_raw["MZ"] / (1 + mz_corrections / 1e6)

    # Create scaled_values dataframe
    scaled_values = pd.DataFrame(
        {
            "Compound_ID": df2_raw["Compound_ID"],
            "RT": scaled_rts,
            "MZ": scaled_mzs,
            "Intensity": df2_raw["Intensity"],
            "Metabolite": df2_raw.get(
                "Annotation_ID", pd.Series(index=df2_raw.index, dtype="object")
            ),
        }
    )

    return scaled_values, smooth_x_rt, smooth_y_rt, smooth_x_mz, smooth_y_mz


def find_all_matches(
    ref: pd.DataFrame, 
    query: pd.DataFrame, 
    rt_threshold: float, 
    mz_threshold: float, 
    alpha_rank: float = 0.5, 
    alpha_rt: float = 0.5, 
    alpha_mz: float = 0.5,
) -> pd.DataFrame:
    """
    Find all potential matches within thresholds and score them.
    """
    # Weight calculation
    denom = np.exp(alpha_rank) + np.exp(alpha_rt) + np.exp(alpha_mz)
    rank_weight, rt_weight, mz_weight = (
        np.exp(alpha_rank) / denom,
        np.exp(alpha_rt) / denom,
        np.exp(alpha_mz) / denom,
    )

    ref = ref.sort_values("RT").assign(RT_rank_1=(np.arange(len(ref))) / (len(ref) - 1))
    query = query.sort_values("RT_adj").assign(
        RT_rank_2=(np.arange(len(query))) / (len(query) - 1)
    )

    ref["key"] = 1
    query["key"] = 1

    potential_matches = pd.merge(ref, query, on="key", suffixes=("_1", "_2")).drop(
        "key", axis=1
    )

    potential_matches = potential_matches[
        (potential_matches["RT_adj_2"] >= potential_matches["RT_1"] - rt_threshold)
        & (potential_matches["RT_adj_2"] <= potential_matches["RT_1"] + rt_threshold)
        & (
            potential_matches["MZ_adj_2"]
            >= potential_matches["MZ_1"]
            - (potential_matches["MZ_1"] * mz_threshold / 1e6)
        )
        & (
            potential_matches["MZ_adj_2"]
            <= potential_matches["MZ_1"]
            + (potential_matches["MZ_1"] * mz_threshold / 1e6)
        )
    ].copy()

    # Score calculation
    rank_diff = abs(potential_matches["RT_rank_2"] - potential_matches["RT_rank_1"])
    rt_diff = abs(potential_matches["RT_adj_2"] - potential_matches["RT_1"])
    mz_diff = abs(
        (potential_matches["MZ_adj_2"] - potential_matches["MZ_1"])
        / potential_matches["MZ_1"]
        * 1e6
    )

    def min_max_scale(s):
        return (s - s.min()) / (s.max() - s.min())

    potential_matches["score"] = (
        rank_weight * min_max_scale(rank_diff)
        + rt_weight * min_max_scale(rt_diff)
        + mz_weight * min_max_scale(mz_diff)
    )
    return potential_matches


def mass_combine(
    ms1_df: pd.DataFrame | pl.DataFrame,
    ms2_df: pd.DataFrame | pl.DataFrame,
    optimize: bool = False,
    rt_delta: float = 0.5,
    mz_delta: float = 15,
    minimum_intensity: float = 10,
    iso_method: str = "manual",
    eps: float = 0.1,
    rt_iso_threshold: float = 0.01,
    mz_iso_threshold: float = 2,
    match_method: str = "unsupervised",
    smooth_method: str = "gam",
    alpha_rank: float = 0.5,
    alpha_rt: float = 0.5,
    alpha_mz: float = 0.5,
    id_name: str = "Compound_ID",
    rt_name: str = "RT",
    mz_name: str = "MZ",
    int_name: str = "Intensity",
    metab_name: str = "Annotation_ID",
) -> pd.DataFrame:
    """
    Python implementation of the R mass_combine function.
    
    Args:
        ms1_df: First mass spectrometry dataset (pandas or polars DataFrame)
        ms2_df: Second mass spectrometry dataset (pandas or polars DataFrame)
        ... (other parameters as documented)
    
    Returns:
        Combined DataFrame with matched features (always pandas DataFrame)
    """
    # Convert polars to pandas if needed
    if isinstance(ms1_df, pl.DataFrame):
        ms1_df = ms1_df.to_pandas()
    if isinstance(ms2_df, pl.DataFrame):
        ms2_df = ms2_df.to_pandas()
    
    # Standardize column names
    ms1 = ms1_df.rename(
        columns={
            id_name: "Compound_ID",
            rt_name: "RT",
            mz_name: "MZ",
            int_name: "Intensity",
            metab_name: "Annotation_ID",
        }
    )
    ms2 = ms2_df.rename(
        columns={
            id_name: "Compound_ID",
            rt_name: "RT",
            mz_name: "MZ",
            int_name: "Intensity",
            metab_name: "Annotation_ID",
        }
    )

    if optimize:
        raise NotImplementedError("Optimization path is not yet implemented.")

    # 1. Isolate compounds
    if match_method == "unsupervised":
        if iso_method == "manual":
            ref_iso_ids = get_vectors(
                ms1, rt_sim=rt_iso_threshold, mz_sim=mz_iso_threshold
            )
            query_iso_ids = get_vectors(
                ms2, rt_sim=rt_iso_threshold, mz_sim=mz_iso_threshold
            )
            ms1_iso = ms1[ms1["Compound_ID"].isin(ref_iso_ids)]
            ms2_iso = ms2[ms2["Compound_ID"].isin(query_iso_ids)]
        elif iso_method == "dbscan":
            ms1_iso = iso_dbscan(ms1, eps=eps)
            ms2_iso = iso_dbscan(ms2, eps=eps)
    else:  # supervised
        ms1_iso = ms1[ms1["Annotation_ID"].notna() & (ms1["Annotation_ID"] != "")]
        ms2_iso = ms2[ms2["Annotation_ID"].notna() & (ms2["Annotation_ID"] != "")]

    # 2. Align isolated compounds
    iso_matched = align_isolated_compounds(
        ms1_iso, ms2_iso, rt_delta=rt_delta, mz_delta=mz_delta
    )

    # 3. Smooth drift
    scaled_values, _, _, _, _ = smooth_drift(
        iso_matched, ms2, smooth_method=smooth_method
    )
    ms2_adj = ms2.copy()
    ms2_adj["RT_adj"] = scaled_values["RT"]
    ms2_adj["MZ_adj"] = scaled_values["MZ"]

    # 4. Final results
    # Use adjusted values for final matching
    potential_matches = find_all_matches(
        ms1,
        ms2_adj,
        rt_threshold=rt_delta,
        mz_threshold=mz_delta,
        alpha_rank=alpha_rank,
        alpha_rt=alpha_rt,
        alpha_mz=alpha_mz,
    )

    # Resolve to 1-to-1 matches by taking the best score
    best_matches = (
        potential_matches.sort_values("score")
        .drop_duplicates("Compound_ID_1")
        .drop_duplicates("Compound_ID_2")
    )

    # Rename columns for final output
    best_matches = best_matches.rename(
        columns={
            "Compound_ID_1": "Compound_ID_DS1",
            "RT_1": "RT_DS1",
            "MZ_1": "MZ_DS1",
            "Intensity_1": "Intensity_DS1",
            "Annotation_ID_1": "Metabolite_DS1",
            "Compound_ID_2": "Compound_ID_DS2",
            "RT_2": "RT_DS2",
            "MZ_2": "MZ_DS2",
            "Intensity_2": "Intensity_DS2",
            "Annotation_ID_2": "Metabolite_DS2",
        }
    )

    return best_matches
