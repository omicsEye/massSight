import numpy as np
import pandas as pd
import time
from dataclasses import dataclass
from .core import KernelDensityMatcher


@dataclass
class MatchingResult:
    """Container for matching results."""

    matches: pd.DataFrame
    method: str
    execution_time: float
    parameters: dict[str, object]
    metrics: dict[str, float]


class MatchingEngine:
    """Interface for kernel density-based feature matching."""

    def match_datasets(
        self, ds1: pd.DataFrame, ds2: pd.DataFrame, **kwargs: object
    ) -> MatchingResult:
        """
        Match two datasets using kernel density-based matching.

        Parameters:
        -----------
        ds1, ds2 : pd.DataFrame
            Input datasets with columns: MZ, RT, and optionally intensity columns
        **kwargs : additional parameters for the matcher

        Returns:
        --------
        MatchingResult containing matches and evaluation metrics
        """
        self._validate_input(ds1, ds2)

        # split kwargs: init only those for constructor, pass rest to match_datasets
        init_keys = ["ppm_bandwidth", "rt_bandwidth", "ppm_weight"]
        init_args = {k: kwargs[k] for k in init_keys if k in kwargs}
        run_args = {k: v for k, v in kwargs.items() if k not in init_keys}
        matcher = KernelDensityMatcher(**init_args)
        start_time = time.time()
        matches = matcher.match_datasets(ds1, ds2, **run_args)
        execution_time = time.time() - start_time
        metrics = self._calculate_metrics(matches, ds1, ds2)

        return MatchingResult(
            matches=matches,
            method="kernel_density",
            execution_time=execution_time,
            parameters=kwargs,
            metrics=metrics,
        )

    def _validate_input(self, ds1: pd.DataFrame, ds2: pd.DataFrame) -> None:
        """Validate input datasets."""
        required_cols = ["MZ", "RT"]

        for i, ds in enumerate([ds1, ds2], 1):
            if ds.empty:
                raise ValueError(f"Dataset {i} is empty")

            missing_cols = [col for col in required_cols if col not in ds.columns]
            if missing_cols:
                raise ValueError(f"Dataset {i} missing columns: {missing_cols}")

            if ds["MZ"].isna().any() or ds["RT"].isna().any():
                raise ValueError(f"Dataset {i} contains NaN values in MZ or RT columns")

    def _calculate_metrics(
        self, matches: pd.DataFrame, ds1: pd.DataFrame, ds2: pd.DataFrame
    ) -> dict[str, float]:
        """Calculate matching quality metrics."""
        if matches.empty:
            return {
                "num_matches": 0,
                "coverage_ds1": 0.0,
                "coverage_ds2": 0.0,
                "precision": 0.0,
                "mean_ppm_error": np.nan,
                "mean_rt_error": np.nan,
                "std_ppm_error": np.nan,
                "std_rt_error": np.nan,
            }

        # Basic counts
        num_matches = len(matches)
        unique_ds1_matches = matches["id1"].nunique()
        unique_ds2_matches = matches["id2"].nunique()

        # Coverage metrics
        coverage_ds1 = unique_ds1_matches / len(ds1) if len(ds1) > 0 else 0
        coverage_ds2 = unique_ds2_matches / len(ds2) if len(ds2) > 0 else 0

        # Quality metrics
        mean_ppm_error = (
            matches["ppm_diff"].mean() if "ppm_diff" in matches.columns else np.nan
        )
        std_ppm_error = (
            matches["ppm_diff"].std() if "ppm_diff" in matches.columns else np.nan
        )
        mean_rt_error = (
            matches["rt_diff"].mean() if "rt_diff" in matches.columns else np.nan
        )
        std_rt_error = (
            matches["rt_diff"].std() if "rt_diff" in matches.columns else np.nan
        )

        # Precision estimate (based on match quality)
        if "ppm_diff" in matches.columns and "rt_diff" in matches.columns:
            # Good matches have low ppm and rt errors
            good_matches = matches[
                (matches["ppm_diff"] < 5.0) & (matches["rt_diff"] < 1.0)
            ]
            precision = len(good_matches) / num_matches if num_matches > 0 else 0
        else:
            precision = np.nan

        return {
            "num_matches": num_matches,
            "unique_ds1_matches": unique_ds1_matches,
            "unique_ds2_matches": unique_ds2_matches,
            "coverage_ds1": coverage_ds1,
            "coverage_ds2": coverage_ds2,
            "precision": precision,
            "mean_ppm_error": mean_ppm_error,
            "mean_rt_error": mean_rt_error,
            "std_ppm_error": std_ppm_error,
            "std_rt_error": std_rt_error,
        }


class MatchResolver:
    """Resolve duplicate matches and optimize final matching."""

    @staticmethod
    def resolve_one_to_one(
        matches: pd.DataFrame, score_column: str = "match_score", ascending: bool = True
    ) -> pd.DataFrame:
        """
        Resolve matches to one-to-one mapping using score-based selection.

        Parameters:
        -----------
        matches : pd.DataFrame
            Matches with id1, id2, and score columns
        score_column : str
            Column to use for scoring (lower is better if ascending=True)
        ascending : bool
            Whether lower scores are better
        """
        if matches.empty:
            return matches

        # Sort by score
        sorted_matches = matches.sort_values(score_column, ascending=ascending)

        # Greedy selection for one-to-one mapping
        used_id1 = set()
        used_id2 = set()
        final_matches = []

        for _, row in sorted_matches.iterrows():
            if row["id1"] not in used_id1 and row["id2"] not in used_id2:
                final_matches.append(row)
                used_id1.add(row["id1"])
                used_id2.add(row["id2"])

        return pd.DataFrame(final_matches)

    @staticmethod
    def resolve_by_intensity(
        matches: pd.DataFrame,
        ds1: pd.DataFrame,
        ds2: pd.DataFrame,
        intensity_col: str = "Intensity",
    ) -> pd.DataFrame:
        """
        Resolve matches by preferring high-intensity features.
        """
        if (
            matches.empty
            or intensity_col not in ds1.columns
            or intensity_col not in ds2.columns
        ):
            return matches

        # Add intensity information
        matches_with_int = matches.merge(
            ds1[[intensity_col]]
            .reset_index()
            .rename(columns={"index": "id1", intensity_col: "int1"}),
            on="id1",
            how="left",
        ).merge(
            ds2[[intensity_col]]
            .reset_index()
            .rename(columns={"index": "id2", intensity_col: "int2"}),
            on="id2",
            how="left",
        )

        # Calculate combined intensity score
        matches_with_int["intensity_score"] = (
            matches_with_int["int1"] * matches_with_int["int2"]
        ) ** 0.5

        # Resolve using intensity score
        return MatchResolver.resolve_one_to_one(
            matches_with_int, score_column="intensity_score", ascending=False
        )


def load_dataset(
    file_path: str,
    mz_col: str = "MZ",
    rt_col: str = "RT",
    intensity_col: str | None = None,
    **read_csv_kwargs: object,
) -> pd.DataFrame:
    """
    Load dataset from file with standardized column names.

    Parameters:
    -----------
    file_path : str
        Path to input file
    mz_col, rt_col : str
        Column names for m/z and retention time
    intensity_col : str, optional
        Column name for intensity
    **read_csv_kwargs : additional arguments for pd.read_csv
    """
    df = pd.read_csv(file_path, **read_csv_kwargs)

    # Standardize column names
    column_mapping = {mz_col: "MZ", rt_col: "RT"}
    if intensity_col and intensity_col in df.columns:
        column_mapping[intensity_col] = "Intensity"

    df = df.rename(columns=column_mapping)

    # Remove rows with missing MZ or RT
    df = df.dropna(subset=["MZ", "RT"])

    return df


# Convenience function for direct usage
def quick_match(
    ds1: pd.DataFrame,
    ds2: pd.DataFrame,
    resolve_duplicates: bool = True,
    **kwargs: object,
) -> pd.DataFrame:
    """
    Quick matching function using kernel density-based matcher.
    """
    engine = MatchingEngine()
    result = engine.match_datasets(ds1, ds2, **kwargs)

    matches = result.matches

    if resolve_duplicates and not matches.empty:
        matches = MatchResolver.resolve_one_to_one(matches)

    return matches

    # Removed benchmark_methods: only one method is available.
