"""
Core matching algorithm for massSight - inter-lab metabolomics feature alignment.

This module provides the main KernelDensityMatcher class for aligning LC-MS features
across different laboratories/studies for meta-analysis.
"""

import warnings
import numpy as np
import pandas as pd
from sklearn.mixture import BayesianGaussianMixture

warnings.filterwarnings("ignore")


def ppm_diff(mz1: float, mz2: float) -> float:
    """Calculate PPM difference between two m/z values."""
    return abs(mz1 - mz2) / mz1 * 1e6


class KernelDensityMatcher:
    """
    Kernel density-based matching for inter-lab LC-MS feature alignment.
    
    This is the core algorithm for massSight, designed for meta-analysis of
    metabolomics data from different laboratories or studies.
    """

    def __init__(
        self,
        ppm_bandwidth: float = 5.0,  # Bandwidth in PPM
        rt_bandwidth: float = 0.1,  # Bandwidth in RT units (minutes)
        ppm_weight: float = 1.0,    # Relative weight of MZ vs RT
    ):
        """
        Initialize the matcher.
        
        Args:
            ppm_bandwidth: MZ bandwidth for kernel density estimation (ppm)
            rt_bandwidth: RT bandwidth for kernel density estimation (minutes)
            ppm_weight: Relative weight of MZ vs RT in distance calculation
        """
        self.ppm_bandwidth = ppm_bandwidth
        self.rt_bandwidth = rt_bandwidth
        self.ppm_weight = ppm_weight

    def prepare_features(self, ds: pd.DataFrame) -> np.ndarray:
        """Prepare MZ-RT features for matching."""
        features = np.column_stack([ds["MZ"].values, ds["RT"].values])
        return features

    def match_datasets(
        self,
        ds1: pd.DataFrame,
        ds2: pd.DataFrame,
        min_score: float = 0.1,
        use_gmm: bool = True,
        resolve_duplicates: bool = False,
    ) -> pd.DataFrame:
        """
        Match datasets using kernel density approach.
        
        Args:
            ds1: First dataset (reference)
            ds2: Second dataset (to be matched)
            min_score: Minimum overlap score threshold
            use_gmm: Whether to use GMM scoring for quality assessment
            resolve_duplicates: Whether to resolve to one-to-one matches
            
        Returns:
            DataFrame with matches including confidence estimates
        """
        # Prepare features for kernel density calculation
        features1 = self.prepare_features(ds1)
        features2 = self.prepare_features(ds2)

        # Vectorized pairwise distance calculation
        mz1 = features1[:, 0][:, np.newaxis]
        mz2 = features2[:, 0][np.newaxis, :]
        rt1 = features1[:, 1][:, np.newaxis]
        rt2 = features2[:, 1][np.newaxis, :]

        # Use scale-invariant PPM for m/z distance
        ppm_dist = np.abs(mz1 - mz2) / mz1 * 1e6
        rt_dist = np.abs(rt1 - rt2)

        # Combined distance using PPM and RT bandwidths
        combined_dist = np.sqrt(
            (self.ppm_weight * ppm_dist / self.ppm_bandwidth) ** 2
            + (rt_dist / self.rt_bandwidth) ** 2
        )
        overlap_score = np.exp(-(combined_dist**2) / 2)

        # Filter by min_score
        idx1, idx2 = np.where(overlap_score >= min_score)
        if len(idx1) == 0:
            return pd.DataFrame()

        # Gather matches efficiently
        mz1_vals = ds1.iloc[idx1]["MZ"].values
        mz2_vals = ds2.iloc[idx2]["MZ"].values
        rt1_vals = ds1.iloc[idx1]["RT"].values
        rt2_vals = ds2.iloc[idx2]["RT"].values
        ppm_diffs = np.abs(mz1_vals - mz2_vals) / mz1_vals * 1e6
        rt_diffs = np.abs(rt1_vals - rt2_vals)
        scores = overlap_score[idx1, idx2]

        matches_df = pd.DataFrame(
            {
                "id1": idx1,
                "id2": idx2,
                "mz1": mz1_vals,
                "mz2": mz2_vals,
                "rt1": rt1_vals,
                "rt2": rt2_vals,
                "overlap_score": scores,
                "ppm_diff": ppm_diffs,
                "rt_diff": rt_diffs,
            }
        )
        
        # Apply hard cutoffs for inter-lab matching (more lenient than intra-lab)
        matches_df = matches_df[
            (matches_df["ppm_diff"] <= 10) & (matches_df["rt_diff"] <= 1.0)
        ].reset_index(drop=True)

        # Initial match score based on kernel density overlap
        matches_df["match_score"] = matches_df["overlap_score"]

        # GMM-based scoring for quality assessment
        if use_gmm and len(matches_df) > 15:
            try:
                # Use 2D GMM on PPM and RT differences
                X = np.column_stack([
                    matches_df["ppm_diff"].values,
                    matches_df["rt_diff"].values,
                ])

                from sklearn.preprocessing import StandardScaler
                scaler = StandardScaler()
                X_scaled = scaler.fit_transform(X)

                # Bayesian GMM with 2 components for robustness
                gmm = BayesianGaussianMixture(
                    n_components=2, random_state=42, max_iter=400
                )
                gmm.fit(X_scaled)

                # Identify the "correct" cluster (smaller variance = tighter cluster)
                cluster_variances = [np.linalg.det(cov) for cov in gmm.covariances_]
                correct_idx = np.argmin(cluster_variances)

                # Get probabilities of belonging to the correct cluster
                probs = gmm.predict_proba(X_scaled)[:, correct_idx]
                matches_df["gmm_prob"] = np.clip(probs, 0, 1)

                # Use GMM probability as match score and confidence
                matches_df["match_score"] = matches_df["gmm_prob"]
                matches_df["match_confidence"] = matches_df["gmm_prob"]
                
                # Create categorical confidence levels
                matches_df["confidence_level"] = pd.cut(
                    matches_df["gmm_prob"], 
                    bins=[0, 0.5, 0.8, 0.95, 1.0],
                    labels=["Low", "Medium", "High", "Very High"],
                    include_lowest=True
                )

            except Exception as e:
                print(f"GMM scoring failed: {e}")
                # If GMM fails, use kernel density overlap as confidence
                matches_df["match_confidence"] = matches_df["overlap_score"]
                matches_df["confidence_level"] = pd.cut(
                    matches_df["overlap_score"], 
                    bins=[0, 0.3, 0.6, 0.8, 1.0],
                    labels=["Low", "Medium", "High", "Very High"],
                    include_lowest=True
                )
        else:
            # No GMM - use kernel density overlap as confidence
            matches_df["match_confidence"] = matches_df["overlap_score"]
            matches_df["confidence_level"] = pd.cut(
                matches_df["overlap_score"], 
                bins=[0, 0.3, 0.6, 0.8, 1.0],
                labels=["Low", "Medium", "High", "Very High"],
                include_lowest=True
            )

        # Sort by score (higher is better)
        result_df = matches_df.sort_values("match_score", ascending=False)
        
        # Optionally resolve to one-to-one matches
        if resolve_duplicates:
            # Keep best match per feature (highest score)
            result_df = result_df.drop_duplicates(subset=['id1'], keep='first')
            result_df = result_df.drop_duplicates(subset=['id2'], keep='first')

        return result_df


def match_datasets(
    ds1: pd.DataFrame,
    ds2: pd.DataFrame,
    ppm_bandwidth: float = 5.0,
    rt_bandwidth: float = 0.1,
    ppm_weight: float = 1.0,
    min_score: float = 0.1,
    use_gmm: bool = True,
    resolve_duplicates: bool = True,
) -> pd.DataFrame:
    """
    Convenience function for matching two datasets.
    
    Args:
        ds1: First dataset with MZ, RT columns
        ds2: Second dataset with MZ, RT columns
        ppm_bandwidth: MZ bandwidth for kernel density (ppm)
        rt_bandwidth: RT bandwidth for kernel density (minutes)
        ppm_weight: Relative weight of MZ vs RT
        min_score: Minimum overlap score threshold
        use_gmm: Whether to use GMM scoring
        resolve_duplicates: Whether to resolve to one-to-one matches
        
    Returns:
        DataFrame with matched features and confidence estimates
    """
    matcher = KernelDensityMatcher(ppm_bandwidth, rt_bandwidth, ppm_weight)
    return matcher.match_datasets(
        ds1, ds2, min_score, use_gmm, resolve_duplicates
    )