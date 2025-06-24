import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from scipy import stats
from sklearn.metrics import precision_recall_curve, roc_curve, auc

from .match_engine import MatchingResult


@dataclass
class EvaluationMetrics:
    """Comprehensive evaluation metrics for matching results."""
    
    # Basic metrics
    num_matches: int
    coverage_ds1: float
    coverage_ds2: float
    
    # Quality metrics
    precision: float
    recall: float
    f1_score: float
    
    # Error metrics
    mean_ppm_error: float
    median_ppm_error: float
    std_ppm_error: float
    ppm_95th_percentile: float
    
    mean_rt_error: float
    median_rt_error: float
    std_rt_error: float
    rt_95th_percentile: float
    
    # Distribution metrics
    ppm_outlier_rate: float
    rt_outlier_rate: float
    
    # Consistency metrics
    ppm_cv: float  # Coefficient of variation
    rt_cv: float
    
    # Performance metrics
    execution_time: float
    matches_per_second: float


class MatchingEvaluator:
    """Comprehensive evaluation of matching results."""
    
    def __init__(self, 
                 ppm_threshold: float = 5.0,
                 rt_threshold: float = 1.0,
                 outlier_factor: float = 3.0):
        self.ppm_threshold = ppm_threshold
        self.rt_threshold = rt_threshold
        self.outlier_factor = outlier_factor
    
    def evaluate_matches(self, 
                        result: MatchingResult,
                        ds1: pd.DataFrame,
                        ds2: pd.DataFrame,
                        ground_truth: Optional[pd.DataFrame] = None) -> EvaluationMetrics:
        """
        Comprehensive evaluation of matching results.
        
        Parameters:
        -----------
        result : MatchingResult
            Matching result to evaluate
        ds1, ds2 : pd.DataFrame
            Original datasets
        ground_truth : pd.DataFrame, optional
            Ground truth matches for precision/recall calculation
        """
        matches = result.matches
        
        if matches.empty:
            return self._empty_metrics(result.execution_time)
        
        # Basic metrics
        num_matches = len(matches)
        unique_ds1_matches = matches["id1"].nunique()
        unique_ds2_matches = matches["id2"].nunique()
        coverage_ds1 = unique_ds1_matches / len(ds1)
        coverage_ds2 = unique_ds2_matches / len(ds2)
        
        # Error metrics
        ppm_errors = matches["ppm_diff"].dropna()
        rt_errors = matches["rt_diff"].dropna()
        
        ppm_metrics = self._calculate_error_metrics(ppm_errors, self.ppm_threshold)
        rt_metrics = self._calculate_error_metrics(rt_errors, self.rt_threshold)
        
        # Precision/Recall if ground truth available
        if ground_truth is not None:
            precision, recall, f1 = self._calculate_precision_recall(matches, ground_truth)
        else:
            # Estimate precision based on error thresholds
            good_matches = matches[
                (ppm_errors < self.ppm_threshold) & 
                (rt_errors < self.rt_threshold)
            ]
            precision = len(good_matches) / num_matches if num_matches > 0 else 0
            recall = np.nan  # Cannot calculate without ground truth
            f1 = np.nan
        
        # Performance metrics
        matches_per_second = num_matches / result.execution_time if result.execution_time > 0 else 0
        
        return EvaluationMetrics(
            num_matches=num_matches,
            coverage_ds1=coverage_ds1,
            coverage_ds2=coverage_ds2,
            precision=precision,
            recall=recall,
            f1_score=f1,
            mean_ppm_error=ppm_metrics["mean"],
            median_ppm_error=ppm_metrics["median"],
            std_ppm_error=ppm_metrics["std"],
            ppm_95th_percentile=ppm_metrics["p95"],
            mean_rt_error=rt_metrics["mean"],
            median_rt_error=rt_metrics["median"],
            std_rt_error=rt_metrics["std"],
            rt_95th_percentile=rt_metrics["p95"],
            ppm_outlier_rate=ppm_metrics["outlier_rate"],
            rt_outlier_rate=rt_metrics["outlier_rate"],
            ppm_cv=ppm_metrics["cv"],
            rt_cv=rt_metrics["cv"],
            execution_time=result.execution_time,
            matches_per_second=matches_per_second
        )
    
    def compare_methods(self, 
                       results: Dict[str, MatchingResult],
                       ds1: pd.DataFrame,
                       ds2: pd.DataFrame,
                       ground_truth: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        Compare multiple matching methods.
        
        Returns a DataFrame with evaluation metrics for each method.
        """
        evaluations = {}
        
        for method, result in results.items():
            try:
                metrics = self.evaluate_matches(result, ds1, ds2, ground_truth)
                evaluations[method] = metrics
            except Exception as e:
                print(f"Warning: Failed to evaluate {method}: {e}")
                continue
        
        if not evaluations:
            return pd.DataFrame()
        
        # Convert to DataFrame
        comparison_data = []
        for method, metrics in evaluations.items():
            row = {
                "Method": method,
                "Matches": metrics.num_matches,
                "Coverage_DS1": round(metrics.coverage_ds1, 3),
                "Coverage_DS2": round(metrics.coverage_ds2, 3),
                "Precision": round(metrics.precision, 3),
                "F1_Score": round(metrics.f1_score, 3) if not np.isnan(metrics.f1_score) else np.nan,
                "Mean_PPM_Error": round(metrics.mean_ppm_error, 2),
                "Median_PPM_Error": round(metrics.median_ppm_error, 2),
                "PPM_95th_Pct": round(metrics.ppm_95th_percentile, 2),
                "Mean_RT_Error": round(metrics.mean_rt_error, 3),
                "Median_RT_Error": round(metrics.median_rt_error, 3),
                "RT_95th_Pct": round(metrics.rt_95th_percentile, 3),
                "PPM_Outlier_Rate": round(metrics.ppm_outlier_rate, 3),
                "RT_Outlier_Rate": round(metrics.rt_outlier_rate, 3),
                "Execution_Time": round(metrics.execution_time, 3),
                "Matches_Per_Sec": round(metrics.matches_per_second, 1)
            }
            comparison_data.append(row)
        
        return pd.DataFrame(comparison_data)
    
    def _calculate_error_metrics(self, errors: pd.Series, threshold: float) -> Dict[str, float]:
        """Calculate comprehensive error metrics."""
        if errors.empty:
            return {
                "mean": np.nan, "median": np.nan, "std": np.nan, "p95": np.nan,
                "outlier_rate": np.nan, "cv": np.nan
            }
        
        mean_error = errors.mean()
        median_error = errors.median()
        std_error = errors.std()
        p95_error = errors.quantile(0.95)
        
        # Outlier detection using IQR method
        q1, q3 = errors.quantile([0.25, 0.75])
        iqr = q3 - q1
        outlier_threshold = q3 + self.outlier_factor * iqr
        outlier_rate = (errors > outlier_threshold).mean()
        
        # Coefficient of variation
        cv = std_error / mean_error if mean_error > 0 else np.nan
        
        return {
            "mean": mean_error,
            "median": median_error,
            "std": std_error,
            "p95": p95_error,
            "outlier_rate": outlier_rate,
            "cv": cv
        }
    
    def _calculate_precision_recall(self, 
                                   matches: pd.DataFrame, 
                                   ground_truth: pd.DataFrame) -> Tuple[float, float, float]:
        """Calculate precision and recall against ground truth."""
        if ground_truth.empty:
            return np.nan, np.nan, np.nan
        
        # Convert matches to set of (id1, id2) tuples
        predicted_pairs = set(zip(matches["id1"], matches["id2"]))
        true_pairs = set(zip(ground_truth["id1"], ground_truth["id2"]))
        
        # Calculate metrics
        true_positives = len(predicted_pairs & true_pairs)
        false_positives = len(predicted_pairs - true_pairs)
        false_negatives = len(true_pairs - predicted_pairs)
        
        precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
        recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        
        return precision, recall, f1
    
    def _empty_metrics(self, execution_time: float) -> EvaluationMetrics:
        """Return metrics for empty results."""
        return EvaluationMetrics(
            num_matches=0,
            coverage_ds1=0.0,
            coverage_ds2=0.0,
            precision=0.0,
            recall=np.nan,
            f1_score=np.nan,
            mean_ppm_error=np.nan,
            median_ppm_error=np.nan,
            std_ppm_error=np.nan,
            ppm_95th_percentile=np.nan,
            mean_rt_error=np.nan,
            median_rt_error=np.nan,
            std_rt_error=np.nan,
            rt_95th_percentile=np.nan,
            ppm_outlier_rate=np.nan,
            rt_outlier_rate=np.nan,
            ppm_cv=np.nan,
            rt_cv=np.nan,
            execution_time=execution_time,
            matches_per_second=0.0
        )


class MatchingVisualizer:
    """Visualization tools for matching results."""
    
    def __init__(self, figsize: Tuple[int, int] = (12, 8)):
        self.figsize = figsize
        plt.style.use('default')
    
    def plot_error_distributions(self, 
                                results: Dict[str, MatchingResult],
                                save_path: Optional[str] = None) -> None:
        """Plot PPM and RT error distributions for all methods."""
        n_methods = len(results)
        fig, axes = plt.subplots(2, n_methods, figsize=(4*n_methods, 8))
        
        if n_methods == 1:
            axes = axes.reshape(2, 1)
        
        for i, (method, result) in enumerate(results.items()):
            matches = result.matches
            
            if matches.empty:
                axes[0, i].text(0.5, 0.5, 'No matches', ha='center', va='center')
                axes[1, i].text(0.5, 0.5, 'No matches', ha='center', va='center')
                continue
            
            # PPM error distribution
            if "ppm_diff" in matches.columns:
                ppm_errors = matches["ppm_diff"].dropna()
                axes[0, i].hist(ppm_errors, bins=30, alpha=0.7, edgecolor='black')
                axes[0, i].axvline(ppm_errors.median(), color='red', linestyle='--', 
                                 label=f'Median: {ppm_errors.median():.2f}')
                axes[0, i].set_xlabel('PPM Error')
                axes[0, i].set_ylabel('Frequency')
                axes[0, i].set_title(f'{method}\nPPM Errors')
                axes[0, i].legend()
            
            # RT error distribution
            if "rt_diff" in matches.columns:
                rt_errors = matches["rt_diff"].dropna()
                axes[1, i].hist(rt_errors, bins=30, alpha=0.7, edgecolor='black')
                axes[1, i].axvline(rt_errors.median(), color='red', linestyle='--',
                                 label=f'Median: {rt_errors.median():.3f}')
                axes[1, i].set_xlabel('RT Error (min)')
                axes[1, i].set_ylabel('Frequency')
                axes[1, i].set_title(f'{method}\nRT Errors')
                axes[1, i].legend()
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def plot_match_scatter(self, 
                          result: MatchingResult,
                          ds1: pd.DataFrame,
                          ds2: pd.DataFrame,
                          title: Optional[str] = None,
                          save_path: Optional[str] = None) -> None:
        """Plot scatter plot of matches in MZ-RT space."""
        matches = result.matches
        
        if matches.empty:
            print("No matches to plot")
            return
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize)
        
        # Dataset 1
        ax1.scatter(ds1["RT"], ds1["MZ"], alpha=0.5, s=20, color='blue', label='All features')
        
        if "id1" in matches.columns:
            matched_ds1 = ds1.iloc[matches["id1"]]
            ax1.scatter(matched_ds1["RT"], matched_ds1["MZ"], 
                       alpha=0.8, s=40, color='red', label='Matched features')
        
        ax1.set_xlabel('Retention Time (min)')
        ax1.set_ylabel('m/z')
        ax1.set_title('Dataset 1')
        ax1.legend()
        
        # Dataset 2
        ax2.scatter(ds2["RT"], ds2["MZ"], alpha=0.5, s=20, color='blue', label='All features')
        
        if "id2" in matches.columns:
            matched_ds2 = ds2.iloc[matches["id2"]]
            ax2.scatter(matched_ds2["RT"], matched_ds2["MZ"], 
                       alpha=0.8, s=40, color='red', label='Matched features')
        
        ax2.set_xlabel('Retention Time (min)')
        ax2.set_ylabel('m/z')
        ax2.set_title('Dataset 2')
        ax2.legend()
        
        if title:
            fig.suptitle(title)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def plot_method_comparison(self, 
                              comparison_df: pd.DataFrame,
                              metrics: List[str] = None,
                              save_path: Optional[str] = None) -> None:
        """Plot comparison of multiple methods."""
        if metrics is None:
            metrics = ["Precision", "Coverage_DS1", "Mean_PPM_Error", "Mean_RT_Error"]
        
        n_metrics = len(metrics)
        fig, axes = plt.subplots(1, n_metrics, figsize=(4*n_metrics, 6))
        
        if n_metrics == 1:
            axes = [axes]
        
        for i, metric in enumerate(metrics):
            if metric in comparison_df.columns:
                bars = axes[i].bar(comparison_df["Method"], comparison_df[metric])
                axes[i].set_title(metric.replace('_', ' '))
                axes[i].set_ylabel(metric)
                
                # Add value labels on bars
                for bar in bars:
                    height = bar.get_height()
                    if not np.isnan(height):
                        axes[i].text(bar.get_x() + bar.get_width()/2., height,
                                   f'{height:.3f}', ha='center', va='bottom')
                
                # Rotate x-axis labels if needed
                if len(comparison_df) > 3:
                    axes[i].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()
    
    def plot_drift_correction(self, 
                             matches_before: pd.DataFrame,
                             matches_after: pd.DataFrame,
                             title: str = "Drift Correction Effect",
                             save_path: Optional[str] = None) -> None:
        """Plot the effect of drift correction."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # PPM errors before/after
        if "ppm_diff" in matches_before.columns and "ppm_diff" in matches_after.columns:
            ax1.hist(matches_before["ppm_diff"], bins=30, alpha=0.7, 
                    label='Before correction', color='red')
            ax1.hist(matches_after["ppm_diff"], bins=30, alpha=0.7, 
                    label='After correction', color='blue')
            ax1.set_xlabel('PPM Error')
            ax1.set_ylabel('Frequency')
            ax1.set_title('PPM Error Distribution')
            ax1.legend()
        
        # RT errors before/after
        if "rt_diff" in matches_before.columns and "rt_diff" in matches_after.columns:
            ax2.hist(matches_before["rt_diff"], bins=30, alpha=0.7, 
                    label='Before correction', color='red')
            ax2.hist(matches_after["rt_diff"], bins=30, alpha=0.7, 
                    label='After correction', color='blue')
            ax2.set_xlabel('RT Error (min)')
            ax2.set_ylabel('Frequency')
            ax2.set_title('RT Error Distribution')
            ax2.legend()
        
        # Scatter plots
        if all(col in matches_before.columns for col in ["mz1", "mz2"]):
            ax3.scatter(matches_before["mz1"], matches_before["mz2"] - matches_before["mz1"], 
                       alpha=0.6, s=20, color='red', label='Before correction')
            ax3.set_xlabel('m/z Dataset 1')
            ax3.set_ylabel('m/z Difference (DS2 - DS1)')
            ax3.set_title('MZ Drift Pattern')
            ax3.legend()
        
        if all(col in matches_before.columns for col in ["rt1", "rt2"]):
            ax4.scatter(matches_before["rt1"], matches_before["rt2"] - matches_before["rt1"], 
                       alpha=0.6, s=20, color='red', label='Before correction')
            ax4.set_xlabel('RT Dataset 1 (min)')
            ax4.set_ylabel('RT Difference (DS2 - DS1)')
            ax4.set_title('RT Drift Pattern')
            ax4.legend()
        
        fig.suptitle(title)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()


# Convenience functions
def evaluate_and_compare(results: Dict[str, MatchingResult],
                        ds1: pd.DataFrame,
                        ds2: pd.DataFrame,
                        ground_truth: Optional[pd.DataFrame] = None,
                        show_plots: bool = True) -> pd.DataFrame:
    """
    Evaluate and compare multiple matching results with optional visualization.
    """
    evaluator = MatchingEvaluator()
    comparison_df = evaluator.compare_methods(results, ds1, ds2, ground_truth)
    
    print("="*80)
    print("MATCHING METHODS EVALUATION")
    print("="*80)
    print(comparison_df.to_string(index=False))
    print("="*80)
    
    if show_plots and len(results) > 0:
        visualizer = MatchingVisualizer()
        visualizer.plot_error_distributions(results)
        visualizer.plot_method_comparison(comparison_df)
    
    return comparison_df