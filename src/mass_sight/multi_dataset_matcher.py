"""
Multi-dataset matching engine for simultaneous alignment of multiple LC-MS datasets.

This module extends the pairwise matching approach to handle multiple datasets
simultaneously, building a unified feature map for meta-analysis.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Optional
from itertools import combinations
import networkx as nx
from dataclasses import dataclass

from .core import KernelDensityMatcher
from .match_engine import MatchResolver


@dataclass
class MultiMatchResult:
    """Results from multi-dataset matching."""
    consensus_features: pd.DataFrame  # Unified feature table
    pairwise_matches: Dict[Tuple[int, int], pd.DataFrame]  # All pairwise matches
    feature_graph: nx.Graph  # Network of feature relationships
    mean_confidence: float  # Overall match confidence score
    dataset_coverage: Dict[int, float]  # Per-dataset coverage


class MultiDatasetMatcher:
    """
    Simultaneous matching of multiple LC-MS datasets.
    
    Strategy:
    1. All-vs-all pairwise matching
    2. Build feature similarity graph
    3. Find connected components (consensus features)
    4. Calculate quality metrics
    5. Output unified feature table
    """
    
    def __init__(self, 
                 ppm_bandwidth: float = 5.0,
                 rt_bandwidth: float = 0.1,
                 ppm_weight: float = 1.0,
                 min_datasets_per_feature: int = 2):
        """
        Initialize multi-dataset matcher.
        
        Args:
            ppm_bandwidth: MZ bandwidth for kernel density (ppm)
            rt_bandwidth: RT bandwidth for kernel density (minutes)
            ppm_weight: Relative weight of MZ vs RT in distance calculation
            min_datasets_per_feature: Minimum datasets a feature must appear in
        """
        self.matcher = KernelDensityMatcher(ppm_bandwidth, rt_bandwidth, ppm_weight)
        self.resolver = MatchResolver()
        self.min_datasets_per_feature = min_datasets_per_feature
    
    def match_multiple_datasets(self, 
                              datasets: List[pd.DataFrame],
                              dataset_names: Optional[List[str]] = None,
                              use_gmm: bool = True) -> MultiMatchResult:
        """
        Match multiple datasets simultaneously.
        
        Args:
            datasets: List of DataFrames with MZ, RT, Intensity columns
            dataset_names: Optional names for datasets
            use_gmm: Whether to use GMM scoring
            
        Returns:
            MultiMatchResult with consensus features and metadata
        """
        n_datasets = len(datasets)
        if n_datasets < 2:
            raise ValueError("Need at least 2 datasets for matching")
        
        if dataset_names is None:
            dataset_names = [f"DS_{i}" for i in range(n_datasets)]
        
        print(f"Matching {n_datasets} datasets: {dataset_names}")
        
        # Step 1: All-vs-all pairwise matching
        print("Step 1: Pairwise matching...")
        pairwise_matches = {}
        all_confidences = []
        
        for i, j in combinations(range(n_datasets), 2):
            print(f"  Matching {dataset_names[i]} vs {dataset_names[j]}")
            
            matches = self.matcher.match_datasets(
                datasets[i], datasets[j],
                use_gmm=use_gmm,
                resolve_duplicates=False
            )
            
            if not matches.empty:
                # Add dataset indices
                matches = matches.copy()
                matches['dataset_i'] = i
                matches['dataset_j'] = j
                matches['pair_name'] = f"{dataset_names[i]}_vs_{dataset_names[j]}"
                
                # Resolve to one-to-one matches
                resolved = self.resolver.resolve_one_to_one(
                    matches, score_column="match_score", ascending=False
                )
                
                pairwise_matches[(i, j)] = resolved
                
                # Collect confidence scores
                if "match_confidence" in resolved.columns:
                    all_confidences.extend(resolved["match_confidence"].tolist())
                
                print(f"    Found {len(resolved)} matches")
            else:
                print(f"    No matches found")
        
        # Step 2: Build feature similarity graph
        print("Step 2: Building feature graph...")
        feature_graph = self._build_feature_graph(datasets, pairwise_matches, dataset_names)
        
        # Step 3: Find consensus features (connected components)
        print("Step 3: Finding consensus features...")
        consensus_features = self._find_consensus_features(
            datasets, feature_graph, dataset_names
        )
        
        # Step 4: Calculate quality metrics
        print("Step 4: Calculating quality metrics...")
        mean_confidence = np.mean(all_confidences) if all_confidences else 0.0
        
        # Step 5: Calculate coverage statistics
        coverage = self._calculate_coverage(datasets, consensus_features, dataset_names)
        
        print(f"Found {len(consensus_features)} consensus features")
        print(f"Mean confidence: {mean_confidence:.3f}")
        print(f"Dataset coverage: {coverage}")
        
        return MultiMatchResult(
            consensus_features=consensus_features,
            pairwise_matches=pairwise_matches,
            feature_graph=feature_graph,
            mean_confidence=mean_confidence,
            dataset_coverage=coverage
        )
    
    def _build_feature_graph(self, 
                           datasets: List[pd.DataFrame],
                           pairwise_matches: Dict[Tuple[int, int], pd.DataFrame],
                           dataset_names: List[str]) -> nx.Graph:
        """Build a graph where nodes are features and edges are matches."""
        G = nx.Graph()
        
        # Add all features as nodes
        for ds_idx, ds in enumerate(datasets):
            for feat_idx in ds.index:
                node_id = f"{ds_idx}_{feat_idx}"
                G.add_node(node_id, 
                          dataset=ds_idx,
                          dataset_name=dataset_names[ds_idx],
                          feature_idx=feat_idx,
                          mz=ds.loc[feat_idx, "MZ"],
                          rt=ds.loc[feat_idx, "RT"],
                          intensity=ds.loc[feat_idx, "Intensity"] if "Intensity" in ds.columns else 0)
        
        # Add edges for matches
        for (i, j), matches in pairwise_matches.items():
            for _, match in matches.iterrows():
                node_i = f"{i}_{match['id1']}"
                node_j = f"{j}_{match['id2']}"
                G.add_edge(node_i, node_j, 
                          match_score=match['match_score'],
                          ppm_diff=match.get('ppm_diff', 0),
                          rt_diff=match.get('rt_diff', 0),
                          confidence=match.get('match_confidence', 0))
        
        return G
    
    def _find_consensus_features(self, 
                               datasets: List[pd.DataFrame],
                               feature_graph: nx.Graph,
                               dataset_names: List[str]) -> pd.DataFrame:
        """Find consensus features from connected components."""
        consensus_features = []
        
        # Get connected components (groups of matched features)
        components = list(nx.connected_components(feature_graph))
        
        for comp_idx, component in enumerate(components):
            # Skip if not enough datasets represented
            datasets_in_component = set()
            for node in component:
                ds_idx = feature_graph.nodes[node]['dataset']
                datasets_in_component.add(ds_idx)
            
            if len(datasets_in_component) < self.min_datasets_per_feature:
                continue
            
            # Calculate consensus MZ, RT (weighted by intensity)
            total_intensity = 0
            weighted_mz = 0
            weighted_rt = 0
            intensities = {}
            confidences = []
            
            for node in component:
                node_data = feature_graph.nodes[node]
                intensity = node_data['intensity']
                
                weighted_mz += node_data['mz'] * intensity
                weighted_rt += node_data['rt'] * intensity
                total_intensity += intensity
                
                ds_name = node_data['dataset_name']
                intensities[f"intensity_{ds_name}"] = intensity
            
            # Get edge confidences
            for edge in feature_graph.subgraph(component).edges(data=True):
                if 'confidence' in edge[2]:
                    confidences.append(edge[2]['confidence'])
            
            if total_intensity > 0:
                consensus_mz = weighted_mz / total_intensity
                consensus_rt = weighted_rt / total_intensity
            else:
                # Fallback to simple average
                consensus_mz = np.mean([feature_graph.nodes[node]['mz'] for node in component])
                consensus_rt = np.mean([feature_graph.nodes[node]['rt'] for node in component])
            
            # Create consensus feature entry
            feature_data = {
                'consensus_id': f"CF_{comp_idx:05d}",
                'consensus_mz': consensus_mz,
                'consensus_rt': consensus_rt,
                'n_datasets': len(datasets_in_component),
                'n_features': len(component),
                'mean_confidence': np.mean(confidences) if confidences else 0.0,
                'datasets': ','.join([dataset_names[i] for i in sorted(datasets_in_component)])
            }
            
            # Add individual intensities
            for ds_name in dataset_names:
                feature_data[f"intensity_{ds_name}"] = intensities.get(f"intensity_{ds_name}", 0)
            
            consensus_features.append(feature_data)
        
        return pd.DataFrame(consensus_features)
    
    def _calculate_coverage(self, 
                          datasets: List[pd.DataFrame],
                          consensus_features: pd.DataFrame,
                          dataset_names: List[str]) -> Dict[int, float]:
        """Calculate what fraction of each dataset's features are in consensus."""
        coverage = {}
        
        for ds_idx, (ds, ds_name) in enumerate(zip(datasets, dataset_names)):
            if len(ds) == 0:
                coverage[ds_idx] = 0.0
                continue
            
            # Count features from this dataset in consensus
            n_in_consensus = 0
            for _, consensus_feat in consensus_features.iterrows():
                if consensus_feat[f"intensity_{ds_name}"] > 0:
                    n_in_consensus += 1
            
            coverage[ds_idx] = n_in_consensus / len(ds)
        
        return coverage


def example_usage():
    """Example of how to use the multi-dataset matcher."""
    # Mock data - replace with real datasets
    datasets = [
        pd.DataFrame({
            "MZ": [100.0, 200.0, 300.0],
            "RT": [1.0, 2.0, 3.0], 
            "Intensity": [1000, 2000, 3000]
        }),
        pd.DataFrame({
            "MZ": [100.1, 200.1, 350.0],
            "RT": [1.1, 2.1, 3.5],
            "Intensity": [1100, 2100, 3500]
        }),
        pd.DataFrame({
            "MZ": [99.9, 250.0, 300.1],
            "RT": [0.9, 2.5, 3.1],
            "Intensity": [900, 2500, 3100]
        })
    ]
    
    matcher = MultiDatasetMatcher()
    result = matcher.match_multiple_datasets(
        datasets, 
        dataset_names=["Study1", "Study2", "Study3"]
    )
    
    print("Consensus Features:")
    print(result.consensus_features)
    print(f"\nMean confidence: {result.mean_confidence:.3f}")
    print(f"Coverage: {result.dataset_coverage}")


if __name__ == "__main__":
    example_usage()