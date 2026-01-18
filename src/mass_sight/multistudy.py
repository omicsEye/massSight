"""Multi-study utilities for massSight.

This module provides lightweight helpers to:
- load untargeted studies from a manifest (multiple input formats)
- run massSight pairwise between studies
- build cross-study clusters (hub-based or symmetric graph)

Normalization and downstream statistics are intentionally out-of-scope for the core
`massSight` library; this module focuses on matching/clustering only.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from .matcher import MassSightConfig, MatchResult, match_features


_FLOAT_PREFIX = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?")


@dataclass(frozen=True)
class StudySpec:
    analysis_id: str
    feature_format: str  # "mw_special_rows" | "id_encodes_mz_rt" | "separate_tables"
    matrix_path: Path
    features_path: Optional[Path] = None
    endpoints_path: Optional[Path] = None
    polarity: str = "positive"
    chromatography: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class LoadedStudy:
    analysis_id: str
    features: pd.DataFrame  # index: feature_id; columns: MZ, RT, Intensity
    expr: Optional[pd.DataFrame] = None  # index: sample_id; columns: feature_id
    spec: Optional[StudySpec] = None


@dataclass
class ClusterResult:
    strategy: str
    cluster_metadata: pd.DataFrame
    cluster_map: pd.DataFrame
    cluster_expr: Dict[str, pd.DataFrame]
    diagnostics: Dict[str, Any]


def _safe_yaml_load(path: Path) -> Any:
    try:
        import yaml  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "YAML manifest requested but PyYAML is not installed. "
            "Install `pyyaml` or provide a JSON manifest."
        ) from exc
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_manifest(path: Path) -> List[StudySpec]:
    """Load a multi-study manifest from YAML or JSON."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    suffix = path.suffix.lower().strip()
    if suffix in {".yaml", ".yml"}:
        obj = _safe_yaml_load(path)
    elif suffix == ".json":
        obj = json.loads(path.read_text(encoding="utf-8"))
    else:
        raise ValueError(f"Unsupported manifest type {suffix!r}; expected .yaml/.yml or .json")

    studies_obj = obj.get("studies") if isinstance(obj, dict) else obj
    if not isinstance(studies_obj, list):
        raise ValueError("Manifest must be a list of studies or a dict with key 'studies'.")

    specs: List[StudySpec] = []
    for item in studies_obj:
        if not isinstance(item, dict):
            raise ValueError("Each study entry must be a mapping/dict.")
        analysis_id = str(item.get("analysis_id") or item.get("study_id") or "").strip()
        if not analysis_id:
            raise ValueError("Study entry missing 'analysis_id'.")
        feature_format = str(item.get("feature_format") or "").strip()
        if not feature_format:
            raise ValueError(f"{analysis_id}: missing feature_format.")
        matrix_path = Path(str(item.get("matrix_path") or "").strip())
        if not matrix_path:
            raise ValueError(f"{analysis_id}: missing matrix_path.")

        features_path = item.get("features_path")
        endpoints_path = item.get("endpoints_path")
        spec = StudySpec(
            analysis_id=analysis_id,
            feature_format=feature_format,
            matrix_path=matrix_path,
            features_path=Path(str(features_path)) if features_path else None,
            endpoints_path=Path(str(endpoints_path)) if endpoints_path else None,
            polarity=str(item.get("polarity") or "positive"),
            chromatography=item.get("chromatography"),
            metadata=dict(item.get("metadata") or {}),
        )
        specs.append(spec)
    return specs


def _parse_feature_name(token: str) -> Tuple[float, float]:
    if "_" not in token:
        raise ValueError(f"Unexpected feature column name (missing '_'): {token!r}")
    first_raw, rest = token.split("_", 1)
    m = _FLOAT_PREFIX.match(rest)
    if m is None:
        raise ValueError(f"Could not parse numeric second token from {token!r}")
    return float(first_raw), float(m.group(0))


def load_mw_untarg_matrix(path: Path, *, na_as_zero: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load a Workbench-style untargeted TSV and return (features, expr).

    Supports two common layouts:
    - special rows `Samples == "m/z"` and `Samples == "rt"`
    - feature IDs that encode `MZ_RT` (order inferred heuristically)
    """
    df = pd.read_csv(path, sep="\t")
    if "Samples" not in df.columns:
        raise ValueError(f"Expected a 'Samples' column in {path}")

    feature_cols = [c for c in df.columns if c not in {"Samples", "group"}]
    if not feature_cols:
        raise ValueError(f"No feature columns found in {path}")

    samples = df["Samples"].astype(str)
    has_mz_row = bool(np.any(samples == "m/z"))
    has_rt_row = bool(np.any(samples == "rt"))

    if has_mz_row and has_rt_row:
        mz = df.loc[samples == "m/z", feature_cols].iloc[0].astype(float)
        rt = df.loc[samples == "rt", feature_cols].iloc[0].astype(float)
        expr = (
            df.loc[~samples.isin(["m/z", "rt"]), ["Samples", *feature_cols]]
            .set_index("Samples")
            .apply(pd.to_numeric, errors="coerce")
        )
        features = pd.DataFrame({"feature_id": feature_cols, "MZ": mz.to_numpy(float), "RT": rt.to_numpy(float)})
    else:
        first = np.zeros(len(feature_cols), dtype=float)
        second = np.zeros(len(feature_cols), dtype=float)
        for idx, name in enumerate(feature_cols):
            first[idx], second[idx] = _parse_feature_name(name)

        def score_assignment(mz_v: np.ndarray, rt_v: np.ndarray) -> float:
            mz_ok = np.mean((mz_v >= 50.0) & (mz_v <= 5000.0))
            rt_ok = np.mean((rt_v >= 0.0) & (rt_v <= 60.0))
            sep = float(np.median(mz_v) - np.median(rt_v))
            return float(mz_ok + rt_ok + sep / 1000.0)

        score_first_second = score_assignment(first, second)
        score_second_first = score_assignment(second, first)
        if abs(score_first_second - score_second_first) < 1e-6:
            raise ValueError(f"Ambiguous feature name parsing for {path}: cannot infer MZ vs RT.")
        if score_first_second > score_second_first:
            mz_v, rt_v = first, second
        else:
            mz_v, rt_v = second, first

        expr = df.loc[:, ["Samples", *feature_cols]].set_index("Samples").apply(pd.to_numeric, errors="coerce")
        features = pd.DataFrame({"feature_id": feature_cols, "MZ": mz_v, "RT": rt_v})

    expr.index = expr.index.astype(str)
    if na_as_zero:
        expr = expr.fillna(0.0)

    features["Intensity"] = expr.mean(axis=0).to_numpy(dtype=float)
    features["feature_id"] = features["feature_id"].astype(str)
    return features, expr


def load_study(spec: StudySpec, *, na_as_zero: bool = True) -> LoadedStudy:
    fmt = str(spec.feature_format).lower().strip()
    matrix_path = Path(spec.matrix_path)
    if not matrix_path.exists():
        raise FileNotFoundError(matrix_path)

    if fmt in {"mw_special_rows", "id_encodes_mz_rt"}:
        features, expr = load_mw_untarg_matrix(matrix_path, na_as_zero=na_as_zero)
    elif fmt == "separate_tables":
        if spec.features_path is None:
            raise ValueError(f"{spec.analysis_id}: features_path is required for feature_format=separate_tables.")
        features_path = Path(spec.features_path)
        if not features_path.exists():
            raise FileNotFoundError(features_path)
        expr = pd.read_csv(matrix_path, sep=None, engine="python")
        if expr.shape[1] < 2:
            raise ValueError(f"{spec.analysis_id}: matrix must have sample_id + feature columns.")
        sample_id = expr.iloc[:, 0].astype(str)
        expr = expr.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        expr.index = sample_id
        if na_as_zero:
            expr = expr.fillna(0.0)

        feat = pd.read_csv(features_path, sep=None, engine="python")
        required = {"feature_id", "mz", "rt"}
        cols = {c.lower().strip() for c in feat.columns}
        if not required.issubset(cols):
            raise ValueError(f"{spec.analysis_id}: features table must contain columns {sorted(required)}.")
        rename = {c: c.lower().strip() for c in feat.columns}
        feat = feat.rename(columns=rename)
        features = feat.loc[:, ["feature_id", "mz", "rt"]].copy()
        features = features.rename(columns={"mz": "MZ", "rt": "RT"})
        features["feature_id"] = features["feature_id"].astype(str)
        if "Intensity" not in features.columns:
            features["Intensity"] = expr.mean(axis=0).reindex(features["feature_id"]).to_numpy(dtype=float)
    else:
        raise ValueError(f"{spec.analysis_id}: unsupported feature_format={spec.feature_format!r}")

    features = features.copy()
    if "feature_id" not in features.columns:
        raise ValueError(f"{spec.analysis_id}: features must include a feature_id column.")
    features["feature_id"] = features["feature_id"].astype(str)
    features = features.set_index("feature_id", drop=False)
    features = features.loc[:, ["feature_id", "MZ", "RT", "Intensity"]].copy()

    expr = expr.loc[:, [c for c in expr.columns if str(c) in set(features.index)]].copy()
    expr.columns = [str(c) for c in expr.columns]
    return LoadedStudy(analysis_id=spec.analysis_id, features=features, expr=expr, spec=spec)


def select_hub(studies: Iterable[LoadedStudy], *, detect_rate_min: float = 0.05) -> str:
    """Choose a hub study for hub-based clustering."""
    best_id = None
    best_score = -float("inf")
    for study in studies:
        if study.expr is None:
            score = float(len(study.features))
        else:
            x = study.expr.to_numpy(dtype=float)
            detect_rate = np.mean(x > 0.0, axis=0)
            n_detect = int(np.sum(detect_rate >= float(detect_rate_min)))
            score = float(n_detect) * math.log1p(float(x.shape[0]))
        if score > best_score:
            best_score = score
            best_id = study.analysis_id
    if best_id is None:
        raise ValueError("No studies provided.")
    return best_id


def _prep_match_df(features: pd.DataFrame) -> Tuple[pd.DataFrame, np.ndarray]:
    feats = features[["feature_id", "MZ", "RT", "Intensity"]].copy().reset_index(drop=True)
    feature_ids = feats["feature_id"].to_numpy(dtype=str)
    ds = feats[["MZ", "RT", "Intensity"]].copy().reset_index(drop=True)
    return ds, feature_ids


def mutual_top1_map(study: LoadedStudy, hub: LoadedStudy, cfg: MassSightConfig) -> pd.DataFrame:
    ds_s, ids_s = _prep_match_df(study.features)
    ds_h, ids_h = _prep_match_df(hub.features)

    res_s2h = match_features(ds_s, ds_h, cfg)
    res_h2s = match_features(ds_h, ds_s, cfg)

    t_s2h = res_s2h.top1.copy()
    t_h2s = res_h2s.top1.copy()

    t_s2h = t_s2h[t_s2h["id2"].to_numpy(dtype=int) >= 0].copy()
    t_h2s = t_h2s[t_h2s["id2"].to_numpy(dtype=int) >= 0].copy()

    map_s2h = dict(zip(t_s2h["id1"].to_numpy(dtype=int), t_s2h["id2"].to_numpy(dtype=int)))
    map_h2s = dict(zip(t_h2s["id1"].to_numpy(dtype=int), t_h2s["id2"].to_numpy(dtype=int)))

    rows: List[Dict[str, object]] = []
    for s_idx, h_idx in map_s2h.items():
        back = map_h2s.get(int(h_idx))
        if back is None or int(back) != int(s_idx):
            continue
        rows.append({"study_idx": int(s_idx), "hub_idx": int(h_idx)})

    if not rows:
        return pd.DataFrame(columns=["cluster_id", "study_id", "feature_id", "support_s2h", "support_h2s"])

    mutual = pd.DataFrame(rows)
    mutual["cluster_id"] = mutual["hub_idx"].map(lambda i: str(ids_h[int(i)]))
    mutual["study_id"] = str(study.analysis_id)
    mutual["feature_id"] = mutual["study_idx"].map(lambda i: str(ids_s[int(i)]))

    support_s2h = (
        dict(zip(res_s2h.top1["id1"].to_numpy(dtype=int), res_s2h.top1["support"].to_numpy(dtype=float)))
        if "support" in res_s2h.top1.columns
        else {}
    )
    support_h2s = (
        dict(zip(res_h2s.top1["id1"].to_numpy(dtype=int), res_h2s.top1["support"].to_numpy(dtype=float)))
        if "support" in res_h2s.top1.columns
        else {}
    )
    mutual["support_s2h"] = mutual["study_idx"].map(support_s2h)
    mutual["support_h2s"] = mutual["hub_idx"].map(support_h2s)

    return mutual.loc[:, ["cluster_id", "study_id", "feature_id", "support_s2h", "support_h2s"]].copy()


def build_cluster_expr(hub: LoadedStudy, other: LoadedStudy, cluster_map: pd.DataFrame) -> pd.DataFrame:
    if other.expr is None:
        raise ValueError("build_cluster_expr requires an expr matrix for the non-hub study.")

    hub_clusters = hub.features["feature_id"].astype(str).tolist()
    # Important: distinguish within-study "0 intensity" (observed / LOD) from clusters that
    # are structurally absent in this study (no matched feature). We encode the latter as NA.
    out = pd.DataFrame(np.nan, index=other.expr.index, columns=hub_clusters, dtype=float)
    if cluster_map.empty:
        return out

    rows = cluster_map.loc[cluster_map["study_id"].astype(str) == str(other.analysis_id)]
    for _, row in rows.iterrows():
        cluster_id = str(row["cluster_id"])
        feat_id = str(row["feature_id"])
        if feat_id not in other.expr.columns or cluster_id not in out.columns:
            continue
        out[cluster_id] = other.expr[feat_id].to_numpy(dtype=float)
    return out


def cluster_hub_mutual_top1(
    studies: List[LoadedStudy],
    cfg: MassSightConfig,
    *,
    hub_id: Optional[str] = None,
    export_expr: bool = True,
    include_unmatched: bool = False,
) -> ClusterResult:
    if not studies:
        raise ValueError("No studies provided.")
    hub_id = hub_id or select_hub(studies)
    hub = next((s for s in studies if s.analysis_id == hub_id), None)
    if hub is None:
        raise ValueError(f"Hub {hub_id!r} not found in studies.")

    hub_clusters = hub.features["feature_id"].astype(str).tolist()
    hub_meta = hub.features.set_index(hub.features["feature_id"].astype(str))[["MZ", "RT", "Intensity"]].copy()
    hub_meta.index.name = "cluster_id"

    mappings: Dict[str, pd.DataFrame] = {}
    unmatched: Dict[str, List[Tuple[str, str]]] = {}
    cluster_meta = hub_meta.copy()

    cluster_map_rows: List[pd.DataFrame] = []
    for study in studies:
        if study.analysis_id == hub.analysis_id:
            continue
        mapping = mutual_top1_map(study, hub, cfg)
        mappings[study.analysis_id] = mapping
        if not mapping.empty:
            cluster_map_rows.append(mapping)

        if not include_unmatched:
            continue

        mapped_ids = set(mapping["feature_id"].astype(str)) if not mapping.empty else set()
        all_ids = study.features["feature_id"].astype(str).tolist()
        # Preserve feature-table order for deterministic output.
        unmapped_ids = [fid for fid in all_ids if fid not in mapped_ids]
        singleton_rows: List[Tuple[str, str]] = []
        for fid in unmapped_ids:
            cid = f"{study.analysis_id}__{fid}"
            singleton_rows.append((cid, fid))
            if cid in cluster_meta.index:
                continue
            row = study.features.loc[str(fid)]
            cluster_meta.loc[cid, ["MZ", "RT", "Intensity"]] = [float(row["MZ"]), float(row["RT"]), float(row["Intensity"])]
        unmatched[study.analysis_id] = singleton_rows

    # Build cluster_map. In full-join mode, include hub self-membership and singleton clusters.
    if include_unmatched:
        records: List[Dict[str, object]] = []
        for fid in hub_clusters:
            records.append(
                {
                    "cluster_id": str(fid),
                    "study_id": str(hub.analysis_id),
                    "feature_id": str(fid),
                    "support_s2h": float("nan"),
                    "support_h2s": float("nan"),
                }
            )
        if cluster_map_rows:
            records.extend(pd.concat(cluster_map_rows, ignore_index=True).to_dict(orient="records"))
        for sid, pairs in unmatched.items():
            for cid, fid in pairs:
                records.append(
                    {
                        "cluster_id": str(cid),
                        "study_id": str(sid),
                        "feature_id": str(fid),
                        "support_s2h": float("nan"),
                        "support_h2s": float("nan"),
                    }
                )
        cluster_map = pd.DataFrame.from_records(
            records,
            columns=["cluster_id", "study_id", "feature_id", "support_s2h", "support_h2s"],
        )
    else:
        cluster_map = pd.concat(cluster_map_rows, ignore_index=True) if cluster_map_rows else pd.DataFrame(
            columns=["cluster_id", "study_id", "feature_id", "support_s2h", "support_h2s"]
        )

    cluster_expr: Dict[str, pd.DataFrame] = {}
    if export_expr:
        all_clusters = cluster_meta.index.astype(str).tolist()

        def _empty_expr(expr: pd.DataFrame) -> pd.DataFrame:
            return pd.DataFrame(np.nan, index=expr.index, columns=all_clusters, dtype=float)

        if hub.expr is not None:
            if include_unmatched:
                mat = _empty_expr(hub.expr)
                for fid in hub_clusters:
                    if fid in hub.expr.columns:
                        mat[fid] = hub.expr[fid].to_numpy(dtype=float)
                cluster_expr[hub.analysis_id] = mat
            else:
                cluster_expr[hub.analysis_id] = hub.expr.loc[:, hub_clusters].copy()

        for study in studies:
            if study.analysis_id == hub.analysis_id or study.expr is None:
                continue
            mapping = mappings.get(study.analysis_id, pd.DataFrame())
            if include_unmatched:
                mat = _empty_expr(study.expr)
                if not mapping.empty:
                    for _, row in mapping.iterrows():
                        cid = str(row["cluster_id"])
                        fid = str(row["feature_id"])
                        if fid in study.expr.columns and cid in mat.columns:
                            mat[cid] = study.expr[fid].to_numpy(dtype=float)
                for cid, fid in unmatched.get(study.analysis_id, []):
                    if fid in study.expr.columns and cid in mat.columns:
                        mat[cid] = study.expr[fid].to_numpy(dtype=float)
                cluster_expr[study.analysis_id] = mat
            else:
                cluster_expr[study.analysis_id] = build_cluster_expr(hub, study, mapping)

    diagnostics: Dict[str, Any] = {
        "hub_id": hub.analysis_id,
        "n_studies": len(studies),
        "include_unmatched": bool(include_unmatched),
        "n_clusters": int(len(cluster_meta)),
        "n_singletons_by_study": {sid: int(len(v)) for sid, v in unmatched.items()} if include_unmatched else {},
    }
    return ClusterResult(
        strategy="hub_mutual_top1_full" if include_unmatched else "hub_mutual_top1",
        cluster_metadata=cluster_meta.reset_index(),
        cluster_map=cluster_map,
        cluster_expr=cluster_expr,
        diagnostics=diagnostics,
    )


def cluster_symmetric_mutual_graph(
    studies: List[LoadedStudy],
    cfg: MassSightConfig,
    *,
    export_expr: bool = True,
) -> ClusterResult:
    """Build clusters from pairwise mutual top-1 edges across all study pairs.

    This is intended as a sensitivity analysis to hub choice; it is not as conservative
    as hub-based mutual matching and uses a simple component + within-study tie-breaking.
    """
    if len(studies) < 2:
        raise ValueError("Need at least two studies.")

    node_index: Dict[Tuple[str, str], int] = {}
    nodes: List[Tuple[str, str]] = []
    for study in studies:
        for feat_id in study.features["feature_id"].astype(str).tolist():
            key = (study.analysis_id, feat_id)
            if key in node_index:
                continue
            node_index[key] = len(nodes)
            nodes.append(key)

    parent = np.arange(len(nodes), dtype=int)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    edges: List[Dict[str, object]] = []
    for i, s1 in enumerate(studies):
        for s2 in studies[i + 1 :]:
            study_a_df, ids_a = _prep_match_df(s1.features)
            study_b_df, ids_b = _prep_match_df(s2.features)
            res_12 = match_features(study_a_df, study_b_df, cfg)
            res_21 = match_features(study_b_df, study_a_df, cfg)

            t12 = res_12.top1.copy()
            t21 = res_21.top1.copy()
            t12 = t12[t12["id2"].to_numpy(int) >= 0].copy()
            t21 = t21[t21["id2"].to_numpy(int) >= 0].copy()
            map_12 = dict(zip(t12["id1"].to_numpy(int), t12["id2"].to_numpy(int)))
            map_21 = dict(zip(t21["id1"].to_numpy(int), t21["id2"].to_numpy(int)))

            support_12 = (
                dict(zip(res_12.top1["id1"].to_numpy(int), res_12.top1["support"].to_numpy(float)))
                if "support" in res_12.top1.columns
                else {}
            )
            support_21 = (
                dict(zip(res_21.top1["id1"].to_numpy(int), res_21.top1["support"].to_numpy(float)))
                if "support" in res_21.top1.columns
                else {}
            )

            for idx1, idx2 in map_12.items():
                back = map_21.get(int(idx2))
                if back is None or int(back) != int(idx1):
                    continue
                f1 = str(ids_a[int(idx1)])
                f2 = str(ids_b[int(idx2)])
                n1 = node_index[(s1.analysis_id, f1)]
                n2 = node_index[(s2.analysis_id, f2)]
                union(n1, n2)
                edges.append(
                    {
                        "study_a": s1.analysis_id,
                        "feature_a": f1,
                        "study_b": s2.analysis_id,
                        "feature_b": f2,
                        "support_ab": float(support_12.get(int(idx1), np.nan)),
                        "support_ba": float(support_21.get(int(idx2), np.nan)),
                    }
                )

    comp_to_nodes: Dict[int, List[int]] = {}
    for idx in range(len(nodes)):
        root = find(idx)
        comp_to_nodes.setdefault(root, []).append(idx)

    cluster_rows: List[Dict[str, object]] = []
    cluster_map_rows: List[Dict[str, object]] = []
    cluster_expr: Dict[str, pd.DataFrame] = {}

    cluster_id_counter = 0
    for comp_nodes in comp_to_nodes.values():
        if len(comp_nodes) == 1:
            continue
        cluster_id_counter += 1
        cluster_id = f"C{cluster_id_counter:05d}"

        chosen: Dict[str, Tuple[str, float]] = {}
        for node_id in comp_nodes:
            study_id, feat_id = nodes[node_id]
            score = 0.0
            for e in edges:
                if e["study_a"] == study_id and e["feature_a"] == feat_id:
                    score += float(e.get("support_ab") or 0.0)
                if e["study_b"] == study_id and e["feature_b"] == feat_id:
                    score += float(e.get("support_ba") or 0.0)
            prev = chosen.get(study_id)
            if prev is None or score > prev[1]:
                chosen[study_id] = (feat_id, score)

        if len(chosen) < 2:
            continue

        mz_vals = []
        rt_vals = []
        intensity_vals = []
        for study in studies:
            pick = chosen.get(study.analysis_id)
            if pick is None:
                continue
            feat_id = pick[0]
            row = study.features.loc[str(feat_id)]
            mz_vals.append(float(row["MZ"]))
            rt_vals.append(float(row["RT"]))
            intensity_vals.append(float(row["Intensity"]))
            cluster_map_rows.append({"cluster_id": cluster_id, "study_id": study.analysis_id, "feature_id": feat_id})

        cluster_rows.append(
            {
                "cluster_id": cluster_id,
                "mz_med": float(np.median(mz_vals)) if mz_vals else np.nan,
                "rt_med": float(np.median(rt_vals)) if rt_vals else np.nan,
                "intensity_med": float(np.median(intensity_vals)) if intensity_vals else np.nan,
                "n_studies": int(len(chosen)),
            }
        )

    cluster_metadata = pd.DataFrame(cluster_rows)
    cluster_map = pd.DataFrame(cluster_map_rows)

    if export_expr:
        for study in studies:
            if study.expr is None:
                continue
            rows = cluster_map.loc[cluster_map["study_id"].astype(str) == str(study.analysis_id)]
            # Encode clusters absent in this study as NA (structural missingness).
            mat = pd.DataFrame(np.nan, index=study.expr.index, columns=cluster_metadata["cluster_id"].astype(str), dtype=float)
            for _, row in rows.iterrows():
                feat_id = str(row["feature_id"])
                cid = str(row["cluster_id"])
                if feat_id in study.expr.columns and cid in mat.columns:
                    mat[cid] = study.expr[feat_id].to_numpy(dtype=float)
            cluster_expr[study.analysis_id] = mat

    diagnostics = {
        "n_edges_mutual": int(len(edges)),
        "n_clusters": int(len(cluster_metadata)),
        "n_studies": int(len(studies)),
    }
    return ClusterResult(
        strategy="symmetric_mutual_graph",
        cluster_metadata=cluster_metadata,
        cluster_map=cluster_map,
        cluster_expr=cluster_expr,
        diagnostics=diagnostics,
    )


def _mutual_top1_pairs(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    cfg: MassSightConfig,
) -> Tuple[pd.DataFrame, MatchResult, MatchResult]:
    """Return mutual top-1 pairs between study_a and study_b plus both match results."""
    res_12 = match_features(study_a, study_b, cfg)
    res_21 = match_features(study_b, study_a, cfg)

    t12 = res_12.top1.copy()
    t21 = res_21.top1.copy()

    t12 = t12[t12["id2"].to_numpy(dtype=int) >= 0].copy()
    t21 = t21[t21["id2"].to_numpy(dtype=int) >= 0].copy()

    map_12 = dict(zip(t12["id1"].to_numpy(dtype=int), t12["id2"].to_numpy(dtype=int)))
    map_21 = dict(zip(t21["id1"].to_numpy(dtype=int), t21["id2"].to_numpy(dtype=int)))

    rows: List[Dict[str, int]] = []
    for i1, i2 in map_12.items():
        back = map_21.get(int(i2))
        if back is None or int(back) != int(i1):
            continue
        rows.append({"id1": int(i1), "id2": int(i2)})

    return pd.DataFrame(rows), res_12, res_21


def _mutual_top1_pairs_one_direction(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    cfg: MassSightConfig,
) -> Tuple[pd.DataFrame, MatchResult]:
    """Return mutual top-1 pairs between study_a and study_b using a single call (study_a->study_b).

    We compute study_a's top-1 choice from `res.top1`, and study_b's top-1 choice by taking the
    argmax over study_a for each study_b column using the same OT coupling (`ot_weight`).
    The intersection enforces ≤1 feature per study per template cluster without running
    a full reverse match.
    """
    res = match_features(study_a, study_b, cfg)
    top = res.top1.copy()
    if top.empty:
        return pd.DataFrame(columns=["id1", "id2"]), res

    top = top[top["id2"].to_numpy(dtype=int) >= 0].copy()
    if top.empty:
        return pd.DataFrame(columns=["id1", "id2"]), res

    cand = res.candidates
    if cand.empty or "ot_weight" not in cand.columns:
        return pd.DataFrame(columns=["id1", "id2"]), res

    id1 = cand["id1"].to_numpy(dtype=int)
    id2 = cand["id2"].to_numpy(dtype=int)
    w = cand["ot_weight"].to_numpy(dtype=float)
    n2 = int(len(study_b))

    ok = (
        (id1 >= 0)
        & (id1 < int(len(study_a)))
        & (id2 >= 0)
        & (id2 < n2)
        & np.isfinite(w)
        & (w > 0)
    )
    id1 = id1[ok]
    id2 = id2[ok]
    w = w[ok]
    if id2.size == 0:
        return pd.DataFrame(columns=["id1", "id2"]), res

    col_sum = np.bincount(id2, weights=w, minlength=n2).astype(float)
    denom = col_sum[id2]
    denom = np.where(denom > 0, denom, 1.0)
    p_col = w / denom

    # For each study_b column, choose the study_a row with max p_col (ties broken deterministically).
    order = np.lexsort((-p_col, id2))
    id2_sorted = id2[order]
    id1_sorted = id1[order]
    uniq_j, first_idx = np.unique(id2_sorted, return_index=True)
    best_i_for_j = np.full(n2, -1, dtype=int)
    best_i_for_j[uniq_j] = id1_sorted[first_idx]

    id1_row = top["id1"].to_numpy(dtype=int)
    id2_row = top["id2"].to_numpy(dtype=int)
    ok2 = (id2_row >= 0) & (id2_row < n2)
    id1_row = id1_row[ok2]
    id2_row = id2_row[ok2]

    mutual = best_i_for_j[id2_row] == id1_row
    out = pd.DataFrame({"id1": id1_row[mutual].astype(int), "id2": id2_row[mutual].astype(int)})
    if not out.empty:
        out = out.sort_values(["id1", "id2"]).reset_index(drop=True)
    return out, res


def _template_features_from_hub(hub: LoadedStudy) -> Tuple[pd.DataFrame, np.ndarray]:
    feats = hub.features[["feature_id", "MZ", "RT", "Intensity"]].copy().reset_index(drop=True)
    ids = feats["feature_id"].to_numpy(dtype=str)
    ds = feats[["MZ", "RT", "Intensity"]].copy().reset_index(drop=True)
    return ds, ids


def _build_cluster_expr_from_map(
    *,
    cluster_ids: List[str],
    study: LoadedStudy,
    cluster_map: pd.DataFrame,
) -> pd.DataFrame:
    if study.expr is None:
        raise ValueError("_build_cluster_expr_from_map requires expr.")
    out = pd.DataFrame(np.nan, index=study.expr.index, columns=cluster_ids, dtype=float)
    if cluster_map.empty:
        return out
    rows = cluster_map.loc[cluster_map["study_id"].astype(str) == str(study.analysis_id)]
    for _, row in rows.iterrows():
        cid = str(row["cluster_id"])
        fid = str(row["feature_id"])
        if cid not in out.columns or fid not in study.expr.columns:
            continue
        out[cid] = study.expr[fid].to_numpy(dtype=float)
    return out


def cluster_hub_consensus_template_ot(
    studies: List[LoadedStudy],
    cfg: MassSightConfig,
    *,
    hub_id: Optional[str] = None,
    n_iters: int = 2,
    use_reverse_call: bool = True,
    export_expr: bool = True,
) -> ClusterResult:
    """
    Multi-study alignment via an iteratively updated hub template (V1).

    This is a lightweight “multi-marginal” extension of hub clustering:
    - Initialize a template from the hub study’s features.
    - Align every other study to the template via mutual top-1 matching.
    - Update template m/z using the matched features’ *effective* m/z (mz + inferred discrete shift).
    - Repeat for a small number of iterations.

    Notes:
    - Template m/z updates use `mz1_eff` from candidate edges (feature m/z adjusted into template space),
      which avoids the failure mode of averaging raw m/z across adduct states.
    - RT/intensity are not updated in V1 (hub values are retained).
    """
    if not studies:
        raise ValueError("No studies provided.")
    hub_id = hub_id or select_hub(studies)
    hub = next((s for s in studies if s.analysis_id == hub_id), None)
    if hub is None:
        raise ValueError(f"Hub {hub_id!r} not found in studies.")

    template_ds, template_ids = _template_features_from_hub(hub)
    template_ids_list = [str(x) for x in template_ids.tolist()]
    template_mz = template_ds["MZ"].to_numpy(dtype=float)

    iter_diag: List[Dict[str, object]] = []
    cluster_map_rows: List[pd.DataFrame] = []

    for it in range(max(int(n_iters), 1)):
        cluster_map_rows = []
        mz_eff_by_cluster: Dict[str, List[float]] = {cid: [] for cid in template_ids_list}

        for study in studies:
            if study.analysis_id == hub.analysis_id:
                continue
            ds_s, ids_s = _prep_match_df(study.features)

            if use_reverse_call:
                mutual_pairs, res_s2t, res_t2s = _mutual_top1_pairs(ds_s, template_ds, cfg)
            else:
                mutual_pairs, res_s2t = _mutual_top1_pairs_one_direction(ds_s, template_ds, cfg)
                res_t2s = None
            if mutual_pairs.empty:
                continue

            # Attach diagnostics from the forward candidate table (mz1_eff, mz_shift_da, support).
            cand = res_s2t.candidates
            cand_sub = cand.merge(mutual_pairs, on=["id1", "id2"], how="inner")
            if not cand_sub.empty and "loglik" in cand_sub.columns:
                cand_sub = (
                    cand_sub.sort_values(["id1", "id2", "loglik"], ascending=[True, True, False])
                    .drop_duplicates(["id1", "id2"], keep="first")
                    .copy()
                )

            # Support maps
            support_s2t = (
                dict(zip(res_s2t.top1["id1"].to_numpy(dtype=int), res_s2t.top1["support"].to_numpy(dtype=float)))
                if "support" in res_s2t.top1.columns
                else {}
            )
            support_t2s = {}
            if res_t2s is not None and "support" in res_t2s.top1.columns:
                support_t2s = dict(zip(res_t2s.top1["id1"].to_numpy(dtype=int), res_t2s.top1["support"].to_numpy(dtype=float)))

            out = mutual_pairs.copy()
            out["cluster_id"] = out["id2"].map(lambda j: str(template_ids[int(j)]))
            out["study_id"] = str(study.analysis_id)
            out["feature_id"] = out["id1"].map(lambda i: str(ids_s[int(i)]))
            out["support_s2t"] = out["id1"].map(support_s2t)
            out["support_t2s"] = out["id2"].map(support_t2s) if support_t2s else np.nan

            if not cand_sub.empty:
                out = out.merge(cand_sub.loc[:, ["id1", "id2", "mz1_eff", "mz_shift_da"]], on=["id1", "id2"], how="left")

            cluster_map_rows.append(out.loc[:, ["cluster_id", "study_id", "feature_id", "support_s2t", "support_t2s"]].copy())

            if "mz1_eff" in out.columns:
                for rec in out.to_dict("records"):
                    cid = str(rec.get("cluster_id", ""))
                    mz_eff = rec.get("mz1_eff", float("nan"))
                    if cid and np.isfinite(float(mz_eff)):
                        mz_eff_by_cluster.setdefault(cid, []).append(float(mz_eff))

        # Update template m/z using the effective m/z values in template space.
        n_updates = 0
        for j, cid in enumerate(template_ids_list):
            vals = mz_eff_by_cluster.get(cid) or []
            if not vals:
                continue
            # Always include the current template mass as a weak prior (stabilizes small groups).
            prior = float(template_mz[j]) if np.isfinite(template_mz[j]) else float("nan")
            use = [v for v in vals if np.isfinite(v) and v > 0]
            if np.isfinite(prior) and prior > 0:
                use.append(prior)
            if not use:
                continue
            template_mz[j] = float(np.median(np.asarray(use, dtype=float)))
            n_updates += 1

        template_ds = template_ds.copy()
        template_ds["MZ"] = template_mz.astype(float)

        iter_diag.append(
            {
                "iter": int(it),
                "n_matches": int(sum(len(df) for df in cluster_map_rows)),
                "n_mz_updates": int(n_updates),
            }
        )

    cluster_map = pd.concat(cluster_map_rows, ignore_index=True) if cluster_map_rows else pd.DataFrame(columns=["cluster_id", "study_id", "feature_id", "support_s2t", "support_t2s"])

    # Cluster metadata = final template values.
    # Note: cluster_id is the hub feature_id (template index) for this strategy.
    counts = cluster_map.groupby("cluster_id")["study_id"].nunique() if not cluster_map.empty else pd.Series(dtype=int)
    n_studies_by_cluster = counts.reindex(template_ids_list).fillna(0).astype(int).to_numpy(dtype=int)
    cluster_metadata = pd.DataFrame(
        {
            "cluster_id": template_ids_list,
            "MZ": template_ds["MZ"].to_numpy(dtype=float),
            "RT": template_ds["RT"].to_numpy(dtype=float),
            "Intensity": template_ds["Intensity"].to_numpy(dtype=float),
            "n_studies": (n_studies_by_cluster + 1).astype(int),  # include hub
        }
    )

    cluster_expr: Dict[str, pd.DataFrame] = {}
    if export_expr:
        for study in studies:
            if study.expr is None:
                continue
            if study.analysis_id == hub.analysis_id:
                # Hub clusters correspond 1–1 to hub features.
                mat = study.expr.reindex(columns=template_ids_list).copy()
                cluster_expr[study.analysis_id] = mat
            else:
                cluster_expr[study.analysis_id] = _build_cluster_expr_from_map(
                    cluster_ids=template_ids_list,
                    study=study,
                    cluster_map=cluster_map,
                )

    diagnostics = {
        "hub_id": str(hub.analysis_id),
        "n_iters": int(max(int(n_iters), 1)),
        "iter_diag": iter_diag,
        "n_clusters": int(len(cluster_metadata)),
        "n_studies": int(len(studies)),
        "n_cluster_assignments": int(len(cluster_map)),
    }
    return ClusterResult(
        strategy="hub_consensus_template_ot",
        cluster_metadata=cluster_metadata,
        cluster_map=cluster_map,
        cluster_expr=cluster_expr,
        diagnostics=diagnostics,
    )


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if values.size == 0:
        return float("nan")
    finite = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values = values[finite]
    weights = weights[finite]
    if values.size == 0:
        return float("nan")
    order = np.argsort(values, kind="mergesort")
    values = values[order]
    weights = weights[order]
    cum = np.cumsum(weights)
    if cum[-1] <= 0:
        return float("nan")
    cutoff = 0.5 * float(cum[-1])
    idx = int(np.searchsorted(cum, cutoff, side="left"))
    idx = min(max(idx, 0), int(values.size - 1))
    return float(values[idx])


def cluster_hub_barycenter_ot(
    studies: List[LoadedStudy],
    cfg: MassSightConfig,
    *,
    hub_id: Optional[str] = None,
    n_iters: int = 2,
    top_k_per_row: int = 3,
    min_p_row_ot: float = 0.0,
    prior_mult: float = 1.0,
    use_reverse_call: bool = False,
    include_unmatched: bool = False,
    export_expr: bool = True,
) -> ClusterResult:
    """
    Multi-study alignment via a hub-anchored OT-barycenter template (V1).

    This differs from `cluster_hub_consensus_template_ot` by updating the template m/z
    using *soft OT couplings* (`ot_weight`) rather than only mutual top‑1 matches.

    Template updates use `mz1_eff` from the study→template direction so discrete shift
    expansion is handled in template space.
    """
    if not studies:
        raise ValueError("No studies provided.")
    hub_id = hub_id or select_hub(studies)
    hub = next((s for s in studies if s.analysis_id == hub_id), None)
    if hub is None:
        raise ValueError(f"Hub {hub_id!r} not found in studies.")

    template_ds, template_ids = _template_features_from_hub(hub)
    template_ids_list = [str(x) for x in template_ids.tolist()]
    template_mz = template_ds["MZ"].to_numpy(dtype=float)

    n_template = int(len(template_ids_list))
    null_mass = float(getattr(cfg, "null_mass", 0.0) or 0.0)
    null_mass = min(max(null_mass, 0.0), 0.9)
    prior_weight = float(prior_mult) * (1.0 - null_mass) / float(max(n_template, 1))

    iter_diag: List[Dict[str, object]] = []
    for it in range(max(int(n_iters), 1)):
        mz_eff: Dict[str, List[float]] = {}
        mz_w: Dict[str, List[float]] = {}
        n_edges_used = 0

        for study in studies:
            if study.analysis_id == hub.analysis_id:
                continue
            ds_s, _ids_s = _prep_match_df(study.features)
            res = match_features(ds_s, template_ds, cfg)
            cand = res.candidates
            if cand.empty:
                continue
            if "mz1_eff" not in cand.columns or "ot_weight" not in cand.columns:
                continue

            use = cand
            if int(top_k_per_row) > 0:
                score_col = "p_row_ot" if "p_row_ot" in use.columns else "ot_weight"
                use = (
                    use.sort_values(["id1", score_col], ascending=[True, False])
                    .groupby("id1", sort=False)
                    .head(int(top_k_per_row))
                    .copy()
                )
            if float(min_p_row_ot) > 0 and "p_row_ot" in use.columns:
                use = use[use["p_row_ot"].to_numpy(dtype=float) >= float(min_p_row_ot)].copy()
            if use.empty:
                continue

            id2 = use["id2"].to_numpy(dtype=int)
            mz1_eff = use["mz1_eff"].to_numpy(dtype=float)
            w = use["ot_weight"].to_numpy(dtype=float)
            ok = np.isfinite(mz1_eff) & (mz1_eff > 0) & np.isfinite(w) & (w > 0) & np.isfinite(id2)
            id2 = id2[ok]
            mz1_eff = mz1_eff[ok]
            w = w[ok]
            if id2.size == 0:
                continue

            n_edges_used += int(id2.size)
            for j, mz_val, w_val in zip(id2.tolist(), mz1_eff.tolist(), w.tolist()):
                if j < 0 or j >= n_template:
                    continue
                cid = template_ids_list[int(j)]
                mz_eff.setdefault(cid, []).append(float(mz_val))
                mz_w.setdefault(cid, []).append(float(w_val))

        n_updates = 0
        for j, cid in enumerate(template_ids_list):
            vals = mz_eff.get(cid)
            if not vals:
                continue
            wts = mz_w.get(cid) or []
            if len(wts) != len(vals):
                continue
            values = np.asarray(vals, dtype=float)
            weights = np.asarray(wts, dtype=float)
            if prior_weight > 0 and np.isfinite(template_mz[j]) and template_mz[j] > 0:
                values = np.concatenate([values, np.array([float(template_mz[j])], dtype=float)])
                weights = np.concatenate([weights, np.array([float(prior_weight)], dtype=float)])

            new_mz = _weighted_median(values, weights)
            if np.isfinite(new_mz) and new_mz > 0:
                template_mz[j] = float(new_mz)
                n_updates += 1

        template_ds = template_ds.copy()
        template_ds["MZ"] = template_mz.astype(float)
        iter_diag.append(
            {
                "iter": int(it),
                "n_edges_used": int(n_edges_used),
                "n_mz_updates": int(n_updates),
                "prior_weight": float(prior_weight),
                "top_k_per_row": int(top_k_per_row),
                "min_p_row_ot": float(min_p_row_ot),
            }
        )

    # Finalize to deterministic clusters via mutual top‑1 against the final template.
    cluster_map_rows: List[pd.DataFrame] = []
    for study in studies:
        if study.analysis_id == hub.analysis_id:
            continue
        ds_s, ids_s = _prep_match_df(study.features)
        if use_reverse_call:
            mutual_pairs, res_s2t, res_t2s = _mutual_top1_pairs(ds_s, template_ds, cfg)
        else:
            mutual_pairs, res_s2t = _mutual_top1_pairs_one_direction(ds_s, template_ds, cfg)
            res_t2s = None
        if mutual_pairs.empty:
            continue

        support_s2t = (
            dict(zip(res_s2t.top1["id1"].to_numpy(dtype=int), res_s2t.top1["support"].to_numpy(dtype=float)))
            if "support" in res_s2t.top1.columns
            else {}
        )
        support_t2s = {}
        if res_t2s is not None and "support" in res_t2s.top1.columns:
            support_t2s = dict(zip(res_t2s.top1["id1"].to_numpy(dtype=int), res_t2s.top1["support"].to_numpy(dtype=float)))

        out = mutual_pairs.copy()
        out["cluster_id"] = out["id2"].map(lambda j: str(template_ids[int(j)]))
        out["study_id"] = str(study.analysis_id)
        out["feature_id"] = out["id1"].map(lambda i: str(ids_s[int(i)]))
        out["support_s2t"] = out["id1"].map(support_s2t)
        out["support_t2s"] = out["id2"].map(support_t2s) if support_t2s else np.nan
        cluster_map_rows.append(out.loc[:, ["cluster_id", "study_id", "feature_id", "support_s2t", "support_t2s"]].copy())

    cluster_map_matches = (
        pd.concat(cluster_map_rows, ignore_index=True)
        if cluster_map_rows
        else pd.DataFrame(columns=["cluster_id", "study_id", "feature_id", "support_s2t", "support_t2s"])
    )

    unmatched: Dict[str, List[Tuple[str, str]]] = {}
    singleton_meta_rows: List[Dict[str, object]] = []
    if include_unmatched:
        # Preserve feature-table order for deterministic output.
        for study in studies:
            if study.analysis_id == hub.analysis_id:
                continue
            mapped_ids = (
                set(
                    cluster_map_matches.loc[
                        cluster_map_matches["study_id"].astype(str) == str(study.analysis_id), "feature_id"
                    ]
                    .astype(str)
                    .tolist()
                )
                if not cluster_map_matches.empty
                else set()
            )
            all_ids = study.features["feature_id"].astype(str).tolist()
            unmapped_ids = [fid for fid in all_ids if fid not in mapped_ids]
            singleton_rows: List[Tuple[str, str]] = []
            for fid in unmapped_ids:
                cid = f"{study.analysis_id}__{fid}"
                singleton_rows.append((cid, fid))
                row = study.features.loc[str(fid)]
                singleton_meta_rows.append(
                    {
                        "cluster_id": str(cid),
                        "MZ": float(row["MZ"]),
                        "RT": float(row["RT"]),
                        "Intensity": float(row["Intensity"]),
                        "n_studies": 1,
                    }
                )
            unmatched[study.analysis_id] = singleton_rows

        template_cluster_ids = template_ids_list
        singleton_cluster_ids = [cid for pairs in unmatched.values() for cid, _ in pairs]
        all_cluster_ids = template_cluster_ids + singleton_cluster_ids

        records: List[Dict[str, object]] = []
        for cid in template_cluster_ids:
            records.append(
                {
                    "cluster_id": str(cid),
                    "study_id": str(hub.analysis_id),
                    "feature_id": str(cid),
                    "support_s2t": float("nan"),
                    "support_t2s": float("nan"),
                }
            )
        if not cluster_map_matches.empty:
            records.extend(cluster_map_matches.to_dict(orient="records"))
        for sid, pairs in unmatched.items():
            for cid, fid in pairs:
                records.append(
                    {
                        "cluster_id": str(cid),
                        "study_id": str(sid),
                        "feature_id": str(fid),
                        "support_s2t": float("nan"),
                        "support_t2s": float("nan"),
                    }
                )
        cluster_map = pd.DataFrame.from_records(
            records,
            columns=["cluster_id", "study_id", "feature_id", "support_s2t", "support_t2s"],
        )

        counts = (
            cluster_map.groupby("cluster_id")["study_id"].nunique()
            if not cluster_map.empty
            else pd.Series(dtype=int)
        )
        n_studies_template = counts.reindex(template_cluster_ids).fillna(0).astype(int).to_numpy(dtype=int)
        cluster_metadata = pd.DataFrame(
            {
                "cluster_id": template_cluster_ids,
                "MZ": template_ds["MZ"].to_numpy(dtype=float),
                "RT": template_ds["RT"].to_numpy(dtype=float),
                "Intensity": template_ds["Intensity"].to_numpy(dtype=float),
                "n_studies": n_studies_template.astype(int),
            }
        )
        if singleton_meta_rows:
            cluster_metadata = pd.concat(
                [cluster_metadata, pd.DataFrame.from_records(singleton_meta_rows)],
                ignore_index=True,
            )
        cluster_metadata = (
            cluster_metadata.set_index(cluster_metadata["cluster_id"].astype(str))
            .reindex([str(c) for c in all_cluster_ids])
            .reset_index(drop=True)
        )

        cluster_expr: Dict[str, pd.DataFrame] = {}
        if export_expr:
            def _empty_expr(expr: pd.DataFrame) -> pd.DataFrame:
                return pd.DataFrame(np.nan, index=expr.index, columns=all_cluster_ids, dtype=float)

            for study in studies:
                if study.expr is None:
                    continue
                mat = _empty_expr(study.expr)
                if study.analysis_id == hub.analysis_id:
                    for fid in template_cluster_ids:
                        if fid in study.expr.columns:
                            mat[fid] = study.expr[fid].to_numpy(dtype=float)
                    cluster_expr[study.analysis_id] = mat
                    continue

                mapping = cluster_map_matches.loc[
                    cluster_map_matches["study_id"].astype(str) == str(study.analysis_id)
                ] if not cluster_map_matches.empty else pd.DataFrame()
                if not mapping.empty:
                    for _, row in mapping.iterrows():
                        cid = str(row["cluster_id"])
                        fid = str(row["feature_id"])
                        if fid in study.expr.columns and cid in mat.columns:
                            mat[cid] = study.expr[fid].to_numpy(dtype=float)
                for cid, fid in unmatched.get(study.analysis_id, []):
                    if fid in study.expr.columns and cid in mat.columns:
                        mat[cid] = study.expr[fid].to_numpy(dtype=float)
                cluster_expr[study.analysis_id] = mat
    else:
        cluster_map = cluster_map_matches

        counts = (
            cluster_map.groupby("cluster_id")["study_id"].nunique()
            if not cluster_map.empty
            else pd.Series(dtype=int)
        )
        n_studies_by_cluster = counts.reindex(template_ids_list).fillna(0).astype(int).to_numpy(dtype=int)
        cluster_metadata = pd.DataFrame(
            {
                "cluster_id": template_ids_list,
                "MZ": template_ds["MZ"].to_numpy(dtype=float),
                "RT": template_ds["RT"].to_numpy(dtype=float),
                "Intensity": template_ds["Intensity"].to_numpy(dtype=float),
                "n_studies": (n_studies_by_cluster + 1).astype(int),  # include hub
            }
        )

        cluster_expr: Dict[str, pd.DataFrame] = {}
        if export_expr:
            for study in studies:
                if study.expr is None:
                    continue
                if study.analysis_id == hub.analysis_id:
                    cluster_expr[study.analysis_id] = study.expr.reindex(columns=template_ids_list).copy()
                else:
                    cluster_expr[study.analysis_id] = _build_cluster_expr_from_map(
                        cluster_ids=template_ids_list,
                        study=study,
                        cluster_map=cluster_map,
                    )

    diagnostics = {
        "hub_id": str(hub.analysis_id),
        "n_iters": int(max(int(n_iters), 1)),
        "iter_diag": iter_diag,
        "prior_mult": float(prior_mult),
        "prior_weight": float(prior_weight),
        "top_k_per_row": int(top_k_per_row),
        "min_p_row_ot": float(min_p_row_ot),
        "use_reverse_call": bool(use_reverse_call),
        "n_clusters": int(len(cluster_metadata)),
        "n_studies": int(len(studies)),
        "n_cluster_assignments": int(len(cluster_map)),
        "include_unmatched": bool(include_unmatched),
        "n_singletons_by_study": {sid: int(len(v)) for sid, v in unmatched.items()} if include_unmatched else {},
    }
    return ClusterResult(
        strategy="hub_barycenter_ot_full" if include_unmatched else "hub_barycenter_ot",
        cluster_metadata=cluster_metadata,
        cluster_map=cluster_map,
        cluster_expr=cluster_expr,
        diagnostics=diagnostics,
    )


__all__ = [
    "StudySpec",
    "LoadedStudy",
    "ClusterResult",
    "load_manifest",
    "load_study",
    "select_hub",
    "cluster_hub_mutual_top1",
    "cluster_hub_consensus_template_ot",
    "cluster_hub_barycenter_ot",
    "cluster_symmetric_mutual_graph",
]
