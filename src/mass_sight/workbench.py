"""Metabolomics Workbench helpers for massSight reuse workflows."""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
from urllib import error, request

from .multistudy import StudySpec


WORKBENCH_BASE = "https://www.metabolomicsworkbench.org/rest/study/analysis_id/{analysis_id}"

_FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
_FEATURELIKE_RE = re.compile(
    r"^\s*[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?[_/@/][-+]?\d*\.?\d+(?:[eE][-+]?\d+)?\s*$"
)

_LOG_TOKENS = (
    "log2",
    "log 2",
    "log10",
    "log 10",
    "log-transformed",
    "log transformed",
    "ln-transformed",
    "ln transformed",
)
_NORM_TOKENS = (
    "normalized",
    "normalised",
    "batch corrected",
    "batch-corrected",
    "loess",
    "quantile normalization",
    "quantile-normalization",
    "median normalization",
    "median-normalized",
    "internal standard normalized",
    "qc normalized",
)
_SCALE_TOKENS = (
    "scaled units",
    "autoscaled",
    "auto-scaled",
    "z-score",
    "z score",
    "pareto",
    "unit variance",
    "uv scaling",
    "uv-scaled",
)
_RAWLIKE_TOKENS = (
    "peak area",
    "abundance",
    "intensity",
    "counts",
    "ms count",
    "ion abundance",
)


def _normalize_str(x: object) -> str:
    return str(x or "").strip()


def _normalize_key(key: str) -> str:
    return "".join(ch for ch in str(key).strip().lower() if ch.isalnum())


def _canonical_analysis_id(raw: object) -> str:
    s = _normalize_str(raw).upper()
    if not s:
        return ""
    if s.startswith("AN"):
        return s
    m = re.search(r"(\d+)", s)
    if not m:
        return s
    return f"AN{int(m.group(1)):06d}"


def _get_ci(row: dict, key: str) -> Any:
    if key in row:
        return row[key]
    want = _normalize_key(key)
    for k, v in row.items():
        if _normalize_key(str(k)) == want:
            return v
    return None


def _canonical_ion_mode(raw: object) -> str:
    s = _normalize_str(raw).lower()
    if "pos" in s:
        return "positive"
    if "neg" in s:
        return "negative"
    return "unknown"


def _canonical_chromatography(raw: object) -> str:
    s = _normalize_str(raw).lower()
    if not s:
        return "unknown"
    if "hilic" in s:
        return "hilic"
    if ("reverse" in s) or ("reversed" in s) or re.search(r"\brp\b", s) or ("c18" in s):
        return "reversed_phase"
    if "gc" in s:
        return "gc"
    return "other"


def _canonical_yes_no(raw: object) -> str:
    s = _normalize_str(raw).lower()
    if not s:
        return "unknown"
    if s.startswith("y"):
        return "yes"
    if s.startswith("n"):
        return "no"
    return "other"


def _canonical_rt_units(raw: object) -> str:
    s = _normalize_str(raw).lower()
    if not s:
        return "unknown"
    if "sec" in s:
        return "seconds"
    if "min" in s:
        return "minutes"
    return "other"


def parse_ms_results_file_meta(ms_results_file: object) -> Dict[str, str]:
    """Parse MW `MS_RESULTS_FILE` inline metadata fields."""
    text = _normalize_str(ms_results_file)
    out = {"raw": text, "units": "", "has_mz": "", "has_rt": "", "rt_units": ""}
    if not text:
        return out

    def _find_field(label: str, stop_labels: List[str]) -> str:
        m = re.search(rf"(?i)\b{re.escape(label)}\s*:\s*", text)
        if not m:
            return ""
        start = m.end()
        stops: List[int] = []
        for stop in stop_labels:
            m2 = re.search(rf"(?i)\b{re.escape(stop)}\s*:\s*", text[start:])
            if m2:
                stops.append(start + m2.start())
        end = min(stops) if stops else len(text)
        return text[start:end].strip()

    out["units"] = _find_field("UNITS", ["Has m/z", "Has RT", "RT units"])
    out["has_mz"] = _find_field("Has m/z", ["Has RT", "RT units"])
    out["has_rt"] = _find_field("Has RT", ["RT units"])
    out["rt_units"] = _find_field("RT units", [])
    return out


def _mz_semantics(has_mz_raw: object) -> str:
    s = _normalize_str(has_mz_raw).lower()
    if not s:
        return "unknown"
    if "neutral" in s:
        return "neutral_masses"
    if s.startswith("y"):
        return "yes"
    if s.startswith("n"):
        return "no"
    return "other"


def _dedupe_nonempty(values: Sequence[str]) -> List[str]:
    out: List[str] = []
    seen: set[str] = set()
    for v in values:
        s = _normalize_str(v)
        if not s:
            continue
        if s in seen:
            continue
        seen.add(s)
        out.append(s)
    return out


def _scale_hints_from_texts(*texts: str) -> Tuple[str, List[str]]:
    corpus = " ".join(_normalize_str(t) for t in texts if _normalize_str(t)).lower()
    if not corpus:
        return "unknown", []

    # Avoid false positives from acquisition parameters.
    corpus = corpus.replace("normalized collision energy", "collision energy")
    corpus = corpus.replace("normalised collision energy", "collision energy")

    hints: List[str] = []
    if any(tok in corpus for tok in _LOG_TOKENS):
        hints.append("log_transformed")
    if any(tok in corpus for tok in _SCALE_TOKENS):
        hints.append("scaled")
    if any(tok in corpus for tok in _NORM_TOKENS):
        hints.append("normalized")
    if "imput" in corpus:
        hints.append("imputed")
    if "batch" in corpus and ("correct" in corpus or "drift" in corpus):
        hints.append("batch_adjusted")

    if "log_transformed" in hints:
        scale = "log_transformed"
    elif "scaled" in hints:
        scale = "scaled"
    elif "normalized" in hints:
        scale = "normalized"
    elif any(tok in corpus for tok in _RAWLIKE_TOKENS):
        scale = "raw_like"
    else:
        scale = "unknown"

    return scale, hints


def extract_workbench_metadata(mwtab: Dict[str, Any], *, analysis_id_hint: str = "") -> Dict[str, Any]:
    """Extract normalized metadata fields used by reuse workflows."""
    wb = mwtab.get("METABOLOMICS WORKBENCH") or {}
    project = mwtab.get("PROJECT") or {}
    study = mwtab.get("STUDY") or {}
    subject = mwtab.get("SUBJECT") or {}
    collection = mwtab.get("COLLECTION") or {}
    chrom = mwtab.get("CHROMATOGRAPHY") or {}
    ms = mwtab.get("MS") or {}
    msmd = mwtab.get("MS_METABOLITE_DATA") or {}
    analysis = mwtab.get("ANALYSIS") or {}

    analysis_id = _canonical_analysis_id(
        _get_ci(wb, "ANALYSIS_ID") or _get_ci(mwtab, "ANALYSIS_ID") or analysis_id_hint
    )
    study_id = _normalize_str(_get_ci(wb, "STUDY_ID") or _get_ci(study, "STUDY_ID"))
    project_id = _normalize_str(_get_ci(wb, "PROJECT_ID") or _get_ci(project, "PROJECT_ID"))

    ion_mode_raw = _normalize_str(_get_ci(ms, "ION_MODE"))
    chrom_raw = _normalize_str(_get_ci(chrom, "CHROMATOGRAPHY_TYPE"))
    instrument_name = _normalize_str(_get_ci(ms, "INSTRUMENT_NAME"))
    instrument_type = _normalize_str(_get_ci(ms, "INSTRUMENT_TYPE"))
    ms_type = _normalize_str(_get_ci(ms, "MS_TYPE"))
    ms_comments = _normalize_str(_get_ci(ms, "MS_COMMENTS"))
    analysis_type = _normalize_str(_get_ci(analysis, "ANALYSIS_TYPE"))

    ms_results_meta = parse_ms_results_file_meta(_get_ci(ms, "MS_RESULTS_FILE"))
    msmd_units = _normalize_str(_get_ci(msmd, "Units") or _get_ci(msmd, "UNITS") or _get_ci(mwtab, "Units"))
    units_candidates = _dedupe_nonempty([ms_results_meta.get("units", ""), msmd_units])
    units_raw = " | ".join(units_candidates)
    value_scale, scale_hints = _scale_hints_from_texts(units_raw, ms_comments)

    mets = msmd.get("Metabolites") if isinstance(msmd, dict) else None
    data_rows = msmd.get("Data") if isinstance(msmd, dict) else None

    return {
        "analysis_id": analysis_id,
        "study_id": study_id,
        "project_id": project_id,
        "study_title": _normalize_str(_get_ci(study, "STUDY_TITLE")),
        "institute": _normalize_str(_get_ci(study, "INSTITUTE") or _get_ci(project, "INSTITUTE")),
        "subject_species": _normalize_str(_get_ci(subject, "SUBJECT_SPECIES")),
        "sample_type": _normalize_str(_get_ci(collection, "SAMPLE_TYPE") or _get_ci(subject, "SUBJECT_TYPE")),
        "analysis_type": analysis_type,
        "ion_mode_raw": ion_mode_raw,
        "ion_mode_norm": _canonical_ion_mode(ion_mode_raw),
        "chromatography_raw": chrom_raw,
        "chromatography_norm": _canonical_chromatography(chrom_raw),
        "instrument_name": instrument_name,
        "instrument_type": instrument_type,
        "ms_type": ms_type,
        "ms_results_file": _normalize_str(ms_results_meta.get("raw", "")),
        "ms_results_units": _normalize_str(ms_results_meta.get("units", "")),
        "ms_results_has_mz_raw": _normalize_str(ms_results_meta.get("has_mz", "")),
        "ms_results_has_mz_semantics": _mz_semantics(ms_results_meta.get("has_mz", "")),
        "ms_results_has_rt": _canonical_yes_no(ms_results_meta.get("has_rt", "")),
        "ms_results_rt_units_raw": _normalize_str(ms_results_meta.get("rt_units", "")),
        "ms_results_rt_units": _canonical_rt_units(ms_results_meta.get("rt_units", "")),
        "value_units_raw": units_raw,
        "value_scale": value_scale,
        "value_scale_hints": "|".join(scale_hints),
        "n_mwtab_metabolites": int(len(mets)) if isinstance(mets, list) else 0,
        "n_mwtab_data_rows": int(len(data_rows)) if isinstance(data_rows, list) else 0,
    }


def _parse_two_floats(token: str) -> Tuple[float, float]:
    s = _normalize_str(token).replace("/", "_")
    nums = _FLOAT_RE.findall(s)
    if not nums:
        return float("nan"), float("nan")
    a = float(nums[0])
    b = float(nums[1]) if len(nums) >= 2 else float("nan")
    return a, b


def estimate_untarg_feature_count(header: str, *, mz_min_plausible: float = 40.0) -> int:
    """Estimate number of feature columns from an untarg header line."""
    h = _normalize_str(header)
    if not h:
        return 0
    delim = "\t" if "\t" in h else ","
    parts = [p.strip() for p in h.split(delim)]
    count = 0
    for token in parts:
        if not _FEATURELIKE_RE.match(token):
            continue
        a, b = _parse_two_floats(token)
        if a <= 0:
            continue
        if max(a, b if b == b else 0.0) < float(mz_min_plausible):
            continue
        count += 1
    return int(count)


def _http_get_text(url: str, *, timeout_s: float, retries: int, delay_s: float) -> Optional[str]:
    headers = {"User-Agent": "mass-sight/0.1"}
    for attempt in range(int(retries)):
        if delay_s:
            time.sleep(float(delay_s))
        try:
            req = request.Request(url, headers=headers)
            with request.urlopen(req, timeout=float(timeout_s)) as resp:  # noqa: S310
                body = resp.read()
            return body.decode("utf-8", errors="replace")
        except error.HTTPError as exc:
            if exc.code == 404:
                return None
            if attempt == int(retries) - 1:
                return None
        except Exception:
            if attempt == int(retries) - 1:
                return None
    return None


def _http_get_json(url: str, *, timeout_s: float, retries: int, delay_s: float) -> Optional[Dict[str, Any]]:
    text = _http_get_text(url, timeout_s=timeout_s, retries=retries, delay_s=delay_s)
    if not text:
        return None
    try:
        obj = json.loads(text)
    except Exception:
        return None
    if isinstance(obj, dict):
        return obj
    return None


@dataclass(frozen=True)
class WorkbenchAnalysisBundle:
    analysis_id: str
    mwtab_path: Path
    untarg_path: Path
    metadata: Dict[str, Any]


def fetch_workbench_analysis(
    analysis_id: str,
    *,
    cache_root: Path,
    timeout_s: float = 60.0,
    retries: int = 3,
    delay_s: float = 0.25,
    refresh: bool = False,
) -> Optional[WorkbenchAnalysisBundle]:
    """Fetch one MW analysis (`mwtab` + `untarg_data`) with on-disk caching."""
    aid = _canonical_analysis_id(analysis_id)
    if not aid:
        return None

    cache_root = Path(cache_root)
    mwtab_dir = cache_root / "mwtab"
    untarg_dir = cache_root / "untarg"
    mwtab_dir.mkdir(parents=True, exist_ok=True)
    untarg_dir.mkdir(parents=True, exist_ok=True)

    mwtab_path = mwtab_dir / f"{aid}_mwtab.json"
    untarg_path = untarg_dir / f"{aid}_untarg.txt"

    mwtab_obj: Optional[Dict[str, Any]] = None
    if mwtab_path.exists() and mwtab_path.stat().st_size > 0 and not refresh:
        try:
            mwtab_obj = json.loads(mwtab_path.read_text(encoding="utf-8"))
        except Exception:
            mwtab_obj = None

    if mwtab_obj is None:
        mwtab_url = f"{WORKBENCH_BASE.format(analysis_id=aid)}/mwtab/json"
        mwtab_obj = _http_get_json(mwtab_url, timeout_s=timeout_s, retries=retries, delay_s=delay_s)
        if not mwtab_obj:
            return None
        mwtab_path.write_text(json.dumps(mwtab_obj), encoding="utf-8")

    untarg_text: Optional[str] = None
    if untarg_path.exists() and untarg_path.stat().st_size > 0 and not refresh:
        try:
            untarg_text = untarg_path.read_text(encoding="utf-8")
        except Exception:
            untarg_text = None

    if untarg_text is None:
        untarg_url = f"{WORKBENCH_BASE.format(analysis_id=aid)}/untarg_data/"
        untarg_text = _http_get_text(untarg_url, timeout_s=timeout_s, retries=retries, delay_s=delay_s)
        if not untarg_text:
            return None
        untarg_path.write_text(untarg_text, encoding="utf-8")

    header = untarg_text.splitlines()[0].strip() if untarg_text.splitlines() else ""
    n_untarg = estimate_untarg_feature_count(header)
    if n_untarg <= 0:
        return None

    metadata = extract_workbench_metadata(mwtab_obj, analysis_id_hint=aid)
    metadata.update(
        {
            "untarg_header_delimiter": "\\t" if "\t" in header else ",",
            "n_untarg_features": int(n_untarg),
            "mwtab_path": str(mwtab_path),
            "untarg_path": str(untarg_path),
        }
    )
    return WorkbenchAnalysisBundle(
        analysis_id=aid,
        mwtab_path=mwtab_path,
        untarg_path=untarg_path,
        metadata=metadata,
    )


def group_workbench_analyses(
    bundles: Sequence[WorkbenchAnalysisBundle],
    *,
    min_group_size: int = 2,
    allow_unknown_strata: bool = False,
) -> Tuple[Dict[str, List[WorkbenchAnalysisBundle]], List[Dict[str, str]]]:
    """Group analyses by compatibility strata for reuse (ion mode + chromatography)."""
    grouped: Dict[str, List[WorkbenchAnalysisBundle]] = {}
    dropped: List[Dict[str, str]] = []

    for bundle in bundles:
        meta = bundle.metadata
        ion = _normalize_str(meta.get("ion_mode_norm", "unknown")).lower() or "unknown"
        chrom = _normalize_str(meta.get("chromatography_norm", "unknown")).lower() or "unknown"
        if not allow_unknown_strata and (ion not in {"positive", "negative"} or chrom in {"unknown", "other"}):
            dropped.append(
                {
                    "analysis_id": bundle.analysis_id,
                    "reason": "unknown_or_unsupported_strata",
                    "ion_mode_norm": ion,
                    "chromatography_norm": chrom,
                }
            )
            continue
        key = f"{ion}__{chrom}"
        grouped.setdefault(key, []).append(bundle)

    for key in list(grouped.keys()):
        if len(grouped[key]) >= int(min_group_size):
            continue
        for bundle in grouped[key]:
            dropped.append(
                {
                    "analysis_id": bundle.analysis_id,
                    "reason": f"group_size_below_min_{int(min_group_size)}",
                    "group_key": key,
                }
            )
        grouped.pop(key, None)

    return grouped, dropped


def resolve_use_intensity_mode(
    metadata_rows: Sequence[Dict[str, Any]],
    *,
    mode: str = "off",
) -> Tuple[bool, str, str]:
    """
    Resolve whether intensity should be used for cross-study matching.

    Returns:
      (use_intensity, intensity_standardize, reason)
    """
    mode_norm = _normalize_str(mode).lower() or "off"
    if mode_norm not in {"off", "on", "auto"}:
        mode_norm = "off"

    if mode_norm == "off":
        return False, "none", "use_intensity_off"
    if mode_norm == "on":
        return True, "robust_zscore", "use_intensity_on"

    scales = {_normalize_str(m.get("value_scale", "")).lower() for m in metadata_rows}
    instr = {_normalize_str(m.get("instrument_type", "")).lower() for m in metadata_rows if _normalize_str(m.get("instrument_type", ""))}
    if any(s in {"log_transformed", "scaled"} for s in scales):
        return False, "none", "auto_disabled_scaled_or_log_values"
    if "normalized" in scales:
        return False, "none", "auto_disabled_normalized_values"
    if "unknown" in scales:
        return False, "none", "auto_disabled_unknown_value_scale"
    if len(instr) > 1:
        return False, "none", "auto_disabled_mixed_instrument_families"
    return True, "robust_zscore", "auto_enabled_raw_like_single_instrument_family"


def bundle_to_study_spec(bundle: WorkbenchAnalysisBundle, *, polarity_override: Optional[str] = None) -> StudySpec:
    """Convert a fetched Workbench analysis bundle into a `StudySpec`."""
    meta = bundle.metadata
    ion = _normalize_str(polarity_override or meta.get("ion_mode_norm", "positive")).lower()
    if ion not in {"positive", "negative"}:
        ion = "positive"
    return StudySpec(
        analysis_id=str(bundle.analysis_id),
        feature_format="mw_special_rows",
        matrix_path=Path(bundle.untarg_path),
        polarity=ion,
        chromatography=_normalize_str(meta.get("chromatography_norm") or meta.get("chromatography_raw")),
        metadata=dict(meta),
    )
