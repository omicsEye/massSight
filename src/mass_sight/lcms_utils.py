from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd


PROTON_MASS = 1.007276466812

ADDUCT_MASS_POS = {
    "H": PROTON_MASS,
    "NH4": 18.033823,
    "NA": 22.989218,
    "K": 38.963158,
}

# Common negative-mode adduct masses (Da) used for discrete shift expansion.
# Values are monoisotopic masses and follow the convention:
#   observed_mz ≈ neutral_mass + adduct_mass
# where for deprotonation [M-H]- we treat the adduct mass as -H.
#
# Note: We keep the conservative NEG default in massSight (±H only) and expose these
# as an optional expansion (see `MassSightConfig.mz_shift_neg_adducts`).
ADDUCT_MASS_NEG = {
    # Chloride adduct: [M+Cl]-
    "CL": 34.968852682,
    # Formate adduct: [M+FA-H]-  (≈ formic_acid - H)
    "FA_H": 44.998202837,
    # Acetate adduct: [M+Ac-H]-  (≈ acetic_acid - H)
    "AC_H": 59.013852901,
}

def get_retention_order(ds: pd.DataFrame, rt_col: str = "RT") -> pd.Series:
    """Normalized retention order in [0,1], robust to monotone drift."""
    return ds[rt_col].rank(pct=True)


@dataclass(frozen=True)
class SchemaConfig:
    mz_col: str = "MZ"
    rt_col: str = "RT"
    intensity_col: Optional[str] = "Intensity"
    annotation_col: str = "Annotation_ID"
    compound_id_col: str = "Compound_ID"


def normalize_schema(
    ds: pd.DataFrame,
    cfg: SchemaConfig,
    *,
    mz_override: Optional[str] = None,
    rt_override: Optional[str] = None,
    annotation_override: Optional[str] = None,
    compound_override: Optional[str] = None,
) -> pd.DataFrame:
    """Return a copy with canonical columns for matching.

    Standardizes to:
      - `MZ`, `RT` (required)
      - `Intensity_log10` (derived when possible; else 0.0)
      - `ro` retention order (derived)
    """
    df = ds.copy()

    mz_src = mz_override or cfg.mz_col
    rt_src = rt_override or cfg.rt_col
    ann_src = annotation_override or cfg.annotation_col
    cmp_src = compound_override or cfg.compound_id_col

    if mz_src != "MZ" and mz_src in df.columns:
        df.rename(columns={mz_src: "MZ"}, inplace=True)
    if rt_src != "RT" and rt_src in df.columns:
        df.rename(columns={rt_src: "RT"}, inplace=True)
    if ann_src != "Annotation_ID" and ann_src in df.columns:
        df.rename(columns={ann_src: "Annotation_ID"}, inplace=True)
    if cmp_src != "Compound_ID" and cmp_src in df.columns:
        df.rename(columns={cmp_src: "Compound_ID"}, inplace=True)

    if "MZ" not in df.columns or "RT" not in df.columns:
        raise ValueError("Input dataset must contain m/z and RT columns (after overrides).")

    eps = 1e-9
    if cfg.intensity_col and cfg.intensity_col in df.columns:
        raw = pd.to_numeric(df[cfg.intensity_col], errors="coerce").astype(float)
        df["Intensity_log10"] = np.log10(np.clip(raw.to_numpy(dtype=float), 0.0, np.inf) + eps)
    else:
        df["Intensity_log10"] = 0.0

    if "ro" not in df.columns:
        df["ro"] = get_retention_order(df, rt_col="RT")

    return df
