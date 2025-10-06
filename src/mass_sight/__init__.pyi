"""Type stubs for mass_sight package (ML matcher + MI fold-change)."""

from typing import Any, Optional
import pandas as pd

class MLMatchConfig:
    ppm: float
    rt: float
    model: str
    seed: int
    calibrate: bool
    neg_row_k: int
    neg_col_k: int
    neg_rand_k: int
    neg_expand: float
    neg_hard_k: int
    mz_col: str
    rt_col: str
    annotation_col: str
    compound_id_col: str
    intensity_col: str | None
    intensity_cols: list[str] | None
    intensity_regex: str | None

class MLMatchResult:
    candidates: pd.DataFrame
    top1: pd.DataFrame
    hungarian: pd.DataFrame

def match_mlnodrift(ds1: pd.DataFrame, ds2: pd.DataFrame, config: Optional[MLMatchConfig] = ...) -> MLMatchResult: ...

def mi_fold_change(
    ds1_expr: pd.DataFrame,
    ds1_meta: pd.DataFrame,
    ds2_expr: pd.DataFrame,
    ds2_meta: pd.DataFrame,
    candidates: pd.DataFrame,
    *,
    M: int = ..., seed: int = ...,
) -> pd.DataFrame: ...

__all__: list[str]
