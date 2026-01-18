"""Type stubs for mass_sight package."""

from typing import Any, Dict, Optional
import pandas as pd

class MassSightConfig: ...

class MatchResult:
    candidates: pd.DataFrame
    top1: pd.DataFrame
    ot_summary: Dict[str, float]
    drift_params: Dict[str, Any]

def match_features(
    study_a: pd.DataFrame,
    study_b: pd.DataFrame,
    config: Optional[MassSightConfig] = ...,
    *,
    expr_a: Optional[object] = ...,
    expr_b: Optional[object] = ...,
) -> MatchResult: ...

__all__: list[str]
