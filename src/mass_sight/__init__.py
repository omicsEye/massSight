"""
mass_sight: ML-based feature matching + MI fold-change for LC–MS across datasets.
"""

from .ml_matcher import MLMatchConfig, MLMatchResult, match_ml
from .multi_match import MultiMatchConfig, MultiMatchResult, align_multi
from .mi_fold_change import mi_fold_change

__version__ = "0.1.0"

__all__ = [
    "MLMatchConfig",
    "MLMatchResult",
    "match_ml",
    "mi_fold_change",
    "MultiMatchConfig",
    "MultiMatchResult",
    "align_multi",
    "__version__",
]
