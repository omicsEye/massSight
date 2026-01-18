"""mass_sight: massSight probabilistic matcher for cross-study LC–MS alignment."""

from .matcher import MassSightConfig, MatchResult, match_features

__version__ = "0.1.0"

__all__ = [
    "MassSightConfig",
    "MatchResult",
    "match_features",
    "__version__",
]
