"""epochmatch: cross-match two catalogues at a common epoch."""

from .core import EpochMatchError, load_catalogue, match_tables, \
    propagate_positions

__version__ = "0.1.0"
__all__ = ["EpochMatchError", "load_catalogue", "match_tables",
           "propagate_positions", "__version__"]
