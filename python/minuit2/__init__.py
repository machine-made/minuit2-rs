"""minuit2 — pure-Rust Minuit-style optimizer with a measured iminuit-compatible subset.

The compiled extension lives in ``minuit2._minuit2``; this package re-exports
its public surface so users can ``from minuit2 import Minuit``.
"""

from ._minuit2 import FMin, MError, Minuit, Param, __version__

__all__ = ["Minuit", "FMin", "Param", "MError", "__version__"]
