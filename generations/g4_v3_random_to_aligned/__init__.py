"""G4 v3: random-to-aligned collagen mechanics built from frozen G4 v2."""

from .model import (
    G4V3Config,
    build_random_network,
    make_random_void_spec,
    run_prescribed_contraction,
)

__all__ = [
    "G4V3Config",
    "build_random_network",
    "make_random_void_spec",
    "run_prescribed_contraction",
]
