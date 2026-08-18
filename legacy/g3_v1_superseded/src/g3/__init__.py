"""G3: minimal 2D emergent cell--collagen guidance model.

G3 is an active mechanism-validation model, not a realistic 3D tumor-migration model.
The public API intentionally lives outside the superseded ``cells`` coupling path.
"""

from .config import G3Config
from .state import ClutchState, G3Snapshot, ProtrusionState, RigidCellState

__all__ = [
    "G3Config",
    "RigidCellState",
    "ClutchState",
    "ProtrusionState",
    "G3Snapshot",
]
