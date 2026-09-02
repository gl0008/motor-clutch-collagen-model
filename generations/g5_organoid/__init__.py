"""Generation 5 - tumor organoid invading a collagen bead-spring network.

G5 scales the validated single-cell motor-clutch + G2 collagen engine up to a
multicellular organoid: N motor-clutch disks held together by a simplified
cell-cell adhesion potential, all coupled to one shared G2 fibre network.

See ``docs/G5_organoid_plan.md`` for the staged plan and provenance.
"""

from .model import (
    OrganoidConfig,
    hex_centers,
    make_organoid,
    cell_cell_forces,
    multi_cell_repulsion,
    organoid_active_forces,
    radial_alignment_profile,
    run_organoid_pull,
)

__all__ = [
    "OrganoidConfig",
    "hex_centers",
    "make_organoid",
    "cell_cell_forces",
    "multi_cell_repulsion",
    "organoid_active_forces",
    "radial_alignment_profile",
    "run_organoid_pull",
]
