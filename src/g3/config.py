"""Configuration and SI-unit provenance for the active G3 model."""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class G3Config:
    """Resolved G3 parameters; all active calculations use metres, newtons, seconds."""

    # Geometry / elastic collagen: Saraswathibhatla et al. 2025 SI Table 2.
    cell_radius: float = 10.0e-6
    fibre_count: int = 8
    fibre_length: float = 40.0e-6
    bead_spacing: float = 1.0e-6
    kappa_s_f: float = 4.0e-3
    kappa_b_f: float = 8.27e-20
    theta0_f: float = 0.0
    bead_drag: float = 3.73e-5

    # Local representation of a point clutch force (D010).
    gaussian_sigma: float = 2.0e-6
    gaussian_support_sigma: float = 3.0
    capture_distance: float = 2.0e-6

    # Adebowale et al. 2021 SI Table 4, converted once to SI.
    n_clutches: int = 200
    n_motors: int = 200
    motor_force: float = 2.0e-12
    bell_force: float = 2.0e-12
    bind_rate: float = 0.2
    unbind_rate: float = 0.02
    clutch_stiffness: float = 5.0e-3
    unloaded_actin_speed: float = 24.0e-9

    # Carey-inspired coarse-grained protrusion feedback; provisional assumptions.
    n_sectors: int = 24
    n_active_protrusions: int = 2
    protrusion_lifetime: float = 120.0
    beta_geometry: float = 2.0
    beta_traction: float = 2.0
    sensing_sigma: float = 3.0e-6
    feedback_time: float = 30.0
    geometry_update_interval: float = 0.1

    # Overdamped integration and calibrated rigid-cell drag.
    dt: float = 0.005
    ecm_substeps: int = 2
    cell_drag: float = 0.3  # 300 nN s / um == 0.3 N s / m
    rotational_drag_factor: float = 1.0
    overlap_tolerance: float = 0.1e-6

    # Recording defaults. Long preregistered runs override duration at the CLI.
    duration_g3a: float = 15.0
    duration_g3b: float = 600.0
    duration_g3c: float = 600.0
    metrics_interval: float = 0.1
    frame_interval: float = 1.0

    @property
    def rotational_drag(self) -> float:
        return self.rotational_drag_factor * self.cell_drag * self.cell_radius**2

    @property
    def beads_per_fibre(self) -> int:
        return int(round(self.fibre_length / self.bead_spacing)) + 1

    def validate(self) -> None:
        positive = (
            "cell_radius", "fibre_length", "bead_spacing", "kappa_s_f", "bead_drag",
            "gaussian_sigma", "capture_distance", "bell_force", "clutch_stiffness",
            "dt", "cell_drag", "protrusion_lifetime", "feedback_time",
        )
        for name in positive:
            if getattr(self, name) <= 0.0:
                raise ValueError(f"{name} must be positive")
        if self.n_clutches < 1 or self.n_motors < 1:
            raise ValueError("n_clutches and n_motors must be positive")
        if self.n_sectors < 2 or self.n_active_protrusions not in (1, 2):
            raise ValueError("G3 supports one or two active protrusions and >=2 sectors")
        if self.n_active_protrusions > self.n_sectors:
            raise ValueError("active protrusions cannot exceed candidate sectors")
        if self.ecm_substeps < 1:
            raise ValueError("ecm_substeps must be a positive integer")
        # A chain bead couples to two neighbours; its largest extensional eigenvalue tends to
        # 4*kappa. Forward Euler therefore requires dt_ecm < zeta/(2*kappa).
        if self.dt / self.ecm_substeps >= self.bead_drag / (2.0 * self.kappa_s_f):
            raise ValueError("ECM substep violates the explicit bead-chain stability bound")

    def to_dict(self) -> dict[str, Any]:
        out = asdict(self)
        out["rotational_drag"] = self.rotational_drag
        out["beads_per_fibre"] = self.beads_per_fibre
        return out

    @classmethod
    def from_dict(cls, values: dict[str, Any]) -> "G3Config":
        allowed = {f.name for f in fields(cls)}
        unknown = set(values) - allowed
        if unknown:
            raise ValueError(f"unknown G3 config keys: {sorted(unknown)}")
        cfg = cls(**values)
        cfg.validate()
        return cfg

    @classmethod
    def from_yaml(cls, path: str | Path) -> "G3Config":
        with open(path, encoding="utf-8") as handle:
            values = yaml.safe_load(handle) or {}
        if "g3" in values:
            values = values["g3"]
        return cls.from_dict(values)


DEFAULT_CONFIG_PATH = Path(__file__).parents[1] / "config" / "params_g3.yaml"
