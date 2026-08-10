"""Standard-linear-solid constitutive updates in force/extension form."""

from dataclasses import dataclass
from math import exp


@dataclass
class StandardLinearSolid:
    """One SLS element: a long-time spring parallel to a Maxwell branch.

    Parameters use a force/length representation. ``k0`` is the instantaneous
    stiffness and ``kinf`` is the long-time stiffness. The Maxwell spring is
    ``k0 - kinf`` and its dashpot is ``eta = tau * (k0 - kinf)``.
    """

    k0: float
    kinf: float
    tau: float
    extension: float = 0.0
    maxwell_force: float = 0.0

    def __post_init__(self) -> None:
        if self.k0 <= 0.0:
            raise ValueError("k0 must be positive")
        if not 0.0 <= self.kinf <= self.k0:
            raise ValueError("kinf must satisfy 0 <= kinf <= k0")
        if self.tau <= 0.0:
            raise ValueError("tau must be positive")

    @property
    def km(self) -> float:
        return self.k0 - self.kinf

    @property
    def eta(self) -> float:
        return self.tau * self.km

    @property
    def force(self) -> float:
        return self.kinf * self.extension + self.maxwell_force

    def set_extension(self, value: float, dt: float = 0.0, instantaneous: bool = False) -> float:
        """Update extension and return force.

        ``instantaneous=True`` applies an ideal step, so the Maxwell spring
        receives ``km * delta_extension`` immediately. Otherwise the extension
        is assumed to change at a constant rate over ``dt`` and the Maxwell
        state is integrated exactly over that interval.
        """

        delta = value - self.extension
        if instantaneous:
            self.maxwell_force += self.km * delta
        else:
            if dt <= 0.0:
                raise ValueError("dt must be positive for a finite-rate update")
            decay = exp(-dt / self.tau)
            rate = delta / dt
            self.maxwell_force = (
                decay * self.maxwell_force
                + self.km * self.tau * (1.0 - decay) * rate
            )
        self.extension = value
        return self.force

