"""One-dimensional single-protrusion motor--clutch model.

The collagen fibre is a serial chain of identical SLS bonds.  Internal beads
are treated as massless points in quasi-static axial equilibrium, so every bond
carries the same force and the end-to-end constitutive response is independent
of the number of beads used to draw the fibre.
"""

from dataclasses import asdict, dataclass, replace
from math import exp
from typing import Dict, List, Tuple

import numpy as np


@dataclass(frozen=True)
class SingleProtrusionConfig:
    """Dimensionless parameters for the one-protrusion experiment."""

    n_beads: int = 10
    fibre_rest_length: float = 9.0
    n_clutches: int = 12
    clutch_stiffness: float = 1.0
    on_rate: float = 1.0
    unloaded_off_rate: float = 0.1
    rupture_force: float = 1.0
    force_dependent_off: bool = True
    unloaded_actin_speed: float = 1.0
    motor_stall_force: float = 12.0
    k0_chain: float = 2.0
    kinf_chain: float = 0.5
    tau_fast: float = 0.1
    tau_slow: float = 10.0
    dt: float = 1.0e-3
    duration: float = 30.0
    sample_interval: float = 0.01

    def __post_init__(self) -> None:
        if self.n_beads < 2:
            raise ValueError("n_beads must be at least 2")
        if self.fibre_rest_length <= 0.0:
            raise ValueError("fibre_rest_length must be positive")
        if self.n_clutches < 1:
            raise ValueError("n_clutches must be positive")
        if self.clutch_stiffness <= 0.0 or self.rupture_force <= 0.0:
            raise ValueError("clutch stiffness and rupture force must be positive")
        if self.on_rate < 0.0 or self.unloaded_off_rate < 0.0:
            raise ValueError("kinetic rates cannot be negative")
        if self.unloaded_actin_speed < 0.0 or self.motor_stall_force <= 0.0:
            raise ValueError("motor parameters must be positive")
        if self.k0_chain <= 0.0 or not 0.0 <= self.kinf_chain <= self.k0_chain:
            raise ValueError("chain stiffnesses must satisfy k0 > 0 and 0 <= kinf <= k0")
        if self.tau_fast <= 0.0 or self.tau_slow <= 0.0:
            raise ValueError("relaxation times must be positive")
        if self.dt <= 0.0 or self.duration <= 0.0 or self.sample_interval <= 0.0:
            raise ValueError("time parameters must be positive")

    @property
    def n_bonds(self) -> int:
        return self.n_beads - 1

    @property
    def k0_bond(self) -> float:
        return self.n_bonds * self.k0_chain

    @property
    def kinf_bond(self) -> float:
        return self.n_bonds * self.kinf_chain


_MASK64 = (1 << 64) - 1


def counter_uniform(seed: int, step: int, clutch: int, channel: int) -> float:
    """Deterministic U[0, 1) variate addressed by event coordinates.

    The counter address, rather than a mutable generator state, makes paired
    fast/slow simulations use the same random draw for the same possible event.
    """

    value = (
        int(seed)
        + 0x9E3779B97F4A7C15 * (int(step) + 1)
        + 0xD1B54A32D192ED03 * (int(clutch) + 1)
        + 0x94D049BB133111EB * (int(channel) + 1)
    ) & _MASK64
    value ^= value >> 30
    value = (value * 0xBF58476D1CE4E5B9) & _MASK64
    value ^= value >> 27
    value = (value * 0x94D049BB133111EB) & _MASK64
    value ^= value >> 31
    return ((value >> 11) & ((1 << 53) - 1)) / float(1 << 53)


def _solve_quasistatic_endpoint(
    config: SingleProtrusionConfig,
    tau: float,
    old_extension: float,
    old_q: float,
    actin_position: float,
    binding_reference: np.ndarray,
    bound: np.ndarray,
) -> Tuple[float, float, np.ndarray]:
    """Solve SLS-chain/clutch force balance over one finite-rate step."""

    decay = exp(-config.dt / tau)
    km = config.k0_chain - config.kinf_chain
    alpha = km * tau * (1.0 - decay) / config.dt
    candidate = bound.copy()
    new_extension = old_extension

    # Tension-only clutches require a tiny active-set solve.  The set can only
    # shrink, so convergence takes at most n_clutches + 1 passes.
    for _ in range(config.n_clutches + 1):
        count = int(np.count_nonzero(candidate))
        if count:
            free_extensions = actin_position - binding_reference[candidate]
            numerator = (
                config.clutch_stiffness * float(np.sum(free_extensions))
                - decay * old_q
                + alpha * old_extension
            )
            denominator = (
                config.clutch_stiffness * count + config.kinf_chain + alpha
            )
        else:
            numerator = -decay * old_q + alpha * old_extension
            denominator = config.kinf_chain + alpha
        new_extension = numerator / denominator
        new_candidate = bound & (
            actin_position - binding_reference > new_extension + 1.0e-14
        )
        if np.array_equal(new_candidate, candidate):
            break
        candidate = new_candidate

    new_q = decay * old_q + alpha * (new_extension - old_extension)
    clutch_force = np.zeros(config.n_clutches, dtype=float)
    clutch_force[bound] = config.clutch_stiffness * np.maximum(
        actin_position - binding_reference[bound] - new_extension,
        0.0,
    )
    return new_extension, new_q, clutch_force


def _phase_label(
    bound_count: int,
    had_binding: bool,
    had_rupture: bool,
    cluster_failed: bool,
    traction: float,
    old_traction: float,
) -> str:
    if cluster_failed:
        return "cluster failure"
    if had_binding:
        return "binding"
    if had_rupture:
        return "rupture"
    if bound_count == 0:
        return "unbound"
    if traction < old_traction - 1.0e-8:
        return "relaxation"
    return "loading"


def _simulate_condition(
    config: SingleProtrusionConfig,
    tau: float,
    seed: int,
    record: bool,
) -> Dict[str, object]:
    n_steps = int(round(config.duration / config.dt))
    sample_every = max(1, int(round(config.sample_interval / config.dt)))

    bound = np.zeros(config.n_clutches, dtype=bool)
    binding_reference = np.zeros(config.n_clutches, dtype=float)
    clutch_force = np.zeros(config.n_clutches, dtype=float)
    actin_position = 0.0
    chain_extension = 0.0
    q = 0.0
    traction = 0.0
    cluster_start = None
    episode_index = -1
    episodes: List[Dict[str, object]] = []
    events: List[Dict[str, object]] = []

    times: List[float] = []
    bead_positions: List[np.ndarray] = []
    sls_q: List[np.ndarray] = []
    clutch_bound: List[np.ndarray] = []
    clutch_forces: List[np.ndarray] = []
    total_traction: List[float] = []
    actin_velocity: List[float] = []
    loading_rate: List[float] = []
    cluster_episode: List[int] = []
    phase: List[str] = []

    latest_velocity = config.unloaded_actin_speed
    latest_rate = 0.0
    latest_phase = "unbound"

    def save(time: float) -> None:
        if not record:
            return
        fraction = np.linspace(0.0, 1.0, config.n_beads)
        positions = fraction * config.fibre_rest_length + fraction * chain_extension
        times.append(time)
        bead_positions.append(positions)
        sls_q.append(np.full(config.n_bonds, q, dtype=float))
        clutch_bound.append(bound.copy())
        clutch_forces.append(clutch_force.copy())
        total_traction.append(traction)
        actin_velocity.append(latest_velocity)
        loading_rate.append(latest_rate)
        cluster_episode.append(episode_index if cluster_start is not None else -1)
        phase.append(latest_phase)

    save(0.0)
    positive_loading_rates: List[float] = []

    for step in range(n_steps):
        time = step * config.dt
        old_traction = traction
        previously_bound = bound.copy()
        had_rupture = False
        had_binding = False
        cluster_failed = False

        for clutch in np.flatnonzero(previously_bound):
            force_ratio = clutch_force[clutch] / config.rupture_force
            multiplier = exp(min(force_ratio, 30.0)) if config.force_dependent_off else 1.0
            off_rate = config.unloaded_off_rate * multiplier
            probability = 1.0 - exp(-off_rate * config.dt)
            if counter_uniform(seed, step, int(clutch), 1) < probability:
                bound[clutch] = False
                clutch_force[clutch] = 0.0
                had_rupture = True
                events.append({"time": time, "type": "rupture", "clutch": int(clutch)})

        if cluster_start is not None and not np.any(bound):
            duration = time - cluster_start
            episodes.append(
                {"start": cluster_start, "end": time, "duration": duration, "observed": True}
            )
            events.append({"time": time, "type": "cluster failure", "clutch": -1})
            cluster_start = None
            cluster_failed = True

        # A clutch that ruptures cannot rebind until the next time step.
        on_probability = 1.0 - exp(-config.on_rate * config.dt)
        for clutch in np.flatnonzero(~previously_bound):
            if counter_uniform(seed, step, int(clutch), 0) < on_probability:
                bound[clutch] = True
                binding_reference[clutch] = actin_position - chain_extension
                had_binding = True
                events.append({"time": time, "type": "binding", "clutch": int(clutch)})

        if cluster_start is None and np.any(bound):
            cluster_start = time
            episode_index += 1

        latest_velocity = config.unloaded_actin_speed * max(
            1.0 - old_traction / config.motor_stall_force,
            0.0,
        )
        actin_position += latest_velocity * config.dt
        chain_extension, q, clutch_force = _solve_quasistatic_endpoint(
            config,
            tau,
            chain_extension,
            q,
            actin_position,
            binding_reference,
            bound,
        )
        traction = float(np.sum(clutch_force))
        latest_rate = (traction - old_traction) / config.dt
        if np.any(bound) and latest_rate > 0.0:
            positive_loading_rates.append(latest_rate)
        latest_phase = _phase_label(
            int(np.count_nonzero(bound)),
            had_binding,
            had_rupture,
            cluster_failed,
            traction,
            old_traction,
        )

        if (step + 1) % sample_every == 0 or step + 1 == n_steps:
            save((step + 1) * config.dt)

    if cluster_start is not None:
        episodes.append(
            {
                "start": cluster_start,
                "end": config.duration,
                "duration": config.duration - cluster_start,
                "observed": False,
            }
        )

    result: Dict[str, object] = {
        "tau": tau,
        "episodes": episodes,
        "events": events,
        "median_positive_loading_rate": (
            float(np.median(positive_loading_rates)) if positive_loading_rates else 0.0
        ),
        "mean_traction": (
            float(np.mean(total_traction))
            if record and total_traction
            else float(traction)
        ),
    }
    if record:
        result.update(
            {
                "time": np.asarray(times, dtype=float),
                "bead_positions": np.asarray(bead_positions, dtype=float),
                "sls_q": np.asarray(sls_q, dtype=float),
                "clutch_bound": np.asarray(clutch_bound, dtype=bool),
                "clutch_force": np.asarray(clutch_forces, dtype=float),
                "total_traction": np.asarray(total_traction, dtype=float),
                "actin_velocity": np.asarray(actin_velocity, dtype=float),
                "loading_rate": np.asarray(loading_rate, dtype=float),
                "cluster_episode": np.asarray(cluster_episode, dtype=int),
                "phase": np.asarray(phase, dtype=str),
            }
        )
    return result


def run_single_protrusion_pair(
    config: SingleProtrusionConfig,
    seed: int,
) -> Dict[str, object]:
    """Run synchronized fast/slow conditions with a common event stream."""

    return {
        "fast": _simulate_condition(config, config.tau_fast, seed, record=True),
        "slow": _simulate_condition(config, config.tau_slow, seed, record=True),
        "config": asdict(config),
        "seed": int(seed),
    }


def _kaplan_meier_median(durations: np.ndarray, observed: np.ndarray) -> float:
    """Return a right-censoring-aware median, or infinity if not reached."""

    if durations.size == 0:
        return float("nan")
    survival = 1.0
    for time in np.unique(durations[observed]):
        at_risk = int(np.count_nonzero(durations >= time))
        events = int(np.count_nonzero((durations == time) & observed))
        if at_risk:
            survival *= 1.0 - events / at_risk
        if survival <= 0.5:
            return float(time)
    return float("inf")


def _paired_trial_summary(args: Tuple[SingleProtrusionConfig, int]) -> Dict[str, Dict[str, object]]:
    config, trial_seed = args
    summary: Dict[str, Dict[str, object]] = {}
    for label, tau in (("fast", config.tau_fast), ("slow", config.tau_slow)):
        result = _simulate_condition(config, tau, trial_seed, record=False)
        trial_episodes = result["episodes"]
        summary[label] = {
            "durations": [float(item["duration"]) for item in trial_episodes],
            "observed": [bool(item["observed"]) for item in trial_episodes],
            "loading": float(result["median_positive_loading_rate"]),
            "completed": sum(bool(item["observed"]) for item in trial_episodes),
        }
    return summary


def _single_trial_summary(
    args: Tuple[SingleProtrusionConfig, float, int]
) -> Dict[str, object]:
    config, tau, trial_seed = args
    result = _simulate_condition(config, tau, trial_seed, record=False)
    episodes = result["episodes"]
    return {
        "durations": [float(item["duration"]) for item in episodes],
        "observed": [bool(item["observed"]) for item in episodes],
        "loading": float(result["median_positive_loading_rate"]),
    }


def run_lifetime_ensemble(
    config: SingleProtrusionConfig,
    trials: int = 200,
    seed: int = 0,
    workers: int = 1,
) -> Dict[str, object]:
    """Aggregate paired cluster lifetimes and loading rates across trials."""

    if trials < 1:
        raise ValueError("trials must be positive")
    if workers < 1:
        raise ValueError("workers must be positive")
    durations = {"fast": [], "slow": []}
    observed = {"fast": [], "slow": []}
    loading = {"fast": [], "slow": []}
    completed_counts = {"fast": [], "slow": []}

    arguments = [(config, int(seed) + trial) for trial in range(trials)]
    if workers == 1:
        summaries = map(_paired_trial_summary, arguments)
        executor = None
    else:
        from concurrent.futures import ProcessPoolExecutor

        try:
            executor = ProcessPoolExecutor(max_workers=min(workers, trials))
            summaries = executor.map(_paired_trial_summary, arguments)
        except (OSError, PermissionError):
            # Restricted runtimes may disallow process semaphores. Correctness
            # does not depend on parallel execution, so use the serial path.
            executor = None
            summaries = map(_paired_trial_summary, arguments)

    try:
        for summary in summaries:
            for label in ("fast", "slow"):
                durations[label].extend(summary[label]["durations"])
                observed[label].extend(summary[label]["observed"])
                loading[label].append(float(summary[label]["loading"]))
                completed_counts[label].append(int(summary[label]["completed"]))
    finally:
        if executor is not None:
            executor.shutdown()

    output: Dict[str, object] = {
        "trials": int(trials),
        "seed": int(seed),
        "config": asdict(config),
    }
    for label in ("fast", "slow"):
        duration_array = np.asarray(durations[label], dtype=float)
        observed_array = np.asarray(observed[label], dtype=bool)
        loading_array = np.asarray(loading[label], dtype=float)
        output[label] = {
            "episode_durations": duration_array,
            "episode_observed": observed_array,
            "trial_loading_rate": loading_array,
            "trial_completed_clusters": np.asarray(completed_counts[label], dtype=int),
            "km_median_lifetime": _kaplan_meier_median(duration_array, observed_array),
            "median_loading_rate": float(np.median(loading_array)),
        }
    return output


def run_mechanism_sweep(
    config: SingleProtrusionConfig,
    de_values: Tuple[float, ...] = (0.01, 0.1, 1.0, 10.0, 100.0),
    stiffness_values: Tuple[float, ...] = (0.5, 2.0, 10.0, 50.0),
    trials: int = 10,
    seed: int = 5000,
    workers: int = 1,
) -> List[Dict[str, object]]:
    """Sweep only De and Kchain/kc after the lifetime check is inconclusive.

    ``stiffness_values`` are instantaneous chain stiffnesses normalized by the
    clutch stiffness. The baseline ``kinf/k0`` ratio is held fixed, so the sweep
    does not silently tune a third mechanical quantity.
    """

    if trials < 1 or workers < 1:
        raise ValueError("trials and workers must be positive")
    if any(value <= 0.0 for value in tuple(de_values) + tuple(stiffness_values)):
        raise ValueError("sweep values must be positive")
    relaxed_ratio = config.kinf_chain / config.k0_chain
    jobs = []
    job_keys = []
    for stiffness in stiffness_values:
        sweep_config = replace(
            config,
            k0_chain=float(stiffness) * config.clutch_stiffness,
            kinf_chain=float(stiffness) * config.clutch_stiffness * relaxed_ratio,
        )
        for de in de_values:
            for trial in range(trials):
                jobs.append((sweep_config, float(de), int(seed) + trial))
                job_keys.append((float(stiffness), float(de)))

    executor = None
    if workers == 1:
        summaries = map(_single_trial_summary, jobs)
    else:
        from concurrent.futures import ProcessPoolExecutor

        try:
            executor = ProcessPoolExecutor(max_workers=min(workers, len(jobs)))
            summaries = executor.map(_single_trial_summary, jobs)
        except (OSError, PermissionError):
            executor = None
            summaries = map(_single_trial_summary, jobs)

    grouped: Dict[Tuple[float, float], Dict[str, list]] = {}
    try:
        for key, item in zip(job_keys, summaries):
            bucket = grouped.setdefault(
                key,
                {"durations": [], "observed": [], "loading": []},
            )
            bucket["durations"].extend(item["durations"])
            bucket["observed"].extend(item["observed"])
            bucket["loading"].append(item["loading"])
    finally:
        if executor is not None:
            executor.shutdown()

    records: List[Dict[str, object]] = []
    for stiffness in stiffness_values:
        for de in de_values:
            bucket = grouped[(float(stiffness), float(de))]
            duration_array = np.asarray(bucket["durations"], dtype=float)
            observed_array = np.asarray(bucket["observed"], dtype=bool)
            records.append(
                {
                    "de": float(de),
                    "k0_chain_over_kc": float(stiffness),
                    "kinf_over_k0": relaxed_ratio,
                    "episodes": int(duration_array.size),
                    "observed_failures": int(np.count_nonzero(observed_array)),
                    "km_median_lifetime": _kaplan_meier_median(
                        duration_array, observed_array
                    ),
                    "median_loading_rate": float(np.median(bucket["loading"])),
                }
            )
    return records


def elastic_limit_config(config: SingleProtrusionConfig) -> SingleProtrusionConfig:
    """Convenience helper used by tests and mechanism controls."""

    return replace(config, kinf_chain=config.k0_chain)
