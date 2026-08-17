# Generation 3 — emergent cell–collagen guidance

**Status:** active 2D mechanism model. G2 remains frozen and is not overwritten.

Generation 3 tests whether spatial motor-clutch attachments, local collagen geometry, and
traction-dependent protrusion persistence can produce a migration axis without the prescribed
left/right bias used in G2 V3.

## Stage sequence

1. **G3A — spatial material-point clutches.** A fixed cell uses a prescribed test protrusion.
   Each clutch binds `(fiber_id, segment_id, alpha)` and remains on that material point until
   Bell-law unbinding. The local same-fibre Gaussian projection conserves total force and first
   moment.
2. **G3B — protrusion feedback.** Twenty-four initially unbiased sectors probe the ECM. Local
   collagen abundance/alignment affects sector selection, and successful traction increases
   protrusion persistence. The cell remains fixed so guidance can be debugged separately from
   motion.
3. **G3C — rigid-body motion.** Equal-and-opposite clutch reactions translate and rotate the
   cell. There is no prescribed velocity, fixed `+x` force, or `polarity_probability=0.65`.

The implementation lives in [`../../src/g3/`](../../src/g3/). Parameters are in
[`../../src/config/params_g3.yaml`](../../src/config/params_g3.yaml). Equations, provenance,
and gates are recorded in [`../../equations_talking_points.md`](../../equations_talking_points.md),
[`../../docs/parameter_provenance.md`](../../docs/parameter_provenance.md), and
[`../../docs/validation_targets.md`](../../docs/validation_targets.md).

## Run

From the repository root in PowerShell:

```powershell
py -3 -m pip install -r requirements.txt
$env:PYTHONPATH = "src"

py -3 -m g3.run --stage g3a --fixture single_fibre `
  --seed 4 --output results/g3/g3a-demo

py -3 -m g3.run --stage g3b --fixture aligned_8 `
  --seed 23 --duration 120 --output results/g3/g3b-demo

py -3 -m g3.run --stage g3c --fixture asymmetric_torque `
  --seed 3 --duration 30 --output results/g3/g3c-demo
```

Each run writes its resolved configuration, seed, Git commit, status, compact numerical traces,
metrics, frames, a summary PNG, and a GIF (or optional MP4).

## Test

```powershell
$env:PYTHONPATH = "src"
py -3 -m pytest tests/unit -q
```

The repository's G2 tests must also remain green:

```powershell
py -3 -m unittest discover -s generations/g2_corrected/v2_crosslink_transmission/tests -v
py -3 -m unittest discover -s generations/g2_corrected/v3_two_sided_migration/tests -v
py -3 -m unittest discover -s generations/g2_corrected/v4_contact_plasticity/tests -v
```

## Evidence boundary

Current saved outputs demonstrate implementation mechanisms, not final emergent-guidance
validation. The 100-seed, 600-s G3B/G3C ensemble is preregistered but has not yet been run.
Absolute speed is a calibrated output because cell drag is provisional.

Every G3 figure must be labeled:

> Personal simulation / mechanism validation; not yet a realistic 3D tumor-migration prediction.
