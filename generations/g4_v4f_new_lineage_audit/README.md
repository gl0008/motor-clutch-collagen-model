# G4 v4F_new — v4E → v4F lineage audit

`F_new` is an audit-only generation. It does **not** introduce a new physical
model, rerun either simulation, or replace the preserved v4E and v4F results.
Its purpose is to answer one question precisely:

> Did v4F rerun v4E while changing only the requested baseline mechanics?

The answer is **no**. Git ancestry is correct and the published scalar model
parameters were inherited, but v4F bundled requested mechanics corrections,
baseline harmonisation, numerical corrections, initial-condition changes and a
representative-run change. Therefore the large v4E → v4F trajectory difference
is a compound result and cannot be assigned to steric reaction alone.

## Read the audit

- [`LINEAGE_AUDIT.md`](LINEAGE_AUDIT.md) is the human-readable record.
- [`lineage_audit.json`](lineage_audit.json) is the machine-readable change
  contract.
- [`validate_audit.py`](validate_audit.py) checks that the audit covers the
  dynamics-changing differences and that the published representative seeds
  remain what the audit states.
- The public-facing explanation is `docs/g4-v4f-new.html`.

## Historical states remain frozen

| State | Implementation commit | Evidence commit |
|---|---:|---:|
| G4 v4E | `d49858b` | `c1248c2` |
| G4 v4F | `921e8fc` | `75d7dd9` |
| v4F merge | `f26751c` | — |

The audit names changes that were not explicitly requested. It does not
retroactively describe those changes as user-authorised one-factor changes.

## Future comparison policy

Direct lineage visualisations use seed 45 and preserve geometry, initial
contacts and stochastic-stream version unless one is the declared variable.
Scientific inference remains a matched 20-seed ensemble (41–60). RNG,
integrator, bias and initial-condition fixes must each be introduced as a
separate microstage before the next biological factor is tested.

