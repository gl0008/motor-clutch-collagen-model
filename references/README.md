# Reference registry by model version

This directory records the literature used to choose model structure,
parameter scales and validation targets. The full evidence notes and explicit
limits of each citation are in [`../CITATIONS.md`](../CITATIONS.md); the
assumption-to-evidence register is in
[`../ASSUMPTIONS.md`](../ASSUMPTIONS.md).

The repository stores bibliographic metadata and stable links, not publisher
PDFs. A linked paper may be copyrighted or require institutional access; it
should only be copied into the repository when redistribution rights are clear.

## Which papers support each version?

| Model | Sources actually used | What they informed |
|---|---|---|
| G1 V0 | R1, R5, R6, R13 | bead framework; Bell rupture; motor–clutch loading; stress-relaxation hypothesis |
| G1 V1 | R1, R2, R8 | connected bead fibres; cell-driven reorientation question; breast-cancer/collagen context |
| G1 V2 | R1–R4 | finite fibre geometry, crosslinks, heterogeneous transmission, alignment and microbuckling follow-up |
| G1 V3 | R5–R7, R9 | Bell clutch kinetics, opposing protrusions, polarity and migration-speed scale |
| G1 V4 | R2, R11, R12 | load–unload remodeling question, crosslink/plasticity candidates and irreversible rearrangement |
| G2 V2 | R1–R4, R10 | physical fibre scale, connected transmission, compression softening and spatial response metrics |
| G2 V3 | R5–R9 | clutch law, motor feedback, polarity, migration speed and future force calibration |
| G2 V4 | R1, R2, R11, R12 | new-contact weak-link hypothesis, load–unload control and alternative sliding/merging mechanism |
| Separate future U937 track | R14 | outward protrusive path opening; deliberately not mixed into G2 pulling |

## Bibliography

- **R1 — Lee et al. (2014).** “A Three-Dimensional Computational Model of
  Collagen Network Mechanics.” *PLOS ONE* 9:e111896.
  [doi:10.1371/journal.pone.0111896](https://doi.org/10.1371/journal.pone.0111896)
- **R2 — Abhilash et al. (2014).** “Remodeling of Fibrous Extracellular
  Matrices by Contractile Cells: Predictions from Discrete Fiber Network
  Simulations.” *Biophysical Journal* 107:1829–1840.
  [doi:10.1016/j.bpj.2014.08.029](https://doi.org/10.1016/j.bpj.2014.08.029)
- **R3 — Wang et al. (2014).** “Long-Range Force Transmission in Fibrous
  Matrices Enabled by Tension-Driven Alignment of Fibers.” *Biophysical
  Journal* 107:2592–2603.
  [doi:10.1016/j.bpj.2014.09.044](https://doi.org/10.1016/j.bpj.2014.09.044)
- **R4 — Notbohm et al. (2015).** “Microbuckling of Fibrin Provides a
  Mechanism for Cell Mechanosensing.” *Journal of the Royal Society Interface*.
  [arXiv:1407.3510](https://arxiv.org/abs/1407.3510)
- **R5 — Bell (1978).** “Models for the Specific Adhesion of Cells to Cells.”
  *Science* 200:618–627.
  [doi:10.1126/science.347575](https://doi.org/10.1126/science.347575)
- **R6 — Bangasser & Odde (2013).** “Master equation-based analysis of a
  motor-clutch model for cell traction force.” *Cellular and Molecular
  Bioengineering* 6:449–459.
  [doi:10.1007/s12195-013-0296-5](https://doi.org/10.1007/s12195-013-0296-5)
- **R7 — Prahl et al. (2020).** “Predicting Confined 1D Cell Migration from
  Parameters Calibrated to a 2D Motor-Clutch Model.” *Biophysical Journal*
  118:1709–1720.
  [doi:10.1016/j.bpj.2020.01.048](https://doi.org/10.1016/j.bpj.2020.01.048)
- **R8 — Steinwachs et al. (2016).** “Three-dimensional force microscopy of
  cells in biopolymer networks.” *Nature Methods* 13:171–176.
  [doi:10.1038/nmeth.3685](https://doi.org/10.1038/nmeth.3685)
- **R9 — Aguilar-Rojas et al. (2024).** “Three-dimensional cell culture
  conditions promoted the Mesenchymal-Amoeboid Transition in the
  Triple-Negative Breast Cancer cell line MDA-MB-231.” *Frontiers in Cell and
  Developmental Biology*.
  [PMCID: PMC11327030](https://pmc.ncbi.nlm.nih.gov/articles/PMC11327030/)
- **R10 — Han et al. (2018).** “Cell contraction induces long-ranged stress
  stiffening in the extracellular matrix.” *PNAS*.
  [PMCID: PMC5910866](https://pmc.ncbi.nlm.nih.gov/articles/PMC5910866/)
- **R11 — Ban et al. (2018).** “Mechanisms of Plastic Deformation in Collagen
  Networks Induced by Cellular Forces.” *Biophysical Journal* 114:450–461.
  [doi:10.1016/j.bpj.2017.11.3739](https://doi.org/10.1016/j.bpj.2017.11.3739)
- **R12 — Kim et al. (2017).** “Stress-induced plasticity of dynamic collagen
  networks.” *Nature Communications* 8:842.
  [doi:10.1038/s41467-017-01011-7](https://doi.org/10.1038/s41467-017-01011-7)
- **R13 — Adebowale et al. (2021).** “Enhanced substrate stress relaxation
  promotes filopodia-mediated cell migration.” *Nature Materials* 20:1290–1299.
  [doi:10.1038/s41563-021-00981-w](https://doi.org/10.1038/s41563-021-00981-w)
- **R14 — Adebowale et al. (2025).** “Monocytes use protrusive forces to
  generate migration paths in viscoelastic collagen-based extracellular
  matrices.” *PNAS* 122:e2309772122.
  [doi:10.1073/pnas.2309772122](https://doi.org/10.1073/pnas.2309772122)

## Evidence boundary

These papers motivate candidate equations and plausible scales. They do not
make a coarse-grained parameter set a fit. Every model notebook therefore says
both what a paper supports and what it does **not** establish.
