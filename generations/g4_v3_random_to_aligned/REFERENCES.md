# Evidence map for G4 v3

Only the claim listed here is imported from each paper. The code does not imply
that one paper validates the entire model.

1. **Abhilash AS, Baker BM, Trappmann B, Chen CS, Shenoy VB (2014).**
   “Remodeling of Fibrous Extracellular Matrices by Contractile Cells:
   Predictions from Discrete Fiber Network Simulations.” *Biophysical Journal*
   107:1829–1840. DOI: [10.1016/j.bpj.2014.08.029](https://doi.org/10.1016/j.bpj.2014.08.029).
   [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4213674/).
   Used for: random discrete fibres, linked intersections, fixed outer boundary,
   prescribed cell-shape contraction, and the bending-to-stretching energy
   transition as a fibre-recruitment mechanism.

2. **Provenzano PP, Eliceiri KW, Campbell JM, Inman DR, White JG, Keely PJ
   (2006).** “Collagen reorganization at the tumor-stromal interface facilitates
   local invasion.” *BMC Medicine* 4:38. DOI:
   [10.1186/1741-7015-4-38](https://doi.org/10.1186/1741-7015-4-38).
   [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC1781458/).
   Used for: the experimental distinction between circumferential collagen near
   non-invasive regions and radially aligned collagen associated with invasion.

3. **Kopanska KS, Alcheikh Y, Staneva R, Vignjevic D, Betz T (2016).**
   “Tensile Forces Originating from Cancer Spheroids Facilitate Tumor Invasion.”
   *PLOS ONE* 11:e0156442. DOI:
   [10.1371/journal.pone.0156442](https://doi.org/10.1371/journal.pone.0156442).
   [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4896628/).
   Used for: initially isotropic collagen around tumor spheroids, later
   surface-parallel and radial bundles, and contraction preceding invasion.

4. **Vader D, Kabla A, Weitz D, Mahadevan L (2009).** “Strain-Induced
   Alignment in Collagen Gels.” *PLOS ONE* 4:e5902. DOI:
   [10.1371/journal.pone.0005902](https://doi.org/10.1371/journal.pone.0005902).
   [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC2691583/).
   Used for: collagen alignment can arise under applied strain and can be
   reversible, so plasticity is not required for the first alignment test.

5. **Wang H, Abhilash AS, Chen CS, Wells RG, Shenoy VB (2014).** “Long-range
   force transmission in fibrous matrices enabled by tension-driven alignment
   of fibers.” *Biophysical Journal* 107:2592–2603. DOI:
   [10.1016/j.bpj.2014.09.044](https://doi.org/10.1016/j.bpj.2014.09.044).
   [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4255175/).
   Used for: tension-driven alignment as a route to long-range force
   transmission in fibrous, rather than continuum, matrices.

6. **Böhringer D et al. (2023).** “Fiber alignment in 3D collagen networks as a
   biophysical marker for cell contractility.” *Matrix Biology* 124:39–48.
   DOI: [10.1016/j.matbio.2023.11.004](https://doi.org/10.1016/j.matbio.2023.11.004).
   [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10872942/).
   Used for: shell-resolved orientation as an observable around contractile cells
   in initially disordered 3D collagen networks.

7. **Lee B, Zhou X, Riching KM, Eliceiri KW, Keely PJ, Guelcher SA, Weaver AM,
   Jiang Y (2014).** “A Three-Dimensional Computational Model of Collagen
   Network Mechanics.” *PLOS ONE* 9:e111896. DOI:
   [10.1371/journal.pone.0111896](https://doi.org/10.1371/journal.pone.0111896).
   [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4227658/).
   Used for: bead-based finite collagen fibres and explicit crosslink control as
   the framework lineage. It is not cited as proof of cell-driven radial
   alignment.

## Evidence boundary

No paper above directly validates the full combination of twelve stochastic
clutch bundles, a two-dimensional random bead-spring collagen network and
repeated nearest-material-point relocation. That combination remains a model
hypothesis. The exact relocation rule will therefore be isolated as the last
cumulative G4 v3 block rather than assumed in the base mechanics.

## Closely related work not imported into this elastic baseline

**Ley AWY, Bersie-Larson LM, Adhikari S, Provenzano PP, Dorfman KD, Barocas VH
(2026).** “Cell-Extracellular Matrix Feedback Results in Spontaneous Cell
Polarization and Heterogeneous Remodeling in 3D Isotropic and Aligned
Discrete-Fiber Models of Cell-Mediated Remodeling.” *Cellular and Molecular
Bioengineering*. DOI:
[10.1007/s12195-026-00922-0](https://doi.org/10.1007/s12195-026-00922-0).

This is highly relevant: it uses a 3D bead-spring fibre model, retracting
cellular tractors, proximity-mediated interactions and cell–ECM feedback, and
reports radial reorientation in initially unaligned networks. It is **not**
used to justify adding those mechanisms to G4 v3, because its reported
measurable remodeling depends on cumulative retraction cycles and plasticity.
Those mechanisms belong after the present elastic traction gate, not inside it.
