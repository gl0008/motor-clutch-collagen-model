# Literature map: which source informed which choice

The sources motivate structure and plausible scales.  They do not turn the
coarse-grained defaults into a fit; assumptions remain explicit in
[`ASSUMPTIONS.md`](ASSUMPTIONS.md).

1. **Lee et al. (2014), “A Three-Dimensional Computational Model of Collagen
   Network Mechanics.”** PLOS ONE 9:e111896.
   [doi:10.1371/journal.pone.0111896](https://doi.org/10.1371/journal.pone.0111896)
   - Explicit bead-and-spring collagen fibres and inter-fibre links.
   - Source for the fibre-level 20–200 µm length, 200–350 nm thickness and
     32 MPa wet-modulus scales.
   - Their ≈300 nm long, ≈1.5 nm diameter object is a molecular fibril-level
     object, not the 20–200 µm fibre represented by one bead chain.

2. **Abhilash et al. (2014), “Remodeling of Fibrous Extracellular Matrices by
   Contractile Cells: Predictions from Discrete Fiber Network Simulations.”**
   Biophysical Journal 107:1829–1840.
   [doi:10.1016/j.bpj.2014.08.029](https://doi.org/10.1016/j.bpj.2014.08.029)
   - Supports discrete-fibre alignment, heterogeneous deformation and
     long-range force transmission through connected fibre tracts.
   - Motivates outer-boundary anchoring and the radial displacement/alignment
     profiles rather than judging only the directly pulled fibre.

3. **Wang et al. (2014), “Long-Range Force Transmission in Fibrous Matrices
   Enabled by Tension-Driven Alignment of Fibers.”** Biophysical Journal
   107:2592–2603.
   [doi:10.1016/j.bpj.2014.09.044](https://doi.org/10.1016/j.bpj.2014.09.044)
   - Supports tension-driven fibre recruitment/alignment as a source of
     long-range nonlinear response; Generation 2 reports strain rather than
     assuming a purely visual reorientation is validated.

4. **Notbohm et al. (2015), “Microbuckling of Fibrin Provides a Mechanism for
   Cell Mechanosensing.”** Journal of the Royal Society Interface.
   [arXiv:1407.3510](https://arxiv.org/abs/1407.3510)
   - Motivates testing compression softening/microbuckling rather than giving
     fibres identical tensile and compressive response.

5. **Bell (1978), “Models for the Specific Adhesion of Cells to Cells.”**
   Science 200:618–627.
   [doi:10.1126/science.347575](https://doi.org/10.1126/science.347575)
   - Basis for the minimal force-accelerated slip-bond off-rate.

6. **Bangasser & Odde (2013), “Master equation-based analysis of a
   motor-clutch model for cell traction force.”** Cellular and Molecular
   Bioengineering 6:449–459.
   [doi:10.1007/s12195-013-0296-5](https://doi.org/10.1007/s12195-013-0296-5)
   - Basis for stochastic effective clutches and motor–clutch loading.

7. **Prahl et al. (2020), “Predicting Confined 1D Cell Migration from
   Parameters Calibrated to a 2D Motor-Clutch Model.”** Biophysical Journal
   118:1709–1720.
   [doi:10.1016/j.bpj.2020.01.048](https://doi.org/10.1016/j.bpj.2020.01.048)
   - Supports opposing protrusion modules and the requirement for persistent
     polarity to obtain directional rather than zero-mean migration.
   - Motivates keeping corrected V3 on one left/right migration axis before
     adding full 2-D translation and rotation.

8. **Steinwachs et al. (2016), “Three-dimensional force microscopy of cells in
   biopolymer networks.”** Nature Methods 13:171–176.
   [doi:10.1038/nmeth.3685](https://doi.org/10.1038/nmeth.3685)
   - Supports using MDA-MB-231 in collagen as the pulling-model calibration
     context and motivates future force fitting.

9. **Aguilar-Rojas et al. (2024), “Three-dimensional cell culture conditions
   promoted the Mesenchymal-Amoeboid Transition in the Triple-Negative Breast
   Cancer cell line MDA-MB-231.”** Frontiers in Cell and Developmental Biology.
   [PMCID: PMC11327030](https://pmc.ncbi.nlm.nih.gov/articles/PMC11327030/)
   - Reports mean instantaneous speeds of about 0.4, 0.3 and 0.07 µm/min in 1,
     3 and 6 mg/mL collagen, respectively.
   - Basis for the mechanism-first 0.2–0.4 µm/min V3 mobility target, not a fit
     to a single concentration-specific trajectory.

10. **Han et al. (2018), “Cell contraction induces long-ranged stress
    stiffening in the extracellular matrix.”** PNAS.
    [PMCID: PMC5910866](https://pmc.ncbi.nlm.nih.gov/articles/PMC5910866/)
    - Demonstrates MDA-MB-231-induced long-range effects in 1.5 mg/mL collagen
      and motivates near/intermediate/far response metrics.

11. **Ban et al. (2018), “Mechanisms of Plastic Deformation in Collagen
    Networks Induced by Cellular Forces.”** Biophysical Journal 114:450–461.
    [doi:10.1016/j.bpj.2017.11.3739](https://doi.org/10.1016/j.bpj.2017.11.3739)
    - Motivates stretch/approach-dependent formation of new weak crosslinks
      between nearby fibres and load–unload comparison.
    - Corrected V4 implements only a minimal candidate of this mechanism and
      currently produces a negative baseline result.

12. **Kim et al. (2017), “Stress-induced plasticity of dynamic collagen
    networks.”** Nature Communications 8:842.
    [doi:10.1038/s41467-017-01011-7](https://doi.org/10.1038/s41467-017-01011-7)
    - Shows irreversible sliding/merging as another plasticity mechanism.
    - This mechanism is not yet implemented; it is a defined alternative if the
      weak-link candidate remains inactive.

13. **Adebowale et al. (2021), “Enhanced substrate stress relaxation promotes
    filopodia-mediated cell migration.”** Nature Materials 20:1290–1299.
    [doi:10.1038/s41563-021-00981-w](https://doi.org/10.1038/s41563-021-00981-w)
    - Motivation for the legacy V0 SLS/clutch-lifetime hypothesis.
    - It does not establish that corrected collagen bonds should contain SLS
      elements.

14. **Adebowale et al. (2025), “Monocytes use protrusive forces to generate
    migration paths in viscoelastic collagen-based extracellular matrices.”**
    PNAS 122:e2309772122.
    [doi:10.1073/pnas.2309772122](https://doi.org/10.1073/pnas.2309772122)
    - U937 monocytes open paths using outward protrusive forces and can migrate
      without the same adhesion-based pulling mechanism.
    - This is reserved for a separate future biology track and is not mixed
      into corrected V2/V3.

15. **Carey et al. (2016), “Local extracellular matrix alignment directs
    cellular protrusion dynamics and migration through Rac1 and FAK.”**
    Integrative Biology 8:821–835.
    [doi:10.1039/C6IB00030D](https://doi.org/10.1039/C6IB00030D)
    - Shows that aligned 3D collagen biases protrusion formation, length and persistence along
      the matrix axis in MDA-MB-231 cells.
    - Motivates the sign of G3's geometry/persistence feedback, but not its coarse-grained
      coefficients or turnover rates.

16. **Nam et al. (2016), “Viscoplasticity Enables Mechanical Remodeling of
    Matrix by Cells.”** Biophysical Journal 111:2296–2308.
    [doi:10.1016/j.bpj.2016.10.002](https://doi.org/10.1016/j.bpj.2016.10.002)
    - Provides the load–unload plasticity metric used as G3's elastic-recovery diagnostic.
    - G3 has no irreversible material mechanism; unresolved recovery is not called plasticity.

17. **Saraswathibhatla et al. (2025), “Swirling motion of breast cancer cells
    radially aligns collagen fibers to enable collective invasion.”** bioRxiv preprint.
    [doi:10.1101/2025.01.31.635980](https://doi.org/10.1101/2025.01.31.635980)
    - Provides G3's bead spacing, elastic stiffness and bending coefficient scales.
    - Its annular spheroid geometry and full-density network are later validation/scale-up
      targets, not features of the current eight-fibre fixtures.

18. **Runser, Vetter & Iber (2024), “SimuCell3D: Three-dimensional simulation
    of tissue mechanics with cell polarization.”** Nature Computational Science 4:299–309.
    [doi:10.1038/s43588-024-00620-9](https://doi.org/10.1038/s43588-024-00620-9)
    - Supports the repulsive signed-distance branch used for G3's one-sided
      conservative bead–cell contact potential.
    - Does not establish the provisional contact stiffness, zero bead radius,
      frictionless interface or rigid-cell assumption.

19. **Mori, Jilkine & Edelstein-Keshet (2008), “Wave-Pinning and Cell
    Polarity from a Bistable Reaction-Diffusion System.”** Biophysical Journal
    94:3684–3697.
    [PMCID: PMC2292363](https://pmc.ncbi.nlm.nih.gov/articles/PMC2292363/)
    - Supports mass conservation, unequal active/inactive diffusivity and
      bistable local kinetics as one minimal route to spontaneous cell polarity.
    - Motivates the proposed G3-R cell-intrinsic activity field, not an
      MDA-MB-231 parameter fit or the claim that wave-pinning is the unique
      polarity mechanism.

20. **Fraley et al. (2010), “A distinctive role for focal adhesion proteins in
    three-dimensional cell motility.”** Nature Cell Biology 12:598–604.
    [doi:10.1038/ncb2062](https://doi.org/10.1038/ncb2062)
    - Supports treating protrusion activity and matrix deformation as a coupled
      three-dimensional motility/traction module.
    - Warns against requiring large, visually discrete 2D-style focal-adhesion
      plaques in an embedded-cell model.

21. **Fraley et al. (2015), “Three-dimensional matrix fiber alignment modulates
    cell migration and MT1-MMP utility by spatially and temporally directing
    protrusions.”** Scientific Reports 5:14580.
    [doi:10.1038/srep14580](https://doi.org/10.1038/srep14580)
    - Supports protrusion orientation/rate and fibre alignment as validation
      observables across collagen conditions.
    - Does not allow collagen concentration, pore size, alignment, crosslinking
      and stiffness to be collapsed into one interchangeable model parameter.

22. **Wolf et al. (2013), “Physical limits of cell migration: Control by ECM
    space and nuclear deformation and tuning by proteolysis and traction
    force.”** Journal of Cell Biology 201:1069–1084.
    [PMCID: PMC3691458](https://pmc.ncbi.nlm.nih.gov/articles/PMC3691458/)
    - Supports pore/nucleus geometry, traction and proteolysis as coupled limits
      on three-dimensional cell passage and arrest.
    - Motivates a later deformable-cell/nucleus generation; it does not justify
      adding an uncalibrated nucleus to the G3-R attachment gate.

23. **Licup et al. (2015), “Stress controls the mechanics of collagen
    networks.”** PNAS 112:9573–9578.
    [doi:10.1073/pnas.1504258112](https://doi.org/10.1073/pnas.1504258112)
    - Supports stress-dependent nonlinear stiffening and the role of network
      geometry/connectivity in its onset.
    - Motivates an independently calibrated nonlinear network branch; it does
      not provide one universal strain threshold for the current synthetic
      network.

24. **Nam et al. (2016), “Strain-enhanced stress relaxation impacts nonlinear
    elasticity in collagen gels.”** PNAS 113:5492–5497.
    [doi:10.1073/pnas.1523906113](https://doi.org/10.1073/pnas.1523906113)
    - Shows faster collagen stress relaxation at larger strain and supports
      force-dependent weak interfibre bond unbinding/rebinding as a candidate
      reversible mechanism.
    - Does not establish that every fibre bond is an SLS element or that
      reversible stress relaxation is permanent remodeling.

25. **Krause et al. (2019), “Fiber stiffness, pore size and adhesion control
    migratory phenotype of MDA-MB-231 cells in collagen gels.”** PLOS ONE.
    [PMCID: PMC6853323](https://pmc.ncbi.nlm.nih.gov/articles/PMC6853323/)
    - Supports treating fibre stiffness, pore size and adhesion as distinct
      controls of MDA-MB-231 migration phenotype.
    - Motivates later local-pore and phenotype validation; it does not map the
      current eight-fibre fixture or 99-fibre planar network to a collagen
      concentration.
