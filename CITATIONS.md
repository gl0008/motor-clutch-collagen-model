# Literature map: which paper informed which model choice

The papers motivate model structure; most default parameter values are still
assumptions and are labeled as such in `ASSUMPTIONS.md`.

1. **Lee et al. (2014), “A Three-Dimensional Computational Model of Collagen
   Network Mechanics.”** PLOS ONE 9:e111896.
   [doi:10.1371/journal.pone.0111896](https://doi.org/10.1371/journal.pone.0111896)
   - Basis for a bead-and-spring collagen network with explicit inter-fiber
     elastic links and local deformation propagation.
   - Source for the reported collagen-fiber scale: 20–200 µm length and
     200–350 nm thickness.  Their tropocollagen/fibril molecular values are nm;
     those are different hierarchical objects, not contradictory measurements.

2. **Vader et al. (2009), “Strain-Induced Alignment in Collagen Gels.”** PLOS
   ONE 4:e5902.
   [doi:10.1371/journal.pone.0005902](https://doi.org/10.1371/journal.pone.0005902)
   - Supports measuring alignment and densification under load.
   - Motivates the elastic permanent-link control: crosslinked gels can align
     while loaded yet recover, so alignment during loading alone is not proof
     of plastic remodeling.

3. **Bangasser & Odde (2013), “Master equation-based analysis of a
   motor-clutch model for cell traction force.”** Cellular and Molecular
   Bioengineering 6:449–459.
   [doi:10.1007/s12195-013-0296-5](https://doi.org/10.1007/s12195-013-0296-5)
   - Basis for stochastic effective clutches and stiffness-dependent traction.

4. **Bell (1978), “Models for the Specific Adhesion of Cells to Cells.”**
   Science 200:618–627.
   [doi:10.1126/science.347575](https://doi.org/10.1126/science.347575)
   - Basis for the minimal force-accelerated slip-bond off-rate in V0/V3.

5. **Prahl et al. (2020), “Predicting Confined 1D Cell Migration from
   Parameters Calibrated to a 2D Motor-Clutch Model.”** Biophysical Journal
   118:1709–1720.
   [doi:10.1016/j.bpj.2020.01.048](https://doi.org/10.1016/j.bpj.2020.01.048)
   - Supports opposing protrusion modules connected to a cell body and motion
     governed by force balance.  V3 is a simplified network-coupled analogue,
     not a reproduction of the full CMS.

6. **Adebowale et al. (2021), “Enhanced substrate stress relaxation promotes
   filopodia-mediated cell migration.”** Nature Materials 20:1290–1299.
   [doi:10.1038/s41563-021-00981-w](https://doi.org/10.1038/s41563-021-00981-w)
   - Motivation for the original V0 SLS hypothesis and protrusion-lifetime
     question.  It does not by itself establish that collagen bonds should be
     SLS elements in V2.

7. **Adebowale et al. (2025), “Monocytes use protrusive forces to generate
   migration paths in viscoelastic collagen-based extracellular matrices.”**
   PNAS 122:e2309772122.
   [doi:10.1073/pnas.2309772122](https://doi.org/10.1073/pnas.2309772122)
   - Biological motivation for eventually calibrating protrusive force and
     collagen-network remodeling in an actual collagen-based matrix.

8. **Kim et al. (2017), “Stress-induced plasticity of dynamic collagen
   networks.”** Nature Communications 8:842.
   [doi:10.1038/s41467-017-01011-7](https://doi.org/10.1038/s41467-017-01011-7)
   - Motivates V4's explicit elastic-versus-plastic and load–unload comparison.
   - The paper models irreversible sliding/merging; V4's rest-state reset is a
     simpler provisional hypothesis and must not be presented as their exact law.
