# Random collagen to tumor-associated alignment

## Executive conclusion

The transition from an initially disordered collagen matrix to a locally
aligned matrix is experimentally established. Tumor explants, isolated cells,
cell pairs and tumor spheroids can contract, displace, densify and orient
collagen. Independent mechanical experiments also show that applied strain is
sufficient to align collagen, including chemically crosslinked collagen that
returns toward its initial organization after unloading.

What is not established is one universal sequence of molecular rules that
quantitatively connects integrin-clutch kinetics to TACS architecture over
hours or days. Existing papers validate different parts of that chain at
different scales. Therefore G4 v4 treats alignment as a measured outcome, uses
papers only for the pieces they directly support, and leaves plasticity,
crosslink rupture, proteolysis, growth and collective-cell feedback outside the
single-cell elastic baseline.

## Evidence matrix

| Study | System and starting ECM | Observation clock | Driver and alignment measure | Supports in G4 v4 | Does **not** support |
|---|---|---:|---|---|---|
| Provenzano et al., 2006 | Mouse mammary tumor tissue and tumor explants embedded in initially random 3D collagen | Explants shown after 8 h | SHG/MPE imaging; noninvasive regions had boundary-wrapped collagen, invasive regions had radial collagen in direct cell contact | Random-to-radial alignment can be produced by tumor-cell activity; boundary-relative orientation matters | A numerical clutch law, a unique force magnitude, or the claim that a cell shrinks |
| Vader et al., 2009 | 1 mg/ml collagen gels; untreated and glutaraldehyde-crosslinked | Applied strain cycles; U87 examples at 10–48 h | Controlled strain, orientation distributions, order parameter and density | Elastic fiber mechanics alone can produce alignment; plasticity is unnecessary for the first alignment test | Tumor migration, motor-clutch parameters, or biochemical crosslink density |
| Conklin et al., 2011 | SHG images from 196 human breast-cancer specimens | Clinical endpoint | TACS-3 scoring relative to tumor boundary | Radial/perpendicular collagen is clinically meaningful and must be measured explicitly | How TACS-3 forms, or that it necessarily causes every observed invasion event |
| Ma et al., 2013 | Fibroblast-seeded type-I collagen; image-derived fibrous finite-element models | 4 and 24 h | Collagen compaction/alignment between cell pairs; comparison with homogeneous materials | Fibrous architecture and aligned paths enable long-range stress transmission | A bead-crosslink probability, single-cell tumor kinetics, or material-point clutch relocation |
| Abhilash et al., 2014 | Random discrete-fiber networks surrounding circular/elliptical contracting inclusions | Mechanical loading, not a biological migration clock | Prescribed contraction, fiber orientation, stretch/bend energy and zone of influence | The contractile-cavity benchmark and bending-to-stretching recruitment diagnostic | Biological shrinkage of a tumor cell or an experimentally measured MDA-MB-231 traction program |
| Wang et al., 2014 | Continuum constitutive model of fibrous matrix with cell-shape/contraction controls | Quasi-static mechanics | Tension-driven reorientation, cell aspect ratio and contraction anisotropy | Dipole as a secondary control; aligned tensile paths transmit force farther | The exact discrete contact rule or a proof that a circular cell must polarize |
| Riching et al., 2014 | MDA-MB-231 cells in 1–4 mg/ml random/wide and aligned/narrow collagen microchannels | Fiber displacement: 2 h at 4-min intervals; migration: 6 h at 10-min intervals; invasion endpoint: 3 d | FITC-collagen displacement, migration tracks, protrusion number/length, alignment coefficient | Primary G4 v4 observation clocks; alignment can increase persistence without increasing speed; protrusions follow aligned fibers | A force-calibrated random-to-aligned model or molecular values for this project’s effective clutches |
| Kopanska et al., 2016 | CT26 spheroids embedded in collagen I | Immediate contraction; invasion around 9–12 h; measurements through 72 h | Collagen flow, laser ablation, asymmetric cuts, spheroid invasion | Spheroid tension can reorganize collagen; spheroid growth and inward ECM motion can coexist | Reusing CT26 collective parameters for one MDA-MB-231 cell, or representing contraction as cell shrinkage |
| Geiger et al., 2022 | Cancer spheroids in microfluidically generated radial or tangential collagen | 1–3 d | Invasion front relative to collagen orientation; experiment plus migration simulation | Radial and tangential orientation must be distinguished; radial fibers bias outward invasion | Formation of alignment from a random network or the current motor-clutch force law |
| Poonja et al., 2023 | Multicellular lattice-free tumor/ECM agent model | Simulated days | Agent interaction rules generate TACS-1/2/3-like patterns | Demonstrates that cumulative rule-based TACS-transition tests are possible | Direct validation of bead-spring mechanics, crosslink force transmission or integrin kinetics |
| Böhringer et al., 2023 | Hepatic stellate cells and glioblastoma lines in 1.0–1.2 mg/ml collagen; simulated dipoles/quadrupoles | Experiments over 1–8 d | Shell-wise collagen orientation and intensity; traction microscopy; 5 nN point-force simulations | Shell-wise orientation is a useful contractility proxy; dipole/quadrupole comparison is literature grounded | MDA-MB-231-specific calibration or uniform isotropic cell shrinkage |
| Ley et al., 2026 | 3D isotropic/aligned discrete-fiber model with cell–ECM feedback | Model-dependent repeated cycles | Tractors, retraction cycles, feedback and plasticity | A future comparison for spontaneous polarization and heterogeneous remodeling | Permission to import its plasticity and feedback into the clean elastic G4 v4 baseline |
| Chan & Odde, 2008 | Chick forebrain growth-cone filopodia on compliant 2D substrates | Nanoscale load-and-fail cycles | Stochastic motor–clutch model and measured actin/traction dynamics | Continued loading after contact, motor stall, stochastic clutch failure | MDA-MB-231 or 3D-collagen-specific clutch constants |
| Bangasser, Rosenfeld & Odde, 2013 | Motor–clutch theory on compliant substrates | Model load-and-fail cycle | Stiffness-dependent traction and retrograde flow | Effective motor–clutch cycling and the need to label inherited constants uncalibrated | A collagen-fiber contact rule or this model's equal-load-sharing site |
| Erdmann & Schwarz, 2004 | Theoretical receptor–ligand adhesion clusters | Stochastic cluster lifetime under shared force/loading | Equal load sharing, Bell rupture, rebinding and complete dissociation | Shared-load cluster as a transparent comparison to independent clutches | Integrin catch bonding, adhesion maturation, or MDA-MB-231 parameter values |

## 1. Tumor-associated collagen and boundary orientation

Provenzano and colleagues established the central biological observation used
here. In tumor tissue, noninvasive interfaces could have collagen wrapped
around the tumor boundary, while locally invasive regions contained radially
oriented collagen. In tumor explants embedded in an initially random collagen
gel, cells pulled collagen toward the explant and generated radial alignment at
invasive regions. This is direct evidence that random-to-aligned collagen is a
valid outcome to test. It does not identify a unique microscopic force law.

Conklin and colleagues subsequently evaluated TACS-3 in 196 human breast-cancer
specimens. TACS-3—straight collagen bundles perpendicular to the tumor
boundary—was associated with poorer disease-specific and disease-free
survival. This makes boundary-relative orientation a relevant observable, but
the clinical association cannot determine which force component created the
fibers.

G4 v4 therefore records two angles rather than using “parallel” ambiguously:

\[
\alpha_r=\cos^{-1}|\mathbf t\cdot\mathbf e_r|,
\qquad
\alpha_b=\cos^{-1}|\mathbf t\cdot\mathbf t_{\mathrm{boundary}}|.
\]

Radial collagen has \(\alpha_r\approx0^\circ\) and
\(\alpha_b\approx90^\circ\). Boundary-parallel collagen has the opposite
values.

## 2. Applied strain can align collagen without plasticity

Vader and colleagues mechanically deformed collagen gels and measured fiber
orientation and densification. Both untreated and glutaraldehyde-crosslinked
networks aligned under applied strain. Crosslinked gels recovered much of the
organization after unloading, whereas untreated networks could retain
irreversible changes. The important modeling consequence is narrow: an elastic
network can be used to test whether force produces alignment. Permanent
remodeling after release is a separate question.

This justifies keeping G4 v4 crosslinks permanent and elastic during the first
mechanics and clutch tests. It does not mean real tumor collagen lacks
viscoelasticity, crosslink turnover or plasticity.

## 3. Why transverse and axial loading have different roles

Abhilash and colleagues modeled random discrete fibers around contracting
cells. Their analysis distinguishes an early bending-dominated regime from a
stretching-dominated regime reached as fibers rotate toward a tensile loading
direction. Wang and colleagues likewise showed that tension-driven alignment,
fiber stiffness, strain hardening, cell shape and contraction anisotropy
control the range of force transmission.

For a contact force \(\mathbf F_s\) and local fiber tangent \(\mathbf t_s\), G4
v4 measures

\[
\mathbf F_{\parallel,s}=(\mathbf F_s\cdot\mathbf t_s)\mathbf t_s,
\qquad
\mathbf F_{\perp,s}=\mathbf F_s-\mathbf F_{\parallel,s}.
\]

The model does not give these components separate biological labels. It tests
the mechanical expectation that a larger transverse component is associated
with bending/rotation, whereas an axial component loads the fiber in tension
and can recruit connected tensile paths. If that association is absent in the
simulated data, the hypothesis fails.

## 4. Long-range force transmission is a network property

Ma and colleagues observed collagen compaction and alignment between nearby
cells over 24 hours. Image-derived fibrous finite-element models transmitted
stress farther than homogeneous linear-elastic and strain-hardening comparison
materials. Abhilash and Wang similarly found that aligned bundles can extend a
cell’s mechanical zone of influence.

This supports indirect motion in G4 v4: a fiber outside the direct cell-contact
band receives no Gaussian cell force. It can move only because permanent
crosslinks transmit force into its stretching and bending degrees of freedom.
The graph-distance classes—direct, one crosslink hop, two or more hops, and
unconnected—make that path testable.

The existing 180 µm square places an outer boundary about nine provisional cell
radii from the center. Because a finite anchor may influence response, G4 v4
does not assume this domain is converged. Larger domains must preserve fiber
density and be compared until near-cell metrics change by less than 5%.

## 5. Single-cell clocks and what Riching et al. actually calibrate

Riching and colleagues provide the primary experimental clock because they
studied MDA-MB-231 cells and explicitly reported imaging intervals. FITC-labeled
collagen and cells were imaged every four minutes for two hours to quantify
fiber displacement. Migration was imaged every ten minutes for six hours.
Aligned and random matrices were made in narrow and wide microchannels over
1–4 mg/ml collagen conditions.

The paper shows that alignment increased directional persistence and net
displacement without necessarily increasing migration speed. Cells in aligned
collagen had fewer, longer and more stable protrusions oriented with the
fibers. These results justify the G4 v4 output clocks and the later requirement
that contact geometry be read from local collagen rather than from a global
polarization value.

They do not provide a direct mapping from collagen concentration to this 2D
bead density, crosslink probability, effective clutch number, or Bell force.
Those remain assumptions.

## 6. Spheroid evidence belongs to a separate scale

Kopanska and colleagues studied CT26 spheroids. Collagen contraction began
immediately after seeding, invasion became evident on roughly 9–12-hour scales,
and experiments extended through 72 hours. Laser cutting and asymmetric matrix
cuts supported a role for collagen tension in invasion. Importantly, ECM
contraction and spheroid growth occurred together. The correct interpretation
is not that the spheroid shrank.

Geiger and colleagues later compared spheroids with radially and tangentially
oriented collagen and found a strong invasion bias toward radial organization.
These papers motivate a future 24–72-hour spheroid branch. They do not justify
reusing its collective-cell force, size or growth parameters in a single-cell
G4 v4 run.

## 7. Later feedback models are informative but not baseline evidence

Poonja and colleagues used a lattice-free multicellular agent model to simulate
transitions among TACS-like patterns. Their result is important because it
frames TACS changes as consequences of explicit cell–ECM interaction rules,
rather than static images. However, an agent-based fibril update rule is not a
validation of the current bead-spring force balance.

Böhringer and colleagues linked shell-wise collagen orientation and density to
cell contractility for hepatic stellate and glioblastoma cells. Their simulated
cells were represented by dipoles or quadrupoles, supporting a dipole control
but not making it the primary MDA-MB-231 condition.

Ley and colleagues introduced repeated tractor/retraction cycles, 3D
discrete-fiber mechanics, feedback, plasticity and spontaneous polarization.
It is the closest conceptual future comparison, but importing all of those
components now would prevent attribution of the first successful alignment
change to a specific added mechanism.

## What is established, unresolved and assumed

### Established strongly enough to test

- Initially random collagen can become locally aligned under strain or
  cell-generated contraction.
- Radial/perpendicular-to-boundary collagen is associated with invasive regions
  and persistent migration.
- Fiber bending/rotation can precede tensile recruitment.
- Aligned fibrous paths can transmit force farther than homogeneous comparison
  materials.
- Single-cell displacement, migration and spheroid remodeling require distinct
  observation clocks.

### Not resolved by one paper

- The quantitative mapping from effective clutches to collagen alignment.
- Which combination of clutch kinetics, adhesion maturation, proteolysis,
  crosslink turnover and cell shape produces TACS-2-to-TACS-3 transitions.
- A universal crosslink density or collagen modulus across preparations.
- Whether an isotropic circular single cell is sufficient for long-range radial
  alignment under experimentally realistic forces.

### Assumptions introduced in G4 v4

- A 10 µm fixed-radius circular MDA-MB-231 proxy. The declared 8–12 µm
  sensitivity range is implemented but has not yet been run as a full
  seed ensemble.
- Two-dimensional geometry with effective probabilistic intersection links.
- A 3 µm contact band and 1.5 µm Gaussian distance width inherited from the
  advisor-approved model design.
- Equal-load-sharing Bell-slip clusters as the primary site model, with an
  independent-clutch redundancy control; all constants are inherited and
  effective rather than molecularly calibrated.
- Contact relocation only after complete site failure. Its direction is read
  from the encountered local collagen tangent; this is a declared
  contact-guidance rule, not a fitted protrusion law.

## References

1. Provenzano PP, Eliceiri KW, Campbell JM, Inman DR, White JG, Keely PJ. “[Collagen reorganization at the tumor-stromal interface facilitates local invasion](https://doi.org/10.1186/1741-7015-4-38).” *BMC Medicine*. 2006;4:38.
2. Vader D, Kabla A, Weitz D, Mahadevan L. “[Strain-induced alignment in collagen gels](https://doi.org/10.1371/journal.pone.0005902).” *PLoS One*. 2009;4:e5902.
3. Conklin MW et al. “[Aligned collagen is a prognostic signature for survival in human breast carcinoma](https://doi.org/10.1016/j.ajpath.2010.11.076).” *American Journal of Pathology*. 2011;178:1221–1232.
4. Ma X et al. “[Fibers in the extracellular matrix enable long-range stress transmission between cells](https://doi.org/10.1016/j.bpj.2013.02.017).” *Biophysical Journal*. 2013;104:1410–1418.
5. Abhilash AS, Baker BM, Trappmann B, Chen CS, Shenoy VB. “[Remodeling of fibrous extracellular matrices by contractile cells](https://doi.org/10.1016/j.bpj.2014.08.029).” *Biophysical Journal*. 2014;107:1829–1840.
6. Wang H, Abhilash AS, Chen CS, Wells RG, Shenoy VB. “[Long-range force transmission in fibrous matrices enabled by tension-driven alignment of fibers](https://doi.org/10.1016/j.bpj.2014.09.044).” *Biophysical Journal*. 2014;107:2592–2603.
7. Riching KM et al. “[3D collagen alignment limits protrusions to enhance breast cancer cell persistence](https://doi.org/10.1016/j.bpj.2014.10.035).” *Biophysical Journal*. 2014;107:2546–2558.
8. Kopanska KS, Alcheikh Y, Staneva R, Vignjevic D, Betz T. “[Tensile forces originating from cancer spheroids facilitate tumor invasion](https://doi.org/10.1371/journal.pone.0156442).” *PLoS One*. 2016;11:e0156442.
9. Geiger F, Schnitzler LG, Brugger MS, Westerhausen C, Engelke H. “[Directed invasion of cancer cell spheroids inside 3D collagen matrices oriented by microfluidic flow in experiment and simulation](https://doi.org/10.1371/journal.pone.0264571).” *PLoS One*. 2022;17:e0264571.
10. Poonja S, Forero Pinto A, Lloyd MC, Damaghi M, Rejniak KA. “[Dynamics of fibril collagen remodeling by tumor cells](https://doi.org/10.3390/cells12232688).” *Cells*. 2023;12:2688.
11. Böhringer D et al. “[Fiber alignment in 3D collagen networks as a biophysical marker for cell contractility](https://doi.org/10.1016/j.matbio.2023.11.004).” *Matrix Biology*. 2023;122:159–177.
12. Ley AWY, Bersie-Larson LM, Adhikari S, Provenzano PP, Dorfman KD, Barocas VH. “[Cell–extracellular matrix feedback results in spontaneous cell polarization and heterogeneous remodeling in 3D isotropic and aligned discrete-fiber models](https://doi.org/10.1007/s12195-026-00922-0).” *Cellular and Molecular Bioengineering*. 2026.
13. Chan CE, Odde DJ. “[Traction dynamics of filopodia on compliant substrates](https://doi.org/10.1126/science.1163595).” *Science*. 2008;322:1687–1691. PMID 19074349. Abstract/metadata used; publisher full text was not treated as open access.
14. Bangasser BL, Rosenfeld SS, Odde DJ. “[Determinants of maximal force transmission in a motor-clutch model of cell traction in a compliant microenvironment](https://doi.org/10.1016/j.bpj.2013.06.027).” *Biophysical Journal*. 2013;105:581–592. PMID 23931306; PMCID PMC3736748.
15. Erdmann T, Schwarz US. “[Adhesion clusters under shared linear loading: a stochastic analysis](https://doi.org/10.1209/epl/i2003-10239-3).” *Europhysics Letters*. 2004;66:603–609. Open manuscript: [arXiv:cond-mat/0403552](https://arxiv.org/abs/cond-mat/0403552).
