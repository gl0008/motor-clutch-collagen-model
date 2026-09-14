# G4 v4 paper-by-paper evidence ledger

This ledger separates observations, mechanics benchmarks and model rules. A
paper supports only the statement listed under **Used for**. Missing details are
marked unknown rather than inferred. No publisher PDF is redistributed.

## TACS observations and tumor-boundary orientation

### Provenzano et al. (2006)

- **Exact citation:** Paolo P. Provenzano, Kevin W. Eliceiri, Jay M. Campbell,
  David R. Inman, John G. White and Patricia J. Keely, “Collagen
  reorganization at the tumor-stromal interface facilitates local invasion,”
  *BMC Medicine* 4:38 (2006). DOI
  [10.1186/1741-7015-4-38](https://doi.org/10.1186/1741-7015-4-38);
  PMID 17190588; PMCID PMC1781458. **Open full text, CC BY.**
- **System / matrix / geometry / clock:** mouse mammary tissue and primary
  tumor explants in 3D collagen; explant example at 8 h; random starting gel in
  the explant experiment.
- **Driver / measurement:** tumor-cell activity; MPE/SHG imaging and collagen
  orientation relative to the tumor boundary.
- **Used for:** random-to-radial alignment is a biologically observed target;
  radial and boundary-parallel orientations must be distinct.
- **Not used for:** force magnitude, clutch kinetics or cell shrinkage.

### Conklin et al. (2011)

- **Exact citation:** Matthew W. Conklin, Jens C. Eickhoff, Kristin M. Riching,
  Carolyn A. Pehlke, Kevin W. Eliceiri, Paolo P. Provenzano, Andreas Friedl and
  Patricia J. Keely, “Aligned collagen is a prognostic signature for survival
  in human breast carcinoma,” *American Journal of Pathology* 178:1221–1232
  (2011). DOI [10.1016/j.ajpath.2010.11.076](https://doi.org/10.1016/j.ajpath.2010.11.076);
  PMID 21356373; PMCID PMC3070581. **Free PMC full text.**
- **System / matrix / geometry / clock:** 196 human breast-cancer specimens;
  native tumor collagen; clinical follow-up rather than a remodeling movie;
  starting orientation is not an experimental initial condition.
- **Driver / measurement:** formation mechanism not tested; SHG-based TACS-3
  scoring relative to the boundary.
- **Used for:** perpendicular-to-boundary collagen is a clinically meaningful
  output.
- **Not used for:** a causal alignment mechanism or numerical parameters.

## Strain-induced collagen alignment

### Vader et al. (2009)

- **Exact citation:** David Vader, Alexandre Kabla, David Weitz and
  Lakshminarayana Mahadevan, “Strain-induced alignment in collagen gels,”
  *PLoS ONE* 4:e5902 (2009). DOI
  [10.1371/journal.pone.0005902](https://doi.org/10.1371/journal.pone.0005902);
  PMID 19529768; PMCID PMC2691583. **Open full text, CC BY.**
- **System / matrix / geometry / clock:** 1 mg/ml collagen gels, untreated or
  glutaraldehyde-crosslinked, mechanically strained in a two-point geometry;
  cell images at 10–48 h are illustrative, while the alignment test is a strain
  cycle; initially disordered collagen.
- **Driver / measurement:** imposed strain; confocal orientation histograms,
  order parameter and density.
- **Used for:** an elastic network can align under load; plasticity is not
  required for the first alignment test.
- **Not used for:** tumor-cell migration, crosslink probability or clutch
  kinetics.

## Discrete-fiber bending-to-stretching transition

### Abhilash et al. (2014)

- **Exact citation:** A. S. Abhilash, Brendon M. Baker, Britta Trappmann,
  Christopher S. Chen and Vivek B. Shenoy, “Remodeling of fibrous extracellular
  matrices by contractile cells: predictions from discrete fiber network
  simulations,” *Biophysical Journal* 107:1829–1840 (2014). DOI
  [10.1016/j.bpj.2014.08.029](https://doi.org/10.1016/j.bpj.2014.08.029);
  PMID 25418164; PMCID PMC4213674. **Free PMC full text.**
- **System / matrix / geometry / clock:** computational random discrete-fiber
  network around circular/elliptical contracting inclusions; mechanical load
  progression, not a biological observation duration.
- **Driver / measurement:** prescribed contraction; fiber orientation,
  stretch/bend energy and zone of influence.
- **Used for:** the separately labelled cavity benchmark and the
  bending-to-stretching recruitment diagnostic.
- **Not used for:** biological cell shrinkage or MDA-MB-231 traction values.

### Wang et al. (2014)

- **Exact citation:** Hailong Wang, A. S. Abhilash, Christopher S. Chen,
  Rebecca G. Wells and Vivek B. Shenoy, “Long-range force transmission in
  fibrous matrices enabled by tension-driven alignment of fibers,”
  *Biophysical Journal* 107:2592–2603 (2014). DOI
  [10.1016/j.bpj.2014.09.044](https://doi.org/10.1016/j.bpj.2014.09.044);
  PMID 25468338; PMCID PMC4255175. **Free PMC full text.**
- **System / matrix / geometry / clock:** constitutive fibrous-matrix model with
  circular/elliptical contraction controls; quasi-static mechanics; model
  orientation prescribed initially through the constitutive distribution.
- **Driver / measurement:** tensile contraction, cell shape and contraction
  anisotropy; alignment range and transmitted stress.
- **Used for:** tension-driven alignment and a force-dipole comparison control.
- **Not used for:** this project’s discrete contact or clutch rules.

## Single-cell traction and long-range force transmission

### Ma et al. (2013)

- **Exact citation:** Xiaoyue Ma, Maureen E. Schickel, Mark D. Stevenson,
  Alisha L. Sarang-Sieminski, Keith J. Gooch, Samir N. Ghadiali and Richard T.
  Hart, “Fibers in the extracellular matrix enable long-range stress
  transmission between cells,” *Biophysical Journal* 104:1410–1418 (2013).
  DOI [10.1016/j.bpj.2013.02.017](https://doi.org/10.1016/j.bpj.2013.02.017);
  PMID 23561517; PMCID PMC3617419. **Free PMC full text.**
- **System / matrix / geometry / clock:** fibroblasts in type-I collagen;
  image-derived fibrous finite-element geometry; observations at 4 and 24 h;
  starting matrix not experimentally prealigned.
- **Driver / measurement:** cell compaction/contraction; collagen alignment
  between cell pairs and computed stress propagation.
- **Used for:** indirect response can arise through a connected fibrous path.
- **Not used for:** MDA-MB-231 kinetics or a 2D intersection-link probability.

### Böhringer et al. (2023)

- **Exact citation:** David Böhringer, Andreas Bauer, Ivana Moravec, Lars
  Bischof, Delf Kah, Christoph Mark, Thomas J. Grundy, Ekkehard Görlach,
  Geraldine M. O'Neill, Silvia Budday, Pamela L. Strissel, Reiner Strick,
  Andrea Malandrino, Richard Gerum, Michael Mak, Martin Rausch and Ben Fabry,
  “Fiber alignment in 3D collagen networks as a biophysical marker for cell
  contractility,” *Matrix Biology* 122:159–177 (2023). DOI
  [10.1016/j.matbio.2023.11.004](https://doi.org/10.1016/j.matbio.2023.11.004);
  PMID 37967726; PMCID PMC10872942. **Free PMC manuscript; publisher copyright
  retained.**
- **System / matrix / geometry / clock:** hepatic stellate cells and
  glioblastoma lines in 1.0–1.2 mg/ml 3D collagen; 1–8 d; random cell-free gel
  is the orientation control.
- **Driver / measurement:** cell traction; shell-wise structure-tensor
  orientation/density, matrix displacement and traction reconstruction;
  dipole/quadrupole simulations used 5 nN point forces.
- **Used for:** shell-wise orientation and the secondary dipole control.
- **Not used for:** MDA-MB-231 calibration or isotropic cell shrinkage.

## Spheroid-scale contraction, growth and invasion

### Kopanska et al. (2016)

- **Exact citation:** Katarzyna S. Kopanska, Yara Alcheikh, Ralitza Staneva,
  Danijela Vignjevic and Timo Betz, “Tensile forces originating from cancer
  spheroids facilitate tumor invasion,” *PLoS ONE* 11:e0156442 (2016). DOI
  [10.1371/journal.pone.0156442](https://doi.org/10.1371/journal.pone.0156442);
  PMID 27271249; PMCID PMC4896628. **Open full text, CC BY.**
- **System / matrix / geometry / clock:** CT26 spheroids in 3D collagen I;
  immediate collagen contraction, invasion around 9–12 h, measurements to 72 h;
  initial matrix described before spheroid-generated reorganization.
- **Driver / measurement:** spheroid traction plus growth; collagen flow,
  orientation, cuts and laser ablation.
- **Used for:** a later separately calibrated spheroid track and the distinction
  between inward matrix motion and cell/spheroid shrinkage.
- **Not used for:** single-cell force, size or six-hour migration parameters.

### Geiger et al. (2022)

- **Exact citation:** Florian Geiger, Lukas G. Schnitzler, Manuel S. Brugger,
  Christoph Westerhausen and Hanna Engelke, “Directed invasion of cancer cell
  spheroids inside 3D collagen matrices oriented by microfluidic flow in
  experiment and simulation,” *PLoS ONE* 17:e0264571 (2022). DOI
  [10.1371/journal.pone.0264571](https://doi.org/10.1371/journal.pone.0264571);
  PMID 35231060; PMCID PMC8887745. **Open full text, CC BY.**
- **System / matrix / geometry / clock:** spheroids in microfluidically
  prealigned radial or tangential 3D collagen; 1–3 d.
- **Driver / measurement:** alignment was imposed by flow, not formed by the
  spheroid; invasion front relative to fiber orientation.
- **Used for:** radial and tangential collagen have different guidance effects.
- **Not used for:** random-to-aligned mechanics or Brownian motion (explicitly
  excluded from G4 v4).

## Contact guidance and migration

### Riching et al. (2014)

- **Exact citation:** Kristin M. Riching, Benjamin L. Cox, Max R. Salick,
  Carolyn Pehlke, Andrew S. Riching, Susan M. Ponik, Benjamin R. Bass, Wendy C.
  Crone, Yi Jiang, Alissa M. Weaver, Kevin W. Eliceiri and Patricia J. Keely,
  “3D collagen alignment limits protrusions to enhance breast cancer cell
  persistence,” *Biophysical Journal* 107:2546–2558 (2014). DOI
  [10.1016/j.bpj.2014.10.035](https://doi.org/10.1016/j.bpj.2014.10.035);
  PMID 25468334; PMCID PMC4255204. **Free PMC full text and supplementary
  movies.**
- **System / matrix / geometry / clock:** MDA-MB-231 single cells in 1–4 mg/ml
  collagen, random/wide and flow-aligned/narrow microchannels; collagen
  displacement every 4 min for 2 h; migration every 10 min for 6 h; invasion at
  3 d.
- **Driver / measurement:** alignment in the comparison gels was prepared
  mechanically; FITC-collagen displacement, cell tracks, speed/persistence and
  protrusion number/length/orientation.
- **Used for:** primary cell type, observation clocks, and local collagen
  contact guidance for new contacts.
- **Not used for:** effective clutch count, Bell force, bead density or
  crosslink probability.

## Motor–clutch loading and shared-load adhesion

### Chan and Odde (2008)

- **Exact citation:** Clarence E. Chan and David J. Odde, “Traction dynamics of
  filopodia on compliant substrates,” *Science* 322:1687–1691 (2008). DOI
  [10.1126/science.1163595](https://doi.org/10.1126/science.1163595);
  PMID 19074349. **Abstract/metadata used; publisher full text not treated as
  open access.**
- **System / matrix / geometry / clock:** chick forebrain growth-cone filopodia
  on compliant 2D substrates; nanoscale actin/traction dynamics; no collagen
  network and no hour-scale tumor observation.
- **Driver / measurement:** myosin-driven retrograde actin flow coupled by
  stochastic clutches; load-and-fail versus frictional slippage.
- **Used for:** continued load after contact comes from relative actin motion,
  and ends at tracking, stall or failure.
- **Not used for:** any MDA-MB-231 or collagen parameter.

### Bangasser, Rosenfeld and Odde (2013)

- **Exact citation:** Benjamin L. Bangasser, Steven S. Rosenfeld and David J.
  Odde, “Determinants of maximal force transmission in a motor-clutch model of
  cell traction in a compliant microenvironment,” *Biophysical Journal*
  105:581–592 (2013). DOI
  [10.1016/j.bpj.2013.06.027](https://doi.org/10.1016/j.bpj.2013.06.027);
  PMID 23931306; PMCID PMC3736748. **Free PMC full text.**
- **System / matrix / geometry / clock:** theoretical motor–clutch model on a
  compliant substrate; model cycle time, not an MDA-MB-231 experiment; no
  collagen orientation state.
- **Driver / measurement:** substrate stiffness, motor/clutch number and
  stochastic loading; traction and retrograde flow.
- **Used for:** load-and-fail cycling as an explicit kinetic block.
- **Not used for:** numerical calibration of this 3D-collagen problem.

### Erdmann and Schwarz (2004)

- **Exact citation:** Thorsten Erdmann and Ulrich S. Schwarz, “Adhesion clusters
  under shared linear loading: a stochastic analysis,” *Europhysics Letters*
  66:603–609 (2004). DOI
  [10.1209/epl/i2003-10239-3](https://doi.org/10.1209/epl/i2003-10239-3);
  arXiv cond-mat/0403552. **Open author manuscript on arXiv; no PMID.**
- **System / matrix / geometry / clock:** theoretical cluster of receptor–ligand
  bonds between surfaces under linearly increasing shared load; no cell type,
  collagen preparation or biological observation duration.
- **Driver / measurement:** equal load sharing, Bell-type dissociation,
  rebinding and stochastic complete cluster failure.
- **Used for:** the shared-load comparison in which decreasing bound count
  increases remaining-bond force and hazard.
- **Not used for:** integrin catch bonds, focal-adhesion maturation or the
  current effective parameter values.

## Integrated feedback and plasticity models

### Poonja et al. (2023)

- **Exact citation:** Sharan Poonja, Ana Forero Pinto, Mark C. Lloyd, Mehdi
  Damaghi and Katarzyna A. Rejniak, “Dynamics of fibril collagen remodeling by
  tumor cells: a model of tumor-associated collagen signatures,” *Cells*
  12:2688 (2023). DOI [10.3390/cells12232688](https://doi.org/10.3390/cells12232688);
  PMID 38067116; PMCID PMC10705683. **Open full text, CC BY.**
- **System / matrix / geometry / clock:** lattice-free multicellular tumor/ECM
  agent model over simulated days; initial patterns are prescribed per TACS
  experiment rather than generated by bead mechanics.
- **Driver / measurement:** agent interaction, movement, compliance and
  feedback rules; TACS-like orientation patterns.
- **Used for:** later cumulative TACS feedback comparisons.
- **Not used for:** validation of bead springs, crosslink forces or clutches.

### Ley et al. (2026)

- **Exact citation:** A. W. Y. Ley, L. M. Bersie-Larson, S. Adhikari, P. P.
  Provenzano, K. D. Dorfman and V. H. Barocas, “Cell–extracellular matrix
  feedback results in spontaneous cell polarization and heterogeneous
  remodeling in 3D isotropic and aligned discrete-fiber models,” *Cellular and
  Molecular Bioengineering* (2026). DOI
  [10.1007/s12195-026-00922-0](https://doi.org/10.1007/s12195-026-00922-0).
  **Metadata/abstract-level access only when reviewed; no quantitative values
  imported.**
- **System / matrix / geometry / clock:** 3D isotropic/aligned discrete-fiber
  model with tractors, retraction, feedback and plasticity; model-dependent
  cycles.
- **Driver / measurement:** repeated cell–ECM feedback; spontaneous
  polarization and heterogeneous remodeling.
- **Used for:** a later comparator.
- **Not used for:** permission to add feedback, plasticity or its parameters to
  this clean elastic baseline.

## One-line evidence boundary

These studies jointly support the *sequence of testable blocks*. They do not
jointly calibrate one unified model. Every G4 v4 stiffness, drag, contact,
crosslink and clutch value without a source-matched measurement remains an
explicit assumption.
