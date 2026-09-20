# G4 v4F evidence ledger

## Riching et al. (2014)

- **Authors:** Kristin M. Riching et al.
- **Title:** 3D Collagen Alignment Limits Protrusions to Enhance Breast Cancer Cell Persistence
- **Identifiers:** DOI 10.1016/j.bpj.2014.10.035; PMID 25468334;
  PMCID PMC4255204
- **Access:** Open full text
- **Used for:** MDA-MB-231 cells migrate more persistently along aligned
  collagen; six-hour migration and two-hour collagen-displacement observation
  clocks. Their computational model also increased the probability that a
  protrusion vector maintained its orientation across timesteps as matrix
  alignment increased.
- **Does not support:** The exact front-memory time, Gaussian probability or
  linear steric penalty.
- **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4255204/

## Carey et al. (2016)

- **Authors:** Shawn P. Carey, Zachary E. Goldblatt, Karen E. Martin,
  Bethsabe Romero, Rebecca M. Williams and Cynthia A. Reinhart-King
- **Title:** Local extracellular matrix alignment directs cellular protrusion
  dynamics and migration through Rac1 and FAK
- **Identifiers:** DOI 10.1039/C6IB00030D; PMID 27384462; PMCID PMC4980151
- **Access:** Open full text
- **Used for:** MDA-MB-231 probing protrusions near aligned collagen persist
  and lengthen along the alignment axis; locally remodeled anisotropy can
  precede directional spreading.
- **Does not support:** A fixed 120 s memory clock, the numerical guidance
  gains or forced rear detachment.
- **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4980151/

## Ray et al. (2017)

- **Authors:** Arja Ray et al.
- **Title:** Anisotropic forces from spatially constrained focal adhesions
  mediate contact guidance directed cell migration
- **Identifiers:** DOI 10.1038/ncomms14923; PMID 28401884; PMCID PMC5394287
- **Access:** Open full text
- **Used for:** Aligned physical cues constrain protrusion orientation and
  increase directional persistence in MDA-MB-231 cells.
- **Does not support:** The exact fiber-search or memory probability used here.
- **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5394287/

## Koch et al. (2019)

- **Authors:** Mar Cóndor, Christoph Mark, Richard C. Gerum, Nadine C. Grummel,
  Andreas Bauer, José M. García-Aznar and Ben Fabry
- **Title:** Breast Cancer Cells Adapt Contractile Forces to Overcome Steric
  Hindrance
- **Identifiers:** DOI 10.1016/j.bpj.2019.02.029; PMID 30902366;
  PMCID PMC6451061
- **Access:** Open full text
- **Used for:** Steric hindrance from collagen pore size is a real mechanical
  constraint for MDA-MB-231 invasion and changes cell contractile behavior.
- **Does not support:** Treating a cell as a rigid circle or calibrating the
  present linear overlap penalty.
- **Link:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6451061/

## Chan and Odde (2008); motor–clutch lineage

- **Authors:** Clarence E. Chan and David J. Odde
- **Title:** Traction dynamics of filopodia on compliant substrates
- **Identifiers:** DOI 10.1126/science.1163595; PMID 19074349
- **Access:** Abstract and bibliographic record; full text was not used here.
- **Used for:** Dynamic clutch binding, load accumulation and force-dependent
  failure.
- **Does not support:** Collagen-specific contact search or protrusion memory.
- **Link:** https://pubmed.ncbi.nlm.nih.gov/19074349/

## Evidence boundary for the new v4F corrections

- **Two-way steric force:** Koch et al. support the biological importance of
  steric hindrance in 3D collagen. Equal-and-opposite force transfer is a
  mechanics requirement of this coupled model. The linear overlap law and its
  stiffness are not calibrated by Koch et al.
- **Dynamic contact search:** no paper above specifies a 5 s search clock or
  the exact candidate-selection algorithm. Updating candidates around the
  current cell and deformed fibers is required for geometric consistency; the
  5 s interval and Gaussian distance rule are explicit model assumptions.
- **E1 matrix cue:** Riching, Carey and Ray directly support aligned ECM as a
  directional cue, but not the exact 60° one-sided tract construction.
- **E2 protrusion memory:** Carey directly supports longer-lived, longer
  protrusions along aligned collagen. Riching supports restricted protrusions,
  increased persistence, and used an alignment-dependent protrusion-direction
  persistence probability in its own model. The present 120 s maturation
  threshold and exponential probability gains remain assumptions.
