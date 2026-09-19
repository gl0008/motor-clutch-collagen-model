# G5 biology audit — leader cells & (partial-)EMT vs our model's assumptions

> Source: deep-research pass (2026-09-18), 21 sources → 19 adversarially-verified claims (3-vote).
> Auditing the G5 R1–R3 modelling choices against the experimental biology. "Confidence" = how
> established the biology is (verify vote + source quality). Our own sim results stay personal testing.

## Our model's 4 assumptions being audited
- **A. Leader = higher per-cell motor traction** (R2 per-site F_stall).
- **B. Partial-EMT = reduced cell–cell adhesion ONLY** (R1; motors/integrins/motility unchanged).
- **C. Leader→follower coupling = a soft cell–cell adhesion spring + shared ECM** (no direct rope; each follower moves under its OWN traction).
- **D. Leader switching = per-cell energy drains with motor work → a fresher follower takes over** (R3).

---

## Q1 — What defines a leader cell; do they exert more traction / remodel collagen?
- Leaders **do exert more traction** than followers (K14+ vs K14− breast organoids, not a size artefact) [PMC10760762, 3-0]. Optical-tweezers tensile force **~0.9 nN** at the invasive front (RhoA-dependent) [science/adz4291, 2-1]. Leaders **actively remodel/align collagen** into parallel bundles/tracks, stiffening it (~3.4 kPa; remodeling reaches ~40 µm in invasive MDA-MB-231) [science/adz4291, 3-0].
- BUT a leader is **not just "more motor"**: it is a **molecular program** — basal-epithelial K14/p63/P-cadherin (K14 marks 94 % of leaders; K14 or p63 knockdown *alone* blocks collective invasion) [Cheung 2013 Cell, 3-0], plus elevated integrin-α2, MMP14, active-YAP, larger focal adhesions [PMC10760762; science/adz4291].
- **Audit of A: ✅ defensible but ⚠️ incomplete.** Leaders really do pull harder, so "leader = higher traction" is directionally right. But real leaders also differ in adhesion, integrins, MMP-proteolysis and active collagen-track building — our scalar per-site-stall captures only the traction axis. **Confidence: HIGH** that leaders pull more; **HIGH** that the definition is broader than traction.

## Q2 — Are followers DRAGGED, or independent? (your key question)
- **Followers are NOT simply roped along.** A single autonomously-activated leader can direct **at most ONE** follower; leading a *group* needs a **cluster-scale asymmetric traction/tension pattern**, i.e. followers contribute their own active traction [biorxiv 2024.01.23.576733, 3-0]. Collective velocity = **integral of active traction over the whole cluster** (supracellular force–velocity law) [same, 3-0].
- The coupling that DOES exist runs through **two channels**, both real:
  1. **Cadherin cell–cell junctions** — leaders **retain E-cadherin** junctions with followers during invasion [Cheung 2013, 2-1]; maintained epithelial junctions (Cdh1/Epcam) enable mechanical coupling [biorxiv 2025.04.04.647177, 3-0].
  2. **ECM tracks** — leaders build/align collagen microtracks that followers **reuse** [science/adz4291; fetch-verified MDA-MB-231/HT-1080]. (The claim that coupling is ECM-track *only* was **refuted** 0-3 — it is not the sole channel.)
- **Audit of C: ✅ this is our model's best-supported choice.** Biology matches our picture almost exactly: **no direct rope; followers move under their own traction; coupling via (i) a cadherin "spring" and (ii) shared/remodelled collagen.** We have both channels. ⚠️ Gaps: real coupling adds (a) a **supracellular coordinated tension pattern** (our cells are only pairwise-coupled), and (b) **active track-building/proteolysis** (our shared-ECM coupling is passive — cells deform but don't cut tracks). **Confidence: HIGH.** → **Direct answer to you: followers are independently motile, coordinated (not dragged) — exactly what we assume.**

## Q3 — Is leader switching real; is the energy-depletion handoff established?
- **Dynamic switching is real & well-established:** leaders are replaced by followers over time in MDA-MB-231 collagen invasion [PMC6475372, 3-0]; positions swap **every ~2–8 h** [Zhang/Konstantopoulos 2019 PNAS, 3-0]. Leader lifetime **~120–480 min**, shortening with collagen density (480→120 min as 1.5→6 mg/mL) [PMC6475372, weak vote 1-1].
- **Energy asymmetry is established:** leaders need **~50 % more energy**, higher glucose uptake/OXPHOS, and forward invasion **depletes** the leader's ATP [PNAS 2019, 3-0; PDAC nanofiber study, 3-0].
- **The CAUSAL "energy-depletion triggers the handoff" is DEBATED, not settled.** The direct causal claim verified only **1-2 (refuted)** / **1-1 (abstain)** here; and a **competing mechanism** exists — **Notch1–Dll4 lateral inhibition** drives leader turnover [ncomms7556] (its votes abstained under the session limit, so unresolved). The AICAR/glucose-starvation causal test [PMC6475372] was left unverified.
- **Audit of D: ✅ phenomenon + energy-asymmetry solid; ❓ causal energy-trigger is a hypothesis.** Our R3 "motor-work drains energy → follower takes over" reproduces the *observed* relay and matches the 120–480 min lifetime we calibrated to. But the *causal* energy mechanism is contested (vs Notch–Dll4 signalling). **We already labelled this a HYPOTHESIS — that was the correct call.** **Confidence: switching HIGH; energy-causality LOW/CONTESTED.**

## Q4 — Does partial-EMT change ONLY cell–cell adhesion?
- **No — partial-EMT changes MANY properties at once.** Hybrid E/M cells upregulate Vimentin + integrin-α2 (via an integrin-α2/TGFβ/Activin axis), gain **actin-rich protrusions + matrix-pulling + MMP-dependent collagen degradation**, WHILE **retaining** E-cadherin/Epcam junctions [biorxiv 2025.04.04.647177, 3-0]. So EMT co-varies **adhesion + traction + integrins + protrusions + proteolysis**.
- **Leader ≈ partial-EMT in some systems** (leaders *are* the partial-EMT cells) [biorxiv 2025, 3-0] — but in others leaders are **basal-EPITHELIAL** with NO classic EMT (no Twist/Slug/vimentin, keep E-cadherin) [Cheung 2013, 2-1]. So leader-role and EMT-state are **linked but system-dependent**, not cleanly independent.
- **Audit of B: ⚠️❌ our biggest oversimplification.** "Partial-EMT = lower cell–cell adhesion, motors unchanged" is biologically too narrow — real partial-EMT would *also* raise traction/integrins/protrusions/MMP. We did adhesion-only **deliberately for clean attribution** (isolate the adhesion effect), which is a defensible *controlled* first step, but the honest caveat is that a real EMT cell is not a pure adhesion knob. **Confidence: HIGH** that EMT is multi-property.

## Q5 — Is lowering cell–cell adhesion SUFFICIENT for single-cell escape?
- **No — not sufficient alone.** E-cadherin⁺ and E-cadherin⁻ tumours **both invade collectively** and metastasize [multiple, fetch-verified]. Switching collective→single-cell needs **cell–cell-adhesion loss + matrix proteolysis (MMP14) + confinement together** [SCC CTNNA1/MMP14 study; Ilina 2020 Nat Cell Biol jamming]. Adhesion mainly tunes **how** cells break off (stronger adhesion → fewer, larger multicellular ruptures), not whether.
- **Audit of our R1 result: ⚠️ oversimplified.** Our "low cc_adhesion → single-cell escape" captures the adhesion axis but **omits proteolysis + confinement**, which biology says are co-required. Our model has **no MMP/track-cutting** (cells push an elastic mesh), so R1 escape is only part of the mechanism — consistent with proteolysis being a DEFERRED item in our plan. **Confidence: HIGH** that adhesion alone is insufficient.

---

## Verdict table
| Our assumption | Biology says | Verdict | Fix for a faithful next version |
|---|---|---|---|
| A. leader = ↑ motor traction | leaders pull ~↑ (0.9 nN) **and** are a K14/p63 program w/ ↑integrin, MMP, YAP, collagen-track building | ✅ right axis, ⚠️ narrow | add integrin/focal-adhesion + MMP/track-building to "leader", not just stall |
| B. partial-EMT = adhesion-only | EMT changes adhesion **+ traction + integrins + protrusions + MMP** together | ⚠️❌ oversimplified (deliberate) | couple EMT to also raise traction/motility (a "coupled EMT-leader" variant) |
| C. coupling = adhesion spring + shared ECM, no rope, own traction | followers independently motile, coordinated by **cadherin junctions + ECM tracks**; not dragged | ✅ **best-supported choice** | add supracellular tension coordination + active track/proteolysis |
| D. switching = motor-work energy drain → handoff | switching + energy-asymmetry **real**; **causal** energy-trigger **contested** (vs Notch-Dll4) | ✅ pheno / ❓ mechanism (already a HYPOTHESIS) | keep labelled hypothesis; optionally add a Notch-Dll4 signalling variant as an alternative |

## Bottom line
- **Best-validated choice:** our **leader→follower coupling** (C) — biology directly supports "followers move under their own traction, coupled by cadherin junctions + shared collagen, not roped." Your intuition to check this paid off: the model is right here.
- **Directionally right, incomplete:** leaders pull more (A); switching + energy-asymmetry are real (D, with the causal trigger correctly flagged as a hypothesis).
- **Biggest simplifications to flag to advisors (both already deferred in our plan):** (1) partial-EMT is multi-property, not adhesion-only (B); (2) collective→single-cell needs **proteolysis + confinement**, not adhesion alone (Q5). Neither invalidates the staged attribution — they are the honest "what's missing" for the next version.

### Key sources
Cheung 2013 Cell (K14/p63 leaders) · Zhang/Konstantopoulos 2019 PNAS (energetic leader-follower, 2–8 h switching) · PMC6475372 (leader lifetime 120–480 min, energy handoff — our R3 calibration target) · science/adz4291 2025 (0.9 nN, collagen tracks, MMP14/YAP) · biorxiv 2025.04.04.647177 (leaders = partial-EMT, retain E-cadherin) · biorxiv 2024.01.23.576733 (supracellular force–velocity, 1 leader → 1 follower) · Ilina 2020 Nat Cell Biol (adhesion + confinement jamming) · ncomms7556 (Notch-Dll4 leader turnover — competing switch mechanism) · PMC10760762 (K14 leader TFM traction).
