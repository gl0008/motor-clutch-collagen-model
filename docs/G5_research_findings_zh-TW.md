# G5 文獻綜整 — 纖維狀 collagen 中的 organoid 侵襲

> Deep-research 產出（2026-09-02），經對抗式驗證（25 項主張中 24 項確認、1 項被駁回）。
> 供 `docs/G5_organoid_plan.md` 使用。每項發現皆記錄 confidence 與 vote。

## 需納入模型的 physics / 性質（依重要性排序）

1. **Force-induced radial（TACS-3-like）alignment 是 collective invasion 的一階（1st-order）
   因果驅動因子**（high, 3-0）。在*同一顆* spheroid 周圍，radially aligned fibers → 比
   circumferential 產生更多擴散出去的 clusters。其產生機制為：tangential swirling → shear →
   negative normal stress → 向遠端傳播的 radial contractile stress。
   *Saraswathibhatla 2025 (bioRxiv 2025.01.31.635980); Su/Kim 2021 (Biomaterials, S0142961221002787);
   Botticelli 2025 (bioRxiv 2025.09.02.673661); Ray 2022 (PMC9033577)。參見 Provenzano TACS-3 / Conklin 2011。*

2. **Nonlinear fiber mechanics（strain-stiffening + buckling + tension-driven alignment）→
   long-range force transmission 約 8–10 個 cell diameters（可達約 200 µm）**，遠超過
   linear-elastic 的預測（high, 3-0 / 其中一項 2-1）。此設定了 remodeling halo 的空間尺度。
   *Shenoy/Abhilash/Wang 2014 (arXiv 1409.6216); Tsingos 2023 (Biophys J S0006-3495(23)00325-9);
   Aghvami 2016 (PMC5018119); Mark 2020 (eLife 51912)。*
   ⚠ 被駁回（1-2）：認為*單靠 discreteness*（與 strain-stiffening 無關）就能造成 long-range
   transmission 的主張。不要在缺乏 nonlinearity 的情況下把它歸因於 fibrous discreteness。

3. **Matrix plasticity / stress-relaxation 是一個獨立、可調控的一階決定因子，主導
   protease-independent 的 confined migration**（high, 3-0）。細胞能在 pore 小於 40 nm 的 matrix 中，
   以 plastic 方式撐開 2–3 µm 寬、超過 50 µm 長的 channels；高 plasticity（100 Pa 下約 30%）→
   migration 約為低 plasticity（約 10%）的 5 倍。
   *Wisdom 2018 (Nat Commun s41467-018-06641-z)。*

4. **Stiffening 具有方向性：LOCAL stiffening 促進 collective invasion；GLOBAL/systemic
   stiffening 抑制之**，且與 fiber thickness / pore size 無關（high, 3-0）。應將 stiffness 建模為
   heterogeneous field，而非單一 global 旋鈕。
   *Ray 2022 (Oncogene, PMC9033577)。*

## 他人如何建模

5. **Paradigm 1 — agent/cell core + DISCRETE bead-spring fiber network**（high, 3-0）。Coupling 方案：
   Reinhart-King「hand-over-hand」dynamic repetitive traction（Reinhardt 2018, ASME JBME）；Tsingos
   2023 hybrid Cellular-Potts + bead-spring，focal-adhesion beads 在 1.25 µm band 內 clamp 到 lattice；
   Zhang/Schwarz 2024 以 deformable polyhedra 的 vertex model + fiber linker springs（arXiv 2403.16784）。
   Fiber network = Hookean（intra-fiber）+ crosslinker + angular bending springs，overdamped Langevin；
   strain-stiffening E~[N_fibers]^1.09。
   → **我們的 g2 engine 正是這個 paradigm**（bead-spring + crosslinks + bending，overdamped）。

6. **Paradigm 2 — continuum/voxel ECM field + agent cells**（high, 3-0）。PhysiCell hybrid：ECM voxels
   儲存 fiber orientation / anisotropy / density + contact-guidance rules。在 organoid 尺度比
   bead-spring 更輕量，但會失去 discrete-fiber 的細節。*Botticelli 2025。*

## 量化驗證目標

7. **Alignment（時間 + 空間）：** radial alignment index 在 Day1→3 上升、Day3→5 趨於 plateau；
   在 0–50 µm border 最高，並持續超過 100 µm；>0.2 = anisotropic；non-invasive organoids 無此現象。
   *Saraswathibhatla 2025; Ray 2022。*

8. **Displacement scaling：** radial vs tangential 指數 **n_r≈1.2 vs n_t≈2.2**；radial 指數隨
   stiffening 由 **α=−2（linear，約 1 Pa）→ −0.2（>1000 Pa，stiffened）** 位移。
   *Saraswathibhatla 2025; Mark 2020。*

9. **Force / deformation / stiffening：** TNBC spheroid（約 4000 cells）→ 24 h 內
   **344±35 µN / 677±68 Pa**；**約 200 µm 的 surface deformation（>50% strain）**；
   **>20× local stiffening**（κ0=1645 Pa）。
   *Mark 2020。*（µN 為 spheroid-aggregate 層級；Pa 為 pressure 而非 modulus；TNBC 非 PDAC。）

10. **最佳 ablation handles：**（a）同一顆 spheroid 周圍的 radial vs circumferential fiber
    orientation（Su/Kim 2021）；（b）high vs low plasticity IPN → 5× migration（Wisdom 2018）。

## 注意事項（Caveats）
- Force/pressure/stiffening 數量級（344 µN、677 Pa、>20×）與 alignment dynamics 皆來自 **breast**
  （TNBC/PyMT），應視為 cross-tissue anchors，而非 PDAC-exact。
- Power-law 指數（n_r, n_t；n=0.35–1.09）是特定 2D geometry 下的 **model outputs** —
  重現它們僅驗證 model-vs-model；請優先 anchor 到實測 metrics。
- Plasticity（Wisdom 2018）與 voxel-continuum（Botticelli 2025）各自僅倚賴單一 primary source
  （兩者皆 3-0 verified）。兩篇關鍵 refs 為 preprints（Saraswathibhatla、Botticelli）。
- Cell–cell adhesion / jamming 與 proliferation-pressure 雖被標記為一階，但 **沒有任何存活下來的
  主張提供 quantitative jamming 或 proliferation 目標** — 屬 open questions。

## 主要來源
Saraswathibhatla 2025 bioRxiv · Su/Kim 2021 Biomaterials · Ray 2022 Oncogene · Wisdom 2018 Nat Commun ·
Mark 2020 eLife · Shenoy/Abhilash 2014 · Tsingos 2023 Biophys J · Zhang/Schwarz 2024 arXiv ·
Reinhardt/Reinhart-King 2018 ASME JBME · Botticelli 2025 bioRxiv (PhysiCell) · Nam 2016/2017 ·
Ilina & Friedl 2020 Nat Cell Biol (jamming)。
