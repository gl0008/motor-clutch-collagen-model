# 第五代 (G5) — 腫瘤 Organoid 侵入纖維狀 Collagen

> **狀態:** 提案計畫(尚未建置)。撰寫於 2026-09-02。
> **目標:** 將已驗證的單細胞 motor-clutch + bead-spring-collagen 核心,擴大為
> 多細胞腫瘤 **organoid**(~150–250 µm),以研究 cell–cell + cell–ECM 交互作用
> 如何決定侵入模式(受限 confined / 集體 collective / 單細胞逃逸 single-cell escape)。
> 以下所有數值皆附來源;自行模擬的輸出屬於*個人測試*,而非已確認的發現
> (CLAUDE.md §7.5)。文獻整合見 `docs/G5_research_findings.md`。
>
> **2026-09-02 鎖定的決策 (Minnie):** (1) 細胞模型 = **Option A**,N 個 motor-clutch disks +
> cell–cell adhesion(§2);(2) 第一個里程碑 = **Stage B** — 收縮性 pull → radial alignment,
> 即優先驗證的目標(§7)。Stages C–F 僅在 B 驗證後才進行。cell–cell adhesion 項的顧問簽核
> (§2、§10 Q1)在建置前仍待確認。

---

## 0. 目標現象(來自 Kolade/Gloria 影像 + 文獻)

實驗影片(4MOSC1 OSCC 與 PDAC organoids 於重構 collagen 中,confocal
reflection + collagen probe;`memory/reference_microscopy_scale_pptx.md`)顯示,於 Day2→Day5:

1. **Organoid** 直徑 ~200–250 µm,個別細胞 ~15–20 µm,嵌於一個無序的
   collagen 網絡中。
2. **放射狀 leader/protrusion 束**向外伸展,抓住 collagen 並將其向內拉。
3. 周圍纖維從無序重組為 **radial alignment(aster / TACS-3-like)**,
   organoid 附近伴隨 **clearance zones(清空區)**。
4. 大多為**集體 (collective)** 推進,偶有**單細胞逃逸 (single-cell escape)**。

G5 必須先重現 (2)+(3)(符合 Hongbo 的最小判準:pull → alignment),然後才是 (4)。
這對應到專案假說 H1(alignment ↔ 侵入模式)與 H3(alignment × EMT →
collective vs escape);CLAUDE.md §1。

---

## 1. 基礎決策 — 重用 g2 核心、擷取 g4_v2 部件、撰寫一個全新的多細胞 driver

程式碼審查(完整報告 `docs/G5_code_assessment.md`)確認**所有世代已經
透過 subclassing 共用同一個物理核心**:

```
CollagenConfig / Network  (generations/g2_corrected/common/model.py, 已驗證, 92 tests)
   ├── SpheroidConfig      → g3   (單一大型 rigid ball)
   └── G4Config → G4V2Config → g4_v2 (Numba integrator + Bell/shared-load clutch + rigid-body)
```

**關鍵阻礙:** 細胞在各處都被寫死為一個*單一 rigid circle*(`center`、
`cell_radius`、`theta` 皆為純量;`repulsion_forces` 假設只有一個 center;clutch state 陣列僅為
一個細胞的量)。**整個 repo 沒有任何多細胞支援。** g3 的「spheroid」仍然是
一顆大球,而非細胞的集合體。

**決策:**

| 選擇 | 結論 |
|--------|--------|
| ECM 核心 | **原封不動重用 g2 `Network` / `elastic_forces`** — 唯一已驗證的 collagen 引擎。 |
| Numba integrator + Bell / shared-load clutch kinetics | **從 `g4_v2_multiscale/model.py` 擷取,作為已驗證的積木**(`_advance_numba`、`bell_off_rate`、shared-load hazard)。 |
| 多細胞 organoid 迴圈 | **全新撰寫**於 `generations/g5_organoid/` — 不要延伸 g4_v2 那 140 行的 `run_clutch` 巨型函式。 |
| 無纖維空隙 (fibre-free-gap) 網絡生成器 + protrusion 概念 | **從 g3 移植**(最接近「大型物體帶 t=0 空隙」的先例)。 |

理由:保留已驗證的 collagen 物理、避免繼承使用者指出的 g4_v2 亂象,
並隔離真正新的工作(M 個細胞、cell–cell 接觸、neighbor lists)。

---

## 2. 細胞表示 — organoid = N 個 motor-clutch disks 的叢集 + cell–cell adhesion

**建議 (Option A):** 將 organoid 表示為 **N 個離散 rigid disks**(半徑 ~8–10 µm),
每個都是一個 motor-clutch 細胞,透過既有的 Gaussian kernel 抓住鄰近纖維,並加上
disks 之間的 **cell–cell adhesion/repulsion** 交互作用(一個軟性的吸引–排斥位能,
即 tissue surface tension 組織表面張力)。細胞在 overdamped rigid-body dynamics 下移動
(重用 g4_v2 的 released-cell 分支,推廣到 M 個 centers)。

為何這是正確的首選:
- **相對現有程式碼改動最小**:disk + Gaussian-kernel + clutch 機制已存在
  且已驗證;我們只是把純量推廣成跨 M 個細胞的陣列。
- **直接產生目標科學**:collective 推進 vs 單細胞逃逸是 cell–cell adhesion(維繫叢集)
  與 cell–ECM traction(把細胞拉出)兩者平衡下的*湧現 (emergent)* 結果 —
  正是 H3 軸線。Ilina & Friedl 2020 (Nat Cell Biol) 將 cadherin 基礎的
  cell–cell adhesion + matrix confinement 列為 jamming/collective ↔ single-cell 轉變的
  兩個一階決定因素 — 兩者皆能自然地由(adhesion 位能)×(fiber network)捕捉。
- **可擴展**:proliferation(增加 disks)、EMT(每細胞較低 adhesion / 較高 motility)、
  以及 leader/follower 異質性,日後都能變成 per-cell 參數。

已考量並延後的替代方案:
- *單一大型可變形物體 / active-gel* — 失去個別細胞逃逸的現象;無法
  表現單細胞散播 (dissemination)。延後。
- *Vertex / Cellular-Potts 細胞*(Tsingos 2023 CPM+bead-spring;Zhang/Schwarz 2024 vertex+fiber
  linkers)— 更豐富的細胞形狀與真實 confluent 組織力學,但需大幅重寫,對第一個問題而言
  太過度。保留為 Stage-F 升級,若 disk-cluster 被證明太粗糙。

> **顧問確認 (Hongbo/Kolade):** cell–cell 交互作用目前依 CLAUDE.md §5 是被*擱置 (disregarded)* 的
> (「open question for advisors」)。加入 cell–cell adhesion 位能是 G5 的核心新物理,
> 建置前應與顧問確認。

---

## 3. 尺度、大小與 domain

| 量 | G4 v2(單細胞) | **G5(organoid)** | 依據 |
|----------|---------------------|-------------------|-------|
| Organoid 直徑 | —(一個 10 µm 細胞) | **150–250 µm**(從 ~180 µm 起) | Kolade/Gloria 影片;~200 µm |
| 細胞直徑 | 20 µm | **16–20 µm** | 影片(細胞核 ~15–20 µm) |
| N 細胞數 (2D disk) | 1 | **~19–37**(hex packing,⌀ ~9 細胞) | 200 µm / 20 µm ≈ 10 across → ~⌊π·5²/…⌋;從小起(19) |
| Domain | 180 µm | **≥ 600 µm**(≈ organoid + 200 µm 重塑 halo ×2) | 力傳播 ~8–10 cell diam ≈ 200 µm(Shenoy 2014;Tsingos 2023) |
| 纖維數 | 99 | **~1000–3000**(依較大 domain 密度匹配) | Lee 2014 密度;依面積縮放 |
| Beads | 幾 ×10³ | **~10⁴–幾 ×10⁴** | bead 間距 1 µm × 總纖維長度 |

**後果:** ~10–30× 更多 beads/纖維以及 M 個細胞。這是效能瓶頸 — 見 §6。
從 2D 起(符合 confocal 平面成像與所有先前世代);3D 屬 Stage-F。

---

## 4. 要加入的物理 — 一階 vs 延後排序(來自研究整合)

依 deep-research 證據排序(`docs/G5_research_findings.md`;每項皆引用其來源):

### 一階(於 Stages A–D 建置)
1. **Cell–cell adhesion + tissue surface tension** — 決定 collective vs single-cell(Ilina & Friedl
   2020)。*新增。* G5 的定義性物理。
2. **力誘導的 radial collagen alignment** — 頭號輸出。由許多細胞同時 pull
   + 纖維非線性力學湧現。已部分存在(10% 壓縮 → microbuckling 於 g2
   `elastic_forces`);**加入 tension-driven fiber reorientation / strain-stiffening**,使 pull
   能傳播 ~200 µm 並讓纖維放射狀對齊(Shenoy 2014;Saraswathibhatla 2025;Su/Kim 2021 —
   同一 spheroid 周圍 radial vs circumferential alignment 會改變侵入)。
3. **非線性 fiber-network 力學(strain-stiffening + buckling)** — 線性彈性 matrix 無法產生的
   長程(~8–10 cell 直徑)傳遞所必需(Shenoy 2014;
   Mark 2020 eLife:radial displacement 指數隨 matrix 變硬由 α=−2 → −0.2;contractile
   spheroid 附近 >20× 局部變硬)。部分屬對既有 stretch 項的*參數/力律*更改,而非新模組。

### 二階(Stage E,於 A–D 驗證後)
4. **Matrix plasticity / 不可逆重塑** — 使 radial aster *持續存在*而非彈性回復。
   Bell's-law crosslinker rupture/reformation(CLAUDE.md §5;Nam 2016
   κ≈0.82;Wisdom 2018:高 vs 低 plasticity → 5× 遷移)。這是 lab meeting 的
   「可呈現結果 (presentable results)」需求。**注意:** SLS viscoelasticity 維持*延後* (CLAUDE.md §3a) —
   此處的 plasticity 是 crosslink 拓撲改變,而非 per-spring SLS。
5. **空間異質變硬 (heterogeneous stiffening)** — 局部(Loxl3-type)變硬促進 collective
   侵入;整體變硬抑制之(Ray 2022)。將 stiffness 模型化為一個場,而非單一旋鈕。

### 延後(Stage F / 後續版本)
- MMP 蛋白水解 matrix 降解 / 開孔(Wisdom 2018)— 真正受限單細胞通道所需;
  日後與 plasticity 耦合。
- Proliferation 驅動的 organoid 內壓(open question:相較於純收縮 pull,它是否改變 alignment?)。
- Durotaxis / 明確接觸引導 (contact-guidance) 轉向;EMT continuum;3D;免疫/IL-6(v3)。

---

## 5. 參數與 provenance

**繼承(保留自 g2/g4,已有來源 — CLAUDE.md §8、`docs/g4_parameter_provenance.md`):**
fiber stretch κ_s=4.0×10⁻³ N/m、bending κ_b=8.27×10⁻²⁰ N·m、crosslink stiffness、bead 間距 1 µm、
Gaussian σ≈1.5–2 µm、bead drag ζ、壓縮比 0.10,以及完整 motor-clutch 組
(N=12/site、k_c=2 nN/µm、k_on=0.055 s⁻¹、k_off0=0.018 s⁻¹、F_b=1.5 nN、v0=0.025 µm/s、
F_stall=8 nN/site;Adebowale 2021 SI Table 4 尺度、Bell 1978)。

**G5 引入的新參數(每項使用前需要來源):**

| 符號 | 意義 | 建議值 | 來源 / 狀態 |
|--------|---------|----------------|-----------------|
| N_cells | organoid 中的細胞數 | 19 → 37 (2D) | 影片;為效能從小起 |
| R_cell | 細胞半徑 | 8–10 µm | 影片 |
| ε_cc, σ_cc | cell–cell adhesion 深度 / 範圍 | **TBD — 擬合** | 調整至靜止時叢集維持凝聚;Ilina 2020 定性 |
| γ_tissue | 有效組織表面張力 | **TBD — 由 ε_cc 推導** | — |
| strain-stiffening 起始 / 指數 | 纖維張力非線性 | **TBD** | Steinwachs 2016 / Mark 2020 collagen 擬合 |
| k_off,xl, F_b,xl | crosslink Bell rupture(plasticity, Stage E) | **TBD** | Nam 2016 κ≈0.82 目標;先檢查 bond strains(CLAUDE.md §5) |

> 紀律:切勿無憑據地憑空發明 ε_cc/σ_cc — 要麼將其擬合到靜止 organoid 目標
> (叢集既不散開也不塌縮),要麼向 Kolade 詢問 cadherin 尺度估計。每個值皆
> 記錄於 `parameter_provenance.md`。

---

## 6. 效能計畫(真正的風險)

引擎成本瓶頸**不是**每步的力律(它可順利 Numba 編譯,且以 O(E+T+links) 縮放);
而是**網絡建置與接觸偵測**,目前為 Python 的 O(F²) 迴圈:
- `build_crosslinks` 對所有纖維對做巢狀迴圈 → O(F²),含完整 segment×segment 測試。
- `make_network_spec` 對整個網絡重試以通過 percolation gate。
- `contact_patches` 每次接觸更新皆遍歷每根纖維 → 對 organoid ×M 個細胞 = M×F。

**擴大規模前必需:**
1. **空間 neighbor grid(cell list)**,用於 (a) crosslink 建置與 (b) per-cell 接觸
   偵測 — 以 O(F) / O(M·local) 取代 O(F²) 與 M×F 的 Python 迴圈。
2. 將 g4_v2 的 `_advance_numba` **repulsion 迴圈推廣到 M 個 cell centers**。
3. 保留 Numba integrator;在 Numba 不適用處重用 g3 的 `np.bincount` scatter。

交付一個 **perf smoke test**(Stage A gate):2000 纖維 + 30 細胞必須以可接受的
速率建置 + 步進,才可加入 kinetics。

---

## 7. 分期實作計畫(對照 CLAUDE.md §3b 紀律 — 依序驗證)

| Stage | 交付物 | Gate 判準 |
|-------|-------------|----------------|
| **G5-A** | 多細胞骨架:N 個 disks、cell–cell adhesion 位能、neighbor grid、大型網絡生成器、跨 M 個 centers 的 Numba repulsion | 靜止 organoid 凝聚且力平衡;perf smoke test 通過(2–3k 纖維,~30 細胞) |
| **G5-B** | 收縮性 organoid 拉扯 collagen(所有細胞經 Gaussian kernel 抓握;尚無細胞運動) | **Radial alignment 湧現** — radial-order 指標上升;displacement halo 延伸 ~200 µm |
| **G5-C** | 加入 strain-stiffening / tension-driven fiber reorientation | displacement 縮放與 alignment 持久性在數量級內符合文獻目標(§8) |
| **G5-D** | 釋放細胞(translation+rotation)配合 cell–cell adhesion → **collective vs single-cell** | 掃描 adhesion 強度 → 模式轉變(凝聚 front ↔ 逃逸細胞) |
| **G5-E** | Matrix plasticity(crosslink Bell rupture/reformation)使 aster **持續存在** | κ 指標 > 0;放鬆 pull 後 aster 仍保留(對照彈性回復控制組) |
| **G5-F** | 延後的豐富度:proteolysis、proliferation-pressure、3D、vertex cells | 僅在 A–E 驗證後 |

每個 stage:保存先前版本(CLAUDE.md §7.2)、加入測試、重新生成 LaTeX equation summary、
更新 `parameter_provenance.md`。

---

## 8. 驗證目標(實驗數字 → 模型指標)

這些是「模型是否愈來愈接近?」的**量化錨點**,也是 ablation 的基礎。
Breast/TNBC/PyMT 數字是跨組織的**數量級錨點**,非 PDAC 精確值(研究提出的警示)。
優先採用實驗量測的指標(alignment index、displacement 大小、contractility、plasticity %),
勝過模型輸出的 power-law 指數。

| 指標(模型) | 實驗目標 | 來源 |
|----------------|--------------------|--------|
| Organoid 附近 radial alignment index | Day1→3 上升、Day3→5 **趨於平穩**;於 0–50 µm 邊界最高;**持續 >100 µm**;>0.2 = 各向異性;~0.5→0.75 上升 | Saraswathibhatla 2025;Ray 2022;Lee 2017(方法) |
| Radial vs tangential displacement 縮放 | 指數 **n_r≈1.2**(radial,較遠)vs **n_t≈2.2**(tangential) | Saraswathibhatla 2025 |
| Radial displacement 指數 vs matrix 非線性 | **α:−2(線性,~1 Pa)→ −0.2(>1000 Pa,已變硬)** | Mark 2020 eLife 51912 |
| 整體 organoid contractility / 壓力 | **344±35 µN / 677±68 Pa**(TNBC ~4000-cell spheroid,24 h,1.2 mg/mL) | Mark 2020 |
| 近表面 strain / 變硬 | ~**200 µm 表面變形(>50% strain)**、**>20× 局部變硬** | Mark 2020 |
| Matrix plasticity index | **κ≈0.82**(collagen);100 Pa 下 10–30% 永久 strain | Nam 2016;Wisdom 2018 |
| Plasticity → 遷移 | 高 vs 低 plasticity → **~5× 遷移** | Wisdom 2018 |
| Radial vs circumferential 纖維取向 → 侵入 | radial 側 → 更多散播叢集(同一 spheroid) | Su/Kim 2021 |

> 缺少的 PDAC/OSCC 專屬數字(traction、alignment 上升率、strand 外伸 µm/h、
> 單細胞速度)是一個 **open question** — 詢問 Kolade/Gloria 他們的影片能否被
> 量化(能:他們有 timelapse + collagen probe)。這將提供*組織匹配*的
> 目標,以取代 breast 代理值。

---

## 9. Ablation 研究設計(機制開/關 → 指標 → 與實驗比較)

最強的實驗 ablations 直接對應到 in-silico 開關:

| Ablation(把 X 關/開) | 模型開關 | 指標 | 預期(實驗) |
|--------------------------|-------------|--------|-----------------------|
| 同一 organoid 周圍 **radial vs circumferential** 初始纖維取向 | seed 網絡各向異性 | 散播 / 逃逸細胞數 | radial → 更多侵入(Su/Kim 2021) |
| **Cell–cell adhesion** 高 vs 低 | ε_cc | 集體 front vs 單細胞逃逸比例 | 低 adhesion → 單細胞(Ilina 2020) |
| **Strain-stiffening** 開/關 | 線性 vs 非線性 stretch 律 | displacement 指數 α、alignment 觸及範圍 | 非線性 → 遠程 alignment(Mark 2020) |
| **Plasticity** 高 vs 低 | crosslink k_off,xl | κ;卸載後 aster 持久性;遷移 | 高 → 持久 aster + ~5× 遷移(Wisdom 2018) |
| **整體 vs 局部變硬** | stiffness 場 | 集體侵入範圍 | 局部促進、整體抑制(Ray 2022) |

每個 ablation 只改**一個**開關並比較**一個**指標對文獻目標 — 乾淨的
Sobol 式歸因,釐清什麼才重要(H1–H4)。這是 Kolade 實驗室演講的交付物。

---

## 10. 給顧問的 open questions

1. **Cell–cell 交互作用**(Hongbo/Kolade):確認現在加入 cell–cell adhesion 位能 —
   它目前依 CLAUDE.md §5 是被*擱置*的。Disk-cluster vs vertex/CPM 細胞?
2. **PDAC/OSCC 專屬目標**(Kolade/Gloria):既有影片能否被量化出
   alignment-index 上升、strand 外伸速度(µm/h)、單細胞速度,以取代 breast
   代理值?
3. **Plasticity 時機**(Kolade):在 G5-E 啟動 crosslink-rupture plasticity,或
   在第一個 organoid 里程碑先讓 aster 保持純彈性?先檢查 bond-strain 大小
   (Gloria 的 V4 負面結果:5 nN 下 weak-bond 形成可忽略)。
4. **Proliferation 壓力**(Hongbo):organoid 內部生長壓力相較於純收縮 pull,是否
   實質改變 collagen alignment,還是可延後?

---

*核准後要生成的伴隨文件:`docs/G5_research_findings.md`(完整引用文獻
整合)、`docs/G5_code_assessment.md`(引擎審查),以及每個 stage 的 LaTeX equation summary。*
