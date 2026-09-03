# ECM 參數溯源 — 目前 Gloria G2 基線模型

本文件清楚區分：校正後的 Generation 2（G2）Gloria 模型**目前實際使用**
的參數，以及 Saraswathibhatla acIPN 實驗所提供、供未來校準使用的目標。
標示為「尚未建模」的數值，不可表述為目前模擬的輸入值。

G2 的基本單位為：長度 µm、力 nN、時間 s。基線模型為平面離散 collagen
網路，受一個 MDA-MB-231 細胞拉動；它尚不是 3-D acIPN 腫瘤球體模型。

英文對等版本：[parameter_provenance.md](parameter_provenance.md)。

## 1. Collagen 濃度與 bulk 剛度

| 項目 | 目前 G2 設定 | Saraswathibhatla acIPN 參考值／未來目標 | 狀態 |
|---|---:|---:|---|
| Type-I collagen 濃度 | 尚未建模 | 1.5–1.6 mg/mL | 不是 G2 輸入參數 |
| Alginate 濃度 | 尚未建模 | 4.8 mg/mL | 不是 G2 輸入參數 |
| Bulk 儲能模數 G′ | 未計算、未校準 | 0.6 / 1.2 / 3 kPa | 需要以虛擬流變測試校準 |
| 純 collagen 比較範圍 | 尚未建模 | 1–4 mg/mL | 僅為文獻背景 |

因此，目前 G2 不可稱為 1.5 mg/mL acIPN 模型，也不代表三個 acIPN 剛度條件中的任一個。必須先以 bulk shear/tensile test 校準其纖維與交聯力學。

## 2. 纖維幾何與離散化

| 參數 | 目前 G2 設定 | Saraswathibhatla／文獻參考值 | 狀態 |
|---|---:|---:|---|
| 有效纖維直徑 | 0.30 µm（300 nm） | Saraswathibhatla 100 nm；量測報告範圍 62–325 nm | 已使用 |
| 纖維輪廓長度 | 20–80 µm | 文獻報告 8–100+ µm；Lee 模型為 8 µm | 已使用 |
| 名義 bead 間距 | 0.75 µm | 1.0 µm | 已使用 |
| 每根纖維 beads 數 | 依輪廓長度重採樣，會變動；儲存 seed 為 7,173 / 99 = 平均 72.5 | Lee: 5；先前專案註記: 3 | 已使用；不是 3 或 5 |
| 纖維中心線 | 隨機直線至彎曲線段；曲率振幅 0.9 µm | 非 Saraswathibhatla SI 輸入 | 已使用 |

0.75 µm 是數值離散的節段間距，不是實體 collagen 纖維的直徑。G2 的每根纖維鏈比五-bead 表示法細得多。

## 3. 纖維彈性

| 項目 | 目前 G2 設定 | Saraswathibhatla／文獻參考值 | 狀態 |
|---|---:|---:|---|
| 單根纖維 Young’s modulus，E | 32 MPa | 32 MPa 模擬值；文獻範圍約 30–800 MPa | 已使用 |
| 軸向力學規則 | `k = EA / l0` | 直接彈簧 `κ_s,f = 4.0×10⁻³ N/m` | 已使用，但參數化方式不同 |
| 名義軸向彈簧（`l0 = 0.75 µm`） | 3,016 nN/µm = 3.02 N/m | 0.004 N/m | 由目前參數推得 |
| 壓縮／拉伸剛度比 | 0.10 | 你提供的 acIPN 表中未指定 | 已使用；microbuckling 近似 |
| Beam 彎曲剛度 EI | 12.72 nN·µm² = `1.27×10⁻²⁰ N·m²` | `κ_b,f = 8.27×10⁻²⁰ N·m` | 已使用；需統一 force law 後才可直接比較 |
| Bulk G′ | 未計算 | 1 mg/mL: ~3 Pa；2 mg/mL: ~44–200 Pa；3 mg/mL: ~97–550 Pa | 尚未建模 |

目前 G2 名義軸向彈簧由 32 MPa、300 nm 有效直徑及 0.75 µm 節段長度推得。其數值約為 Saraswathibhatla 的 0.004 N/m 的 754 倍；這是材料模型差異，不能只交換單一常數，必須校準。

## 4. Bead–spring 與交聯參數

| 參數 | 目前 G2 設定 | Saraswathibhatla SI Table 2 參考值 | 狀態 |
|---|---:|---:|---|
| 纖維節段平衡長度 | 每段由重採樣決定，名義值 0.75 µm | 1.0 µm | 已使用 |
| 纖維直徑 | 300 nm | 100 nm | 已使用 |
| 纖維伸縮常數 | `EA/l0` 推得，名義值 3.02 N/m | `4.0×10⁻³ N/m` | 已使用；形式與數值不同 |
| 纖維彎曲律 | 離散 beam `EI/l0³` | `8.27×10⁻²⁰ N·m` | 已使用；慣例不同 |
| 永久交聯位置 | 兩纖維節段的精確交點 | 20 nm arm | 已使用；臂長等效為零 |
| 交聯劑直徑 | 尚未建模 | 10 nm | 尚未建模 |
| 永久交聯剛度 | 75 nN/µm = 0.075 N/m | `2.0×10⁻³ N/m` | 已使用 |
| 交聯 binding-site 間距 | 尚未建模 | 每 100 nm 一個 | 尚未建模 |
| 交聯型態 | 永久、自由鉸接，初始交點處無應力 | 具實體 arm 的 crosslinker | 模型不同 |
| V4 新形成 weak link 剛度 | 15 nN/µm = 0.015 N/m | 不屬於所提供的 Saraswathibhatla 基線參數 | 僅 V4 |

`crosslink_modulus_mpa = 0.40` 仍保留在 G2 組態作為溯源資訊，但未進入力學計算；永久交聯實際使用的是 `crosslink_stiffness = 75 nN/µm`。

## 5. 網路密度與幾何

| 參數 | 目前 G2 設定 | Saraswathibhatla／Lee 參考值 | 狀態 |
|---|---:|---:|---|
| 模擬區域 | 180 × 180 µm 的平面正方形 | 內徑 50 µm、外徑 200 µm、高 1 µm 的環形 | 幾何不同 |
| 纖維數 | 99 | 3-D 密度 ~2.85 fibres/µm³ | 已使用 |
| 2-D 面密度 | 0.00306 fibres/µm²（99 / 180²） | 無直接對應值 | 推得；未定義厚度 |
| 儲存 seed 17 的交聯數 | 383 個永久交點交聯 | crosslinker/fibre ratio ~5.72 | 已使用的快照 |
| 儲存 seed 17 的每纖維交聯數 | 3.87（383 / 99） | ~5.72 | 推得的快照 |
| 交聯數規則 | 由纖維交點生成；會隨網路幾何／seed 改變 | 以 crosslinker ratio 控制的網路 | 規則不同 |
| 邊界條件 | 僅固定外側 2.5 µm 邊帶內的 beads | 非 acIPN 環形邊界 | 已使用 |

目前模型為 2-D；未指定有效厚度，故不可直接將 3-D fibre density 或 acIPN 的 crosslinker/fibre ratio 移植進來。

## 6. 細胞或腫瘤球體尺寸

| 項目 | 目前 G2 設定 | acIPN 腫瘤球體參考值 | 狀態 |
|---|---:|---:|---|
| 生物對象 | 單一 MDA-MB-231 細胞 | 250 cells 聚集的 tumour spheroid | 生物系統不同 |
| 半徑 | 9 µm | 約 50 µm（day-0 直徑約 100 µm） | 已使用；非等效 spheroid |
| 初始內邊界 | 無；僅有圓形 cell exclusion/pull surface | 環形內徑 50 µm | 尚未建模 |
| 主動拉力 | V2/V4 總力 5 nN（另有 2.5/5/10 nN sensitivity）；V3 使用 motor–clutch | 不同於你提供的 spheroid force protocol | 已使用 |

## 7. 位移場驗證

| 驗證項目 | 目前 G2 狀態 | Saraswathibhatla 目標 | 狀態 |
|---|---|---|---|
| Radial displacement 擬合 | 尚未實作 | `u_radial ∼ r⁻¹·²` | 待完成 |
| Tangential displacement 擬合 | 尚未實作 | `u_tangential ∼ r⁻²·²` | 待完成 |
| 目前空間輸出 | near/intermediate/far shells 的平均位移與 alignment | 完整 radial/tangential displacement profiles | 僅部分對應 |

上述冪次律應作為未來驗證目標，但目前 G2 尚未計算或擬合它們，因此不能標示為通過。

## 程式溯源

- 基礎幾何、材料參數、載入與積分預設值：`generations/g2_corrected/common/model.py` 的 `CollagenConfig`。
- 軸向與彎曲換算：同檔案的 `CollagenConfig.axial_rigidity` 與 `CollagenConfig.bending_rigidity`。
- 僅限交點的永久交聯建立：同檔案的 `build_crosslinks`。
- seed 17 儲存網路計數：`generations/g2_corrected/validation_summary.json`。
- V3 motor–clutch 覆寫參數：`generations/g2_corrected/v3_two_sided_migration/model.py` 的 `MigrationConfig`。
- V4 weak-link 覆寫參數：`generations/g2_corrected/v4_contact_plasticity/model.py` 的 `PlasticityConfig`。

本文提到的 Saraswathibhatla、Lee、Picu、Stein 與 Licup 數值，是使用者提供、用於模型比較的研究參考；除非上表明列為「目前 G2 設定」，否則不會自動成為 G2 的輸入參數。
