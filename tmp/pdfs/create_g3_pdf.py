from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.platypus import (
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "pdf" / "G3_model_overview_zh-TW.pdf"
OUT.parent.mkdir(parents=True, exist_ok=True)

pdfmetrics.registerFont(UnicodeCIDFont("MSung-Light"))
FONT = "MSung-Light"

styles = getSampleStyleSheet()
title = ParagraphStyle("TitleZH", parent=styles["Title"], fontName=FONT, fontSize=23, leading=31,
                       textColor=colors.HexColor("#18364A"), alignment=TA_CENTER, spaceAfter=10 * mm)
subtitle = ParagraphStyle("SubtitleZH", parent=styles["Normal"], fontName=FONT, fontSize=10.5, leading=16,
                          textColor=colors.HexColor("#4C6675"), alignment=TA_CENTER, spaceAfter=10 * mm)
h1 = ParagraphStyle("H1ZH", parent=styles["Heading1"], fontName=FONT, fontSize=16, leading=24,
                    textColor=colors.HexColor("#0D596B"), spaceBefore=7 * mm, spaceAfter=4 * mm)
h2 = ParagraphStyle("H2ZH", parent=styles["Heading2"], fontName=FONT, fontSize=12.5, leading=19,
                    textColor=colors.HexColor("#174E63"), spaceBefore=4 * mm, spaceAfter=2.5 * mm)
body = ParagraphStyle("BodyZH", parent=styles["BodyText"], fontName=FONT, fontSize=9.8, leading=16,
                      textColor=colors.HexColor("#172B36"), spaceAfter=2.5 * mm, wordWrap="CJK")
small = ParagraphStyle("SmallZH", parent=body, fontSize=8.6, leading=13)
callout = ParagraphStyle("CalloutZH", parent=body, fontSize=11, leading=18, textColor=colors.HexColor("#173E51"),
                         leftIndent=4 * mm, rightIndent=4 * mm, spaceAfter=3 * mm)
table_head = ParagraphStyle("TableHead", parent=small, textColor=colors.white, alignment=TA_CENTER)
table_body = ParagraphStyle("TableBody", parent=small, alignment=TA_LEFT)


def p(text, style=body):
    return Paragraph(text, style)


def bullets(items):
    return [p(f"• {item}") for item in items]


def table(rows, widths):
    cooked = []
    for index, row in enumerate(rows):
        style = table_head if index == 0 else table_body
        cooked.append([p(str(cell), style) for cell in row])
    result = Table(cooked, colWidths=widths, repeatRows=1, hAlign="LEFT")
    result.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1B667B")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#B8CBD2")),
        ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#F3F7F8")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.HexColor("#F3F7F8"), colors.white]),
        ("LEFTPADDING", (0, 0), (-1, -1), 2.6 * mm),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2.6 * mm),
        ("TOPPADDING", (0, 0), (-1, -1), 2.2 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2 * mm),
    ]))
    return result


def section_heading(text):
    return p(text, h1)


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#B7CCD3"))
    canvas.line(doc.leftMargin, 12 * mm, A4[0] - doc.rightMargin, 12 * mm)
    canvas.setFillColor(colors.HexColor("#55717D"))
    canvas.setFont(FONT, 8)
    canvas.drawString(doc.leftMargin, 7.5 * mm, "G3 Motor-Clutch-Collagen Model")
    canvas.drawRightString(A4[0] - doc.rightMargin, 7.5 * mm, f"{doc.page}")
    canvas.restoreState()


story = []
story.extend([
    Spacer(1, 30 * mm),
    p("G3 Motor-Clutch-Collagen", title),
    p("模型機制、目前證據與尚待驗證事項", subtitle),
    Table([[p("一句話先說：G3 把 G2 裡「先指定細胞往右」拿掉，改成讓細胞向四周試探，看看 collagen 能不能自己引導細胞選方向。", callout)]], colWidths=[160 * mm], style=TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#E7F2F4")),
        ("BOX", (0, 0), (-1, -1), 0.8, colors.HexColor("#7DAEBC")),
        ("LEFTPADDING", (0, 0), (-1, -1), 4 * mm),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4 * mm),
        ("TOPPADDING", (0, 0), (-1, -1), 3 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3 * mm),
    ])),
    Spacer(1, 9 * mm),
    p("它目前仍是 2D minimal mechanism model，不是真實的 tumor migration simulation。", body),
    Spacer(1, 12 * mm),
    p("本文整理 G3 三個 stages（G3A、G3B、G3C）的機制、現有 smoke-run 結果，以及可以與不可以宣稱的結論。", body),
    PageBreak(),
    section_heading("整體架構"),
    p("可以把它想成一顆圓球，周圍有很多小手，在幾條繩子之間試著抓住、拉動。"),
    table([
        ["流程", "作用"],
        ["Rigid cell", "圓形細胞，表面配置 protrusions 與 clutches。"],
        ["24 個 protrusion sectors", "細胞向四周試探；最多兩個同時 active。"],
        ["200 個 motor-clutches", "抓住特定 fibre 位置，並施加局部力量。"],
        ["Local Gaussian projection", "把 point force 分配到附近 beads。"],
        ["Elastic bead-spring collagen", "纖維的 stretching 與 bending deformation。"],
        ["Geometry + traction feedback", "回饋哪些 protrusion directions 較成功。"],
        ["Reaction force + torque", "collagen 對 cell 的反作用力，驅動平移或旋轉。"],
    ], [54 * mm, 106 * mm]),
    p("<b>核心三層</b>"),
    *bullets([
        "<b>Cell</b>：圓形 rigid cell，表面有 protrusions 和 clutches。",
        "<b>Interface</b>：每個 clutch 抓住特定 collagen material point，再用 Gaussian kernel 把力量傳給附近 beads。",
        "<b>Collagen</b>：由 beads、stretching springs、bending elasticity 組成的纖維。",
    ]),
    p("<b>Baseline</b>"),
    *bullets([
        "8 條 fibres；每條長 40 µm；bead spacing 1 µm；cell radius 10 µm。",
        "200 clutches、24 個可能的 protrusion directions，同時最多 2 個 active protrusions。",
        "純 elastic ECM，沒有 SLS、plasticity 或 transient crosslinks。",
    ]),
    section_heading("G3A：Clutch 到底抓在哪裡？"),
    p("G2 比較像細胞在一整片區域灑一股平均力量。G3A 則讓第 17 個 clutch 明確抓住第 3 條 fibre、第 12 段、距離左端 35% 的位置。"),
    p("每個 attachment 儲存 <b>fiber_id</b>、<b>segment_id</b> 與 <b>alpha</b>。其中 alpha 代表 clutch 位於該 segment 的哪一點。當 fibre 變形時，clutch 會跟著同一個 material point 移動，不會每個 timestep 重新跳到最近的 bead；只有 unbind 後才能重新找位置。"),
    p("<b>力怎麼傳到 collagen？</b>Clutch 產生 point force，但電腦中的 fibre 是離散 beads。因此利用 Gaussian kernel，把力量分給同一條 fibre 上鄰近的 beads。"),
    *bullets([
        "bead forces 的總和等於 clutch force。",
        "投影後的 torque 等於原本 point force 的 torque。",
        "cell 受到大小相同、方向相反的 reaction force。",
    ]),
    table([
        ["G3A：15 秒 single-fibre smoke run", "結果"],
        ["FOI", "0.6046 → 0.6239"],
        ["ΔFOI", "+0.0193"],
        ["最大 bead displacement", "0.395 µm"],
        ["最多 bound clutches", "173 / 200"],
        ["Force conservation error", "小於約 4×10⁻¹⁶"],
        ["Torque error", "小於約 9×10⁻¹⁸"],
    ], [75 * mm, 85 * mm]),
    p("結論：clutch 已經能穩定抓住 collagen 的特定位置，而且力量傳遞與 fibre deformation 正常。這是目前三個 stages 中證據最完整的一層。", callout),
    PageBreak(),
    section_heading("G3B：細胞怎麼自己選方向？"),
    p("<b>G2 的問題。</b>G2 V3 有 polarity_probability = 0.65，等於事先告訴細胞有 65% 的傾向支持某一側。因此持續方向主要是輸入，不是 collagen 自己產生的。"),
    p("<b>G3B 的做法。</b>G3B 把 cell surface 分成 24 個 sectors，每格 15°。一開始沒有固定 +x preference，也沒有 polarity_probability = 0.65；最多隨機開啟兩個 protrusions，並將 200 個 clutches 分配給 active protrusions。"),
    p("每個方向會看兩件事："),
    *bullets([
        "<b>附近 collagen geometry</b>：附近有沒有 fibre、fibre 是否沿著這個 protrusion axis、是否在 clutch capture distance 內。",
        "<b>Traction success</b>：有多少 clutches 成功 bound、能不能建立 traction、clutches 是否很快 rupture。",
    ]),
    p("若某個 protrusion 附近 collagen 較多、fibre orientation 比較合適、clutches 綁得住而且可以維持 traction，它就會活得比較久。反之，黏不住的方向會消失，cell 再嘗試其他 sector。"),
    p("<b>目前結果：120 秒 aligned-fibre smoke run</b>"),
    *bullets([
        "active protrusion sectors 確實發生 turnover。",
        "clutches 能在不同方向形成 spatial attachments。",
        "最多約 11 個 clutches 同時 bound。",
        "程式沒有固定向右的偏好。",
    ]),
    p("但目前還不能說 collagen 已經成功讓 cell 自發選出 migration direction，因為只跑了少數 demonstration trajectories。"),
    p("<b>要證明 emergent guidance，仍需要：</b>"),
    *bullets([
        "100 個 independent seeds，以及 isotropic ECM control、aligned ECM、rotated-aligned ECM、feedback-off control。",
        "確認正、反方向是否接近 50 / 50。",
        "確認 ECM 旋轉 30° 時，選擇方向是否也旋轉 30°。",
    ]),
    p("結論：<b>protrusion feedback mechanism 已經能運作，但 emergent guidance 尚未完成 statistical validation。</b>", callout),
    section_heading("G3C：力量能不能讓 Cell 移動和旋轉？"),
    p("G3B 的 cell 仍然固定，方便單獨檢查方向選擇；G3C 才解除固定。每個 clutch 拉 collagen 時，collagen 也會反過來拉 cell。所有 reaction forces 的總和提供 translation force（推動 cell center）和 torque（轉動 cell body）。"),
    p("G3C 沒有 prescribed speed、self-propulsion v0、固定的 +x force，也沒有 polarity_probability = 0.65。因此 cell 只有在 clutches 產生不對稱 reaction force 時才會動。"),
    PageBreak(),
    section_heading("G3C：反作用力、平移與旋轉"),
    p("如果 clutch force 永遠指向圓心，moment arm 與 force 平行，torque 永遠是零。G3 使用的設定是：clutch 剛 bind 時沿 cell surface normal；之後 force direction 跟著實際 clutch vector 演化。當 fibre 和 cell 的相對位置改變，force 就可能產生 tangential component，因而產生 torque。"),
    table([
        ["G3C：asymmetric-torque fixture，30 秒 smoke run", "結果"],
        ["Cell displacement", "0.00187 µm"],
        ["Cell rotation", "6.68×10⁻⁵ rad"],
        ["運動來源", "完全來自 clutch reaction"],
        ["Empty-ECM control", "displacement 和 rotation 均為零"],
        ["力學測試", "Mirror / rotation mechanics tests 通過"],
    ], [92 * mm, 68 * mm]),
    p("結論：<b>不對稱的 spatial clutch forces 確實可以讓 rigid cell 平移和旋轉。</b>但位移非常小，因此不能說已經模擬出 realistic migration。", callout),
    section_heading("FOI 與 κ 在看什麼？"),
    p("<b>FOI（fibre orientation index）</b>在問 collagen fibres 有沒有因為 cell pulling 而變得更一致、更 aligned。會比較 pulling 前、pulling 中、關閉 clutches 後與 recovery 後。"),
    p("<b>κ</b>則在問 cell 停止拉動後，alignment 有多少留下來。G3 是 permanent elastic network，沒有真正 plasticity mechanism，所以理論上應該大部分恢復。"),
    table([
        ["Load-unload 指標", "結果"],
        ["Initial FOI", "0.6046"],
        ["Pull 後 FOI", "0.6239"],
        ["Recovery 後 FOI", "0.6141"],
        ["κ", "0.489"],
        ["600 秒後 elastic energy", "約為 peak 的 1.27%"],
    ], [75 * mm, 85 * mm]),
    p("這沒有通過設定的 κ &lt; 0.1 recovery gate。但不能解釋成 plasticity，因為模型裡根本沒有 irreversible mechanism。正確解釋是：600 秒內仍未完全 relaxation，目前標記為 <b>unresolved_recovery</b>。", callout),
    section_heading("G2 與 G3 的最重要差別"),
    table([
        ["問題", "G2 V3", "G3"],
        ["方向怎麼來？", "0.65 預先指定", "collagen geometry + traction feedback"],
        ["Clutch 黏在哪？", "左右兩側 averaged coupling", "特定 fibre segment material point"],
        ["力怎麼傳？", "一側總力再 Gaussian 分配", "每個 clutch 各自局部投影"],
        ["Cell 可以怎麼動？", "只沿 x 軸", "2D translation + rotation"],
        ["Torque", "沒有完整 spatial torque", "由 evolving clutch vector 產生"],
        ["Cell speed", "主要由 drag calibration 決定", "仍受 drag 影響，不作為主要 prediction"],
        ["現在能宣稱什麼？", "Prescribed imbalance 可推動 cell", "Spatial clutch、feedback、reaction mechanics 已建立"],
    ], [35 * mm, 58 * mm, 67 * mm]),
    PageBreak(),
    section_heading("最誠實的總結"),
    p("目前 G3 已經成功完成「模型零件」："),
    *bullets([
        "Clutch 可以抓住特定 collagen material point。",
        "Cell 可以無偏向地向四周試探。",
        "Collagen geometry 和 traction 可以影響 protrusion persistence。",
        "Clutch reaction 可以讓 cell 平移與旋轉。",
        "沒有偷偷加入 0.65 polarity 或 self-propulsion。",
    ]),
    Spacer(1, 6 * mm),
    Table([[p("尚未完成的核心統計問題：在不指定方向時，aligned collagen 是否真的能可靠地引導 cell axis 和 trajectory？", callout)]], colWidths=[160 * mm], style=TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#FCEEDC")),
        ("BOX", (0, 0), (-1, -1), 0.8, colors.HexColor("#D6A766")),
        ("LEFTPADDING", (0, 0), (-1, -1), 4 * mm),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4 * mm),
        ("TOPPADDING", (0, 0), (-1, -1), 3 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3 * mm),
    ])),
    Spacer(1, 7 * mm),
    p("這要等完整 100-seed G3B / G3C validation 才能回答。", body),
])

doc = SimpleDocTemplate(str(OUT), pagesize=A4, leftMargin=25 * mm, rightMargin=25 * mm,
                        topMargin=20 * mm, bottomMargin=22 * mm, title="G3 Motor-Clutch-Collagen 模型說明")
doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(OUT)
