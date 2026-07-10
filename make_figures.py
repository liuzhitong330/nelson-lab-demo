"""
make_figures.py — Nelson Lab Demo
Figure 1: Scatter plot of Drd1 vs Grin2b, dot size = Drd2
           highlights CP and ACB as the striatal zones with combined dopamine
           receptor + glutamate receptor expression relevant to LID
Figure 2: Schematic of the TRAPed dMSN mechanism (Ryan et al. 2024, Cell Reports)
           showing how dMSN subpopulation heterogeneity drives levodopa-induced dyskinesia
"""

import csv, os, math

OUT = os.path.dirname(os.path.abspath(__file__))

rows = []
with open(os.path.join(OUT, "striatum_genes.tsv")) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append(row)

# ── Figure 1: Scatter Drd1 vs Grin2b, size = Drd2 ───────────────────────────
FW, FH = 680, 460
PAD_L, PAD_R, PAD_T, PAD_B = 78, 30, 60, 68
AW = FW - PAD_L - PAD_R
AH = FH - PAD_T - PAD_B

x_gene = "Drd1"
y_gene = "Grin2b"
s_gene = "Drd2"

xmax = max(float(r[x_gene]) for r in rows) * 1.15
ymax = max(float(r[y_gene]) for r in rows) * 1.12
smax = max(float(r[s_gene]) for r in rows)

def px(v): return PAD_L + v / xmax * AW
def py(v): return PAD_T + AH - v / ymax * AH
def sz(v): return 6 + v / smax * 18

COLORS = {
    "CP":   "#c0392b",  # dorsal striatum — red (primary LID site)
    "ACB":  "#e74c3c",  # nucleus accumbens — lighter red
    "OT":   "#e67e22",  # olfactory tubercle — orange
    "MOp":  "#2980b9",  # motor cortex — blue
    "Thal": "#8e44ad",  # thalamus — purple
    "SNr":  "#7f8c8d",  # SNr — grey
    "Cereb":"#95a5a6",  # cerebellum — light grey
    "Hipp": "#27ae60",  # hippocampus — green
}

STRIATAL = {"CP", "ACB", "OT"}

# Absolute label positions (computed to avoid overlap): {region: (lx, ly, anchor, dx, dy)}
# dx, dy = offset from dot center for leader line endpoint
LABELS = {
    "CP":    (px(0.0138)+14,  py(0.1439)-18,  "start",   0, -8),
    "ACB":   (px(0.1051)+12,  py(0.1138)+5,   "start",   4,  0),
    "OT":    (px(0.0220)+12,  py(0.1316)-14,  "start",   4, -6),
    "MOp":   (px(0.0045)-12,  py(0.1393)-14,  "end",    -4, -6),
    "Thal":  (px(0.0049)-12,  py(0.0964)+16,  "end",    -4,  6),
    "SNr":   (px(0.0008)-12,  py(0.0426)+16,  "end",    -4,  6),
    "Cereb": (px(0.0006)-12,  py(0.0115)+14,  "end",    -4,  6),
    "Hipp":  (px(0.0037)+12,  py(0.0972)-14,  "start",   4, -6),
}

dots = ""
labels = ""

for row in rows:
    reg = row["region"]
    xv = float(row[x_gene])
    yv = float(row[y_gene])
    sv = float(row[s_gene])
    col = COLORS.get(reg, "#aaa")
    r = sz(sv)
    cx = px(xv); cy = py(yv)
    border = "2.5" if reg in STRIATAL else "1.5"
    ring = "" if reg not in STRIATAL else (
        f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r+4:.1f}" '
        f'fill="none" stroke="{col}" stroke-width="1" opacity="0.35"/>'
    )
    dots += ring
    dots += (f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
             f'fill="{col}" opacity="0.82" stroke="white" stroke-width="{border}"/>')

    lx, ly, anchor, ddx, ddy = LABELS[reg]
    # Leader line from dot to label
    labels += (f'<line x1="{cx+ddx:.1f}" y1="{cy+ddy:.1f}" x2="{lx-4 if anchor=="start" else lx+4:.1f}" '
               f'y2="{ly:.1f}" stroke="{col}" stroke-width="0.8" opacity="0.5"/>')
    fw = "700" if reg in STRIATAL else "400"
    labels += (f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="{anchor}" '
               f'font-size="10" fill="{col}" font-weight="{fw}">{reg}</text>')

# Gridlines
grids = ""
for xv in [0.02, 0.04, 0.06, 0.08, 0.10]:
    gx = px(xv)
    grids += (f'<line x1="{gx:.1f}" y1="{PAD_T}" x2="{gx:.1f}" y2="{PAD_T+AH}" '
              f'stroke="#eee" stroke-width="1"/>'
              f'<text x="{gx:.1f}" y="{PAD_T+AH+18}" text-anchor="middle" '
              f'font-size="9" fill="#bbb">{xv:.2f}</text>')
for yv in [0.04, 0.08, 0.12]:
    gy = py(yv)
    grids += (f'<line x1="{PAD_L}" y1="{gy:.1f}" x2="{PAD_L+AW}" y2="{gy:.1f}" '
              f'stroke="#eee" stroke-width="1"/>'
              f'<text x="{PAD_L-8}" y="{gy+4:.1f}" text-anchor="end" '
              f'font-size="9" fill="#bbb">{yv:.2f}</text>')

# Axes
axes = (f'<line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{PAD_T+AH}" '
        f'stroke="#ccc" stroke-width="1.2"/>'
        f'<line x1="{PAD_L}" y1="{PAD_T+AH}" x2="{PAD_L+AW}" y2="{PAD_T+AH}" '
        f'stroke="#ccc" stroke-width="1.2"/>'
        f'<text x="{PAD_L+AW/2:.0f}" y="{PAD_T+AH+36}" text-anchor="middle" '
        f'font-size="10.5" fill="#555">Drd1 expression density (D1 receptor, Allen ISH)</text>'
        f'<text x="{PAD_L-52}" y="{PAD_T+AH//2}" text-anchor="middle" '
        f'font-size="10.5" fill="#555" '
        f'transform="rotate(-90,{PAD_L-52},{PAD_T+AH//2})">Grin2b expression (GluN2B)</text>')

# Size legend
leg_x = PAD_L + AW - 110; leg_y = PAD_T + AH - 60
slegend = (
    f'<text x="{leg_x+30}" y="{leg_y-2}" text-anchor="middle" font-size="9" fill="#888">bubble = Drd2</text>'
)
for sv, lab in [(0.096*0.25,"small"),(0.096*0.75,"large")]:
    r2 = sz(sv)
    cx2 = leg_x + (10 if lab=="small" else 50)
    slegend += (f'<circle cx="{cx2}" cy="{leg_y+14}" r="{r2:.1f}" '
                f'fill="#aaa" opacity="0.5" stroke="white" stroke-width="1.5"/>')

svg1 = f"""<svg viewBox="0 0 {FW} {FH}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW//2}" y="20" text-anchor="middle" font-size="13" font-weight="700" fill="#222">
    Dorsal Striatum: High Glutamate Drive Meets D1/D2 Dopamine Balance
  </text>
  <text x="{FW//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Allen Mouse Brain Atlas ISH · Drd1 vs Grin2b (GluN2B) · dot size = Drd2 expression
  </text>
  <text x="{FW//2}" y="54" text-anchor="middle" font-size="10" fill="#444">
    CP (dorsal striatum) has among the highest Grin2b of any region — peak glutamatergic excitability in the LID circuit
  </text>
  {grids}{axes}{dots}{labels}{slegend}
</svg>"""

with open(os.path.join(OUT, "striatum_scatter.svg"), "w") as f:
    f.write(svg1)
print("Wrote striatum_scatter.svg")


# ── Figure 2: TRAPed dMSN mechanism schematic ────────────────────────────────
FW2, FH2 = 700, 400

# Three columns: normal dMSN | TRAPed dMSN | outcome
# Show: Motor cortex → Striatum (dMSN) → SNr → Thalamus/Motor cortex
# Key: TRAPed dMSN has MORE Drd1 + MORE cortical glutamate input

svg2_body = ""

# Background panels
PANEL_W = 190; PANEL_H = 200; PANEL_Y = 80
cx_normal = 90; cx_trap = 350; cx_outcome = 600

def panel(cx, label, col, sublabel=""):
    px0 = cx - PANEL_W//2; py0 = PANEL_Y
    out = (f'<rect x="{px0}" y="{py0}" width="{PANEL_W}" height="{PANEL_H}" '
           f'rx="8" fill="{col}" opacity="0.07"/>'
           f'<rect x="{px0}" y="{py0}" width="{PANEL_W}" height="{PANEL_H}" '
           f'rx="8" fill="none" stroke="{col}" stroke-width="1.4" opacity="0.4"/>'
           f'<text x="{cx}" y="{py0-10}" text-anchor="middle" font-size="10.5" '
           f'font-weight="700" fill="{col}">{label}</text>')
    if sublabel:
        out += (f'<text x="{cx}" y="{py0-0}" text-anchor="middle" font-size="9" fill="{col}">'
                f'{sublabel}</text>')
    return out

svg2_body += panel(cx_normal, "Non-TRAPed dMSN", "#1a5c8a", "(most dorsal striatum dMSNs)")
svg2_body += panel(cx_trap,   "TRAPed dMSN", "#c0392b", "(dyskinesia-associated subpopulation)")
svg2_body += panel(cx_outcome,"LID outcome", "#922b21", "")

# Normal dMSN contents
def neuron_icon(cx, cy, col, r=22):
    return (f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{col}" opacity="0.15"/>'
            f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="{col}" stroke-width="1.5"/>')

CELL_Y = PANEL_Y + 80

# Normal dMSN
svg2_body += neuron_icon(cx_normal, CELL_Y, "#1a5c8a")
svg2_body += (f'<text x="{cx_normal}" y="{CELL_Y+5}" text-anchor="middle" '
              f'font-size="9" fill="#1a5c8a" font-weight="700">dMSN</text>')
# D1 receptor (small, few)
for dx in [-12, 0, 12]:
    svg2_body += (f'<circle cx="{cx_normal+dx}" cy="{CELL_Y-26}" r="3.5" '
                  f'fill="#1a5c8a" opacity="0.4"/>')
svg2_body += (f'<text x="{cx_normal}" y="{CELL_Y-35}" text-anchor="middle" '
              f'font-size="8.5" fill="#1a5c8a" opacity="0.8">Drd1 (baseline)</text>')
# Few glutamate arrows
for dx in [-8, 8]:
    svg2_body += (f'<line x1="{cx_normal+dx}" y1="{CELL_Y-52}" x2="{cx_normal+dx}" '
                  f'y2="{CELL_Y-28}" stroke="#2980b9" stroke-width="2"/>'
                  f'<polygon points="{cx_normal+dx},{CELL_Y-28} {cx_normal+dx-3},{CELL_Y-34} '
                  f'{cx_normal+dx+3},{CELL_Y-34}" fill="#2980b9"/>')
svg2_body += (f'<text x="{cx_normal}" y="{CELL_Y-56}" text-anchor="middle" '
              f'font-size="8.5" fill="#2980b9">Glut. input (basal)</text>')

# Output arrow (thin)
svg2_body += (f'<line x1="{cx_normal}" y1="{CELL_Y+24}" x2="{cx_normal}" y2="{CELL_Y+54}" '
              f'stroke="#1a5c8a" stroke-width="2"/>'
              f'<polygon points="{cx_normal},{CELL_Y+54} {cx_normal-4},{CELL_Y+48} '
              f'{cx_normal+4},{CELL_Y+48}" fill="#1a5c8a"/>'
              f'<text x="{cx_normal}" y="{CELL_Y+64}" text-anchor="middle" '
              f'font-size="8.5" fill="#1a5c8a">→ SNr (normal)</text>')

# Levodopa arrow
svg2_body += (f'<text x="{cx_normal}" y="{CELL_Y+82}" text-anchor="middle" '
              f'font-size="8.5" fill="#7f8c8d">Levodopa → moderate D1 activation</text>')
svg2_body += (f'<text x="{cx_normal}" y="{CELL_Y+94}" text-anchor="middle" '
              f'font-size="8.5" fill="#7f8c8d">→ normal motor output</text>')

# TRAPed dMSN
svg2_body += neuron_icon(cx_trap, CELL_Y, "#c0392b", r=28)
svg2_body += (f'<text x="{cx_trap}" y="{CELL_Y+5}" text-anchor="middle" '
              f'font-size="9" fill="#c0392b" font-weight="700">dMSN</text>'
              f'<text x="{cx_trap}" y="{CELL_Y+17}" text-anchor="middle" '
              f'font-size="8" fill="#c0392b">(TRAPed)</text>')
# More D1 receptors
for dx in [-16, -8, 0, 8, 16]:
    svg2_body += (f'<circle cx="{cx_trap+dx}" cy="{CELL_Y-31}" r="3.5" '
                  f'fill="#c0392b" opacity="0.7"/>')
svg2_body += (f'<text x="{cx_trap}" y="{CELL_Y-40}" text-anchor="middle" '
              f'font-size="8.5" fill="#c0392b" font-weight="700">Drd1 ↑↑ (higher expression)</text>')
# More glutamate arrows
for dx in [-20, -10, 0, 10, 20]:
    svg2_body += (f'<line x1="{cx_trap+dx}" y1="{CELL_Y-60}" x2="{cx_trap+dx}" '
                  f'y2="{CELL_Y-32}" stroke="#c0392b" stroke-width="2.2"/>'
                  f'<polygon points="{cx_trap+dx},{CELL_Y-32} {cx_trap+dx-3},{CELL_Y-38} '
                  f'{cx_trap+dx+3},{CELL_Y-38}" fill="#c0392b"/>')
svg2_body += (f'<text x="{cx_trap}" y="{CELL_Y-65}" text-anchor="middle" '
              f'font-size="8.5" fill="#c0392b" font-weight="700">Glut. input ↑↑ (MOp + Thal)</text>')

svg2_body += (f'<line x1="{cx_trap}" y1="{CELL_Y+30}" x2="{cx_trap}" y2="{CELL_Y+54}" '
              f'stroke="#c0392b" stroke-width="3.5"/>'
              f'<polygon points="{cx_trap},{CELL_Y+54} {cx_trap-5},{CELL_Y+46} '
              f'{cx_trap+5},{CELL_Y+46}" fill="#c0392b"/>'
              f'<text x="{cx_trap}" y="{CELL_Y+64}" text-anchor="middle" '
              f'font-size="8.5" fill="#c0392b" font-weight="700">→ SNr (hyperactive)</text>')
svg2_body += (f'<text x="{cx_trap}" y="{CELL_Y+82}" text-anchor="middle" '
              f'font-size="8.5" fill="#c0392b">Levodopa → excessive D1 activation</text>'
              f'<text x="{cx_trap}" y="{CELL_Y+94}" text-anchor="middle" '
              f'font-size="8.5" fill="#c0392b">→ disinhibited thalamocortical circuit</text>')

# Outcome box
OB_Y = PANEL_Y + 20; OB_H = 160
svg2_body += (
    f'<text x="{cx_outcome}" y="{PANEL_Y+30}" text-anchor="middle" font-size="10.5" '
    f'fill="#922b21" font-weight="700">Levodopa-induced</text>'
    f'<text x="{cx_outcome}" y="{PANEL_Y+44}" text-anchor="middle" font-size="10.5" '
    f'fill="#922b21" font-weight="700">Dyskinesia (LID)</text>'
    f'<text x="{cx_outcome}" y="{PANEL_Y+66}" text-anchor="middle" font-size="9" fill="#555">'
    f'Involuntary movements:</text>'
    f'<text x="{cx_outcome}" y="{PANEL_Y+80}" text-anchor="middle" font-size="9" fill="#555">'
    f'abnormal involuntary</text>'
    f'<text x="{cx_outcome}" y="{PANEL_Y+94}" text-anchor="middle" font-size="9" fill="#555">'
    f'movements (AIMs)</text>'
    f'<text x="{cx_outcome}" y="{PANEL_Y+118}" text-anchor="middle" font-size="9" fill="#888">'
    f'FosTRAP identifies the</text>'
    f'<text x="{cx_outcome}" y="{PANEL_Y+131}" text-anchor="middle" font-size="9" fill="#888">'
    f'specific dMSNs activated</text>'
    f'<text x="{cx_outcome}" y="{PANEL_Y+144}" text-anchor="middle" font-size="9" fill="#888">'
    f'during LID episodes</text>'
    f'<text x="{cx_outcome}" y="{PANEL_Y+158}" text-anchor="middle" font-size="8.5" fill="#c0392b">'
    f'Drd1 ↑, cortical syn. ↑</text>'
    f'<text x="{cx_outcome}" y="{PANEL_Y+170}" text-anchor="middle" font-size="8.5" fill="#c0392b">'
    f'= TRAPed cell signature</text>'
)

# Arrow between normal and TRAPed (showing levodopa transforms)
mid_x = (cx_normal + cx_trap) // 2
svg2_body += (
    f'<text x="{mid_x}" y="{CELL_Y-10}" text-anchor="middle" font-size="9" fill="#e67e22">'
    f'Levodopa</text>'
    f'<text x="{mid_x}" y="{CELL_Y+2}" text-anchor="middle" font-size="9" fill="#e67e22">'
    f'reveals</text>'
    f'<line x1="{cx_normal+50}" y1="{CELL_Y}" x2="{cx_trap-60}" y2="{CELL_Y}" '
    f'stroke="#e67e22" stroke-width="1.5" stroke-dasharray="4,3"/>'
    f'<polygon points="{cx_trap-60},{CELL_Y} {cx_trap-68},{CELL_Y-4} {cx_trap-68},{CELL_Y+4}" '
    f'fill="#e67e22"/>'
)

# Key finding callout
KF_Y = PANEL_Y + PANEL_H + 20
svg2_body += (
    f'<rect x="30" y="{KF_Y}" width="{FW2-60}" height="40" rx="4" '
    f'fill="#fff5f5" stroke="#c0392b" stroke-width="1.2"/>'
    f'<text x="{FW2//2}" y="{KF_Y+15}" text-anchor="middle" font-size="9.5" '
    f'fill="#c0392b" font-weight="700">Ryan et al. 2024 (Cell Reports) key finding:</text>'
    f'<text x="{FW2//2}" y="{KF_Y+30}" text-anchor="middle" font-size="8.5" fill="#444">'
    f'TRAPed dMSNs have higher Drd1 expression, greater dopamine-dependent excitability, and more excitatory input from '
    f'motor cortex and thalamus</text>'
)

svg2 = f"""<svg viewBox="0 0 {FW2} {FH2}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW2//2}" y="22" text-anchor="middle" font-size="13" font-weight="700" fill="#222">
    dMSN Subpopulation Heterogeneity Underlies Levodopa-Induced Dyskinesia
  </text>
  <text x="{FW2//2}" y="40" text-anchor="middle" font-size="10" fill="#666">
    Schematic based on Ryan et al. (2024) Cell Reports · FosTRAP identifies dyskinesia-activated dMSNs in dorsal striatum
  </text>
  {svg2_body}
  <text x="{FW2//2}" y="{FH2-6}" text-anchor="middle" font-size="8.5" fill="#aaa">
    Ryan MB et al. (2024) Cell Rep 43:114483 · Nelson lab, UCSF
  </text>
</svg>"""

with open(os.path.join(OUT, "dmSN_mechanism.svg"), "w") as f:
    f.write(svg2)
print("Wrote dmSN_mechanism.svg")
