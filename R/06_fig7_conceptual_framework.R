# ══════════════════════════════════════════════════════════════════════════════
# 06_fig7_conceptual_framework.R
# ══════════════════════════════════════════════════════════════════════════════
# Figure 7: Conceptual framework linking landscape controls to hysteresis
# regime and management responses.
#
# Output: figures/Fig7_Conceptual_Framework.png
# Dependencies: DiagrammeR, DiagrammeRsvg, rsvg
# ══════════════════════════════════════════════════════════════════════════════

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

fig7 <- grViz("
digraph conceptual_framework {

  graph [
    layout = dot, rankdir = TB, fontname = 'Arial',
    bgcolor = '#f7f7f7', pad = 0.5, nodesep = 0.6, ranksep = 0.8
  ]

  node [
    fontname = 'Arial', fontsize = 11, shape = rectangle,
    style = 'filled,rounded', penwidth = 1.5, margin = '0.2,0.15'
  ]

  edge [
    fontname = 'Arial', fontsize = 9, color = '#555555',
    penwidth = 1.5, arrowsize = 0.8, arrowhead = vee
  ]

  { rank = same; soil; landuse; topo }

  soil [
    label = 'SOIL PROPERTIES\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 Erodibility\\n\\u2022 Drainage class\\n\\u2022 Organic content',
    fillcolor = '#d8b48a', fontcolor = '#3d2b1f', color = '#a0784a'
  ]

  landuse [
    label = 'LAND USE\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 Disturbance intensity\\n\\u2022 Vegetation cover\\n\\u2022 Drainage modification',
    fillcolor = '#d8b48a', fontcolor = '#3d2b1f', color = '#a0784a'
  ]

  topo [
    label = 'TOPOGRAPHY &\\nCONNECTIVITY\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 Slope gradient\\n\\u2022 Flow pathway dominance\\n\\u2022 Drainage network density',
    fillcolor = '#d8b48a', fontcolor = '#3d2b1f', color = '#a0784a'
  ]

  hydro [
    label = 'EVENT HYDROLOGY\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 Runoff magnitude and timing  \\u2022  Hydrograph flashiness\\n\\u2022 Baseflow contribution  \\u2022  Antecedent conditions',
    fillcolor = '#c6dbef', fontcolor = '#08306b', color = '#4393c3'
  ]

  sediment [
    label = 'SEDIMENT SUPPLY & DELIVERY\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 Near-channel source availability  \\u2022  Hillslope connectivity\\n\\u2022 Transport capacity',
    fillcolor = '#f0c08a', fontcolor = '#3d1f00', color = '#b35806'
  ]

  { rank = same; clockwise; figureeight; anticlockwise }

  clockwise [
    label = 'CLOCKWISE\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 TSS peaks before Q\\n\\u2022 Near-channel source\\n   dominated\\n\\u2022 Efficient delivery',
    fillcolor = '#f4a582', fontcolor = '#67000d', color = '#d73027'
  ]

  figureeight [
    label = 'FIGURE-EIGHT\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 Multiple source activation\\n\\u2022 Direction shifts within event\\n\\u2022 Not predictable from\\n   static predictors',
    fillcolor = '#d9d9d9', fontcolor = '#333333', color = '#878787', margin = '0.3,0.2'
  ]

  anticlockwise [
    label = 'ANTICLOCKWISE\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 TSS peaks after Q\\n\\u2022 Distal or stored\\n   sediment sources\\n\\u2022 Transport-limited',
    fillcolor = '#92c5de', fontcolor = '#08306b', color = '#2166ac'
  ]

  { rank = same; mgmt_cw; mgmt_mixed; mgmt_acw }

  mgmt_cw [
    label = 'NEAR-CHANNEL\\nINTERVENTIONS\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 Riparian buffers\\n\\u2022 Stock exclusion\\n\\u2022 Bank stabilisation\\n\\u2022 Drain management',
    fillcolor = '#c7e9c0', fontcolor = '#00441b', color = '#41ab5d'
  ]

  mgmt_mixed [
    label = 'COMBINED\\nINTERVENTIONS\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 Retention & attenuation\\n   infrastructure\\n\\u2022 Adaptive management\\n\\u2022 Broad-spectrum measures',
    fillcolor = '#e0e0e0', fontcolor = '#333333', color = '#878787', margin = '0.3,0.2'
  ]

  mgmt_acw [
    label = 'HILLSLOPE &\\nCATCHMENT\\nINTERVENTIONS\\n\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\u2500\\n\\u2022 Hillslope erosion control\\n\\u2022 Forestry management\\n\\u2022 Detention structures\\n\\u2022 Upstream stabilisation',
    fillcolor = '#c7e9c0', fontcolor = '#00441b', color = '#41ab5d'
  ]

  soil -> hydro; landuse -> hydro; topo -> hydro
  hydro -> sediment
  sediment -> clockwise; sediment -> figureeight; sediment -> anticlockwise
  clockwise -> mgmt_cw; figureeight -> mgmt_mixed; anticlockwise -> mgmt_acw
}
")

dir.create("figures", showWarnings = FALSE)
svg_content <- export_svg(fig7)
writeLines(svg_content, "figures/Fig7_Conceptual_Framework.svg")
rsvg_png(svg = "figures/Fig7_Conceptual_Framework.svg",
         file = "figures/Fig7_Conceptual_Framework.png",
         width = 3000, height = 4000)
cat("Figure 7 saved\n")
