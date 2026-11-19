# ══════════════════════════════════════════════════════════════════════════════
# 07_table4_load_partitioning.R
# ══════════════════════════════════════════════════════════════════════════════
# Table 4: Suspended sediment load partitioned by hysteresis direction and
# monitoring site for enrichment events (classes I–IV).
#
# Network totals use outlet sites only (the most downstream monitored site on
# each river branch) to avoid double-counting at nested sites (Reviewer 2,
# Comment 6). PKO-UW and PKO-UE are nested within PKO-D, and PNG-TR1 within
# PNG-M; these are reported as a share of their parent outlet site's load.
#
# Input:  data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv
# Output: tables/Table4_Load_By_Direction.tex
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)

df <- read_csv("data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv",
               show_col_types = FALSE)

# ── 1. Site lookup, nesting structure and display order ─────────────────────
site_lookup <- tribble(
  ~Site,    ~Site_Code, ~Sub_Catchment, ~Parent,
  "Tuangi", "PKO-D",    "Pokopoko",     NA,
  "Inanga", "PKO-UW",   "Pokopoko",     "PKO-D",
  "Kokopu", "PKO-UE",   "Pokopoko",     "PKO-D",
  "Matata", "PNG-M",    "Pongakawa",    NA,
  "Tuna",   "PNG-TR1",  "Pongakawa",    "PNG-M",
  "Kokako", "MGT",      "Mangatoetoe",  NA,
  "Pipi",   "PUN",      "Puanene",      NA,
  "Kotuku", "WHR-D",    "Wharere",      NA
)
site_order <- site_lookup$Site_Code

# ── 2. Enrichment events with direction ─────────────────────────────────────
df <- df |>
  inner_join(site_lookup, by = "Site") |>
  filter(Zuec_class %in% 1:4) |>
  mutate(Direction = if_else(Zuec_class %in% c(1, 2), "CW", "ACW"))

# ── 3. Site x direction summary ─────────────────────────────────────────────
site_dir <- df |>
  group_by(Sub_Catchment, Site_Code, Parent) |>
  summarise(
    n        = n(),
    load_cw  = sum(Load_kg[Direction == "CW"],  na.rm = TRUE),
    load_acw = sum(Load_kg[Direction == "ACW"], na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(
    load_total      = load_cw + load_acw,
    mean_event_load = load_total / n,
    pct_cw          = 100 * load_cw  / load_total,
    pct_acw         = 100 * load_acw / load_total,
    is_outlet       = is.na(Parent)
  )

# Network total from outlet sites only
outlets       <- filter(site_dir, is_outlet)
network_total <- sum(outlets$load_total)

# Outlet sites: % of network; nested sites: % of parent outlet's load
site_dir <- site_dir |>
  left_join(select(outlets, Parent = Site_Code, parent_load = load_total),
            by = "Parent") |>
  mutate(
    pct_network = if_else(is_outlet, 100 * load_total / network_total,
                          100 * load_total / parent_load),
    Site_Code   = factor(Site_Code, levels = site_order)
  ) |>
  arrange(Site_Code)

# ── 4. Console summary ──────────────────────────────────────────────────────
cat("Network total (outlet sites):", round(network_total / 1000), "t\n")
cat("CW share:",  round(100 * sum(outlets$load_cw)  / network_total, 1), "%\n")
cat("ACW share:", round(100 * sum(outlets$load_acw) / network_total, 1), "%\n\n")
print(site_dir |> select(Site_Code, n, load_total, pct_cw, pct_acw, pct_network),
      width = Inf)

# ── 5. LaTeX table ──────────────────────────────────────────────────────────
t  <- function(kg) formatC(kg / 1000, format = "f", digits = 0, big.mark = ",")
p1 <- function(x)  formatC(x, format = "f", digits = 1)

body <- character(0)
prev_sc <- ""
for (i in seq_len(nrow(site_dir))) {
  r <- site_dir[i, ]
  sc_label  <- if (r$Sub_Catchment != prev_sc) r$Sub_Catchment else ""
  prev_sc   <- r$Sub_Catchment
  site_lab  <- if (r$is_outlet) as.character(r$Site_Code) else paste0(r$Site_Code, "$^{b}$")
  net_lab   <- if (r$is_outlet) p1(r$pct_network) else paste0("(", p1(r$pct_network), "\\%$^{c}$)")
  body <- c(body, paste(
    sc_label, site_lab, r$n, t(r$load_cw), t(r$load_acw), t(r$load_total),
    t(r$mean_event_load), p1(r$pct_cw), p1(r$pct_acw), net_lab,
    sep = " & "
  ) |> paste0(" \\\\"))
}

totals <- paste0(
  "\\textbf{Total (outlet sites)} & & ",
  "\\textbf{", sum(outlets$n), "} & ",
  "\\textbf{", t(sum(outlets$load_cw)), "} & ",
  "\\textbf{", t(sum(outlets$load_acw)), "} & ",
  "\\textbf{", t(network_total), "} & ",
  "\\textbf{", t(network_total / sum(outlets$n)), "} & ",
  "\\textbf{", p1(100 * sum(outlets$load_cw)  / network_total), "} & ",
  "\\textbf{", p1(100 * sum(outlets$load_acw) / network_total), "} & ",
  "\\textbf{100.0} \\\\"
)

caption <- paste0(
  "\\caption{Suspended sediment load partitioned by hysteresis direction and monitoring site ",
  "for enrichment events (classes I--IV). CW = clockwise (classes I+II); ACW = anticlockwise ",
  "(classes III+IV). Network (\\%) is calculated from outlet sites only (the most downstream ",
  "monitored site on each river branch) to avoid double-counting sediment at nested upstream ",
  "sites. PKO-UW and PKO-UE are nested within PKO-D, and PNG-TR1 is nested within PNG-M; loads ",
  "at these sites are reported as a share of their respective outlet site's total rather than ",
  "as an independent network contribution. Absolute values reflect site-specific monitoring ",
  "coverage and should be interpreted as indicative of relative contributions rather than ",
  "precise load budgets.}"
)

latex <- c(
  "\\begin{sidewaystable}", "\\centering", caption,
  "\\label{tab:load_by_direction}", "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{llrrrrrrrr}", "\\toprule",
  paste0("Sub-catchment & Site & $n$ & CW load (t) & ACW load (t) & Total load (t) & ",
         "Mean event load (t) & CW (\\%) & ACW (\\%) & Network (\\%)$^{a}$ \\\\"),
  "\\midrule", body, "\\midrule", totals, "\\bottomrule", "\\end{tabular}", "",
  "\\vspace{0.5em}", "\\raggedright", "\\footnotesize",
  paste0("$^{a}$ Calculated from outlet sites only (PKO-D, PNG-M, MGT, PUN, WHR-D); nested ",
         "upstream sites are excluded from this total to avoid double-counting. \\\\"),
  paste0("$^{b}$ Nested within the downstream outlet site for its sub-catchment (PKO-UW and ",
         "PKO-UE within PKO-D; PNG-TR1 within PNG-M). \\\\"),
  paste0("$^{c}$ Expressed as a percentage of the parent outlet site's total load, not as an ",
         "independent network contribution."),
  "\\end{sidewaystable}"
)

dir.create("tables", showWarnings = FALSE)
writeLines(latex, "tables/Table4_Load_By_Direction.tex")