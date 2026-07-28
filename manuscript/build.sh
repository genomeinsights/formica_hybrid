#!/usr/bin/env bash
# Build the manuscript documents in dependency order.
#
# supplementary_methods_ld_modules_ABC.tex cross-references the separate
# Supplementary Figures document via \externaldocument{supplementary_figures}
# (xr package), so the figures document MUST be compiled first to produce its
# .aux; otherwise the (Fig. S#) references come out as "??".
#
# Requires: latexmk (TeX Live). Run from anywhere.
set -e
cd "$(dirname "$0")"

order=(
  supplementary_figures            # produces the .aux that xr reads
  supplementary_tables
  main_text_results_discussion_and_figure_plan
  supplementary_methods_ld_modules_ABC   # last: resolves the xr cross-refs
)
for f in "${order[@]}"; do
  [ -f "$f.tex" ] && latexmk -pdf -interaction=nonstopmode "$f.tex"
done

latexmk -c >/dev/null 2>&1        # remove aux/log/etc., keep the PDFs
echo "Built: $(printf '%s.pdf ' "${order[@]}")"
