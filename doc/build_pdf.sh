#!/usr/bin/env bash
# Build a single combined PDF from the doc/*.md sources.
# Concatenates README + the four stage documents (in order), embeds every
# referenced ../Figures/*.png right after its caption line, and renders with
# pandoc + XeLaTeX (needed for the Unicode maths/glyphs the text uses).
#
# Requires: pandoc, xelatex (TeX Live).  Run from anywhere: it cd's to doc/.
set -euo pipefail
cd "$(dirname "$0")"

out="pipeline_documentation.pdf"
tmp="$(mktemp -t formica_doc).md"
trap 'rm -f "$tmp"' EXIT

# --- document metadata / title page ---
cat > "$tmp" <<'YAML'
---
title: "Formica hybrid pipeline"
subtitle: "LD complexity reduction and the descriptive analysis of ancestry sorting (Stage 0 through Modules A, B, C)"
date: 2026-07-27
toc: true
toc-depth: 2
geometry: margin=1in
mainfont: "Helvetica Neue"
monofont: "Menlo"
colorlinks: true
linkcolor: RoyalBlue
urlcolor: RoyalBlue
---
YAML

# --- body: each source file, with figures embedded after their caption line ---
for f in README.md 01_*.md 02_*.md 03_*.md 04_*.md 05_*.md; do
  awk '{
    print
    # match only concrete single-file paths: excludes {..} grids and Chr*.. shorthands
    if (match($0, /\.\.\/Figures\/[A-Za-z0-9_]+\.png/)) {
      p = substr($0, RSTART, RLENGTH)
      printf "\n![](%s){width=85%%}\n\n", p
    }
  }' "$f" >> "$tmp"
  printf '\n\n\\newpage\n\n' >> "$tmp"
done

# Arrows that the text font lacks -> render them via the math font instead.
# (These glyphs appear only in prose/tables here, never inside code blocks.)
perl -i -CSD -pe 's/\x{2192}/\$\\rightarrow\$/g; s/\x{2194}/\$\\leftrightarrow\$/g;' "$tmp"

pandoc "$tmp" -o "$out" \
  --pdf-engine=xelatex \
  --toc --toc-depth=2 \
  --include-in-header=pdf-header.tex \
  -V geometry:margin=1in

echo "wrote doc/$out"
