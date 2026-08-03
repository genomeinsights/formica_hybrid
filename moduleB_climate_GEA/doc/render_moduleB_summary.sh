#!/usr/bin/env bash
# Render moduleB_summary.md -> styled PDF (pandoc + xelatex). Run from repo root.
set -euo pipefail
cd "$(dirname "$0")"
pandoc moduleB_summary.md -o moduleB_summary.pdf \
  --pdf-engine=xelatex --resource-path=. \
  --include-in-header=moduleB_summary_header.tex \
  -V geometry:a4paper,margin=1.6cm -V fontsize=10pt \
  -V colorlinks=true -V linkcolor=Accent -V urlcolor=Accent \
  --toc -V toc-title="Contents"
echo "wrote moduleB_summary.pdf"
