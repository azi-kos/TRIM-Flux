#!/usr/bin/env bash
# Zgradi diploma/izhod/diploma.docx iz poglavij v Markdownu.
# Uporaba (iz korena repozitorija ali iz mape diploma):  bash diploma/zgradi.sh
# Če obstaja uradni vzorec FKKT (pravila/originali/10_FKKT_Vzorec_diploma_magisterij_2024.docx),
# ga lahko uporabiš kot referenco:  REF=pravila/originali/10_FKKT_Vzorec_diploma_magisterij_2024.docx bash diploma/zgradi.sh
set -euo pipefail
cd "$(dirname "$0")"
REF="${REF:-predloga/fkkt_referenca.docx}"
[ -f "$REF" ] || python3 predloga/naredi_referenco.py
mkdir -p izhod
pandoc metadata.yaml poglavja/*.md \
  --from markdown \
  --citeproc --bibliography literatura/literatura.bib --csl pravila/fkkt-ul.csl \
  --reference-doc "$REF" \
  --toc --toc-depth=3 -M toc-title="Kazalo" \
  --number-sections \
  --resource-path=.:slike \
  -o izhod/diploma.docx
echo "Zgrajeno: $(pwd)/izhod/diploma.docx"
