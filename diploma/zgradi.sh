#!/usr/bin/env bash
# Zgradi diploma/izhod/diploma.docx iz poglavij v Markdownu.
# Uporaba (iz korena repozitorija ali iz mape diploma):  bash diploma/zgradi.sh
# Oblika sledi Vzorcu UL FKKT 2026 (predloga/naredi_referenco.py); naslovnici, izjavo o avtorstvu,
# lihe strani za vsako poglavje in številčenje strani od Uvoda doda predloga/razdelki.py.
set -euo pipefail
cd "$(dirname "$0")"
REF="${REF:-predloga/fkkt_referenca.docx}"
[ -f "$REF" ] || python3 predloga/naredi_referenco.py
mkdir -p izhod
pandoc metadata.yaml poglavja/*.md \
  --from markdown \
  --citeproc --bibliography literatura/literatura.bib --csl pravila/fkkt-ul.csl \
  --reference-doc "$REF" \
  --number-sections \
  --resource-path=.:slike \
  -o izhod/diploma.docx
python3 predloga/razdelki.py izhod/diploma.docx metadata.yaml
echo "Zgrajeno: $(pwd)/izhod/diploma.docx"
