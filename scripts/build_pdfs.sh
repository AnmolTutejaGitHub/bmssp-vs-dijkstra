#!/bin/sh
# Regenerate figures/tables from results/ and typeset both documents.
# Needs: .venv with matplotlib + numpy (see README), and tectonic.
set -e
cd "$(dirname "$0")/.."
PY=.venv/bin/python
[ -x "$PY" ] || PY=python3
"$PY" scripts/make_figures.py
(cd paper && tectonic main.tex)
(cd explanation && tectonic explanation.tex)
cp paper/main.pdf research-paper.pdf
cp explanation/explanation.pdf explanation.pdf
echo "wrote research-paper.pdf and explanation.pdf"
