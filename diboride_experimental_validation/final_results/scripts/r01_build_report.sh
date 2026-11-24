#!/bin/bash
# FINAL_REPORT.md -> FINAL_REPORT.docx (pandoc 3.9 from the isolated report venv, figures embedded, citations kept as
# text) -> FINAL_REPORT.pdf (LibreOffice headless). Run from anywhere.
set -e
FR=/home/ichasekatz/BLADE/BLADE/diboride_experimental_validation/final_results
PANDOC=$(/home/ichasekatz/BLADE/mlip_envs/report/.venv/bin/python -c "import pypandoc; print(pypandoc.get_pandoc_path())")
cd $FR
$PANDOC FINAL_REPORT.md -f gfm+tex_math_dollars -t docx --resource-path=$FR:$FR/figures --toc --toc-depth=2 \
  -o FINAL_REPORT.docx
libreoffice --headless --convert-to pdf --outdir $FR FINAL_REPORT.docx > /dev/null 2>&1
ls -la FINAL_REPORT.*
