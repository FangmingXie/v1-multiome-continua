"""Warn if README.md does not link every figure in figures/, or links a figure that does not exist.

The README's script tables are the report: each figure script's row links the PDFs it draws. run_all.sh
runs this check last and prints a WARNING for each mismatch (a renamed, new or removed figure, or a stale
file); it never fails the run. Reads README.md and the file list of figures/; writes nothing.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
README = REPO_ROOT / 'README.md'
FIGURES = REPO_ROOT / 'figures'

links = set(re.findall(r'\]\((figures/[^)]+\.pdf)\)', README.read_text()))
pdfs = {p.relative_to(REPO_ROOT).as_posix() for p in FIGURES.rglob('*.pdf')}
unlinked = sorted(pdfs - links)
broken = sorted(links - pdfs)
for f in unlinked:
    print(f'WARNING: {f} is not linked from README.md (new or renamed figure, or a stale file)')
for f in broken:
    print(f'WARNING: README.md links {f}, which does not exist')
if not (unlinked or broken):
    print(f'README.md links all {len(pdfs)} figures in figures/')
