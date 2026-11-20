"""Check that README.md links every figure in figures/, and only figures that exist.

The README's script tables are the report: each figure script's row links the PDFs it draws. run_all.sh
runs this check last, so a renamed, new or removed figure fails the run until the README is updated.
Reads README.md and the file list of figures/; writes nothing.
"""

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
README = REPO_ROOT / 'README.md'
FIGURES = REPO_ROOT / 'figures'

links = set(re.findall(r'\]\((figures/[^)]+\.pdf)\)', README.read_text()))
pdfs = {p.relative_to(REPO_ROOT).as_posix() for p in FIGURES.rglob('*.pdf')}
unlinked = sorted(pdfs - links)
broken = sorted(links - pdfs)
for f in unlinked:
    print(f'ERROR: {f} is not linked from README.md (new or renamed figure, or a stale file)', file=sys.stderr)
for f in broken:
    print(f'ERROR: README.md links {f}, which does not exist', file=sys.stderr)
if unlinked or broken:
    sys.exit(1)
print(f'README.md links all {len(pdfs)} figures in figures/')
