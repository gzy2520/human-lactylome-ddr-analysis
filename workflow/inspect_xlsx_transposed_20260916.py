#!/usr/bin/env python3
"""Read-only transposed dump of an xlsx sheet (stdlib only; the server has no pandas/openpyxl)."""
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

M = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
NS = {'m': M[1:-1]}
path = Path(sys.argv[1])


def shared_strings(zf):
    if 'xl/sharedStrings.xml' not in zf.namelist():
        return []
    root = ET.fromstring(zf.read('xl/sharedStrings.xml'))
    return [''.join(t.text or '' for t in si.iter(M + 't')) for si in root.findall('m:si', NS)]


def col_index(ref):
    n = 0
    for ch in ''.join(c for c in ref if c.isalpha()):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


with zipfile.ZipFile(path) as zf:
    strings = shared_strings(zf)
    for sheet in sorted(n for n in zf.namelist() if n.startswith('xl/worksheets/sheet')):
        root = ET.fromstring(zf.read(sheet))
        rows = []
        for row in root.iter(M + 'row'):
            cells = {}
            for c in row.findall('m:c', NS):
                v = c.find('m:v', NS)
                if v is None or v.text is None:
                    continue
                cells[col_index(c.get('r') or '')] = strings[int(v.text)] if c.get('t') == 's' else v.text
            if cells:
                rows.append([cells.get(i, '') for i in range(max(cells) + 1)])
        ncol = max(len(r) for r in rows)
        rows = [r + [''] * (ncol - len(r)) for r in rows]
        print(f'### {sheet}: {len(rows)} rows x {ncol} cols')
        for i in range(ncol):
            print(f'{i+1}|{rows[0][i]} :: ' + ' | '.join(r[i] for r in rows[1:]))
