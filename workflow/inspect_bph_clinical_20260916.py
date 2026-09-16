#!/usr/bin/env python3
"""Read-only dump of the GSE132714 BPH clinical annotation table (stdlib only; server has no pandas)."""
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
path = Path(sys.argv[1])


def shared_strings(zf):
    if 'xl/sharedStrings.xml' not in zf.namelist():
        return []
    root = ET.fromstring(zf.read('xl/sharedStrings.xml'))
    out = []
    for si in root.findall('m:si', NS):
        out.append(''.join(t.text or '' for t in si.iter(
            '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t')))
    return out


def col_index(ref):
    letters = ''.join(c for c in ref if c.isalpha())
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


with zipfile.ZipFile(path) as zf:
    strings = shared_strings(zf)
    sheets = [n for n in zf.namelist() if n.startswith('xl/worksheets/sheet')]
    for sheet in sorted(sheets):
        root = ET.fromstring(zf.read(sheet))
        rows = []
        for row in root.iter('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}row'):
            cells = {}
            for c in row.findall('m:c', NS):
                ref = c.get('r') or ''
                v = c.find('m:v', NS)
                if v is None or v.text is None:
                    continue
                val = strings[int(v.text)] if c.get('t') == 's' else v.text
                cells[col_index(ref)] = val
            if cells:
                width = max(cells) + 1
                rows.append([cells.get(i, '') for i in range(width)])
        print(f'### SHEET {sheet}: {len(rows)} rows x {max(len(r) for r in rows)} cols')
        header = rows[0]
        print('\t'.join(header))
        for r in rows[1:]:
            print('\t'.join(r))
