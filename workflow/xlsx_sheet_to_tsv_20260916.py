#!/usr/bin/env python3
"""Convert the first worksheet of an xlsx workbook to TSV (stdlib only).

The server has no pandas/openpyxl, so the workbook is read through zipfile +
ElementTree. All rows are emitted verbatim, including the leading note rows, so
that the R side can locate the real header row itself.

Usage: xlsx_sheet_to_tsv_20260916.py <input.xlsx> <output.tsv>
"""
import sys
import zipfile
import xml.etree.ElementTree as ET

M = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
NS = {'m': M[1:-1]}


def col_index(ref):
    n = 0
    for ch in ''.join(c for c in ref if c.isalpha()):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def shared_strings(zf):
    if 'xl/sharedStrings.xml' not in zf.namelist():
        return []
    root = ET.fromstring(zf.read('xl/sharedStrings.xml'))
    return [''.join(t.text or '' for t in si.iter(M + 't')) for si in root.findall('m:si', NS)]


def clean(text):
    # a cell value containing a newline or tab would break the row/column grid
    return text.replace('\r\n', ' ').replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')


def main(src, dst):
    with zipfile.ZipFile(src) as zf:
        strings = shared_strings(zf)
        sheets = sorted(n for n in zf.namelist() if n.startswith('xl/worksheets/sheet'))
        if not sheets:
            raise SystemExit('no worksheet found in %s' % src)
        root = ET.fromstring(zf.read(sheets[0]))
    grid = []
    for row in root.iter(M + 'row'):
        cells = {}
        for c in row.findall('m:c', NS):
            v = c.find('m:v', NS)
            if v is None or v.text is None:
                continue
            cells[col_index(c.get('r') or '')] = strings[int(v.text)] if c.get('t') == 's' else v.text
        grid.append(cells)
    width = max((max(c) + 1 for c in grid if c), default=0)
    # every row is padded to the full sheet width so that the column count is
    # unambiguous for the reader
    with open(dst, 'w', encoding='utf-8') as fh:
        for cells in grid:
            fh.write('\t'.join(clean(cells.get(i, '')) for i in range(width)) + '\n')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
