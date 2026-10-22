"""Write the PatternedHexMeshGenerator pattern of the GFMR core from the shared channel map.

The pattern is a 23-row hexagon (rings 0-11 around the central cell). Row r (top to bottom) is
y = +11 ... -11 lattice rows; within a row the entries run in +x. Unit-cell indices:
  0 inner fuel channel   1 outer fuel channel   2 control-rod channel   3 solid graphite
(3 is used for the 18 ring-11 positions that hold no channel).

Usage: python make_pattern.py > core_pattern.txt  (the text is pasted into gfmr_mesh.i)
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'common'))
import gfmr  # noqa: E402

CODE = {'inner': 0, 'outer': 1, 'cr': 2, None: 3}


def rows():
    out = []
    for j in range(11, -12, -1):
        row = []
        for i in range(-11, 12):
            if gfmr.ring(i, j) <= 11:
                row.append(CODE[gfmr.channel_type(i, j)])
        out.append(row)
    return out


if __name__ == '__main__':
    r = rows()
    n = len(r)
    lines = []
    for k, row in enumerate(r):
        pad = ' ' * abs(k - n // 2)
        end = ';' if k < n - 1 else ''
        lines.append('    ' + pad + ' '.join(str(c) for c in row) + end)
    print('\n'.join(lines))
    flat = [c for row in r for c in row]
    sys.stderr.write('inner %d outer %d cr %d graphite %d\n' % tuple(flat.count(c) for c in range(4)))
