#!/usr/bin/env python3
"""Convertit P1.json / P2.json en entrée exacte pour cgal_check.cpp (sommets rationnels recalculés
à partir des plans par la règle de Cramer, sans passer par les sommets stockés dans le JSON)."""
import json
import sys
from fractions import Fraction as Fr


def det3(a, b, c):
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])


d = json.load(open(sys.argv[1]))
planes = [[Fr(x) for x in p] for p in d["planes"]]
out = [f"{len(d['tri'])} {len(d['faces'])}"]
for fs in d["tri"]:
    rows = [planes[f] for f in fs]
    D = det3(*[r[:3] for r in rows])
    pt = [det3(*[[r[3] if j == k else r[j] for j in range(3)] for r in rows]) / D for k in range(3)]
    out.append(" ".join(f"{x.numerator}/{x.denominator}" for x in pt))
for c in d["faces"]:
    out.append(f"{len(c)} " + " ".join(map(str, c)))
print("\n".join(out))
