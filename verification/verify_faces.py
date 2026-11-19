#!/usr/bin/env python3
"""Vérification exacte (arithmétique rationnelle) d'une solution .sol de search.cpp.

Les plans flottants sont arrondis en rationnels ; les sommets sont recalculés
exactement comme intersections de 3 plans (donc les faces sont planes par construction).
On vérifie : surface fermée, genre, adjacence de toutes les paires de faces,
polygones simples, et absence de toute intersection entre faces hors arêtes communes.
Usage : verify_faces.py fichier.(sol|json) [denominateur_max] [--json sortie.json]
"""
import sys, json
from fractions import Fraction as Fr
from itertools import combinations


def load(path, den):
    if path.endswith(".json"):
        d = json.load(open(path))
        planes = [[Fr(x) for x in p] for p in d["planes"]]
        return len(planes), len(d["tri"]), planes, d["tri"], None, d["faces"]
    L = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    F, T = int(L[0][0]), int(L[0][1])
    planes = [[Fr(float(x)).limit_denominator(den) for x in L[1 + f]] for f in range(F)]
    tri = [list(map(int, L[1 + F + t][:3])) for t in range(T)]
    adj = [list(map(int, L[1 + F + t][3:])) for t in range(T)]
    cyc = [list(map(int, L[1 + F + T + f][1:])) for f in range(F)]
    return F, T, planes, tri, adj, cyc


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))


def meet(p, q, r):
    d = det3(p[:3], q[:3], r[:3])
    if d == 0:
        return None
    cols = lambda k: det3(*[[row[3] if j == k else row[j] for j in range(3)] for row in (p, q, r)])
    return tuple(cols(k) / d for k in range(3))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def proj2(n):
    """Projection sur 2 coordonnées (on retire l'axe dominant de la normale) : affine, exacte."""
    k = max(range(3), key=lambda i: abs(n[i]))
    idx = [i for i in range(3) if i != k]
    return lambda p: (p[idx[0]], p[idx[1]])


def orient(a, b, c):
    v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (v > 0) - (v < 0)


def on_seg(a, b, p):
    return orient(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def segs_meet(a, b, c, d):
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    if o1 != o2 and o3 != o4:
        return True
    return on_seg(a, b, c) or on_seg(a, b, d) or on_seg(c, d, a) or on_seg(c, d, b)


def in_closed_poly(poly, p):
    n = len(poly)
    for i in range(n):
        if on_seg(poly[i], poly[(i + 1) % n], p):
            return True
    inside = False
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        if (a[1] > p[1]) != (b[1] > p[1]):
            x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if p[0] < x:
                inside = not inside
    return inside


def verify(path, den=10**6, verbose=True):
    F, T, planes, tri, adj, cyc = load(path, den)
    err = []
    X = [meet(*[planes[f] for f in tri[t]]) for t in range(T)]
    if any(x is None for x in X):
        return ["sommet dégénéré"], None
    if len(set(X)) != T:
        err.append("sommets confondus")
    # combinatoire
    edges = {}
    for f in range(F):
        c = cyc[f]
        if any(f not in tri[t] for t in c):
            err.append(f"face {f}: sommet hors face")
        for k in range(len(c)):
            a, b = c[k], c[(k + 1) % len(c)]
            edges.setdefault((min(a, b), max(a, b)), []).append(f)
    if any(len(v) != 2 for v in edges.values()):
        err.append("arête non partagée par exactement 2 faces")
    E = len(edges)
    chi = T - E + F
    mult = {}
    for (a, b), fs in edges.items():
        key = tuple(sorted(fs))
        mult.setdefault(key, []).append((a, b))
    missing = [p for p in combinations(range(F), 2) if p not in mult]
    if missing:
        err.append(f"paires non adjacentes : {missing}")
    # orientabilité : tri/adj orientés de façon cohérente
    for t in range(T if adj else 0):
        for i in range(3):
            u = adj[t][i]
            p, q = tri[t][i], tri[t][(i + 1) % 3]
            if not any(tri[u][j] == q and tri[u][(j + 1) % 3] == p and adj[u][j] == t for j in range(3)):
                err.append("orientation incohérente")
    # polygones simples
    P2 = []
    for f in range(F):
        pr = proj2(planes[f][:3])
        poly = [pr(X[t]) for t in cyc[f]]
        P2.append((pr, poly))
        n = len(poly)
        for k in range(n):
            if orient(poly[k - 1], poly[k], poly[(k + 1) % n]) == 0:
                err.append(f"face {f}: sommets alignés")
            for l in range(k + 1, n):
                if l == k + 1 or (k == 0 and l == n - 1):
                    continue
                if segs_meet(poly[k], poly[(k + 1) % n], poly[l], poly[(l + 1) % n]):
                    err.append(f"face {f}: polygone non simple ({k},{l})")
    # intersections entre faces : exactement les arêtes communes
    for i, j in combinations(range(F), 2):
        ni, nj = planes[i][:3], planes[j][:3]
        d = cross(ni, nj)
        if d == (0, 0, 0):
            if all(a * planes[j][3] == b * planes[i][3] for a, b in zip(ni, nj)):
                err.append(f"faces {i},{j} coplanaires")
            continue
        # point de base sur la droite
        k = max(range(3), key=lambda a: abs(d[a]))
        aux = [Fr(0)] * 4
        aux[k] = Fr(1)
        base = meet(planes[i], planes[j], aux)
        par = lambda p: dot(sub(p, base), d) / dot(d, d)
        pt = lambda s: tuple(base[a] + s * d[a] for a in range(3))
        ev = set()
        for f, g in ((i, j), (j, i)):
            c = cyc[f]
            for kk in range(len(c)):
                A, B = X[c[kk]], X[c[(kk + 1) % len(c)]]
                sa = dot(planes[g][:3], A) - planes[g][3]
                sb = dot(planes[g][:3], B) - planes[g][3]
                if sa == 0:
                    ev.add(par(A))
                if sb == 0:
                    ev.add(par(B))
                if (sa > 0 and sb < 0) or (sa < 0 and sb > 0):
                    tt = sa / (sa - sb)
                    ev.add(par(tuple(A[a] + tt * (B[a] - A[a]) for a in range(3))))
        ev = sorted(ev)
        shared = [tuple(sorted((par(X[a]), par(X[b])))) for a, b in mult.get((i, j), [])]
        both = lambda s: in_closed_poly(P2[i][1], P2[i][0](pt(s))) and in_closed_poly(P2[j][1], P2[j][0](pt(s)))
        on_shared = lambda s: any(lo <= s <= hi for lo, hi in shared)
        for s in ev:
            if both(s) and not on_shared(s):
                err.append(f"faces {i},{j}: contact parasite")
        for s0, s1 in zip(ev, ev[1:]):
            mid = (s0 + s1) / 2
            if both(mid) != on_shared(mid):
                err.append(f"faces {i},{j}: intersection parasite")
        for a, b in combinations(shared, 2):
            if max(a[0], b[0]) <= min(a[1], b[1]):
                err.append(f"faces {i},{j}: arêtes communes chevauchantes")
    info = {
        "F": F, "V": T, "E": E, "chi": chi, "genus": (2 - chi) // 2,
        "degrees": [len(c) for c in cyc],
        "pairs_double": sum(1 for v in mult.values() if len(v) == 2),
        "pairs_triple": sum(1 for v in mult.values() if len(v) >= 3),
        "planes": [[str(x) for x in p] for p in planes],
        "vertices": [[str(x) for x in v] for v in X],
        "vertices_float": [[float(x) for x in v] for v in X],
        "faces": [[t + 1 for t in c] for c in cyc],
    }
    if verbose:
        print({k: info[k] for k in ("F", "V", "E", "chi", "genus", "degrees", "pairs_double", "pairs_triple")})
        print("VALIDE" if not err else "INVALIDE : " + "; ".join(sorted(set(err))[:8]))
    return err, info


if __name__ == "__main__":
    den = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 10**6
    err, info = verify(sys.argv[1], den)
    if "--json" in sys.argv and info and not err:
        json.dump(info, open(sys.argv[sys.argv.index("--json") + 1], "w"))
    sys.exit(1 if err else 0)
