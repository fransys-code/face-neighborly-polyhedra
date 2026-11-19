#!/usr/bin/env python3
"""Seconde vérification, indépendante de verify_faces.py : triangulation + tests exacts triangle/triangle.

Méthode (différente de verify_faces.py, qui raisonne face par face sur les droites d'intersection) :
  1. sommets recalculés par élimination de Gauss exacte sur les 3 plans de chaque sommet ;
     on vérifie aussi qu'aucun sommet ne tombe sur un 4e plan (position générale) ;
  2. chaque face est triangulée par « oreilles » en arithmétique exacte ;
  3. la surface triangulée est plongée si et seulement si deux triangles quelconques ne se
     rencontrent que selon leur « simplexe commun » (rien, un sommet commun, ou une arête commune).
     On le teste exactement : pour chaque arête de chaque triangle, l'ensemble de ses points qui
     sont dans l'autre triangle (intervalle de paramètres calculé en fractions) doit être contenu
     dans ce simplexe commun ;
  4. combinatoire recalculée indépendamment : variété (lien de chaque sommet = un cycle),
     orientabilité, caractéristique d'Euler, et pour chaque paire de faces le nombre d'arêtes communes.

Bibliothèque standard uniquement. Usage : verify_triangles.py fichier.json [...]
Format JSON : "planes" (a, b, c, d pour ax+by+cz=d, en chaînes rationnelles), "faces" (cycles
d'indices de sommets, à partir de 0) et "tri" (pour chaque sommet, les 3 faces qui s'y rencontrent).
"""
import json
import sys
from fractions import Fraction as Fr
from itertools import combinations


# ---------- algèbre exacte ----------

def solve3(rows):
    """Résout un système 3x3 (lignes [a, b, c, d]) par Gauss exact ; None si singulier."""
    m = [list(r) for r in rows]
    for col in range(3):
        piv = next((r for r in range(col, 3) if m[r][col] != 0), None)
        if piv is None:
            return None
        m[col], m[piv] = m[piv], m[col]
        for r in range(3):
            if r != col and m[r][col] != 0:
                f = m[r][col] / m[col][col]
                m[r] = [x - f * y for x, y in zip(m[r], m[col])]
    return tuple(m[i][3] / m[i][i] for i in range(3))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def sgn(x):
    return (x > 0) - (x < 0)


def orient2(a, b, c):
    return sgn((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


# ---------- triangulation par oreilles ----------

def drop_axis(n):
    k = max(range(3), key=lambda i: abs(n[i]))
    return [i for i in range(3) if i != k]


def ear_clip(idx, P2):
    """idx : indices du polygone ; P2 : coordonnées 2D exactes. Renvoie une liste de triangles."""
    poly = list(idx)
    area2 = sum(P2[poly[i]][0] * P2[poly[(i + 1) % len(poly)]][1] - P2[poly[(i + 1) % len(poly)]][0] * P2[poly[i]][1] for i in range(len(poly)))
    s = sgn(area2)
    assert s != 0, "polygone d'aire nulle"
    tris = []
    while len(poly) > 3:
        n = len(poly)
        for i in range(n):
            a, b, c = poly[i - 1], poly[i], poly[(i + 1) % n]
            if orient2(P2[a], P2[b], P2[c]) != s:
                continue  # coin rentrant ou plat
            ok = True
            for q in poly:
                if q in (a, b, c):
                    continue
                o1, o2, o3 = orient2(P2[a], P2[b], P2[q]), orient2(P2[b], P2[c], P2[q]), orient2(P2[c], P2[a], P2[q])
                if (o1 * s >= 0) and (o2 * s >= 0) and (o3 * s >= 0):
                    ok = False  # un autre sommet dans l'oreille (fermée)
                    break
            if ok:
                tris.append((a, b, c))
                poly.pop(i)
                break
        else:
            raise AssertionError("pas d'oreille : polygone non simple")
    tris.append(tuple(poly))
    return tris


# ---------- segment ∩ triangle fermé, exact ----------

def seg_tri(P, Q, A, B, C):
    """Ensemble des t ∈ [0,1] tels que P + t(Q-P) ∈ triangle fermé ABC : None ou (t0, t1)."""
    n = cross(sub(B, A), sub(C, A))
    dP, dQ = dot(n, sub(P, A)), dot(n, sub(Q, A))
    ax = drop_axis(n)
    pr = lambda X: (X[ax[0]], X[ax[1]])
    a2, b2, c2 = pr(A), pr(B), pr(C)
    s = orient2(a2, b2, c2)
    if dP == 0 and dQ == 0:
        # segment dans le plan du triangle : on le découpe par les trois demi-plans (convexe)
        p2, q2 = pr(P), pr(Q)
        lo, hi = Fr(0), Fr(1)
        for u, v in ((a2, b2), (b2, c2), (c2, a2)):
            f = lambda X: s * ((v[0] - u[0]) * (X[1] - u[1]) - (v[1] - u[1]) * (X[0] - u[0]))  # >= 0 à l'intérieur
            fp, fq = f(p2), f(q2)
            if fp < 0 and fq < 0:
                return None
            if fp < 0:
                lo = max(lo, fp / (fp - fq))
            elif fq < 0:
                hi = min(hi, fp / (fp - fq))
            if lo > hi:
                return None
        return (lo, hi)
    if (dP > 0 and dQ > 0) or (dP < 0 and dQ < 0):
        return None
    t = dP / (dP - dQ)
    X = tuple(p + t * (q - p) for p, q in zip(P, Q))
    x2 = pr(X)
    o1, o2, o3 = orient2(a2, b2, x2), orient2(b2, c2, x2), orient2(c2, a2, x2)
    if o1 * s >= 0 and o2 * s >= 0 and o3 * s >= 0:
        return (t, t)
    return None


def edge_inside_allowed(e, iv, shared):
    """e = (u, v) indices ; iv = intervalle trouvé ; shared = ensemble des sommets communs aux deux triangles."""
    if iv is None:
        return True
    u, v = e
    lo, hi = iv
    if u in shared and v in shared:
        return True  # l'arête entière est commune
    if u in shared and lo == 0 and hi == 0:
        return True  # ne touche qu'en u
    if v in shared and lo == 1 and hi == 1:
        return True  # ne touche qu'en v
    return False


def tri_pair_ok(t1, t2, X):
    shared = set(t1) & set(t2)
    for ta, tb in ((t1, t2), (t2, t1)):
        A, B, C = (X[i] for i in tb)
        for e in ((ta[0], ta[1]), (ta[1], ta[2]), (ta[2], ta[0])):
            if not edge_inside_allowed(e, seg_tri(X[e[0]], X[e[1]], A, B, C), shared):
                return False
    return True


def bbox(pts):
    return [min(p[i] for p in pts) for i in range(3)], [max(p[i] for p in pts) for i in range(3)]


# ---------- vérification complète ----------

def verify(path):
    d = json.load(open(path))
    planes = [[Fr(x) for x in p] for p in d["planes"]]
    faces = d["faces"]
    trip = d["tri"]
    F, V = len(planes), len(trip)
    err = []
    X = []
    for v, fs in enumerate(trip):
        x = solve3([planes[f] for f in fs])
        if x is None:
            return [f"sommet {v} : plans dépendants"], None
        X.append(x)
        others = [f for f in range(F) if f not in fs and dot(planes[f][:3], x) == planes[f][3]]
        if others:
            err.append(f"sommet {v} sur un 4e plan {others}")
    if len(set(X)) != V:
        err.append("sommets confondus")
    for f, cyc in enumerate(faces):
        for v in cyc:
            if f not in trip[v]:
                err.append(f"face {f} : sommet {v} hors de son plan")
    # combinatoire : arêtes, variété, orientabilité
    edge_dir = {}
    for f, cyc in enumerate(faces):
        for k in range(len(cyc)):
            a, b = cyc[k], cyc[(k + 1) % len(cyc)]
            edge_dir.setdefault(frozenset((a, b)), []).append((f, a, b))
    E = len(edge_dir)
    if any(len(l) != 2 for l in edge_dir.values()):
        err.append("arête non partagée par exactement deux faces")
    # orientabilité : on cherche un sens de parcours par face (±1) tel que chaque arête soit
    # parcourue dans des sens opposés par ses deux faces ; échec = surface non orientable
    sens = {0: 1}
    pile = [0]
    orientable = True
    while pile and orientable:
        f = pile.pop()
        for l in edge_dir.values():
            if len(l) != 2 or f not in (l[0][0], l[1][0]):
                continue
            (f1, a1, b1), (f2, a2, b2) = l if l[0][0] == f else (l[1], l[0])
            meme_sens = (a1, b1) == (a2, b2)
            voulu = -sens[f1] if meme_sens else sens[f1]
            if f2 not in sens:
                sens[f2] = voulu
                pile.append(f2)
            elif sens[f2] != voulu:
                orientable = False
    if not orientable:
        err.append("orientations incohérentes (ou surface non orientable)")
    for v in range(V):  # lien du sommet : les faces autour de v forment un seul cycle
        inc = [f for f, cyc in enumerate(faces) if v in cyc]
        adj = {f: set() for f in inc}
        for e, l in edge_dir.items():
            if v in e and len(l) == 2:
                adj[l[0][0]].add(l[1][0])
                adj[l[1][0]].add(l[0][0])
        seen, st = {inc[0]}, [inc[0]]
        while st:
            for g in adj[st.pop()]:
                if g not in seen:
                    seen.add(g)
                    st.append(g)
        if len(seen) != len(inc) or any(len(a) != 2 for a in adj.values()):
            err.append(f"sommet {v} : voisinage non disque")
    mult = {}
    for l in edge_dir.values():
        if len(l) == 2:
            k = tuple(sorted((l[0][0], l[1][0])))
            mult[k] = mult.get(k, 0) + 1
    missing = [p for p in combinations(range(F), 2) if p not in mult]
    if missing:
        err.append(f"{len(missing)} paires de faces sans arête commune")
    # triangulation exacte
    tris = []
    for f, cyc in enumerate(faces):
        ax = drop_axis(planes[f][:3])
        P2 = {v: (X[v][ax[0]], X[v][ax[1]]) for v in cyc}
        try:
            tris += [(f, t) for t in ear_clip(cyc, P2)]
        except AssertionError as e:
            err.append(f"face {f} : {e}")
    # tests triangle/triangle exacts (avec pré-filtre par boîtes englobantes exactes)
    boxes = [bbox([X[i] for i in t]) for _, t in tris]
    bad = 0
    for i, j in combinations(range(len(tris)), 2):
        (lo1, hi1), (lo2, hi2) = boxes[i], boxes[j]
        if any(hi1[k] < lo2[k] or hi2[k] < lo1[k] for k in range(3)):
            continue
        if not tri_pair_ok(tris[i][1], tris[j][1], X):
            bad += 1
            if bad <= 3:
                err.append(f"intersection parasite entre triangles des faces {tris[i][0]} et {tris[j][0]}")
    if bad > 3:
        err.append(f"... {bad} paires de triangles fautives au total")
    chi = V - E + F
    info = {"F": F, "V": V, "E": E, "chi": chi, "genre": (2 - chi) // 2 if orientable else None,
            "triangles": len(tris), "paires_par_multiplicite": {m: sum(1 for x in mult.values() if x == m) for m in sorted(set(mult.values()))}}
    return err, info


if __name__ == "__main__":
    code = 0
    for p in sys.argv[1:]:
        err, info = verify(p)
        print(p, info)
        print("  VALIDE" if not err else "  INVALIDE : " + " ; ".join(err[:6]))
        code |= bool(err)
    sys.exit(code)
