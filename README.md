# Two face-neighbourly polyhedra with 12 faces and genus 11

A polyhedron is *face-neighbourly* when every two of its faces share at least one edge
(« totocaèdres » in the Micmaths video that started this project). This repository contains
two such polyhedra found by computer search, the exact data and three independent verifiers.

Interactive 3D figures: https://totocaedres.fransys.io (French) · https://totocaedres.fransys.io/en

> **Status.** P1 and P2 pass three independent exact verifiers written for this project, but
> they have **not yet been checked by a third party**, and their novelty is **not confirmed**.
> The classical open problem (12 faces, genus 6, exactly one edge per pair of faces) is **not**
> solved here. Independent checks and pointers to prior work are very welcome: please open an issue.

*Résumé en français : deux polyèdres à 12 faces et 11 trous dont toutes les faces se touchent
deux à deux, avec leurs coordonnées exactes, et trois vérificateurs indépendants. Pas encore vérifiés par un tiers ; nouveauté non confirmée.*

## The two polyhedra

| | Faces | Vertices | Edges | Genus | Sides per face | Pairs sharing 1 / 2 / 3 / 4 edges |
|---|---|---|---|---|---|---|
| P1 | 12 | 64 | 96 | 11 | 4×14, 4×15, 4×19 | 44 / 16 / 4 / 2 |
| P2 | 12 | 64 | 96 | 11 | 4×13, 4×15, 4×20 | 48 / 12 / 0 / 6 |

Both have the symmetry group D2 (three half-turns about the coordinate axes) and are
"at least one edge per pair" neighbourly, like the genus-3 octahedra of Mizhaev and of
Röst–Vígh (2026). All vertices are trivalent and have rational coordinates.

Each polyhedron is bounded by 12 planes: the images, under the three half-turns, of three
generating planes `a x + b y + c z = d`.

| | Generating planes `[a, b, c, d]` |
|---|---|
| P1 | `[-4, 3, -4, 3]`, `[-3, 5, -4, -1]`, `[-4, 1, -3, 5]` |
| P2 | `[-5, 4, -6, 5]`, `[-5, 7, -5, -1]`, `[-6, 2, -4, 7]` |

Two observations from the search:

- the plane arrangement of P1 carries exactly two D2-symmetric face-neighbourly surfaces, P1
  and a polyhedron of the combinatorial type of P2 (and conversely); in a given arrangement,
  the P1-type solid is contained in the P2-type one, which has 8 more cells;
- in 175 nearby arrangement types that carry a solution, the solver proves that no
  face-neighbourly surface of lower genus exists in the same arrangement.

## Data format

`results/P1.json` and `results/P2.json`:

- `planes` — the 12 planes `[a, b, c, d]` (integers);
- `vertices` — 64 points with exact rational coordinates, as strings such as `"-7/3"`;
- `faces` — for each face, the cyclic list of its vertices (0-indexed);
- `tri` — for each vertex, the three planes meeting there;
- `gens` — the generating planes.

## Verify

```bash
sudo apt-get install libcgal-dev     # once, for the CGAL verifier
verification/run_all.sh
```

The script runs three independent exact checks on `results/*.json` and on control cases:

1. `verify_faces.py` — face by face, intersections along the lines common to two planes (rational arithmetic);
2. `verify_triangles.py` — triangulation, then exact triangle/triangle tests (rational arithmetic);
3. `cgal/check.cpp` — CGAL with the exact kernel (Epeck): triangulation and self-intersection test by CGAL.

Together they check that the surface is closed and orientable, that faces are planar simple
polygons, that the surface does not self-intersect, and that all 66 pairs of faces share an
edge. Expected output: P1, P2, Mizhaev, Röst–Vígh and Szilassi `VALIDE`; Razorcross and a
corrupted copy of P1 `INVALIDE`. The two Python scripts only use the standard library.

Control cases in `verification/controls/` and their sources:

- `szilassi.json` — Szilassi polyhedron, coordinates from Wikimedia Commons (CMG Lee);
- `mizhaev.json` — R. Mizhaev, arXiv:2609.17700;
- `rostvigh.json` — G. Röst and V. Vígh, arXiv:2609.32998;
- `neg_razorcross.json` — self-intersecting candidate from HackerPoet's NeighborlyPolyhedra;
- `neg_P1_sommets_permutes.json` — P1 with vertices permuted (must be rejected).

## Search method

The surface is sought as the boundary of a union of bounded cells of a generic arrangement of
planes, one face per plane. For a fixed arrangement this is a SAT problem (allowed vertex
patterns, edge alternation, one edge at least per pair of planes, and lazily added clauses
forcing each face to be a single simply connected polygon). Symmetry ties the cell variables
together. The answer depends only on the combinatorial type of the arrangement (signs of all
3×3 and 4×4 minors), so a GPU enumerates millions of integer arrangements exactly and keeps
one representative per type.

The search code is not published yet. The data and the verifiers in this repository are
sufficient to check both polyhedra independently of it.

Negative results so far (each one is a proof for its arrangement type only): no solution in
about 660,000 arrangement types of three other symmetric families with 10 or 12 faces, and no
genus-6 solution was found.

## Acknowledgements

The search and verification code was written with the help of Claude (Anthropic). The project
was prompted by Mickaël Launay's Micmaths video on « totocaèdres ».

## Licence

Code: MIT (`LICENSE`). Data in `results/`: CC BY 4.0 (`LICENSE-DATA.md`).
