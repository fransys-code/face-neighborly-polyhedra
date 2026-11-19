#!/bin/bash
# Lance les trois vérifications indépendantes sur les résultats et sur les cas de contrôle.
#   1. verify_faces.py      : face par face, droites d'intersection, fractions exactes
#   2. verify_triangles.py  : triangulation + tests triangle/triangle, fractions exactes
#   3. cgal/check.cpp       : CGAL (noyau exact Epeck), triangulation et auto-intersections par CGAL
# Prérequis pour CGAL : sudo apt-get install libcgal-dev (tire GMP, MPFR et Boost).
# Usage : verification/run_all.sh            (code de sortie 0 si tout est conforme)
set -u
V="$(dirname "$(readlink -f "$0")")"
ROOT="$V/.."
BUILD="$V/cgal/build"; mkdir -p "$BUILD"
if [ ! -x "$BUILD/check" ] || [ "$V/cgal/check.cpp" -nt "$BUILD/check" ]; then
  g++ -O2 -std=c++17 -o "$BUILD/check" "$V/cgal/check.cpp" -lgmp -lmpfr || { echo "compilation CGAL impossible"; exit 2; }
fi
echec=0
verifie() {  # fichier, attendu (VALIDE|INVALIDE)
  local f="$1" attendu="$2" r1 r2 r3
  r1=$(python3 -I "$V/verify_faces.py" "$f" 2>&1 | grep -o "^VALIDE\|^INVALIDE" | tail -1)
  r2=$(python3 -I "$V/verify_triangles.py" "$f" 2>&1 | grep -o "VALIDE\|INVALIDE" | tail -1)
  python3 -I "$V/cgal/to_input.py" "$f" > "$BUILD/in.txt"
  r3=$("$BUILD/check" "$BUILD/in.txt" | tail -1)
  local ok="ok"
  for r in "$r1" "$r2" "$r3"; do [ "$r" = "$attendu" ] || ok="ÉCART"; done
  [ "$ok" = "ok" ] || echec=1
  printf '%-34s attendu %-8s | faces %-8s | triangles %-8s | CGAL %-8s | %s\n' "$(basename "$f")" "$attendu" "$r1" "$r2" "$r3" "$ok"
}
for f in "$ROOT"/results/*.json; do verifie "$f" VALIDE; done
for f in "$V"/controls/{mizhaev,rostvigh,szilassi}.json; do verifie "$f" VALIDE; done
for f in "$V"/controls/neg_*.json; do verifie "$f" INVALIDE; done
exit $echec
