// Troisième vérification, par CGAL en arithmétique exacte (noyau Epeck, rationnels GMP).
// Entrée (générée par cgal_input.py) :
//   V F
//   V lignes : x y z  (rationnels exacts, forme "p/q")
//   F lignes : n i1 ... in  (cycle de sommets, indices à partir de 0)
// Contrôles réalisés par CGAL :
//   - construction d'un Surface_mesh valide (is_valid_polygon_mesh), fermé (is_closed) ;
//   - planéité exacte de chaque face ;
//   - triangulation des faces par CGAL (PMP::triangulate_faces, CDT dans le plan de la face) ;
//   - absence d'auto-intersection (PMP::does_self_intersect, prédicats et constructions exacts) ;
//   - V - E + F, orientabilité implicite (Surface_mesh orienté et fermé), nombre de composantes ;
//   - adjacence : pour chaque paire de faces d'origine, nombre d'arêtes communes.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/helpers.h>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv) {
  std::ifstream in(argv[1]);
  int V, F;
  in >> V >> F;
  Mesh m;
  std::vector<Mesh::Vertex_index> vid;
  for (int i = 0; i < V; i++) {
    std::string a, b, c;
    in >> a >> b >> c;
    typedef K::FT::ET ET;  // rationnel GMP exact
    K::FT x{ET(a)}, y{ET(b)}, z{ET(c)};
    vid.push_back(m.add_vertex(Point(x, y, z)));
  }
  std::vector<std::vector<int>> faces(F);
  auto label = m.add_property_map<Mesh::Face_index, int>("f:orig", -1).first;
  bool ok = true;
  for (int f = 0; f < F; f++) {
    int n;
    in >> n;
    faces[f].resize(n);
    std::vector<Mesh::Vertex_index> cyc;
    for (int k = 0; k < n; k++) in >> faces[f][k], cyc.push_back(vid[faces[f][k]]);
    Mesh::Face_index fi = m.add_face(cyc);
    if (fi == Mesh::null_face()) {  // orientation incompatible avec les faces déjà posées : on inverse
      std::reverse(cyc.begin(), cyc.end());
      fi = m.add_face(cyc);
    }
    if (fi == Mesh::null_face()) {
      std::cout << "ECHEC : face " << f << " impossible à insérer (non-variété ou non orientable)\nINVALIDE\n";
      return 1;
    }
    label[fi] = f;
    // planéité exacte
    const Point &p0 = m.point(cyc[0]);
    for (int k = 1; k + 1 < n; k++) {
      if (!CGAL::collinear(p0, m.point(cyc[k]), m.point(cyc[k + 1]))) {
        K::Plane_3 pl(p0, m.point(cyc[k]), m.point(cyc[k + 1]));
        for (auto v : cyc)
          if (!pl.has_on(m.point(v))) ok = false, std::cout << "ECHEC : face " << f << " non plane\n";
        break;
      }
    }
  }
  bool valid = CGAL::is_valid_polygon_mesh(m), closed = CGAL::is_closed(m);
  std::size_t nv = m.number_of_vertices(), ne = m.number_of_edges(), nf = m.number_of_faces();
  long chi = (long)nv - (long)ne + (long)nf;
  // adjacence des faces d'origine (avant triangulation)
  std::map<std::pair<int, int>, int> mult;
  for (auto e : m.edges()) {
    Mesh::Face_index f0 = m.face(m.halfedge(e, 0)), f1 = m.face(m.halfedge(e, 1));
    if (f0 == Mesh::null_face() || f1 == Mesh::null_face()) continue;  // arête de bord : surface non fermée
    int a = label[f0], b = label[f1];
    if (a > b) std::swap(a, b);
    mult[{a, b}]++;
  }
  int pairs = 0;
  std::map<int, int> hist;
  for (int a = 0; a < F; a++)
    for (int b = a + 1; b < F; b++)
      if (mult.count({a, b})) pairs++, hist[mult[{a, b}]]++;
  auto fcc = m.add_property_map<Mesh::Face_index, std::size_t>("f:cc").first;
  std::size_t ncc = PMP::connected_components(m, fcc);
  if (!valid || !closed) {
    std::cout << "sommets " << nv << ", aretes " << ne << ", faces " << nf << " : maillage " << (valid ? "" : "invalide ")
              << (closed ? "" : "non fermé") << "\nINVALIDE\n";
    return 1;
  }
  // triangulation par CGAL puis test exact d'auto-intersection
  bool tri_ok = PMP::triangulate_faces(m);
  std::vector<std::pair<Mesh::Face_index, Mesh::Face_index>> bad;
  PMP::self_intersections(m, std::back_inserter(bad));
  std::cout << "sommets " << nv << ", aretes " << ne << ", faces " << nf << ", V-E+F = " << chi
            << ", genre " << (2 - chi) / 2 << ", composantes " << ncc << "\n";
  std::cout << "maillage valide " << valid << ", ferme " << closed << ", triangulation " << tri_ok
            << " (" << m.number_of_faces() << " triangles)\n";
  std::cout << "paires de faces adjacentes " << pairs << "/" << F * (F - 1) / 2 << ", multiplicites :";
  for (auto& [k, v] : hist) std::cout << " " << v << "x" << k;
  std::cout << "\npaires de triangles qui s'intersectent : " << bad.size() << "\n";
  ok = ok && valid && closed && tri_ok && ncc == 1 && pairs == F * (F - 1) / 2 && bad.empty();
  std::cout << (ok ? "VALIDE" : "INVALIDE") << "\n";
  return ok ? 0 : 1;
}
