#pragma once

#include "const.h"
#include "eigen.h"

#include <vector>
#include <iostream>
#include <regex>
#include <memory>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/Point_3.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Unique_hash_map.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

const num EP = std::numeric_limits<num>::epsilon();
const num NEP = -1.0 * EP;
const num OPEP = 1.0 + EP;
typedef std::string string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

int factorial(int);
int kChoosei(int n, int k);
num bound(num x);
num bound(num x, num minVal, num maxVal);
num pow2(num x);
num scaleneArea(Vector3d p1, Vector3d p2, Vector3d p3);
template <typename DerivedA, typename DerivedB>
num distance(const Eigen::MatrixBase<DerivedA>& one, const Eigen::MatrixBase<DerivedB>& two) {
    return (one - two).norm(); 
}
num tetrahedronVolume(Vector3d, Vector3d, Vector3d, Vector3d);
Vector3d triangleNormal(MatrixXd &points, Vector3i &face);
using Vector3dIn = const Eigen::Ref<const Eigen::Vector3d>&;
bool barycentric(Vector3dIn p, Vector3dIn a, Vector3dIn b, Vector3dIn c, num &u, num &v, num &w, num &d, bool, num);
bool barycentric(Vector3dIn p, Vector3dIn a, Vector3dIn b, Vector3dIn c, num &u, num &v, num &w, num &d, bool);
bool barycentric(Vector3dIn p, Vector3dIn a, Vector3dIn b, Vector3dIn c, num &u, num &v, num &w, num &d, num);
bool barycentric(Vector3dIn p, Vector3dIn a, Vector3dIn b, Vector3dIn c, num &u, num &v, num &w);
