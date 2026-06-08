#pragma once

#include "util.h"
#include "opencl.hpp"
#include "raytrace.h"
#include "bezier.h"

typedef CGAL::Simple_cartesian<double>                           PCKernel;
typedef PCKernel::Point_2                                        Point_2;
typedef PCKernel::Point_3                                        Point_3;
typedef CGAL::Polyhedron_3<PCKernel>                             PolyMesh;
typedef boost::graph_traits<PolyMesh>::halfedge_descriptor       halfedge_descriptor;
typedef boost::graph_traits<PolyMesh>::vertex_descriptor         vertex_descriptor;
typedef boost::graph_traits<PolyMesh>::face_descriptor           face_descriptor;
typedef boost::graph_traits<PolyMesh>::vertex_iterator           vertex_iterator;
typedef std::array<vertex_descriptor, 4>                       Vd_array;
typedef CGAL::Unique_hash_map<vertex_descriptor, Point_2>        UV_uhm;
typedef boost::associative_property_map<UV_uhm>                  UV_pmap;
namespace SMP = CGAL::Surface_mesh_parameterization;

class Face {
    public:
        Face();
        ~Face() = default;
        Array3i vertices;
        std::set<Face*> neighbors;
        std::set<Face*> pointNeighbors;
        bool isEdge;
        int index;
};

class Mesh {
    public:
        Mesh(string inputFile, string inputType, const Vector4i &paramCorners, double scale);
        Mesh(string inputFile, string inputType, double scale);
        Mesh(string inputFile, string inputType, const Vector4i &paramCorners);
        Mesh(string inputFile, string inputType);
        ~Mesh() = default;
        void initializeMesh();
        Vector3d barycentricLoopingF(const Vector2d &point);
        MatrixXd barycentricLoopingF(const MatrixXd &points);
        Vector3d barycentricLoopingF(const Vector2d &point, int &prevTri);
        Vector2d barycentricLoopingR(const Vector3d &point);
        MatrixXd barycentricLoopingR(const MatrixXd &points);
        Vector2d barycentricLoopingR(const Vector3d &point, Vector3d &normal);
        Vector2d barycentricLoopingR(const Vector3d &point, int &prevTri);
        Vector2d barycentricLoopingR(const Vector3d &point, int &prevTri, Vector3d &normal);

        std::unique_ptr<RTSettings> rtsettings = nullptr;
        vector<Face*> faces;
        MatrixXd paramPoints;
        MatrixXd realPoints;
        int prevTriF, prevTriR, triRows, triCols;
        double scale;
        bool resetPrevTri;
        Vector4i paramCorners;
        Vector4d rayBox;
        MatrixXd pointMap;
    
    private:
        void printModelAndPoint(Vector3d point);
        void printModelAndPoint(MatrixXd &points);
        void buildFaceMap();
        void preprocessBezier();
        void preprocessRayTrace();
        MatrixXd barycentricLoopingOCL(const MatrixXd &points, MatrixXd &normals, bool reverse);
        MatrixXd barycentricLoopingOCL(const MatrixXd &points, bool reverse);
        Vector3d barycentricLoopingOCL(const Vector3d &point, Vector3d &normal, bool reverse);
        Vector3d barycentricLoopingOCL(const Vector3d &point, bool reverse);
        Vector3d barycentricLooping(const Vector3d &point, Vector3d &normal, bool reverse);
        bool read_vertices(const PolyMesh& mesh, Vd_array& fixed_vertices);
        Memory<double> *AF, *BF, *CF, *AR, *BR, *CR, *NP, *N, *V, *PTF, *PTR, *FC;
        Memory<uint> *FN;
        Memory<int> *H, *F;
        Device *device;
        uint face_N, model_point_N;
        string modelFile, inputType, inputFile;
};
