#include "mesh.h"

Face::Face() {}

double error;
vector<double> errors;

bool Mesh::read_vertices(const PolyMesh& mesh, Vd_array& fixed_vertices) {
    // The selection file is a list of integers, so we must build a correspondence
    // between vertices and the integers.
    std::vector<vertex_descriptor> vds;
    vds.reserve(CGAL::num_vertices(mesh));
    vertex_iterator vi = CGAL::vertices(mesh).begin(), vi_end = CGAL::vertices(mesh).end();
    CGAL_For_all(vi, vi_end) {
        vds.push_back(*vi);
    }

    // Get the first line and read the fixed vertex indices
    //std::size_t counter = 0;
    int counter = 0;
    std::size_t s;
    std::unordered_set<std::size_t> indices;
    for (auto i : paramCorners) {
        s = static_cast<size_t>(i);
        if(s >= vds.size())
        {
            cerr << "Error: Vertex index too large" << std::endl;
            return false;
        }

        vertex_descriptor vd = vds[s];
        if(!CGAL::is_border(vd, mesh)) { // must be on the border
            cerr << "Error: vertex is not on the border of the mesh" << std::endl;
            return false;
        }

        if(counter >= 4) { // too many border vertices
            cerr << "Error: Too many vertices are fixed" << std::endl;
            return false;
        }

        fixed_vertices[counter++] = vd;
        indices.insert(s);
    }

    if(indices.size() < 4) {
        cerr << "Error: at least four unique vertices must be provided" << std::endl;
        return false;
    }

    return true;
}

Mesh::Mesh(string inputFile, string inputType) : Mesh(inputFile, inputType,
        Vector4i{-1, -1, -1, -1}, 1.0) {}

Mesh::Mesh(string inputFile, string inputType, double scale) : Mesh(inputFile,
        inputType, Vector4i{-1, -1, -1, -1}, scale) {}

Mesh::Mesh(string inputFile, string inputType, const Vector4i &paramCorners) :
    Mesh(inputFile, inputType, paramCorners, 1.0) {}

Mesh::Mesh(string inputFile, string inputType, const Vector4i &paramCorners, double scale) {
    this->paramCorners = paramCorners;
    this->inputType = inputType;
    this->inputFile = inputFile;
    if (!inputType.compare("raytrace")) {
        rtsettings = std::make_unique<RTSettings>(rayBox);
        rtsettings->inputFileName = inputFile;
    }
    this->scale = scale;
    prevTriF = 0;
    prevTriR = 0;
    triRows = 20;
    triCols = 20;
    resetPrevTri = true;
}

void Mesh::initializeMesh() {
    if (!inputType.compare("bezier")) {
        preprocessBezier();
        modelFile = "temp.obj";
    } else if (!inputType.compare("raytrace")) {
        preprocessRayTrace();
        modelFile = "temp.obj";
    } else {
        modelFile = inputFile;
    }
    // CGAL
    std::ifstream in(modelFile);
    std::vector<Point_3> pointsCGAL;
    std::vector<std::vector<std::size_t>> polygons;
    CGAL::IO::read_OBJ(in, pointsCGAL, polygons);

    for (auto& p : pointsCGAL) {
        p = Point_3(p.x() * scale, p.y() * scale, p.z() * scale);
    }

    std::stringstream ss;
    CGAL::IO::write_OFF(ss, pointsCGAL, polygons);
    //ss.seekg(0, std::ios::beg);
    PolyMesh sm;
    ss >> sm;

    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

    UV_uhm uv_uhm;
    UV_pmap uv_map(uv_uhm);

    typedef SMP::Square_border_arc_length_parameterizer_3<PolyMesh> Border_parameterizer;
    typedef SMP::Mean_value_coordinates_parameterizer_3<PolyMesh, Border_parameterizer> Parameterizer;

    Border_parameterizer *border_param;

    if (paramCorners[0] > -1) {
        Vd_array vda;
        if(!read_vertices(sm, vda)) {
            cerr << "Error: problem loading the square corners" << std::endl;
            exit(1);
        }
        border_param = new Border_parameterizer(vda[0], vda[1], vda[2], vda[3]);
    } else {
        border_param = new Border_parameterizer(); // the border parameterizer will compute the corner vertices
    }

    SMP::Error_code err = SMP::parameterize(sm, Parameterizer(*border_param), bhd, uv_map);
    delete border_param;
    boost::graph_traits<PolyMesh>::vertex_iterator vit, vend;
    boost::tie(vit, vend) = vertices(sm);
    pointMap = MatrixXd::Zero(vertices(sm).size(), 5);
    int g = 0;
    while(vit!=vend) {
        vertex_descriptor vd = *vit++;
        auto par = get(uv_map, vd);
        pointMap.row(g) << par.x(), par.y(), vd->point().x(), vd->point().y(), vd->point().z();
        g++;
    }

    if(err != SMP::OK) {
        cerr << "Error: " << SMP::get_error_message(err) << std::endl;
        exit(1);
    }

    // Generate and read in the results from .off format
    // Not the cleanest way to handle the CGAL data but it is fairly straightforward
    std::stringstream ss1;
    SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, ss1);

    int count = 0, vertCount = 0, faceCount = 0, one, two, three;
    num x, y, z;
    for (std::string line; std::getline(ss1, line); ) {
        std::stringstream lineStream(line);
        std::string temp;
        if (!count) {
        } else if (count == 1) {
            std::getline(lineStream, temp, ' ');
            vertCount = stoi(temp);
            std::getline(lineStream, temp, ' ');
            faceCount = stoi(temp);
            paramPoints = MatrixXd::Zero(vertCount, 3);
            for (int i=0; i<faceCount; i++) {
                auto face = new Face();
                face->index = i;
                faces.push_back(face);
            }
        } else if (count < vertCount + 2) {
            std::getline(lineStream, temp, ' ');
            x = stod(temp);
            std::getline(lineStream, temp, ' ');
            y = stod(temp);
            std::getline(lineStream, temp, ' ');
            z = stod(temp);
            paramPoints.row(count - 2) << x, y, z;
        } else {
            std::getline(lineStream, temp, ' ');
            std::getline(lineStream, temp, ' ');
            one = stoi(temp);
            std::getline(lineStream, temp, ' ');
            two = stoi(temp);
            std::getline(lineStream, temp, ' ');
            three = stoi(temp);
            faces[count - vertCount - 2]->vertices << one, two, three;
        }
        count++;
    }
    realPoints = MatrixXd::Zero(vertCount, 3);
    for (int s=0; s<paramPoints.rows(); s++) {
        auto row = paramPoints.row(s);
        for (auto mapRow : pointMap.rowwise()) {
            if (abs(mapRow(0) - row(0)) < 1e-3 && abs(mapRow(1) - row(1)) < 1e-3) {
                realPoints.row(s) = mapRow.tail(3);
                break;
            }
        }
    }
    buildFaceMap();
    /*surfaceArea = 0.0;
    for (auto f : faces) {
        Vector3d p1 = realPoints.row(f->vertices(0));
        Vector3d p2 = realPoints.row(f->vertices(1));
        Vector3d p3 = realPoints.row(f->vertices(2));
        surfaceArea += scaleneArea(p1, p2, p3);
    }*/
}

Vector3d Mesh::barycentricLoopingF(const Vector2d &point) {
    return barycentricLoopingF(point, prevTriF);
}

MatrixXd Mesh::barycentricLoopingF(const MatrixXd &points) {
    return barycentricLoopingOCL(points, false);
}

MatrixXd Mesh::barycentricLoopingR(const MatrixXd &points) {
    return barycentricLoopingOCL(points, true);
}

Vector2d Mesh::barycentricLoopingR(const Vector3d &point) {
    return barycentricLoopingR(point, prevTriR);
}

Vector2d Mesh::barycentricLoopingR(const Vector3d &point, Vector3d &normal) {
    return barycentricLoopingR(point, prevTriR, normal);
}

Vector3d Mesh::barycentricLoopingF(const Vector2d &point, int &prevTri) {
    Vector3d p;
    p << point, 0.0;
    return barycentricLoopingOCL(p, false);
}

Vector2d Mesh::barycentricLoopingR(const Vector3d &point, int &prevTri) {
    Vector3d normal;
    return barycentricLoopingR(point, prevTri, normal);
}

Vector2d Mesh::barycentricLoopingR(const Vector3d &point, int &prevTri, Vector3d &normal) {
    Vector2d v;
    Vector3d v3;
    v3 << barycentricLoopingOCL(point, normal, true);
    v << v3(0), v3(1);
    return v;
}

std::chrono::duration<double> reverseTime;
std::chrono::duration<double> forwardTime;

Vector3d Mesh::barycentricLoopingOCL(const Vector3d &point, bool reverse) {
    Vector3d normal;
    return barycentricLoopingOCL(point, normal, reverse);
}

Vector3d Mesh::barycentricLoopingOCL(const Vector3d &point, Vector3d &normal, bool reverse) {
    MatrixXd points = MatrixXd::Zero(1,3);
    MatrixXd normals;
    points.row(0) << point;
    MatrixXd outputMat = barycentricLoopingOCL(points, normals, reverse);
    Vector3d output = outputMat.row(0);
    normal << normals.row(0);
    return output;
}

MatrixXd Mesh::barycentricLoopingOCL(const MatrixXd &points, bool reverse) {
    MatrixXd normals;
    return barycentricLoopingOCL(points, normals, reverse);
}

void Mesh::buildFaceMap() {
    Array3i fv, ov;
    int c;
    for (auto f : faces) {
        fv = f->vertices;
        for (auto of : faces) {
            ov = of->vertices;
            if (f != of) {
                c = 0;
                for (int i=0; i<3; i++)
                    c += (fv(i) == ov).any();
                if (c == 2)
                    f->neighbors.insert(of);
                else if (c)
                    f->pointNeighbors.insert(of);
            }
        }
        f->isEdge = f->neighbors.size() != 3;
    }

    // Silence the OpenCLWrapper call
    std::ofstream null_file("/dev/null");
    std::streambuf* original_cout = std::cout.rdbuf();
    std::cout.rdbuf(null_file.rdbuf());
    // compile OpenCL C code for the fastest available device
    device = new Device(select_device_with_most_flops());
    // Restore cout
    std::cout.rdbuf(original_cout);

    face_N = static_cast<uint>(faces.size());
    model_point_N = static_cast<uint>(realPoints.rows());
    FN = new Memory<uint>(*device, 2);
    NP = new Memory<double>(*device, MAX_OCL_PT_COUNT*face_N*3);
    N = new Memory<double>(*device, MAX_OCL_PT_COUNT*face_N*3);
    V = new Memory<double>(*device, MAX_OCL_PT_COUNT*face_N);
    FC = new Memory<double>(*device, MAX_OCL_PT_COUNT*face_N);
    H = new Memory<int>(*device, MAX_OCL_PT_COUNT*face_N);
    
    PTF = new Memory<double>(*device, model_point_N*3);
    PTR = new Memory<double>(*device, model_point_N*3);
    F = new Memory<int>(*device, face_N*3);

    (*FN)[0] = face_N;
    uint i, j, i3, i3j;
    Vector3i inds;
    Vector3d rowF, rowR;
    for (i=0; i<face_N; i++) {
        inds = faces[i]->vertices;
        for (j=0; j<3; j++) {
            (*F)[i*3+j] = inds(j) * 3;
        }
    }
    for (i=0; i<model_point_N; i++) {
        i3 = i*3;
        rowF = realPoints.row(i);
        rowR = paramPoints.row(i);
        for (j=0; j<3; j++) {
            i3j = i3 + j;
            (*PTF)[i3j] = rowF(j);
            (*PTR)[i3j] = rowR(j);
        }
    }

    FN->write_to_device();
    PTF->write_to_device();
    PTR->write_to_device();
    F->write_to_device();
}

int counterB = 0;

void Mesh::printModelAndPoint(Vector3d point) {
    MatrixXd m(1, 3);
    m.row(0) << point;
    printModelAndPoint(m);
}

void Mesh::printModelAndPoint(MatrixXd &points) {
    std::ofstream outMP("out_model_point.txt");
    outMP << points << endl;
    // NOTE: This creates duplicates, but it doesn't really matter
    for (auto f : faces) {
        Array3i fv = f->vertices;
        Vector3d p0 = realPoints.row(fv(0));
        Vector3d p1 = realPoints.row(fv(1));
        Vector3d p2 = realPoints.row(fv(2));
        outMP << p0 << " " << p1 << endl;
        outMP << p0 << " " << p2 << endl;
        outMP << p2 << " " << p1 << endl;
    }
}

MatrixXd Mesh::barycentricLoopingOCL(const MatrixXd &points, MatrixXd &normals, bool reverse) {
    counterB++;
    if (points.rows() > MAX_OCL_PT_COUNT) {
        cout << "ERROR: Cannot pass more than " << MAX_OCL_PT_COUNT << " (MAX_OCL_PT_COUNT) points to bLOCL\n";
        exit(1);
    }
    normals = MatrixXd::Zero(points.rows(), 3);
    Memory<double> *PT, *PF;
    if (reverse) {
        PT = PTR;
        PF = PTF;
    } else {
        PT = PTF;
        PF = PTR;
    }
    uint i, j;
    const uint point_N = static_cast<uint>(points.rows());
    const uint point_N3 = point_N * 3;
    const uint ocl_N = point_N * face_N;

    Kernel *barycentric_loop, *find_best_points = nullptr;
    (*FN)[1] = point_N;
    FN->write_to_device();
    Memory<double> P(*device, point_N3);
    Memory<double> PO(*device, point_N3);
    Memory<double> NO, VO, HO;
    if (reverse) {
        NO = Memory<double>(*device, point_N3);
        VO = Memory<double>(*device, point_N);
    } else {
        HO = Memory<double>(*device, point_N);
    }
    if (reverse) {
        barycentric_loop = new Kernel(*device, ocl_N, "barycentric_loop_R", *H,
                *NP, *V, *N, P, *PF, *F, *FN, *FC);
        find_best_points = new Kernel(*device, point_N, "find_best_points", *H,
                *PT, *F, *NP, *V, *N, PO, VO, NO, *FN, *FC);
    } else {
        barycentric_loop = new Kernel(*device, ocl_N, "barycentric_loop_F", HO,
                PO, P, *PF, *F, *PT, *FN);
    }

    for (i=0; i<point_N; i++) {
        Vector3d point;
        if (points.cols() == 2)
            point << points.row(i), 0;
        else
            point << points.row(i);
        for (j=0; j<3; j++) {
            P[i*3+j] = point[j];
        }
    }

    P.write_to_device();
    barycentric_loop->run(); // run kernel on the device
    if (reverse)
        find_best_points->run();
    delete barycentric_loop;
    if (reverse) {
        delete find_best_points;
        VO.read_from_device();
        NO.read_from_device();
    } else {
        HO.read_from_device();
    }
    PO.read_from_device();
    MatrixXd output = MatrixXd::Zero(point_N, 3);
    bool quit = false;

    std::vector<Vector3d> failedPointsVec;
    for (i=0; i<point_N; i++) {
        if ((reverse && VO[i] >= 0.0) || (!reverse && HO[i] != -1)) {
            int i3 = i * 3;
            output.row(i) << PO[i3], PO[i3+1], PO[i3+2];
            if (reverse)
                normals.row(i) << NO[i3], NO[i3+1], NO[i3+2];
        } else {
            quit = true;
            failedPointsVec.push_back(points.row(i));
            output.row(i) << 0, 0, 0;
            normals.row(i) << 1, 0, 0;
        }
    }
    if (quit) {
        MatrixXd failedPoints(failedPointsVec.size(), 3);
        for (i=0; i<failedPointsVec.size(); i++)
            failedPoints.row(i) = failedPointsVec[i];
        printModelAndPoint(failedPoints);
        exit(0);
    }
    return output;
}

Vector3d Mesh::barycentricLooping(const Vector3d &point, Vector3d &normal, bool reverse) {
    MatrixXd pointsFrom, pointsTo;
    int *prevTri;
    if (reverse) {
        pointsFrom = realPoints;
        pointsTo = paramPoints;
        prevTri = &prevTriR;
    } else {
        pointsFrom = paramPoints;
        pointsTo = realPoints;
        prevTri = &prevTriF;
    }
    Vector3d pointOut;
    vector<Face*> orFaces;
    auto currTri = faces[*prevTri];
    bool found = false;
    num u, v, w;
    if (resetPrevTri) {
        orFaces = faces;
    } else {
        orFaces.push_back(currTri);
        for (auto f : currTri->neighbors) {
            orFaces.push_back(f);
        }
        for (auto f : currTri->pointNeighbors) {
            orFaces.push_back(f);
        }
        for (auto f : faces) {
            if (f != currTri && !currTri->neighbors.count(f) && !currTri->pointNeighbors.count(f)) {
                orFaces.push_back(f);
            }
        }
    }
    unsigned int i=0;
    for (; !found && i<orFaces.size(); i++) {
        auto face = orFaces[i];
        Vector3i inds = face->vertices;
        Vector3d a = pointsFrom.row(inds(0));
        Vector3d b = pointsFrom.row(inds(1));
        Vector3d c = pointsFrom.row(inds(2));
        if (scaleneArea(a, b, c) > 0.0) {
            auto check = barycentric(point, a, b, c, u, v, w);//, dd, false, .001);
            if (check) {
                found = true;
                resetPrevTri = false;
                *prevTri = face->index;
                MatrixXd p(3,3);
                p.row(0) = pointsTo.row(inds(0));
                p.row(1) = pointsTo.row(inds(1));
                p.row(2) = pointsTo.row(inds(2));
                pointOut = p.row(0) * u +
                            p.row(1) * v +
                            p.row(2) * w;
            }
        }
    }
    if (!found) {
        cerr << "Error: Point " << point << fmt::format(" not found in {} space\n", reverse ? "parameter" : "model");
        abort();
        exit(1);
    }
    return pointOut;
}

// TODO: rows and cols need to be part of the config files
void Mesh::preprocessBezier() {
    MatrixXd controlPoints = MatrixXd::Zero(0, 3),
             points(triRows*triCols, 3),
             normals(triRows*triCols, 3);
    vector<Vector3i> faces;
    initialize(&controlPoints, inputFile);
    const double vStep = 1.0 / (triRows - 1);
    const double uStep = 1.0 / (triCols - 1);
    num u, v;
    for (int i = 0; i < triCols; i++) {
        for (int j = 0; j < triRows; j++) {
            u = uStep * i;
            v = vStep * j;
            points.row(i*triRows+j) = p(&controlPoints, u, v);
            normals.row(i*triRows+j) = normal(&controlPoints, u, v);
            if (i < triCols - 1 && j < triRows - 1) {
                Vector3i f = {i * triRows + j,
                              (i + 1) * triRows + j,
                              i * triRows + j + 1};
                faces.push_back(f);
                Vector3i f2 = {i * triRows + j + 1,
                              (i + 1) * triRows + j,
                              (i + 1) * triRows + j + 1};
                faces.push_back(f2);
            }
        }
    }
    std::ofstream out("temp.obj");
    for (Vector3d row : points.rowwise())
        out << fmt::format("v {} {} {}\n", row(0), row(1), row(2));
    for (Vector3d row : normals.rowwise())
        out << fmt::format("vn {} {} {}\n", row(0), row(1), row(2));
    for (auto r : faces)
        out << fmt::format("f {}//{} {}//{} {}//{}\n", r(0)+1, r(0)+1, r(1)+1, r(1)+1, r(2)+1, r(2)+1);
    out.close();
    paramCorners << 0, triRows-1, triCols*(triRows-1), triCols*triRows-1;
}

// TODO: rows and cols need to be part of the config files
void Mesh::preprocessRayTrace() {
    rtsettings->x = triRows;
    rtsettings->y = triCols;
    MatrixXd points, normals(triRows*triCols, 3);
    vector<Vector3i> faces;
    vector<Vector3d> triNorms;
    MatrixXi facesDummy;
    /*num uStart = setts["model"]["ray_start"][0].value_or(0.0);
    num vStart = setts["model"]["ray_start"][1].value_or(0.0);
    num uLength = setts["model"]["ray_box_size"][0].value_or(0.0);
    num vLength = setts["model"]["ray_box_size"][1].value_or(0.0);

    double rayBox[4] = {uStart, vStart, uLength, vLength};*/
    // Get ray tracing data with default u x v sizing for analysis
    //modelFile, points, facesDummy, rayBox, triRows, triCols
    rtsettings->updateRayBox(rayBox);
    auto rayTraceResults = rayTrace(rtsettings.get());
    points = std::get<2>(rayTraceResults);
    for (int i = 0; i < triCols; i++) {
        for (int j = 0; j < triRows; j++) {
            //normals.row(i*triRows+j) = normal(&controlPoints, u, v);
            if (i < triCols - 1 && j < triRows - 1) {
                Vector3i f = {i * triRows + j,
                              (i + 1) * triRows + j,
                              i * triRows + j + 1};
                faces.push_back(f);
                Vector3i f2 = {i * triRows + j + 1,
                              (i + 1) * triRows + j,
                              (i + 1) * triRows + j + 1};
                faces.push_back(f2);
            }
        }
    }

    std::ofstream out("temp.obj");
    for (Vector3d row : points.rowwise())
        out << fmt::format("v {} {} {}\n", row(0), row(1), row(2));
    for (auto r : faces) {
        out << fmt::format("f {}//{} {}//{} {}//{}\n", r(0)+1, r(0)+1, r(1)+1, r(1)+1, r(2)+1, r(2)+1);
        triNorms.push_back(triangleNormal(points, r));
    }
    for (long unsigned int i=0; i<triNorms.size(); i++) {
        Vector3i f = faces[i];
        Vector3d n = triNorms[i];
        for (auto v : f) {
            normals.row(v) += n;
        }
    }
    for (int i=0; i<normals.rows(); i++) {
        Vector3d row = normals.row(i).normalized();
        out << fmt::format("vn {} {} {}\n", row(0), row(1), row(2));
    }
    out.close();
    paramCorners << 0, triRows-1, triCols*(triRows-1), triCols*triRows-1;
}
