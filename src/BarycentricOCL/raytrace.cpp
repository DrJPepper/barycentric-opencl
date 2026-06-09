#include "./raytrace.h"

std::tuple<int, int, MatrixXd&> rayTrace(RTSettings *settings) {
    //string inputFileName, MatrixXd &vertices, MatrixXi &faces, Vector4d rayBox, int xPoints, int yPoints

    //if (xPoints && xPoints) {
        //settings->x = xPoints;
        //settings->y = yPoints;
    //}

    // Run primary functions
    MatrixXd pixelMap = MatrixXd::Zero(settings->x * settings->y, 3);
    std::ifstream currMesh(settings->inputFileName);
    if (!currMesh.good()) {
        cout << "ERROR: input mesh file " << settings->inputFileName
             << " invalid or doesn't exist\n";
        exit(1);
    }
    Intersectable *temp = new SMFModel(settings, &currMesh);
    ACGBVH *newObj = new ACGBVH(temp, settings->maxObjCount);
    objects.push_back(newObj);
    intersectForAllPixels(objects, settings, &pixelMap);
    settings->vertices = reinterpret_cast<SMFModel*>(objects[0]->obj)->vertices;
    settings->faces = reinterpret_cast<SMFModel*>(objects[0]->obj)->facesStored;
    MatrixXd *intPoints = new MatrixXd;
    (*intPoints) = tBuffer.block(0, 1, tBuffer.rows(), 3);
    return {settings->x, settings->y, std::ref(*intPoints)};
}

// Loops over all pixels
void intersectForAllPixels(vector<ACGBVH *> objects, RTSettings *settings,
                           MatrixXd *result) {
    Vector3d djk = settings->Zv;
    int j, k, index;
    int size = static_cast<int>(objects.size());
    double uStep, vStep;
    if (!settings->uLength || !settings->vLength) {
        uStep = objects[0]->obj->maxVec(1) / settings->x;
        vStep = objects[0]->obj->maxVec(2) / settings->y;
    } else {
        uStep = settings->uLength / settings->x;
        vStep = settings->vLength / settings->y;
    }

    // t, x, y, z (of intersection), x, y, z (of normal)
    tBuffer = MatrixXd::Zero(result->rows(), 7);
    // Stores the front object for each pixel since tBuffer is a double matrix
    // and it seems like a bad idea to store indices as doubles
    tObjsBuffer = VectorXi::Zero(result->rows());
    for (int i = 0; i < tBuffer.rows(); i++) tObjsBuffer(i) = -1;
    for (int i = 0; i < size; i++) {
        for (j = 0; j < settings->x; j++) {
            for (k = 0; k < settings->y; k++) {
                index = k * settings->x + j;
                settings->cameraLoc << -DBL_MAX / 500.0, settings->uStart + uStep * j, settings->vStart + vStep * k;
                objects[i]->intersectPixel(djk.normalized(), settings, index, i);
            }
        }
    }
}

// Parse the input file and generate resulting settings and objects
RTSettings::RTSettings(Vector4d &rayBox) {
    // Default values
    maxObjCount = 200;
    d = 3.0;
    x = 20;
    y = 20;
    theta = 0;
    scale = 1.0;
    cameraLoc << 3.0, 0.0, 0.0;
    Zv << -1.0, 0.0, 0.0;
    Vup << 0.0, 1.0, 0.0;
    rotate << 0.0, 0.0, 0.0;

    updateRayBox(rayBox);

    Xv = Zv.cross(Vup);
    Yv = Xv.cross(Zv);
    Xv.normalize();
    Yv.normalize();
    Zv.normalize();
    h = d * tan(theta / 2.0);
    Sj = 2 * h;
    Sk = Sj * (static_cast<double>(y) / static_cast<double>(x));
    P00 = cameraLoc + d * Zv - (Sj / 2) * Xv + (Sk / 2) * Yv;
}

void RTSettings::updateRayBox(Vector4d &rayBox) {
    uStart = rayBox[0];
    vStart = rayBox[1];
    uLength = rayBox[2];
    vLength = rayBox[3];
}

// smf files are parsed in this constructor
SMFModel::SMFModel(RTSettings *settings, std::ifstream *inSMF) {
    const std::regex vertexR("^v\\s*(.*)$");
    const std::regex faceR("^f\\s*(.*)$");
    std::smatch sm;
    string line;
    double temp1, temp2, temp3;
    MatrixXd tempRot = MatrixXd::Zero(4, 4);
    Eigen::Vector4d tempPt;
    Vector3d a, b, c, tempNorm;

    type = SMF;

    scale = settings->scale;
    rotate = settings->rotate * PI / 180.0;
    vertices = MatrixXd::Zero(0, 3);
    facesStored = MatrixXi::Zero(0, 3);

    minVec << DBL_MAX, DBL_MAX, DBL_MAX;
    maxVec << DBL_MIN, DBL_MIN, DBL_MIN;

    getline(*inSMF, line);

    while (!inSMF->eof()) {
        if (line.front() == '#' || !line.compare("")) {
        } else if (std::regex_search(line, sm, vertexR)) {
            vertices.conservativeResize(vertices.rows() + 1, Eigen::NoChange);
            std::stringstream ss(sm[1]);
            ss >> temp1;
            ss >> temp2;
            ss >> temp3;
            vertices.row(vertices.rows() - 1) << temp1, temp2, temp3;
        } else if (std::regex_search(line, sm, faceR)) {
            facesStored.conservativeResize(facesStored.rows() + 1, Eigen::NoChange);
            std::stringstream ss(sm[1]);
            ss >> temp1;
            ss >> temp2;
            ss >> temp3;
            facesStored.row(facesStored.rows() - 1) << temp1, temp2, temp3;
        }
        getline(*inSMF, line);
    }
    int i;
    // Transformations
    for (i = 0; i < vertices.rows(); i++) {
        tempPt << vertices(i, 0), vertices(i, 1), vertices(i, 2), 1;
        tempRot << 1, 0, 0, 0, 0, cos(rotate(0)), -sin(rotate(0)), 0, 0,
            sin(rotate(0)), cos(rotate(0)), 0, 0, 0, 0, 1;
        tempPt = tempRot * tempPt;
        tempRot << cos(rotate(1)), 0, sin(rotate(1)), 0, 0, 1, 0, 0,
            -sin(rotate(1)), 0, cos(rotate(1)), 0, 0, 0, 0, 1;
        tempPt = tempRot * tempPt;
        tempRot << cos(rotate(2)), -sin(rotate(2)), 0, 0, sin(rotate(2)),
            cos(rotate(2)), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
        tempPt = tempRot * tempPt;
        for (int l = 0; l < 3; l++) {
            if (tempPt(l) < minVec(l)) minVec(l) = tempPt(l);
            if (tempPt(l) > maxVec(l)) maxVec(l) = tempPt(l);
        }
        vertices.row(i) = tempPt.head(3);
    }
    vertices = vertices.rowwise() - minVec;//.transpose();
    maxVec -= minVec;
    minVec -= minVec;
    // Normal averaging calculation
    normals = MatrixXd::Zero(vertices.rows(), 3);
    for (i = 0; i < facesStored.rows(); i++) {
        a = vertices.row(facesStored(i, 0) - 1);
        b = vertices.row(facesStored(i, 1) - 1);
        c = vertices.row(facesStored(i, 2) - 1);
        tempNorm = (b - a).cross(c - a).normalized();
        for (int j = 0; j < 3; j++)
            normals.row(facesStored(i, j) - 1) += tempNorm;
    }
    for (i = 0; i < normals.rows(); i++)
        normals.row(i) = normals.row(i).normalized();
}

bool SMFModel::intersectPixel(Vector3d ray, RTSettings *settings, int pixelNum,
                              int objNum) {
    return intersectPixel(ray, settings, &facesStored, pixelNum, objNum);
}

// Intersects pixel with triangles of an SMF model and updates data stores
// accordingly
bool SMFModel::intersectPixel(Vector3d ray, RTSettings *settings, MatrixXi *facesIn,
                              int pixelNum, int objNum) {
    MatrixXi faces = *facesIn;
    Matrix3d temp;
    Vector3d tempVec, tempNorm;
    MatrixXd calc = MatrixXd::Zero(3, 4);
    Vector3d R = settings->cameraLoc;
    Vector3d a, b, c;
    double Adet, beta, gamma, t;
    bool updated = false;
    for (int i = 0; i < faces.rows(); i++) {
        a = vertices.row(faces.row(i)(0) - 1);
        b = vertices.row(faces.row(i)(1) - 1);
        c = vertices.row(faces.row(i)(2) - 1);
        // Cache these calculations for small performance gain
        calc << a(0) - b(0), a(0) - c(0), a(0) - R(0), ray(0), a(1) - b(1),
            a(1) - c(1), a(1) - R(1), ray(1), a(2) - b(2), a(2) - c(2),
            a(2) - R(2), ray(2);
        temp = calc(all, {0, 1, 3});
        Adet = temp.determinant();
        temp = calc(all, {2, 1, 3});
        beta = temp.determinant() / Adet;
        if (beta >= 0) {
            temp = calc(all, {0, 2, 3});
            gamma = temp.determinant() / Adet;
            if (gamma >= 0 && beta + gamma <= 1) {
                temp = calc(all, {0, 1, 2});
                t = temp.determinant() / Adet;
                if (t >= 0) {
                    if (tObjsBuffer(pixelNum) == -1 ||
                        t <= tBuffer(pixelNum, 0)) {
                        updated = true;
                        tempVec = (1 - beta - gamma) * a + beta * b + gamma * c;
                        tempNorm = tempVec;
                        tempNorm =
                            ((1 - beta - gamma) * normals.row(faces(i, 0) - 1) +
                             beta * normals.row(faces(i, 1) - 1) +
                             gamma * normals.row(faces(i, 2) - 1))
                                .normalized();
                        tBuffer.row(pixelNum) << t, tempVec(0), tempVec(1),
                            tempVec(2), tempNorm(0), tempNorm(1),
                            tempNorm(2);
                        tObjsBuffer(pixelNum) = objNum;
                    }
                }
            }
        }
    }
    return updated;
}

Box::Box(SMFModel *smf) {
    minVec = Vector3d(smf->vertices.col(0).minCoeff(),
                      smf->vertices.col(1).minCoeff(),
                      smf->vertices.col(2).minCoeff());
    maxVec = Vector3d(smf->vertices.col(0).maxCoeff(),
                      smf->vertices.col(1).maxCoeff(),
                      smf->vertices.col(2).maxCoeff());
    type = BOX;
}

Box::Box(SubSMF *subsmf) {
    MatrixXi faces = subsmf->faces;
    // Just defaulting these to something
    Vector3d currRow = subsmf->smf->vertices.row(faces(0, 0) - 1);
    minVec = currRow.replicate(1, 1);
    maxVec = currRow.replicate(1, 1);
    // Loop through faces in sub smf only
    for (int i = 0; i < faces.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            currRow = subsmf->smf->vertices.row(faces(i, j) - 1);
            for (int k = 0; k < 3; k++) {
                if (currRow(k) < minVec(k))
                    minVec(k) = currRow(k);
                else if (currRow(k) > maxVec(k))
                    maxVec(k) = currRow(k);
            }
        }
    }
    type = BOX;
}

// Perform box intersection algorithm from the slides
bool Box::intersectPixel(Vector3d Rd, RTSettings *settings, int pixelNum,
                         int objNum) {
    double temp;
    Vector3d t1 = (minVec - settings->cameraLoc).array() / Rd.array();
    Vector3d t2 = (maxVec - settings->cameraLoc).array() / Rd.array();
    for (int i = 0; i < 3; i++) {
        if (t1(i) > t2(i)) {
            temp = t1(i);
            t1(i) = t2(i);
            t2(i) = temp;
        }
    }
    double tnear = t1.maxCoeff();
    double tfar = t2.minCoeff();
    return (tnear <= tfar) && (tfar >= 0);
}

// Sets up the tree data structure for an object's bounding volume hierarchy
ACGBVH::ACGBVH(Intersectable *inObj, int inMaxObjCount) {
    type = HIERARCH;
    obj = inObj;
    maxObjCount = inMaxObjCount;
    if (obj->type == SMF) {
        boundingVol = new Box(reinterpret_cast<SMFModel *>(inObj));
        if ((reinterpret_cast<SMFModel *>(inObj))->facesStored.rows() > maxObjCount) {
            isLeaf = false;
            createChildren();
        } else {
            isLeaf = true;
        }
    } else if (obj->type == SUBSMF) {
        boundingVol = new Box(reinterpret_cast<SubSMF *>(inObj));
        if ((reinterpret_cast<SubSMF *>(inObj))->faces.rows() > maxObjCount) {
            isLeaf = false;
            createChildren();
        } else {
            isLeaf = true;
        }
    } else {
        cout << "Error: invalid shape type '" << obj->type
             << "' passed to ACGBVH\n";
        exit(1);
    }
}

// Defers to its children or object stored
bool ACGBVH::intersectPixel(Vector3d ray, RTSettings *settings, int pixelNum,
                            int objNum) {
    if (boundingVol->intersectPixel(ray, settings, pixelNum, objNum)) {
        if (isLeaf) {
            return obj->intersectPixel(ray, settings, pixelNum, objNum);
        } else {
            children[0]->intersectPixel(ray, settings, pixelNum, objNum);
            children[1]->intersectPixel(ray, settings, pixelNum, objNum);
            return true;
        }
    }
    return false;
}

// Does the volume division
void ACGBVH::createChildren() {
    Eigen::Array3d diff = boundingVol->maxVec - boundingVol->minVec;
    Eigen::Array3d currRow;
    double split;
    int axis;
    SMFModel *smf;
    MatrixXi newFaces1 = MatrixXi::Zero(0, 3);
    MatrixXi newFaces2 = MatrixXi::Zero(0, 3);
    MatrixXi currFaces;
    if (diff(0) >= diff(1) && diff(0) >= diff(2))
        axis = 0;
    else if (diff(1) > diff(0) && diff(1) >= diff(2))
        axis = 1;
    else
        axis = 2;
    split = boundingVol->minVec(axis) + (diff(axis) / 2.0);
    if (obj->type == SUBSMF) {
        smf = (reinterpret_cast<SubSMF *>(obj))->smf;
        currFaces = (reinterpret_cast<SubSMF *>(obj))->faces.replicate(1, 1);
    } else {
        smf = reinterpret_cast<SMFModel *>(obj);
        currFaces = smf->facesStored.replicate(1, 1);
    }
    for (int i = 0; i < currFaces.rows(); i++) {
        currRow << smf->vertices(currFaces(i, 0) - 1, axis),
            smf->vertices(currFaces(i, 1) - 1, axis),
            smf->vertices(currFaces(i, 2) - 1, axis);
        if ((currRow > split).any()) {
            newFaces1.conservativeResize(newFaces1.rows() + 1,
                                         newFaces1.cols());
            newFaces1.row(newFaces1.rows() - 1) = currFaces.row(i);
        } else {
            newFaces2.conservativeResize(newFaces2.rows() + 1,
                                         newFaces2.cols());
            newFaces2.row(newFaces2.rows() - 1) = currFaces.row(i);
        }
    }
    SubSMF *subsmf1 = new SubSMF(smf, newFaces1);
    SubSMF *subsmf2 = new SubSMF(smf, newFaces2);
    ACGBVH *temp1 = new ACGBVH(subsmf1, maxObjCount);
    children.push_back(temp1);
    ACGBVH *temp2 = new ACGBVH(subsmf2, maxObjCount);
    children.push_back(temp2);
}

SubSMF::SubSMF(SMFModel *smfIn, MatrixXi facesIn) {
    type = SUBSMF;
    smf = smfIn;
    faces = facesIn;
}

bool SubSMF::intersectPixel(Vector3d ray, RTSettings *settings, int pixelNum,
                            int objNum) {
    return smf->intersectPixel(ray, settings, &faces, pixelNum, objNum);
}
