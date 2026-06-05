#include "./bezier.h"

// Populates matrix of points
void initialize(MatrixXd *points, string inFileName) {
    std::ifstream inFile(inFileName);
    // Check the file
    if (!inFile.good()) {
        cerr << "ERROR: input file invalid or doesn't exist" << endl;
        exit(1);
    }

    float x, y, z;

    for (string line; getline(inFile, line);) {
        points->conservativeResize(points->rows() + 1, points->cols());
        std::istringstream in(line);
        in >> x >> y >> z;
        points->row(points->rows() - 1) << x, y, z;
    }
}

// Runs p equation
MatrixXd p(MatrixXd *points, float u, float v) {
    MatrixXd acc = MatrixXd::Zero(1, 3);
    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3; j++) {
            acc += B(u, i) * B(v, j) * points->row(i + j * 4);
        }
    }
    return acc;
}

MatrixXd Q(MatrixXd *points, float u) {
    MatrixXd acc = MatrixXd::Zero(1, 3);
    for (int i = 0; i < 4; i++) {
        acc += points->row(i) * kChoosei(3, i) * pow(1 - u, 3 - i) * pow(u, i);
    }
    return acc;
}

// Runs B equation
float B(float u, int i) {
    // n is always 3
    return kChoosei(3, i) * pow(u, i) * pow(1 - u, 3 - i);
}

// Calculates the partial derivative with respect to u at a point
Vector3d dQdu(MatrixXd *points, float u, float v) {
    MatrixXd curveV = MatrixXd::Zero(4, 3);
    MatrixXd temp = MatrixXd::Zero(4, 3);
    for (int i = 0; i < 4; i++) {
        temp.row(0) = points->row(i);
        temp.row(1) = points->row(i + 4);
        temp.row(2) = points->row(i + 8);
        temp.row(3) = points->row(i + 12);
        curveV.row(i) = Q(&temp, v);
    }
    return -3 * (1 - u) * (1 - u) * (Vector3d)curveV.row(0) +
           (3 * (1 - u) * (1 - u) - 6 * u * (1 - u)) * (Vector3d)curveV.row(1) +
           (6 * u * (1 - u) - 3 * u * u) * (Vector3d)curveV.row(2) +
           3 * u * u * (Vector3d)curveV.row(3);
}

// Calculates the partial derivative with respect to v at a point
Vector3d dQdv(MatrixXd *points, float u, float v) {
    MatrixXd curveU = MatrixXd::Zero(4, 3);
    for (int i = 0; i < 4; i++) {
        MatrixXd temp = points->block(i * 4, 0, 4, 3);
        curveU.row(i) = Q(&temp, u);
    }
    return -3 * (1 - v) * (1 - v) * (Vector3d)curveU.row(0) +
           (3 * (1 - v) * (1 - v) - 6 * v * (1 - v)) * (Vector3d)curveU.row(1) +
           (6 * v * (1 - v) - 3 * v * v) * (Vector3d)curveU.row(2) +
           3 * v * v * (Vector3d)curveU.row(3);
}

// Calculates the normal vector at a point
MatrixXd normal(MatrixXd *points, float u, float v) {
    Vector3d tanU = dQdu(points, u, v);
    Vector3d tanV = dQdv(points, u, v);
    return tanV.cross(tanU).normalized();
}

int kChoosei(int k, int i) {
    return factorial(k) / (factorial(i) * factorial(k - i));
}

int factorial(int x) {
    if (x == 0) return 1;
    return x * factorial(x - 1);
}
