#include "util.h"

Vector3d triangleNormal(MatrixXd &points, Vector3i &face) {
    Vector3d normal, p1 = points.row(face(0)), p2 = points.row(face(1)), p3 = points.row(face(2)), v1, v2;
    
    v1 = p2 - p1;
    v2 = p3 - p1;
    
    normal = v1.cross(v2).normalized();
    
    return normal;
}

num pow2(num x) {
    return x * x;
}

num scaleneArea(Vector3d p1, Vector3d p2, Vector3d p3) {
    num e1 = distance(p1, p2);
    num e2 = distance(p1, p3);
    num e3 = distance(p2, p3);
    num s = (e1 + e2 + e3) / 2;

    return sqrt(s * (s - e1) * (s - e2) * (s - e3));
}

num tetrahedronVolume(Vector3d p1, Vector3d p2, Vector3d p3, Vector3d p4) {
    num U = distance(p1, p3);
    num V = distance(p2, p3);
    num W = distance(p1, p2);
    num u = distance(p2, p4);
    num v = distance(p1, p4);
    num w = distance(p3, p4);

    num X = (w - U + v) * (U + v + w);
    num Y = (u - V + w) * (V + w + u);
    num Z = (v - W + u) * (W + u + v);
    num x = (U - v + w) * (v - w + U);
    num y = (V - w + u) * (w - u + V);
    num z = (W - u + v) * (u - v + W);

    num p = sqrt(x * Y * Z);
    num q = sqrt(y * Z * X);
    num r = sqrt(z * X * Y);
    num s = sqrt(x * y * z);

    num result =  sqrt((-p + q + r + s) * (p - q + r + s) *
            (p + q - r + s) * (p + q + r - s)) / (192 * u * v * w);
    if (isnan(result))
        result = 0.0;
    return result;
}

bool barycentric(const Vector3d &p, Vector3d &a, Vector3d &b, Vector3d &c, num &u, num &v, num &w, num &d, num tol) {
    return barycentric(p, a, b, c, u, v, w, d, true, tol);
}

bool barycentric(const Vector3d &p, Vector3d &a, Vector3d &b, Vector3d &c, num &u, num &v, num &w, num &d) {
    return barycentric(p, a, b, c, u, v, w, d, true, BARY_TOL);
}

bool barycentric(const Vector3d &p, Vector3d &a, Vector3d &b, Vector3d &c, num &u, num &v, num &w) {
    num d;
    return barycentric(p, a, b, c, u, v, w, d, false, BARY_TOL);
}

bool barycentric(const Vector3d &p, Vector3d &a, Vector3d &b, Vector3d &c, num &u, num &v, num &w, num &d, bool doTetrahedron, num tol) {
    Vector3d v0 = b - a, v1 = c - a, v2 = p - a;
    num d00 = v0.dot(v0);
    num d01 = v0.dot(v1);
    num d11 = v1.dot(v1);
    num d20 = v2.dot(v0);
    num d21 = v2.dot(v1);
    num denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0 - v - w;
    int checkSum = (u >= NEP) + (u <= OPEP) + (v >= NEP) + (v <= OPEP) + (w >= NEP) + (w <= OPEP);
    bool check = checkSum == 6;
    //bool check = (u >= -1e-3 && u <= 1.001 && v >= -1e-3 && v <= 1.001 && w >= -1e-3 && w <= 1.001);
    if (doTetrahedron && checkSum >= 4)
        d = tetrahedronVolume(p, a, b, c);
    else
        d = -1.0;
    return check;
}

num bound(num x, num minVal, num maxVal) {
    return std::min(maxVal, std::max(minVal, x));
}


num bound(num x) {
    return bound(x, 0.00001, 0.99999);
}

int kChoosei(int k, int i) {
    return factorial(k) / (factorial(i) * factorial(k - i));
}

int factorial(int x) {
    if (x == 0) return 1;
    return x * factorial(x - 1);
}
