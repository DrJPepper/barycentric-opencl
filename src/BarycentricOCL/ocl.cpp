#include "kernel.hpp"

string opencl_c_container() { return R( // begin of OpenCL C code

inline double distance3(double3 one, double3 two) {
    double a = two[0] - one[0];
    double b = two[1] - one[1];
    double c = two[2] - one[2];
    return sqrt(a*a + b*b + c*c);
}

kernel void barycentric_loop_F(global int* check, global double* new_point,
        global double* point, global double* points_from, global int* faces,
        global double* points_to, global uint* Ns) {
    const uint n = get_global_id(0);
    const uint n3 = n * 3;
    const uint pn = n / Ns[0];
    const uint fn = n % Ns[0];
    const uint pn3 = pn * 3;
    const uint fn3 = fn * 3;
    double3 an, bn, cn;
    int3 vs = (int3)(faces[fn3], faces[fn3+1], faces[fn3+2]);
    uint linear_tri = 0;
    for (int j=0; j<3; j++) {
        an[j] = points_from[vs[0]+j];
        bn[j] = points_from[vs[1]+j];
        cn[j] = points_from[vs[2]+j];
        if (j < 2)
            linear_tri |= an[j] == bn[j] && an[j] == cn[j] && bn[j] == cn[j];
    }
    double3 p = (double3)(point[pn3], point[pn3+1], point[pn3+2]);
    double3 v0 = bn - an;
    double3 v1 = cn - an;
    double3 v2 = p - an;
    double d00 = 0, d01 = 0, d11 = 0, d20 = 0, d21 = 0;
    for (uint i=0; i<3; i++) {
        d00 += v0[i] * v0[i];
        d01 += v0[i] * v1[i];
        d11 += v1[i] * v1[i];
        d20 += v2[i] * v0[i];
        d21 += v2[i] * v1[i];
    }
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;

    const double EP = 1e-8;
    const double NM = 0.0 - EP;
    const double PM = 1.0 + EP;
    const uint curr_check = u >= NM && u <= PM && v >= NM && v <= PM && w >= NM && w <= PM;
    //if (pn == 49 && curr_check)
        //printf("%d,%d,%d,%d,%d,%d\n", u >= NM, u <= PM, v >= NM, v <= PM, w >= NM, w <= PM);
        //printf("%f,%f,%f\n", v, w, u);
        //printf("%d|%d,%d,%d|%f,%f,%f|%f,%f,%f|%f,%f,%f|%f,%f,%f\n", linear_tri, vs[0]/3, vs[1]/3, vs[2]/3, p[0], p[1], p[2], an[0], an[1], an[2], bn[0], bn[1], bn[2], cn[0], cn[1], cn[2]);
    if (curr_check && !linear_tri) {
        for (int j=0; j<3; j++) {
            an[j] = points_to[vs[0]+j];
            bn[j] = points_to[vs[1]+j];
            cn[j] = points_to[vs[2]+j];
            new_point[pn3+j] = u * an[j] + v * bn[j] + w * cn[j];
        }
        //new_point[pn3] = u * (an[0] + bn[0] + cn[0]);
        //new_point[pn3+1] = v * (an[1] + bn[1] + cn[1]);
        //new_point[pn3+2] = w * (an[2] + bn[2] + cn[2]);

        check[pn] = fn3;
        //if (pn == 57)
            //printf("%d,%d|%f,%f,%f\n", pn, fn, new_point[pn3], new_point[pn3+1], new_point[pn3+1]);
    }
}

kernel void barycentric_loop_R(global int* check, global double* new_point,
        global double* volume, global double* normal, global double* point,
        global double* points_to, global int* faces, global uint* Ns,
        global double* float_check) {
    uint setVals = 0;
    const uint n = get_global_id(0);
    const uint n3 = n * 3;
    const uint pn = n / Ns[0];
    const uint fn = n % Ns[0];
    const uint pn3 = pn * 3;
    const uint fn3 = fn * 3;
    double3 an, bn, cn;
    int3 vs = (int3)(faces[fn3], faces[fn3+1], faces[fn3+2]);
    for (int j=0; j<3; j++) {
        an[j] = points_to[vs[0]+j];
        bn[j] = points_to[vs[1]+j];
        cn[j] = points_to[vs[2]+j];
    }
    double3 p = (double3)(point[pn3], point[pn3+1], point[pn3+2]);
    double3 v0 = bn - an;
    double3 v1 = cn - an;
    double3 v2 = p - an;
    double d00 = 0, d01 = 0, d11 = 0, d20 = 0, d21 = 0;
    for (uint i=0; i<3; i++) {
        d00 += v0[i] * v0[i];
        d01 += v0[i] * v1[i];
        d11 += v1[i] * v1[i];
        d20 += v2[i] * v0[i];
        d21 += v2[i] * v1[i];
    }
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;

    const double EP = 1e-2;
    const double NM = 0.0 - EP;
    const double PM = 1.0 + EP;
    const uint curr_check = u >= NM && u <= PM && v >= NM && v <= PM && w >= NM && w <= PM;
    if (curr_check) {
        double e1 = distance3(an, bn);
        double e2 = distance3(an, cn);
        double e3 = distance3(bn, cn);
        double ss = (e1 + e2 + e3) / 2;
        double sa = sqrt(ss * (ss - e1) * (ss - e2) * (ss - e3));
        if (sa > 0.0) {
            // TODO: Don't recalculate these over and over
            //double U = distance3(an, cn);
            double U = e2;
            //double V = distance3(bn, cn);
            double V = e3;
            //double W = distance3(an, bn);
            double W = e1;

            double uu = distance3(bn, p);
            double vv = distance3(an, p);
            double ww = distance3(cn, p);

            double X = (ww - U + vv) * (U + vv + ww);
            double Y = (uu - V + ww) * (V + ww + uu);
            double Z = (vv - W + uu) * (W + uu + vv);
            double x = (U - vv + ww) * (vv - ww + U);
            double y = (V - ww + uu) * (ww - uu + V);
            double z = (W - uu + vv) * (uu - vv + W);

            double pp = sqrt(x * Y * Z);
            double q = sqrt(y * Z * X);
            double r = sqrt(z * X * Y);
            double s = sqrt(x * y * z);

            double numerator = sqrt((-pp + q + r + s) * (pp - q + r + s) *
                    (pp + q - r + s) * (pp + q + r - s));
            double denominator = 192 * uu * vv * ww;

            double volval;
            if (!numerator || !denominator)
                volval = 0.0;
            else
                volval =  numerator / denominator;
            setVals = 1;
            volume[n] = volval;

            new_point[n3] = u;
            new_point[n3+1] = v;
            new_point[n3+2] = w;

            check[n] = fn3;

            double3 UU = bn - an;
            double3 VV = cn - an;
            normal[n3] = UU[1] * VV[2] - UU[2] * VV[1];
            normal[n3+1] = UU[2] * VV[0] - UU[0] * VV[2];
            normal[n3+2] = UU[0] * VV[1] - UU[1] * VV[0];
            //printf("%d\n", pn);
        }
    }
    if (!setVals) {
        volume[n] = -1.0;

        new_point[n3] = 0.0;
        new_point[n3+1] = 0.0;
        new_point[n3+2] = 0.0;

        check[n] = -1;

        normal[n3] = 0.0;
        normal[n3+1] = 0.0;
        normal[n3+2] = 0.0;
        float_check[n] = fmin(u-NM,0) + fmin(PM-u,0) + fmin(v-NM,0) +
            fmin(PM-u,0) + fmin(w-NM,0) + fmin(PM-w,0);
    }
}

kernel void find_best_points(global int* indices, global double* points_to, global int* faces,
        global double* uvw_in, global double* volumes_in, global double* normals_in,
        global double* points_out, global double* volumes_out, global double* normals_out,
        global uint* Ns, global double* float_check) {
    const int pn = get_global_id(0);
    if (pn < Ns[1]) {
        const int pn3 = pn * 3;
        const int start_ind = pn * Ns[0];
        int c, ind, bestInd = -1, bestFace = -1, bestIndBackup = -1, bestFaceBackup = -1;
        double minVol = DBL_MAX, vol, minCheck = DBL_MAX, check;
        for (int fn=0; fn<Ns[0]; fn++) {
            ind = start_ind + fn;
            c = indices[ind];
            vol = volumes_in[ind];
            check = float_check[ind];
            if (c != -1 && vol < minVol) {
                minVol = vol;
                bestInd = ind * 3;
                bestFace = fn*3;
            } else if (check < minCheck) {
                minCheck = check;
                bestIndBackup = ind * 3;
                bestFaceBackup = fn*3;
            }
        }
        if (bestInd == -1) {
            //printf("%d\n", bestIndBackup);
            bestInd = bestIndBackup;
            bestFace = bestFaceBackup;
        }
        int3 vs = (int3)(faces[bestFace], faces[bestFace+1], faces[bestFace+2]);
        double3 p1 = (double3)(points_to[vs[0]], points_to[vs[0]+1], points_to[vs[0]+2]);
        double3 p2 = (double3)(points_to[vs[1]], points_to[vs[1]+1], points_to[vs[1]+2]);
        double3 p3 = (double3)(points_to[vs[2]], points_to[vs[2]+1], points_to[vs[2]+2]);
        double3 p_out = p1 * uvw_in[bestInd] + p2 * uvw_in[bestInd+1] + p3 * uvw_in[bestInd+2];
        volumes_out[pn] = minVol;
        for (int j=0; j<3; j++) {
            points_out[pn3+j] = p_out[j];
            normals_out[pn3+j] = normals_in[bestInd+j];
        }
        //} else {
            //volumes_out[pn] = -1.0;
        //}
    }
}

);} // end of OpenCL C code
