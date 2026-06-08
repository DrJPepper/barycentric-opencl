#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Geometry>

typedef double num;

typedef Eigen::Translation<num, 2> Translation2d;
typedef Eigen::Matrix<num, 1, Eigen::Dynamic> VectorXd;
typedef Eigen::Matrix<num, 2, 1> ColVector2d;
typedef Eigen::Matrix<num, 1, 2> Vector2d;
typedef Eigen::Matrix<num, 1, 3> Vector3d;
typedef Eigen::Matrix<num, 1, 4> Vector4d;
typedef Eigen::DenseBase<Eigen::Matrix<num, -1, -1, 1> >::RowXpr Row3d;
typedef Eigen::Matrix<num, 3, 1> ColVector3d;
typedef Eigen::Matrix<int, 1, 2> Vector2i;
typedef Eigen::Matrix<int, 1, 3> Vector3i;
typedef Eigen::Matrix<int, 1, 4> Vector4i;
typedef Eigen::Matrix<int, 1, Eigen::Dynamic> VectorXi;
typedef Eigen::Array<int, 1, 3> Array3i;
typedef Eigen::Matrix3d Matrix3d;
typedef Eigen::Matrix<num,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatrixXd;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatrixXi;
#if EIGEN_VERSION_AT_LEAST(3, 4, 0)
    constexpr auto all = Eigen::placeholders::all;
    constexpr auto last = Eigen::placeholders::all;
#else
    constexpr auto all = Eigen::all;
    constexpr auto last = Eigen::all;
#endif
