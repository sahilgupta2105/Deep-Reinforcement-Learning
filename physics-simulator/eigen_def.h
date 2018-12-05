#ifndef EIGEN_DEF_H
#define EIGEN_DEF_H

#include <cassert>
#include "Eigen/Dense"
#include "Eigen/Sparse"

// used
typedef Eigen::Vector2d VEC2F;
typedef Eigen::Vector3d VEC3F;
typedef Eigen::Rotation2Dd ROTVEC2F;
typedef Eigen::Triplet<double> TRIPLET;
typedef Eigen::VectorXd VECTOR;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SimplicialLDLT< SpMat > PSolver;
typedef Eigen::Matrix3d MATRIX3;

#endif
