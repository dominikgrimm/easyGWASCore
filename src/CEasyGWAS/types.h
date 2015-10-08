#ifndef TYPES
#define TYPES

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "time.h"


typedef double float64;
typedef float float32;
typedef long long int64;
typedef unsigned long long uint64;
typedef unsigned int uint;

typedef Eigen::Matrix<float64,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> MatrixXd;
typedef Eigen::Matrix<float64,Eigen::Dynamic,1,Eigen::ColMajor> VectorXd;
typedef Eigen::DiagonalMatrix<float64,Eigen::Dynamic> DiagXd;
typedef Eigen::SparseMatrix<float64,Eigen::ColMajor> SparseMatrixXd;

//For wrapping python with c++
typedef Eigen::Matrix<float64,Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdNumpy;
typedef Eigen::Matrix<float64,Eigen::Dynamic, 1> VectorXdNumpy;

typedef Eigen::Triplet<float64> eigen_triplet;
typedef std::vector<eigen_triplet> sparse_triplet_vector;

#endif //TYPES
