#ifndef CMATRIX_HELPER
#define CMATRIX_HELPER

#include <fstream>

#include "CEasyGWAS/globals.h"

/*
*Slice Matrix: Remove all rows in indices from a matrix 
*/
inline MatrixXd sliceRowsMatrix(MatrixXd const& m, VectorXd const& indices) {
	MatrixXd out(indices.rows(),m.cols());
	for(int64 i=0; i<indices.rows(); i++) {
		out.block(i,0,1,m.cols()) = m.block(indices(i),0,1,m.cols());
	}
	return out;
}

/*
*Slice Matrix: Remove all cols in indices from a matrix 
*/
inline MatrixXd sliceColsMatrix(MatrixXd const& m, VectorXd const& indices) {
	MatrixXd out(m.rows(),indices.rows());
	for(int64 i=0; i<indices.rows(); i++) {
		out.col(i) = m.col(indices(i)); 
	}
	return out;
}

/*
*Insert data at indices ind 
*/
inline void insertColumnVectorAtIndices(VectorXd* result, VectorXd const& in, VectorXd const& indices) {
	for(int64 i=0; i<in.rows(); i++) {
		(*result)(indices(i)) = in(i);
	}
}

/*
*Save matrix into a binary file 
*/
inline void dumpMatrix(MatrixXd const& m, char const* filename) {
	std::ofstream f(filename, std::ios::binary);
	int64 cols = m.cols();
	int64 rows = m.rows();
	f.write((char*)&cols, sizeof(cols));
	f.write((char*)&rows, sizeof(rows));
	f.write((char*)(m.data()), sizeof(MatrixXd::Scalar)*cols*rows);
	f.close();
}

/*
*Load matrix from a binary file 
*/
inline void loadMatrix(char const* filename, MatrixXd* m) {
	std::ifstream f(filename, std::ios::binary);
	MatrixXd::Index rows,cols;
	f.read((char*)&cols,sizeof(cols));
	f.read((char*)&rows,sizeof(rows));
	m->resize(rows,cols);
	f.read((char*)(m->data()),sizeof(MatrixXd::Scalar)*rows*cols);
	if(f.bad()) logging(ERROR, "Loading matrix file");
	f.close();
}

/*
*Permute Vector <permutations> times and return permutation matrix
*/
inline MatrixXd permuteVector(VectorXd const& y, uint const& permutations) {
	uint n_samples = y.rows();
	MatrixXd Y = MatrixXd::Zero(n_samples,permutations);
	Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(n_samples);
	perm.setIdentity();
	for(uint i=0; i<permutations;i++) {
		std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
		Y.block(0,i,n_samples,1) = perm*y;
	}
	return Y;
}

#endif //CMatrixHelper
