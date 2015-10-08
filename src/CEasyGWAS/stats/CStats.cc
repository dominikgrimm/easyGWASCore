/*
*CStats is a class containing several small statistical helper functions
*This helper functions is using the Eigen library. More information about Eigen
*can be found at http://eigen.tuxfamily.org
*
*@author: Dominik Gerhard Grimm
*@year: 2014
*/
#include <cmath>

#include "CEasyGWAS/stats/CStudentT.h"
#include "CEasyGWAS/utils/CMisc.h"
#include "CEasyGWAS/utils/CMatrixHelper.h"
#include "CStats.h"

extern VectorXd* vpointer;

/*Pearson_corr is a method computing Pearsons correlation coefficent between two
 * vectors v1 and v2 
 *@param: v1 first dynamic Eigen VectorXd
 *@param: v2 second dynamic Eigen VectorXd
 *@return: double correlation coefficent
 */
float64 CStats::pearson_corr(VectorXd const& v1, VectorXd const& v2) {
	VectorXd v1_ = v1.array() - v1.mean();
	VectorXd v2_ = v2.array() - v2.mean();
	return (v1_.cwiseProduct(v2_)).sum()/(sqrt(v1_.array().pow(2).sum()*v2_.array().pow(2).sum()));
}

float64 CStats::pearson_pval(float64 const& r, uint64 const& n) {
	return pearson_pval(r,n,"both");
}

float64 CStats::pearson_pval(float64 const& r, uint64 const& n,std::string const& tail) {
	float64 t = r/(sqrt((1.0-pow(r,2))/(n-2.0)));
	if(tail.compare("right")==0) return 1.0-CStudentT::cdf(t,n-2);
	else if(tail.compare("left")==0) return CStudentT::cdf(t,n-2);
	else return 2.0*(1.0-CStudentT::cdf(t,n-2));
}

/*
*Various helper methods
*/
float64 CStats::varf(Eigen::VectorXd const& x) {
	VectorXd x_ = x.array() - x.mean();
    return 1.0/x.rows() * x_.array().pow(2).sum();
}

float64 CStats::stdf(Eigen::VectorXd const& v) {
	return CStats::stdf(v,false);
}

float64 CStats::stdf(Eigen::VectorXd const& v, bool const& flag) {
	float64 n = 0.0;
	if (flag==true) {
		n = 1.0/(v.rows()-1.0); 
	} else {
		n = 1.0/(v.rows()); 
	}
	return sqrt(n*((v.array()-v.mean()).array().pow(2)).sum());
}

VectorXd CStats::std(MatrixXd const& m) {
	return CStats::std(m,1,false);
}

VectorXd CStats::std(MatrixXd const& m, uint const& dim) {
	return CStats::std(m,dim,false);
}

VectorXd CStats::std(MatrixXd const& m, uint const& dim, bool const& flag) {
	VectorXd result;
	if(dim==0) {
		result.resize(m.rows());
		for(uint i=0; i<m.rows(); i++) {
			result(i) = CStats::stdf(m.row(i),flag);
		}
	} else if(dim==1) {
		result.resize(m.cols());
		for(uint i=0; i<m.cols(); i++) {
			result(i) = CStats::stdf(m.col(i),flag);
		}
	} else {
		logging(ERROR, "Dim has to be 0 or 1");
		throw -1;
	}
	return result;
}

VectorXd CStats::mean(MatrixXd const& m) {
	return CStats::mean(m,1);
}

VectorXd CStats::mean(MatrixXd const& m, uint const& dim) {
	VectorXd result;
	if(dim==0) {
		result.resize(m.rows());
		for(uint i=0; i<m.rows(); i++) {
			result(i) = m.row(i).mean();
		}
	} else if(dim==1) {
		result.resize(m.cols());
		for(uint i=0; i<m.cols(); i++) {
			result(i) = m.col(i).mean();
		}
	} else {
		logging(ERROR, "Dim has to be 0 or 1");
		throw -1;
	}
	return result;
}

MatrixXd CStats::principle_components(MatrixXd const& X) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(X);
    VectorXd evalues = solver.eigenvalues();
    MatrixXd pcs = solver.eigenvectors();
    //Sort eigenvalues and principle components in decreasing order
    ArgSort argsort(-evalues);
    VectorXd ids = argsort.getIndices();
    //Return sorted Eigenvectors
    return sliceColsMatrix(pcs,ids);
}
