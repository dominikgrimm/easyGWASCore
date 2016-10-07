#include "CKernels.h"
#include "CEasyGWAS/stats/CStats.h"

/*
*Compute realized Relationship Kernel
*/
MatrixXd CKernels::realizedRelationshipKernel(MatrixXd const& X) {
	MatrixXd Xn = (X-CStats::mean(X,1).transpose().replicate(X.rows(),1)).array()/CStats::std(X,1).transpose().replicate(X.rows(),1).array();
	return 1.0/Xn.cols() * Xn * Xn.transpose();
}

/*
*Compute Linear Kernel
*/
MatrixXd CKernels::linearKernel(MatrixXd const& X) {
	MatrixXd Xn = (X-CStats::mean(X,1).transpose().replicate(X.rows(),1)).array()/CStats::std(X,1).transpose().replicate(X.rows(),1).array();
	return 1.0/Xn.cols() * Xn * Xn.transpose();
}

/*
*Center a kernel X
*/
MatrixXd CKernels::centerKernel(MatrixXd const& X) {
    MatrixXd identity = MatrixXd::Identity(X.rows(),X.cols());
    MatrixXd H = identity - (1.0/X.rows())*MatrixXd::Ones(X.rows(),X.cols());
    return (H*(X*H));
}
