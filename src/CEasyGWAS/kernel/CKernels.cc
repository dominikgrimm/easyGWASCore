#include "CKernels.h"
#include "CEasyGWAS/stats/CStats.h"

MatrixXd CKernels::realizedRelationshipKernel(MatrixXd const& X) {
	MatrixXd Xn = (X-CStats::mean(X,1).transpose().replicate(X.rows(),1)).array()/CStats::std(X,1).transpose().replicate(X.rows(),1).array();
	return 1.0/Xn.cols() * Xn * Xn.transpose();
}
