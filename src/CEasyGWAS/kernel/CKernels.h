#ifndef CKERNEL_CLASS
#define CKERNEL_CLASS

#include "CEasyGWAS/globals.h"

class CKernels {
	
	public:
		static MatrixXd realizedRelationshipKernel(MatrixXd const&);
        static MatrixXd linearKernel(MatrixXd const&);
        static MatrixXd centerKernel(MatrixXd const&);
};

#endif //CKERNEL_CLASS
