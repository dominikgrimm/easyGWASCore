#ifndef CStatsClass
#define CStatsClass

#include <Eigen/Dense>

#include <string>

#include "CEasyGWAS/globals.h"

class CStats {
	public:
		//Pearsons correlation and Pvalue calculations
		static float64 pearson_corr(VectorXd const&, VectorXd const&);
		static float64 pearson_pval(float64 const&, uint64 const&);
		static float64 pearson_pval(float64 const&, uint64 const&,std::string const&);
	    
        //Fisher exact test
        //static float64 fisher_exact(MatrixXd const&)

		//Various helper methods
		
		//Compute VARIANCE of a Vector X and return a float
		static float64 varf(Eigen::VectorXd const&);
		
        //Compute STD and return a float
		static float64 stdf(Eigen::VectorXd const&);
		static float64 stdf(Eigen::VectorXd const&, bool const&);
		//Compute STD on an Matrix array and return a Vector
		static VectorXd std(MatrixXd const&);
		static VectorXd std(MatrixXd const&, uint const&);
		static VectorXd std(MatrixXd const&, uint const&, bool const&);

		static VectorXd mean(MatrixXd const&);
		static VectorXd mean(MatrixXd const&, uint const&);

        static MatrixXd principle_components(MatrixXd const&);
};

#endif //CStatsClass
