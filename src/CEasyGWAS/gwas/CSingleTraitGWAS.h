#ifndef C_SINGLETRAITGWAS_CLASS
#define C_SINGLETRAITGWAS_CLASS

#include "CEasyGWAS/globals.h"

#include "CEasyGWAS/regression/CRegression.h"
#include "CEasyGWAS/gwas/CGWASData.h"

/*
*CGWAS Exception Class
*/
class CGWASException {
	private:
		std::string __error_msg;
	public:
		CGWASException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CGWAS Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

namespace CSingleTraitGWAS {

	class LinearRegression {
		private:
			bool __covs_set;
			bool __intercept;
			uint __statistical_test;
			uint __n_features;
			uint __n_permutations;
			VectorXd __y;
			MatrixXd __X;
			MatrixXd __covs;
			CLinearRegression __null_model;
			VectorXd __p_values;
			VectorXd __p_values_permutations;
			VectorXd __permutation_pool;
			MatrixXd __betas;
			MatrixXd __se_betas;
			VectorXd __test_statistics;
			VectorXd __logLikelihoods;

			void __checkdata() throw (CGWASException);

		public:
			static const uint LRTEST = 0;
			static const uint FTEST = 1;
			
			LinearRegression();
			LinearRegression(VectorXd const&, MatrixXd const&) throw (CGWASException);
			LinearRegression(VectorXd const&, MatrixXd const&, MatrixXd const&) throw (CGWASException);

			void test_associations();
			void permutations();
			void permutations(uint const&);
			//VectorXd predict(MatrixXd const&) throw (CGWASException);
			//VectorXd predict(MatrixXd const&, MatrixXd const&) throw (CGWASException);
			void setTestStatistic(uint const&) throw (CGWASException);
			void setPhenotype(VectorXd const&) throw (CGWASException);
			void setGenotype(MatrixXd const&) throw (CGWASException);
			void setCovariates(MatrixXd const&) throw (CGWASException);
			void setIntercept(bool const&);

			float64 getLogLikelihoodNullModel();
			VectorXd getLogLikelihoodAlternativeModels();
			
			VectorXd getPValues();
			VectorXd getPermutationPValues();
			MatrixXd getBetas();
			MatrixXd getSEBetas();
			VectorXd getTestStatistics();
			float64 getAIC();
			float64 getBIC();
			float64 getAICc();

			GWASResults getResults();
	};
	
	class LogisticRegression {
		private:
			bool __covs_set;
			bool __intercept;
			uint __n_features;
			uint __n_permutations;
			VectorXd __y;
			MatrixXd __X;
			MatrixXd __covs;
			CLogisticRegression __null_model;
			VectorXd __p_values;
			VectorXd __p_values_permutations;
			VectorXd __permutation_pool;
			MatrixXd __betas;
			MatrixXd __se_betas;
			VectorXd __test_statistics;
			VectorXd __logLikelihoods;

			void __checkdata() throw (CGWASException);

		public:
			LogisticRegression();
			LogisticRegression(VectorXd const&, MatrixXd const&) throw (CGWASException);
			LogisticRegression(VectorXd const&, MatrixXd const&, MatrixXd const&) throw (CGWASException);

			void test_associations();
			void permutations();
			void permutations(uint const&);
			//VectorXd predict(MatrixXd const&) throw (CGWASException);
			//VectorXd predict(MatrixXd const&, MatrixXd const&) throw (CGWASException);
			void setPhenotype(VectorXd const&) throw (CGWASException);
			void setGenotype(MatrixXd const&) throw (CGWASException);
			void setCovariates(MatrixXd const&) throw (CGWASException);
			void setIntercept(bool const&);

			float64 getLogLikelihoodNullModel();
			VectorXd getLogLikelihoodAlternativeModels();
			
			VectorXd getPValues();
			VectorXd getPermutationPValues();
			MatrixXd getBetas();
			MatrixXd getSEBetas();
			VectorXd getTestStatistics();
			float64 getAIC();
			float64 getBIC();
			float64 getAICc();
			
			GWASResults getResults();
	};
	
	class EMMAX {
		private:
			bool __covs_set;
			bool __intercept;
			bool __REML;
			bool __use_brent;
			uint __n_features;
			uint __n_permutations;
			VectorXd __permutation_pool;
			VectorXd __y;
			MatrixXd __X;
			MatrixXd __K;
			MatrixXd __covs;
			CLinearMixedRegression __null_model;
			VectorXd __p_values;
			VectorXd __p_values_permutations;
			MatrixXd __betas;
			MatrixXd __se_betas;
			VectorXd __test_statistics;
			VectorXd __logLikelihoods;

			void __checkdata() throw (CGWASException);

		public:
			EMMAX();
			EMMAX(VectorXd const&, MatrixXd const&, MatrixXd const&) throw (CGWASException);
			EMMAX(VectorXd const&, MatrixXd const&, MatrixXd const&, MatrixXd const&) throw (CGWASException);

			void test_associations();
			void permutations();
			void permutations(uint const&);
			//VectorXd predict(MatrixXd const&, MatrixXd const&) throw (CGWASException);
			//VectorXd predict(MatrixXd const&, MatrixXd const&, MatrixXd const&) throw (CGWASException);
			void setPhenotype(VectorXd const&) throw (CGWASException);
			void setGenotype(MatrixXd const&) throw (CGWASException);
			void setCovariates(MatrixXd const&) throw (CGWASException);
			void setK(MatrixXd const&) throw (CGWASException);
			void setIntercept(bool const&);
			void setREML(bool const&);
			void setBrent(bool const&);

            float64 computeVarianceExplainedNullModel(uint const&);
			float64 getHeritabilityEstimate();
			float64 getGeneticVariance();
            float64 getNoiseVariance();
            float64 getLogLikelihoodNullModel();
			VectorXd getLogLikelihoodAlternativeModels();
			
			float64 getLogDelta0();
			VectorXd getPValues();
			VectorXd getPermutationPValues();
			MatrixXd getBetas();
			MatrixXd getSEBetas();
			VectorXd getTestStatistics();
			float64 getAIC();
			float64 getBIC();
			float64 getAICc();
			
			GWASResults getResults();
	};
	
	class FaSTLMM {
		private:
			bool __covs_set;
			bool __intercept;
			bool __REML;
			bool __use_brent;
			uint __n_features;
			uint __n_permutations;
			VectorXd __y;
			MatrixXd __X;
			MatrixXd __K;
			MatrixXd __covs;
			CLinearMixedRegression __null_model;
			VectorXd __p_values;
			VectorXd __p_values_permutations;
			VectorXd __permutation_pool;
			MatrixXd __betas;
			MatrixXd __se_betas;
			VectorXd __test_statistics;
			VectorXd __logLikelihoods;

			void __checkdata() throw (CGWASException);

		public:
			FaSTLMM();
			FaSTLMM(VectorXd const&, MatrixXd const&, MatrixXd const&) throw (CGWASException);
			FaSTLMM(VectorXd const&, MatrixXd const&, MatrixXd const&, MatrixXd const&) throw (CGWASException);

			void test_associations();
			void permutations();
			void permutations(uint const&);
			//VectorXd predict(MatrixXd const&, MatrixXd const&) throw (CGWASException);
			//VectorXd predict(MatrixXd const&, MatrixXd const&, MatrixXd const&) throw (CGWASException);
			void setPhenotype(VectorXd const&) throw (CGWASException);
			void setGenotype(MatrixXd const&) throw (CGWASException);
			void setCovariates(MatrixXd const&) throw (CGWASException);
			void setK(MatrixXd const&) throw (CGWASException);
			void setIntercept(bool const&);
			void setREML(bool const&);
			void setBrent(bool const&);

            VectorXd computeVarianceExplainedNullModel(uint const&);
			float64 getHeritabilityEstimate();
			float64 getGeneticVariance();
            float64 getNoiseVariance();
			float64 getLogLikelihoodNullModel();
			VectorXd getLogLikelihoodAlternativeModels();
			
			float64 getLogDelta0();
			VectorXd getPValues();
			VectorXd getPermutationPValues();
			MatrixXd getBetas();
			MatrixXd getSEBetas();
			VectorXd getTestStatistics();
			float64 getAIC();
			float64 getBIC();
			float64 getAICc();
			
			GWASResults getResults();
	};

}; //namespace CSingleTraitGWAS

#endif //C_SINGLETRAITGWAS_CLASS
