#ifndef CSCONES_CLASS
#define CSCONES_CLASS

#include "CEasyGWAS/globals.h"

#include "CEasyGWAS/regression/CRegression.h"

#define ROBUSTNESS 0
#define CONSISTENCY 1

#define SKAT 0

/*
 *CSconesException Class
 */
class CSconesException {
	private:
		std::string __error_msg;
	public:
		CSconesException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CScones Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

typedef class CSconesSettings {
	public:
		uint folds; //Number of folds for cross-validation
		float64 seed; //Seed for cross-validation 
		bool dump_intermediate_results; //results can easily get very large in memory space, if true results are dumped
		std::string dump_path; //path to folder for dumping results 
		uint selection_criterion; //Selection criterion
		float64 selection_ratio; //ONLY for ROBUSTNESS selection: ratio of hits e.g. in more than 0.8 of all folds 
		uint test_statistic; //Test statistic for computation: SKAT
		bool autoParameters; //set parameters (lambda and eta) values automatically
		//if autoParameters == True
		uint nParameters; //number of parameters 
		VectorXd lambdas; //vector of lambda values
		VectorXd etas; //vector of eta values
		bool evaluateObjective; //flag if objective function should be evaluated

		CSconesSettings();
} scones_settings;

class CScones {
	private:
		CSconesSettings __settings; //Scones settings struct
		MatrixXd __X; //Genotype Matrix
		VectorXd __y; //Phenotype Vector
		SparseMatrixXd __L; //Sparse Network Adjacency Matrix
		DiagXd __W; //Diagonal Matrix with weights for SCAT statistic
		MatrixXd __covs; //Covariate MatrixXd

		//Store results in datastructure
		std::vector<std::vector<std::vector<SparseMatrixXd> > > __result_stack;
		//std::vector<std::vector<MatrixXd>> results_stack;
		//std::vector<MatrixXd> lambda_stack;

		bool __covs_set; //flag if covariates are set 
		uint64 __n_samples; //number of samples
		uint64 __n_features; //number of features
		float64 __best_c; //best CONSISTENCY
		MatrixXd __cMat; //matrix with all consistency values for eta x lambda
		VectorXd __indicator_vector; //indicator vector
		float64 __objective_score; //objective value 
		float64 __best_eta;
		float64 __best_lambda;

        bool __binary_y;
		CLinearRegression __linear_regression; //Regression model, either LinearRegression for continuous phenotypes or LogisticRegression for binary
		CLogisticRegression __logistic_regression; //Regression model, either LinearRegression for continuous phenotypes or LogisticRegression for binary

		void __checkdata() throw (CSconesException);
		void __autoParameters();
		void __selectRegressionModel();
		void __optimize_objective(VectorXd const&, float64 const&, VectorXd*, float64*);
		void __maxflow(SparseMatrixXd const&, MatrixXd const&, VectorXd*);
		void __gridsearch(VectorXd const&, MatrixXd const&, MatrixXd const&) throw (CSconesException);
	
		VectorXd __computeScoreStatistic(MatrixXd const&, VectorXd const&);
		VectorXd __computeSKATScore(MatrixXd const&, VectorXd const&);
	public:
		CScones();
		CScones(VectorXd const&, MatrixXd const&, SparseMatrixXd const&) throw (CSconesException);
		CScones(VectorXd const&, MatrixXd const&, SparseMatrixXd const&, MatrixXd const&) throw (CSconesException);
		
		CScones(CSconesSettings const&);
		CScones(VectorXd const&, MatrixXd const&, SparseMatrixXd const&, CSconesSettings const&) throw (CSconesException);
		CScones(VectorXd const&, MatrixXd const&, SparseMatrixXd const&, MatrixXd const&, CSconesSettings const&) throw (CSconesException);
	
		void test_associations() throw (CSconesException);
		void test_associations(float64 const&, float64 const&);

		//Setter and Getter
		void setSKATWeights(VectorXd const&);
		CSconesSettings getSettings();
        VectorXd getIndicatorVector();
		float64 getObjectiveScore();
		float64 getBestLambda();
		float64 getBestEta();
		MatrixXd getCMatrix(); //Matrix with all consistency/stability values for all etas x lambdas 
		std::vector<std::vector<std::vector<SparseMatrixXd> > > getResultStack(); //get sparse output of all indicator vectors evaluated in cross-validation and gridsearch
};

#endif //CSCONES_CLASS
