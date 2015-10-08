#include "CScones.h"
#include "math.h"
#include "CEasyGWAS/utils/CCrossValidation.h"
#include "CEasyGWAS/utils/CMatrixHelper.h"
#include "CEasyGWAS/utils/StringHelper.h"
#include "maxflow/maxflow.h"
/*
*Initialize default scones settings
*/
CSconesSettings::CSconesSettings() {
	folds = 10;
	seed = 0;
	selection_criterion = CONSISTENCY;
	selection_ratio = 0.8;
	test_statistic = SKAT;
	autoParameters = true;
	nParameters = 10;
	lambdas = VectorXd::Zero(nParameters);
	etas = VectorXd::Zero(nParameters);
	evaluateObjective = false;
	dump_intermediate_results = true;
	dump_path = "tmp/";
}

/*
* Default constructor for Scones
*/
CScones::CScones() {
	CSconesSettings settings;
	__settings = settings;
	__covs_set = false;
    __binary_y = false;
}

CScones::CScones(CSconesSettings const& settings) {
	__settings = settings;
	__covs_set = false;
    __binary_y = false;
}

CScones::CScones(VectorXd const& y, MatrixXd const& X, SparseMatrixXd const& L) throw (CSconesException) {
	CSconesSettings settings;
	__settings = settings;
	__covs_set = false;
    __binary_y = false;
	
	__y = y;
	__X = X;
	__L = L;

	__checkdata();

	__n_samples = X.rows();
	__n_features = X.cols();

	__W = DiagXd(__n_features);
	__W.diagonal() = VectorXd::Ones(__n_features);
}

CScones::CScones(VectorXd const& y, MatrixXd const& X, SparseMatrixXd const& L, CSconesSettings const& settings) throw (CSconesException) {
	__settings = settings;
	__covs_set = false;
    __binary_y = false;
	
	__y = y;
	__X = X;
	__L = L;

	__checkdata();

	__n_samples = X.rows();
	__n_features = X.cols();

	__W = DiagXd(__n_features);
	__W.diagonal() = VectorXd::Ones(__n_features);
}

CScones::CScones(VectorXd const& y, MatrixXd const& X, SparseMatrixXd const& L, MatrixXd const& covs) throw (CSconesException) {
	CSconesSettings settings;
	__settings = settings;
	__covs_set = true;
    __binary_y = false;
	
	__covs = covs;
	__y = y;
	__X = X;
	__L = L;
	
	__checkdata();

	__n_samples = X.rows();
	__n_features = X.cols();

	__W = DiagXd(__n_features);
	__W.diagonal() = VectorXd::Ones(__n_features);
}

CScones::CScones(VectorXd const& y, MatrixXd const& X, SparseMatrixXd const& L, MatrixXd const& covs, CSconesSettings const& settings) throw (CSconesException) {
	__settings = settings;
	__covs_set = true;
    __binary_y = false;
	
	__covs = covs;
	__y = y;
	__X = X;
	__L = L;
	
	__checkdata();

	__n_samples = X.rows();
	__n_features = X.cols();

	__W = DiagXd(__n_features);
	__W.diagonal() = VectorXd::Ones(__n_features);
}

void CScones::__checkdata() throw (CSconesException) {
	if(__y.cols()>1) throw CSconesException("Phenotype y has wrong dimensions! (n x 1)!");
	if(__X.rows() != __y.rows()) throw CSconesException("Genotype X and Phenotype y must have the same number of samples n!");
	if(__L.rows() != __L.cols()) throw CSconesException("Network L must be a squared matrix with size features m x m!");
	//if(__L.rows() != __X.cols()) throw CSconesException("Network L and Genotype X must have the same number of features m!");
	if(__covs_set==true) {
		if (__covs.rows()!=__X.rows()) throw CSconesException("Covariates and genotype must have the same number of samples n!");
	}
}

/*
 *Getter and Setter Methods
 */

//Set weight vector for SKAT association test
void CScones::setSKATWeights(VectorXd const& w) {
	if(w.rows()!=__X.cols()) throw CSconesException("Weight vector w has wrong dimensions");
	__W = DiagXd(__n_features);
	__W.diagonal() = VectorXd::Ones(__n_features).cwiseProduct(w);
}

VectorXd CScones::getIndicatorVector() {
	return __indicator_vector;
}

float64 CScones::getBestLambda() {
	return __best_lambda;
}

float64 CScones::getBestEta() {
	return __best_eta;
}

MatrixXd CScones::getCMatrix() {
	return __cMat;
}
CSconesSettings CScones::getSettings() {
    return __settings;
}

std::vector<std::vector<std::vector<SparseMatrixXd> > > CScones::getResultStack() {
	return __result_stack;
}

/*
* Main algorithm
*/
void CScones::__selectRegressionModel() {
	bool binary = true;
	for(int64 i=0; i<__y.rows();i++) {
		if(!(__y(i)==0 || __y(i)==1)) {
			binary = false;
			break;
		}
	}
	if(binary) { //if phenotype is binary select LogisticRegression
        __logistic_regression = CLogisticRegression();
        __binary_y = true;
	} else { //if phenotype is continuous select LinearRegression
		__linear_regression = CLinearRegression();
        __binary_y = false;
	}
}

VectorXd CScones::__computeSKATScore(MatrixXd const& X, VectorXd const& r) {
	VectorXd nonweighted_skat = ((X.transpose()*r).array().pow(2));
	return __W*nonweighted_skat;
}

VectorXd CScones::__computeScoreStatistic(MatrixXd const& X, VectorXd const& r) {
	VectorXd score;
	if(__settings.test_statistic==SKAT) {
		score = __computeSKATScore(X,r);
	}
	return score;
}

void CScones::__autoParameters() {
	VectorXd r;
	//fit linear regression without evaluating the likelihood function
	if(__covs_set) {
		if(__binary_y==true) {
            __logistic_regression.fit(__y,__covs);
            r = __logistic_regression.getResiduals();
        } else {
            __linear_regression.fit(__y,__covs);
            r = __linear_regression.getResiduals();
	    }
    } else { //if no covariates are set use the zero meaned y ase residuals
		r = __y.array() - __y.mean();
	}
	//set parameters vectors
    VectorXd c = __computeScoreStatistic(__X,r);
	VectorXd tmpc = (c.array()==0).select(c.maxCoeff(),c);
	float64 minc = floor(log10(tmpc.minCoeff()));
	float64 maxc = ceil(log10(c.maxCoeff()));
	__settings.etas = VectorXd::LinSpaced(__settings.nParameters,minc,maxc);
	for(int i=0; i<__settings.etas.rows(); i++)
		__settings.etas(i) = pow(10,__settings.etas(i));
	__settings.lambdas = __settings.etas;
}

//void CScones::__maxflow(SparseMatrixXd const& A, MatrixXd const& T, VectorXd* indicator_vector) {
void CScones::__maxflow(SparseMatrixXd const& A, MatrixXd const& T, VectorXd* indicator_vector) {
	indicator_vector->resize(__n_features);
	//create graph out of ajdencency matrix
	typedef Graph<float64, float64, float64> MaxGraph;
	MaxGraph *g = new MaxGraph(A.rows(),A.nonZeros());
	//Initialize nodes
	g->add_node(indicator_vector->rows());
	//traverse the sparse adjacency matrix A
	for(int64 k=0; k<A.outerSize(); k++) {
		for(Eigen::SparseMatrix<float64>::InnerIterator it(A,k); it; ++it) {
			g->add_edge(it.row(),it.col(),it.value(),0.0);
		}
	}
	//traverse the T matrix
	for(int64 i=0;i<2;i++) {
		for(int64 k=0; k<T.rows(); k++) {
			if(i==0) {
				if(T(k,i)!=0) {
					g->add_tweights(k,T(k,i),0.0);
				}
			} else {
				if(T(k,i)!=0) {
					g->add_tweights(k,0.0,T(k,i));
				}
			}	
		}
	}
	//Run maxflow algorithm
	g->maxflow();

	//create indicator_vector
	for(int64 i=0; i<indicator_vector->rows();i++)
		(*indicator_vector)(i) = g->what_segment(i);
	//delete graph
	delete g;
}

void CScones::__optimize_objective(VectorXd const& c, float64 const& lambda, VectorXd* indicator_vector, float64* objective_score) {
	//Main Graph
	SparseMatrixXd A = lambda * __L;
	//Add source and sink
	//if((c.array()>0).count()==0 || (c.array()<=0).count()==0)
	//	logging(WARNING,"WARNING: Trivial solution ahead!");
	
	MatrixXd T(__n_features,2);
	//connect positive c values to sink
	VectorXd posc = (c.array()<=0).select(0,c);
	//connect nevative c values to source 
	VectorXd negc = -c;
	negc = (negc.array()<=0).select(0,negc);
	//Store data
	T.col(0) = negc;
	T.col(1) = posc;
	
	//compute maxflow
	__maxflow(A,T,indicator_vector);	
	
	//compute objective score
	if(__settings.evaluateObjective) {
		(*objective_score) = 0.0;
		for(uint i=0; i<c.rows(); i++) {
			if((*indicator_vector)(i)==0) {
				(*objective_score) += c(i); 
			}
		}
	
		for(int64 k=0; k<A.outerSize(); k++) {
			for(Eigen::SparseMatrix<float64>::InnerIterator it(A,k); it; ++it) {
				if((*indicator_vector)(it.row())==1)
					if((*indicator_vector)(it.col())==0)
						(*objective_score) += it.value();
			}
		}
		(*objective_score) -= c.sum();
	} else {
		(*objective_score) = NAN;
	}
}

void CScones::__gridsearch(VectorXd const& y, MatrixXd const& x, MatrixXd const& covs) throw (CSconesException) {
	VectorXd r;
	if(__covs_set) {
		if(__binary_y==true) {
            __logistic_regression.fit(y,covs);
            r = __logistic_regression.getResiduals();
        } else {
            __linear_regression.fit(y,covs);
            r = __linear_regression.getResiduals();
	    }
	} else {
		r = y.array() - y.mean();
	}
	VectorXd c = __computeScoreStatistic(x,r);
	VectorXd c_transformed;
	std::vector<std::vector<SparseMatrixXd> > lambda_stack;
	//Perform gridsearch to find optimal lambda and eta values 
	for(uint i=0; i<__settings.lambdas.rows(); i++) {
		std::vector<SparseMatrixXd> eta_stack;
		for(uint j=0;j<__settings.etas.rows(); j++) {
			//Substract eta value from score-statistic c
			c_transformed = c.array()-__settings.etas(j);
			//Optimize objective
			VectorXd indicator_vector;
			float64 objective_score;
			__optimize_objective(c_transformed, __settings.lambdas(i), &indicator_vector, &objective_score);
			eta_stack.push_back(indicator_vector.sparseView());
		}
		lambda_stack.push_back(eta_stack);
	}
	__result_stack.push_back(lambda_stack);
}

void CScones::test_associations() throw (CSconesException) {
	//Select the correct Regression Model
	__selectRegressionModel();
	if(__settings.autoParameters) {
		__autoParameters();
	} else {
		if(__settings.etas.rows() == 0 || __settings.lambdas.rows() == 0)
			throw CSconesException("Array of lambda or eta values cannot be empty!");
	}
	//Perform crossvalidated gridsearch to find parameters 
	CCrossValidation cv(__settings.seed);
	cv.kFold(__settings.folds,__n_samples);
	
	for(uint k=0;k<__settings.folds;k++) {
		VectorXd tr_indices = cv.getTrainingIndices(k);
		//Slice data according to the indices
		MatrixXd x_train = sliceRowsMatrix(__X,tr_indices);
		VectorXd y_train = sliceRowsMatrix(__y,tr_indices);
		MatrixXd cov_train;
		if(__covs_set) {
			cov_train = sliceRowsMatrix(__covs,tr_indices);
		}
		__gridsearch(y_train,x_train,cov_train);
	}
	//Evaluate all solutions and select the best eta and lambda according to the selection criterion
	if(__settings.selection_criterion==CONSISTENCY) {
		float64 n = __L.outerSize();
        float64 n1 = ceil(n*0.01);
		__cMat = MatrixXd::Zero(__settings.etas.rows(),__settings.lambdas.rows());
		for(int i=0; i<__settings.etas.rows();i++){
			for(int j=0; j<__settings.lambdas.rows();j++) {
				float64 cindex = 0.0;
				for(uint k=0; k<__settings.folds; k++) {
					//inds.push_back(__result_stack[k][j][i]);
					float64 si = __result_stack[k][j][i].nonZeros();
                    if (si==0 || si == n || si > n1) continue;
					for(uint h=k+1; h<__settings.folds; h++) {
						float64 sj = __result_stack[h][j][i].nonZeros();
						if (sj==0 || sj == n || sj > n1) continue;
						cindex += (n*(__result_stack[k][j][i].cwiseProduct(__result_stack[h][j][i])).sum() - si * sj)/(n*fmin(si,sj)-si*sj);
					}
				}
				cindex = 2.0 * cindex / (__settings.folds * (__settings.folds - 1));
				__cMat(i,j) = cindex;
			}
		}
		MatrixXd::Index best_eta_index, best_lambda_index;
		__best_c = __cMat.maxCoeff(&best_eta_index,&best_lambda_index);
		__best_lambda = __settings.lambdas[best_lambda_index];
		__best_eta = __settings.etas[best_eta_index];
		//Selected Features for the best eta and lambda are those selected in all folds
        __indicator_vector = VectorXd::Zero(__n_features);
        if(__best_lambda>0.0 && __best_eta>0) { 
            for(uint k=0; k<__settings.folds; k++) {
                __indicator_vector += __result_stack[k][best_lambda_index][best_eta_index];
            }
            //Select only those features selected in all folds
            __indicator_vector /= __settings.folds;
            for(uint64 i=0; i<__indicator_vector.rows();i++) {
                if(__indicator_vector(i)!=1) __indicator_vector(i)=0;
            }
        }
        //ALTERNATIVE TO GET THE FINAL INDICATOR VECTOR, retrain model using the optimized paramters: HOWEVER this sometimes leads to different selections.
        //test_associations(__best_lambda,__best_eta);
	}
}

void CScones::test_associations(float64 const& lambda, float64 const& eta) {
	//Select the correct Regression Model
	__selectRegressionModel();
	//Test associations for a fixed lambda and eta
	VectorXd r;
	if(__covs_set) {
		if(__binary_y==true) {
            __logistic_regression.fit(__y,__covs);
            r = __logistic_regression.getResiduals();
        } else {
            __linear_regression.fit(__y,__covs);
            r = __linear_regression.getResiduals();
	    }
	} else {
		r = __y.array() - __y.mean();
	}
	VectorXd c = __computeScoreStatistic(__X,r).array() - eta;
	//solve objective function
	__optimize_objective(c, lambda, &__indicator_vector, &__objective_score);
}
