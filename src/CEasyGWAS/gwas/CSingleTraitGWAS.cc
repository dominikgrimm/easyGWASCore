#include "CSingleTraitGWAS.h"
#include "CEasyGWAS/stats/CChi2.h"
#include "CEasyGWAS/stats/CFisherF.h"
#include "CEasyGWAS/utils/CMatrixHelper.h"
#include "CEasyGWAS/utils/CCrossValidation.h"
#include "CEasyGWAS/stats/CStats.h"

#include <stdlib.h>

namespace CSingleTraitGWAS {

LinearRegression::LinearRegression() {
	__intercept = true;
	__covs_set = false;
	__statistical_test = LRTEST;
	__n_permutations = 1000;
}

LinearRegression::LinearRegression(VectorXd const& y, MatrixXd const& X) throw (CGWASException) {
	__intercept = true;
	__covs_set = false;
	__y = y;
	__X = X;
	__statistical_test = LRTEST;
	__n_permutations = 1000;
	__checkdata();
}

LinearRegression::LinearRegression(VectorXd const& y, MatrixXd const& X, MatrixXd const& covs) throw (CGWASException) {
	__intercept = true;
	__covs_set = true;
	__y = y;
	__X = X;
	__covs = covs;
	__statistical_test = LRTEST;
	__n_permutations = 1000;
	__checkdata();
}

void LinearRegression::__checkdata() throw (CGWASException) {
	if(__y.cols()>1) throw CGWASException("Phenotype y has wrong dimensions! (n x 1)!");
	if(__X.rows() != __y.rows()) throw CGWASException("Genotype X and Phenotype y must have the same number of samples n!");
	if(__covs_set==true) {
		if (__covs.rows()!=__X.rows()) throw CGWASException("Covariates and genotype must have the same number of samples n!");
	}
}

void LinearRegression::setTestStatistic(uint const& statistical_test) throw (CGWASException) {
	if(!(statistical_test == LRTEST or statistical_test== FTEST)) 
		throw CGWASException("Statistical test must be one of the following ones [LRTest, FTEST]");
	__statistical_test = statistical_test;
}

void LinearRegression::setPhenotype(VectorXd const& y) throw (CGWASException) {
	__y = y;
	if(__y.cols()>1) throw CGWASException("Phenotype y has wrong dimensions! (n x 1)!");
}


void LinearRegression::setGenotype(MatrixXd const& X) throw (CGWASException) {
	__X = X;
	if(__X.rows() != __y.rows()) throw CGWASException("Genotype X and Phenotype y must have the same number of samples n!");
}

void LinearRegression::setCovariates(MatrixXd const& covs) throw (CGWASException) {
	__covs = covs;
	if (__covs.rows()!=__X.rows()) throw CGWASException("Covariates and genotype must have the same number of samples n!");
	__covs_set = true;
}

void LinearRegression::setIntercept(bool const& intercept) {
	__intercept = intercept;
}

float64 LinearRegression::getLogLikelihoodNullModel() {
	return __null_model.getLogLikelihood();
}

VectorXd LinearRegression::getLogLikelihoodAlternativeModels() {
	return __logLikelihoods;
}

VectorXd LinearRegression::getPValues() {
	return __p_values;
}

MatrixXd LinearRegression::getBetas() {
	return __betas;
}

MatrixXd LinearRegression::getSEBetas() {
	return __se_betas;
}

VectorXd LinearRegression::getTestStatistics() {
	return __test_statistics;
}

GWASResults LinearRegression::getResults() {
	GWASResults results;
	results.p_values = __p_values;
	results.betas = __betas;
	results.se_betas = __se_betas;
	results.test_statistics = __test_statistics;
	results.alternative_loglikelihoods = __logLikelihoods;
	results.null_loglikelihood = __null_model.getLogLikelihood();
	return results;
}

float64 LinearRegression::getAIC() {
	return __null_model.getAIC();
}

float64 LinearRegression::getBIC() {
	return __null_model.getBIC();
}

float64 LinearRegression::getAICc() {
	return __null_model.getAICc();
}

VectorXd LinearRegression::getPermutationPValues() {
	return __p_values_permutations;
}

void LinearRegression::permutations() {
	permutations(1000);
}

void LinearRegression::permutations(uint const& n_perm) {
	__n_permutations = n_perm;
	VectorXd original_y = __y;
	__p_values_permutations = VectorXd::Zero(__X.cols());
	__permutation_pool = VectorXd::Zero(__X.cols()*__n_permutations);
	//permute phenotypes
	MatrixXd Y = permuteVector(__y,__n_permutations-1);
	for(uint i=0;i<__n_permutations-1;i++) {
		__y = Y.col(i);
		test_associations();
		__permutation_pool.block(__X.cols()*i,0,__X.cols(),1) = __test_statistics;
    }
	__y = original_y;
	test_associations();
    __permutation_pool.block(__X.cols()*(__n_permutations-1),0,__X.cols(),1) = __test_statistics;
	
    //sort permutated p_values
	std::sort(__permutation_pool.data(),__permutation_pool.data()+__permutation_pool.rows(),std::greater<float64>());
	
    //Compute Permutation P-Values
    for (int64 i=0;i<__test_statistics.rows();i++) {
        uint64 n_elements = 0;
        for(int64 j=0; j<__permutation_pool.rows();j++) {
            if(__test_statistics(i)<=__permutation_pool(j)) {
                n_elements++;
            } else {
                break;
            }
        }
        float64 p_val=(float64)n_elements/((float64)__permutation_pool.rows());
        __p_values_permutations(i) = p_val;
    }
}

void LinearRegression::test_associations() {
	//Init data strucutres
	__p_values = VectorXd::Zero(__X.cols());
	__test_statistics = VectorXd::Zero(__X.cols());
	__logLikelihoods = VectorXd::Zero(__X.cols());
	if (__intercept)
		__n_features = 2;
	else
		__n_features = 1;
	if (__covs_set)
		__n_features += __covs.cols();
	__betas = MatrixXd::Zero(__X.cols(),__n_features);
	__se_betas = MatrixXd::Zero(__X.cols(),__n_features);
	//Compute NUll Model
	__null_model = CLinearRegression(false);
	if (__covs_set) {
		MatrixXd tmp(__y.rows(),1+__covs.cols());
	        tmp << VectorXd::Ones(__y.rows()), __covs;
		__null_model.fit(__y,tmp);
	} else {
		__null_model.fit(__y,VectorXd::Ones(__y.rows()));
	}
	//Compute Alternative Models
	CLinearRegression a_model;
	for(uint i=0; i < __X.cols(); i++) {
		if (__covs_set) {
			MatrixXd tmp(__y.rows(),1+__covs.cols());
			tmp << __X.col(i), __covs;
			a_model.fit(__y,tmp);
		} else {
			a_model.fit(__y,__X.col(i));
		}
		//store estimated data
		__betas.block(i,0,1,__n_features) = a_model.getBetas().transpose();
		__se_betas.block(i,0,1,__n_features) = a_model.getStdBetas().transpose();
		__logLikelihoods(i) = a_model.getLogLikelihood();
		
		//perform a statistical_test
		if(__statistical_test==LRTEST) {
			__test_statistics(i) = (a_model.getLogLikelihood()-__null_model.getLogLikelihood());
			__p_values(i) = CChi2::sf(2.0*__test_statistics(i),1);
		} else if(__statistical_test==FTEST) {
			__test_statistics(i) = (__X.rows()-__n_features) * a_model.getRSquared()/(1.0-a_model.getRSquared());
			__p_values(i) = CFisherF::sf(__test_statistics(i),1,__X.rows()-__n_features);
		}
	}
}

/*
*Logistic Regression
*/
LogisticRegression::LogisticRegression() {
	__intercept = true;
	__covs_set = false;
}

LogisticRegression::LogisticRegression(VectorXd const& y, MatrixXd const& X) throw (CGWASException) {
	__intercept = true;
	__covs_set = false;
	__y = y;
	__X = X;
	__checkdata();
}

LogisticRegression::LogisticRegression(VectorXd const& y, MatrixXd const& X, MatrixXd const& covs) throw (CGWASException) {
	__intercept = true;
	__covs_set = true;
	__y = y;
	__X = X;
	__covs = covs;
	__checkdata();
}

void LogisticRegression::__checkdata() throw (CGWASException) {
	if(__y.cols()>1) throw CGWASException("Phenotype y has wrong dimensions! (n x 1)!");
	if(__X.rows() != __y.rows()) throw CGWASException("Genotype X and Phenotype y must have the same number of samples n!");
	if(__covs_set==true) {
		if (__covs.rows()!=__X.rows()) throw CGWASException("Covariates and genotype must have the same number of samples n!");
	}
}

void LogisticRegression::setPhenotype(VectorXd const& y) throw (CGWASException) {
	__y = y;
	if(__y.cols()>1) throw CGWASException("Phenotype y has wrong dimensions! (n x 1)!");
}


void LogisticRegression::setGenotype(MatrixXd const& X) throw (CGWASException) {
	__X = X;
	if(__X.rows() != __y.rows()) throw CGWASException("Genotype X and Phenotype y must have the same number of samples n!");
}

void LogisticRegression::setCovariates(MatrixXd const& covs) throw (CGWASException) {
	__covs = covs;
	if (__covs.rows()!=__X.rows()) throw CGWASException("Covariates and genotype must have the same number of samples n!");
	__covs_set = true;
}

void LogisticRegression::setIntercept(bool const& intercept) {
	__intercept = intercept;
}

float64 LogisticRegression::getLogLikelihoodNullModel() {
	return __null_model.getLogLikelihood();
}

VectorXd LogisticRegression::getLogLikelihoodAlternativeModels() {
	return __logLikelihoods;
}

VectorXd LogisticRegression::getPValues() {
	return __p_values;
}

MatrixXd LogisticRegression::getBetas() {
	return __betas;
}

MatrixXd LogisticRegression::getSEBetas() {
	return __se_betas;
}

VectorXd LogisticRegression::getTestStatistics() {
	return __test_statistics;
}

float64 LogisticRegression::getAIC() {
	return __null_model.getAIC();
}

float64 LogisticRegression::getBIC() {
	return __null_model.getBIC();
}

float64 LogisticRegression::getAICc() {
	return __null_model.getAICc();
}

GWASResults LogisticRegression::getResults() {
	GWASResults results;
	results.p_values = __p_values;
	results.betas = __betas;
	results.se_betas = __se_betas;
	results.test_statistics = __test_statistics;
	results.alternative_loglikelihoods = __logLikelihoods;
	results.null_loglikelihood = __null_model.getLogLikelihood();
	return results;
}

VectorXd LogisticRegression::getPermutationPValues() {
	return __p_values_permutations;
}

void LogisticRegression::permutations() {
	permutations(1000);
}

void LogisticRegression::permutations(uint const& n_perm) {
	__n_permutations = n_perm;
	VectorXd original_y = __y;
	__p_values_permutations = VectorXd::Zero(__X.cols());
	__permutation_pool = VectorXd::Zero(__X.cols()*__n_permutations);
	//permute phenotypes
	MatrixXd Y = permuteVector(__y,__n_permutations-1);
	for(uint i=0;i<__n_permutations-1;i++) {
		__y = Y.col(i);
		test_associations();
		__permutation_pool.block(__X.cols()*i,0,__X.cols(),1) = __test_statistics;
    }
	__y = original_y;
	test_associations();
    __permutation_pool.block(__X.cols()*(__n_permutations-1),0,__X.cols(),1) = __test_statistics;
	
    //sort permutated p_values
	std::sort(__permutation_pool.data(),__permutation_pool.data()+__permutation_pool.rows(),std::greater<float64>());
	
    //Compute Permutation P-Values
    for (int64 i=0;i<__test_statistics.rows();i++) {
        uint64 n_elements = 0;
        for(int64 j=0; j<__permutation_pool.rows();j++) {
            if(__test_statistics(i)<=__permutation_pool(j)) {
                n_elements++;
            } else {
                break;
            }
        }
        float64 p_val=(float64)n_elements/((float64)__permutation_pool.rows());
        __p_values_permutations(i) = p_val;
    }
}

void LogisticRegression::test_associations() {
	//Init data strucutres
	__p_values = VectorXd::Zero(__X.cols());
	__test_statistics = VectorXd::Zero(__X.cols());
	__logLikelihoods = VectorXd::Zero(__X.cols());
	if (__intercept)
		__n_features = 2;
	else
		__n_features = 1;
	if (__covs_set)
		__n_features += __covs.cols();
	__betas = MatrixXd::Zero(__X.cols(),__n_features);
	__se_betas = MatrixXd::Zero(__X.cols(),__n_features);
	//Compute NUll Model
	__null_model = CLogisticRegression(false);
	if (__covs_set) {
		MatrixXd tmp(__y.rows(),1+__covs.cols());
	        tmp << VectorXd::Ones(__y.rows()), __covs;
		__null_model.fit(__y,tmp);
	} else {
		__null_model.fit(__y,VectorXd::Ones(__y.rows()));
	}
	//Compute Alternative Models
	CLogisticRegression a_model;
	for(uint i=0; i < __X.cols(); i++) {
		if (__covs_set) {
			MatrixXd tmp(__y.rows(),1+__covs.cols());
			tmp << __X.col(i), __covs;
			a_model.fit(__y,tmp);
		} else {
			a_model.fit(__y,__X.col(i));
		}
		//store estimated data
		__betas.block(i,0,1,__n_features) = a_model.getBetas().transpose();
		__se_betas.block(i,0,1,__n_features) = a_model.getStdBetas().transpose();
		__logLikelihoods(i) = a_model.getLogLikelihood();
		
		//perform a statistical_test
		__test_statistics(i) = (a_model.getLogLikelihood()-__null_model.getLogLikelihood());
		__p_values(i) = CChi2::sf(2.0*__test_statistics(i),1);
	}
}

/*
* EMMAX Model
*/
EMMAX::EMMAX() {
	__intercept = true;
	__covs_set = false;
	__REML = false;
	__use_brent = true;
}

EMMAX::EMMAX(VectorXd const& y, MatrixXd const& X, MatrixXd const& K) throw (CGWASException) {
	__intercept = true;
	__covs_set = false;
	__REML = false;
	__use_brent = true;
	__y = y;
	__X = X;
	__K = K;
	__checkdata();
}

EMMAX::EMMAX(VectorXd const& y, MatrixXd const& X, MatrixXd const& K, MatrixXd const& covs) throw (CGWASException) {
	__intercept = true;
	__covs_set = true;
	__REML = false;
	__use_brent = true;
	__y = y;
	__X = X;
	__K = K;
	__covs = covs;
	__checkdata();
}

void EMMAX::__checkdata() throw (CGWASException) {
	if(__y.cols()>1) throw CGWASException("Phenotype y has wrong dimensions! (n x 1)!");
	if(__X.rows() != __y.rows()) throw CGWASException("Genotype X and Phenotype y must have the same number of samples n!");
	if (__K.rows()!=__X.rows()) throw CGWASException("Kernel and genotype must have the same number of samples n!");
	if (__K.rows()!=__K.cols()) throw CGWASException("Kernel K must be a samples n x samples n matrix!");
	if(__covs_set==true) {
		if (__covs.rows()!=__X.rows()) throw CGWASException("Covariates and genotype must have the same number of samples n!");
	}
}

void EMMAX::setPhenotype(VectorXd const& y) throw (CGWASException) {
	__y = y;
	if(__y.cols()>1) throw CGWASException("Phenotype y has wrong dimensions! (n x 1)!");
}


void EMMAX::setGenotype(MatrixXd const& X) throw (CGWASException) {
	__X = X;
	if(__X.rows() != __y.rows()) throw CGWASException("Genotype X and Phenotype y must have the same number of samples n!");
}

void EMMAX::setCovariates(MatrixXd const& covs) throw (CGWASException) {
	__covs = covs;
	if (__covs.rows()!=__X.rows()) throw CGWASException("Covariates and genotype must have the same number of samples n!");
	__covs_set = true;
}

void EMMAX::setK(MatrixXd const& K) throw (CGWASException) {
	__K = K;
	//if (__K.rows()!=__X.rows()) throw CGWASException("Kernel and genotype must have the same number of samples n!");
	//if (__K.rows()!=__K.cols()) throw CGWASException("Kernel K must be a samples n x samples n matrix!");
}

void EMMAX::setIntercept(bool const& intercept) {
	__intercept = intercept;
}

void EMMAX::setREML(bool const& reml) {
	__REML = reml;
}

void EMMAX::setBrent(bool const& brent) {
	__use_brent = brent;
}

float64 EMMAX::getLogLikelihoodNullModel() {
	return __null_model.getLogLikelihood();
}

VectorXd EMMAX::getLogLikelihoodAlternativeModels() {
	return __logLikelihoods;
}

VectorXd EMMAX::getPValues() {
	return __p_values;
}

MatrixXd EMMAX::getBetas() {
	return __betas;
}

MatrixXd EMMAX::getSEBetas() {
	return __se_betas;
}

float64 EMMAX::getLogDelta0() {
	return __null_model.getLogDelta();
}

VectorXd EMMAX::getTestStatistics() {
	return __test_statistics;
}

float64 EMMAX::getAIC() {
	return __null_model.getAIC();
}

float64 EMMAX::getBIC() {
	return __null_model.getBIC();
}

float64 EMMAX::getAICc() {
	return __null_model.getAICc();
}

GWASResults EMMAX::getResults() {
	GWASResults results;
	results.p_values = __p_values;
	results.betas = __betas;
	results.se_betas = __se_betas;
	results.test_statistics = __test_statistics;
	results.alternative_loglikelihoods = __logLikelihoods;
	results.null_loglikelihood = __null_model.getLogLikelihood();
	return results;
}

VectorXd EMMAX::getPermutationPValues() {
	return __p_values_permutations;
}

void EMMAX::permutations() {
	permutations(1000);
}

void EMMAX::permutations(uint const& n_perm) {
	__n_permutations = n_perm;
	VectorXd original_y = __y;
	__p_values_permutations = VectorXd::Zero(__X.cols());
	__permutation_pool = VectorXd::Zero(__X.cols()*__n_permutations);
	//permute phenotypes
	MatrixXd Y = permuteVector(__y,__n_permutations-1);
	for(uint i=0;i<__n_permutations-1;i++) {
		__y = Y.col(i);
		test_associations();
		__permutation_pool.block(__X.cols()*i,0,__X.cols(),1) = __test_statistics;
    }
	__y = original_y;
	test_associations();
    __permutation_pool.block(__X.cols()*(__n_permutations-1),0,__X.cols(),1) = __test_statistics;
	
    //sort permutated p_values
	std::sort(__permutation_pool.data(),__permutation_pool.data()+__permutation_pool.rows(),std::greater<float64>());
	
    //Compute Permutation P-Values
    for (int64 i=0;i<__test_statistics.rows();i++) {
        uint64 n_elements = 0;
        for(int64 j=0; j<__permutation_pool.rows();j++) {
            if(__test_statistics(i)<=__permutation_pool(j)) {
                n_elements++;
            } else {
                break;
            }
        }
        float64 p_val=(float64)n_elements/((float64)__permutation_pool.rows());
        __p_values_permutations(i) = p_val;
    }
}

float64 EMMAX::getHeritabilityEstimate() {
    float64 sigma_g = exp(__null_model.getLogSigma());
    float64 sigma_e = exp(__null_model.getLogDelta())*sigma_g;
    return sigma_g/(sigma_e+sigma_g);
}

float64 EMMAX::getGeneticVariance() {
    return exp(__null_model.getLogSigma());
}

float64 EMMAX::getNoiseVariance() {
    float64 sigma_g = exp(__null_model.getLogSigma());
    return exp(__null_model.getLogDelta())*sigma_g;
}

float64 EMMAX::computeVarianceExplainedNullModel(uint const& folds) {
	MatrixXd covariates = __covs;
    if (__covs_set) {
        covariates = MatrixXd(__y.rows(),1+__covs.cols());
        covariates << VectorXd::Ones(__y.rows()), __covs;
    } else {
		covariates = VectorXd::Ones(__y.rows());
    }

    CCrossValidation cv(rand()%1000);
    cv.kFold(folds,__y.rows());

    VectorXd y_estimated = VectorXd::Zero(__y.rows());
    
    for(uint k=0; k<folds;k++) {
        VectorXd tr_indices = cv.getTrainingIndices(k);
        VectorXd te_indices = cv.getTestingIndices(k);
        VectorXd y_train = sliceRowsMatrix(__y,tr_indices);
        VectorXd y_test = sliceRowsMatrix(__y,te_indices);
        MatrixXd cov_train = sliceRowsMatrix(covariates,tr_indices);
        MatrixXd cov_test = sliceRowsMatrix(covariates,te_indices);
        MatrixXd K_train = sliceRowsMatrix(__K,tr_indices);
        K_train = sliceColsMatrix(K_train,tr_indices);
        MatrixXd K_test = sliceRowsMatrix(__K,te_indices);
        K_test = sliceColsMatrix(K_test,tr_indices);
	
        CLinearMixedRegression null_model(false);
	    null_model.setInterval(100);
	    null_model.setLogDeltaMin(-5.0);
	    null_model.setLogDeltaMax(5.0);
	    null_model.setREML(__REML);
    	null_model.setBrent(__use_brent);
        null_model.fit(y_train,cov_train,K_train);

        VectorXd y_tmp;
        null_model.predict(&y_tmp,cov_test,K_test);
        insertColumnVectorAtIndices(&y_estimated,y_tmp,te_indices);       
    }
    float64 variance_explained = 1.0 - CStats::varf(__y.array() - y_estimated.array())/CStats::varf(__y);
    return variance_explained;
}

void EMMAX::test_associations() {
	//Init data strucutres
	__p_values = VectorXd::Zero(__X.cols());
	__test_statistics = VectorXd::Zero(__X.cols());
	__logLikelihoods = VectorXd::Zero(__X.cols());
	if (__intercept)
		__n_features = 2;
	else
		__n_features = 1;
	Eigen::SelfAdjointEigenSolver<MatrixXd> solver(__K);
	MatrixXd S = solver.eigenvalues();
	MatrixXd U = solver.eigenvectors();
	MatrixXd Uy = U.transpose()*__y;
	MatrixXd Ucovs;
	MatrixXd Ux = U.transpose()*__X;
	if (__covs_set) {
        __n_features += __covs.cols();
		MatrixXd tmp(__y.rows(),1+__covs.cols());
	        tmp << VectorXd::Ones(__y.rows()), __covs;
		Ucovs = U.transpose()*tmp;
	} else {
		Ucovs = U.transpose()*VectorXd::Ones(__y.rows());
	}
	__betas = MatrixXd::Zero(__X.cols(),__n_features);
	__se_betas = MatrixXd::Zero(__X.cols(),__n_features);
	//Compute NUll Model
	
	__null_model = CLinearMixedRegression(Ucovs,Uy,S,false);
	__null_model.setInterval(100);
	__null_model.setLogDeltaMin(-5.0);
	__null_model.setLogDeltaMax(5.0);
	__null_model.setREML(__REML);
	__null_model.setBrent(__use_brent);
	if (__covs_set) {
		MatrixXd tmp(__y.rows(),1+__covs.cols());
	        tmp << VectorXd::Ones(__y.rows()), __covs;
		__null_model.fit(__y,tmp,__K);
	} else {
		__null_model.fit(__y,VectorXd::Ones(__y.rows()),__K);
	}

	for(uint i=0; i < __X.cols(); i++) {
		//Compute Alternative Models
		CLinearMixedRegression a_model;
		if (__covs_set) {
			MatrixXd Utmp(__y.rows(),2+__covs.cols());
			Utmp << VectorXd::Ones(__y.rows()), __X.col(i), __covs;
			Utmp = U.transpose()*Utmp;
			a_model = CLinearMixedRegression(Utmp,Uy,S,false);
			a_model.setInterval(0);
			a_model.setLogDeltaMin(-5.0);
			a_model.setLogDeltaMax(5.0);
			a_model.setREML(__REML);
			a_model.setBrent(__use_brent);
			a_model.setLogDelta(__null_model.getLogDelta());
			MatrixXd tmp(__y.rows(),1+__covs.cols());
			tmp << __X.col(i), __covs;
			a_model.fit(__y,tmp,__K);
		} else {
			MatrixXd Utmp(__y.rows(),2);
			//Utmp << VectorXd::Ones(__y.rows()), Ux.col(i);
			Utmp << VectorXd::Ones(__y.rows()), __X.col(i);
			Utmp = U.transpose()*Utmp;
			a_model = CLinearMixedRegression(Utmp,Uy,S,false);
			a_model.setInterval(0);
			a_model.setLogDeltaMin(-5.0);
			a_model.setLogDeltaMax(5.0);
			a_model.setREML(__REML);
			a_model.setBrent(__use_brent);
			a_model.setLogDelta(__null_model.getLogDelta());
			a_model.fit(__y,__X.col(i),__K);
		}
		//store estimated data
		__betas.block(i,0,1,__n_features) = a_model.getBetas().transpose();
		__se_betas.block(i,0,1,__n_features) = a_model.getStdBetas().transpose();
		__logLikelihoods(i) = a_model.getLogLikelihood();
		
		//perform a statistical_test
		__test_statistics(i) = (a_model.getLogLikelihood()-__null_model.getLogLikelihood());
		__p_values(i) = CChi2::sf(2.0*__test_statistics(i),1);
	}
}

/*
* FaSTLMM Model
*/
FaSTLMM::FaSTLMM() {
	__intercept = true;
	__covs_set = false;
	__REML = false;
	__use_brent = true;
}

FaSTLMM::FaSTLMM(VectorXd const& y, MatrixXd const& X, MatrixXd const& K) throw (CGWASException) {
	__intercept = true;
	__covs_set = false;
	__REML = false;
	__use_brent = true;
	__y = y;
	__X = X;
	__K = K;
	__checkdata();
}

FaSTLMM::FaSTLMM(VectorXd const& y, MatrixXd const& X, MatrixXd const& K, MatrixXd const& covs) throw (CGWASException) {
	__intercept = true;
	__covs_set = true;
	__REML = false;
	__use_brent = true;
	__y = y;
	__X = X;
	__K = K;
	__covs = covs;
	__checkdata();
}

void FaSTLMM::__checkdata() throw (CGWASException) {
	if(__y.cols()>1) throw CGWASException("Phenotype y has wrong dimensions! (n x 1)!");
	if(__X.rows() != __y.rows()) throw CGWASException("Genotype X and Phenotype y must have the same number of samples n!");
	if (__K.rows()!=__X.rows()) throw CGWASException("Kernel and genotype must have the same number of samples n!");
	if (__K.rows()!=__K.cols()) throw CGWASException("Kernel K must be a samples n x samples n matrix!");
	if(__covs_set==true) {
		if (__covs.rows()!=__X.rows()) throw CGWASException("Covariates and genotype must have the same number of samples n!");
	}
}

void FaSTLMM::setPhenotype(VectorXd const& y) throw (CGWASException) {
	__y = y;
	if(__y.cols()>1) throw CGWASException("Phenotype y has wrong dimensions! (n x 1)!");
}


void FaSTLMM::setGenotype(MatrixXd const& X) throw (CGWASException) {
	__X = X;
	if(__X.rows() != __y.rows()) throw CGWASException("Genotype X and Phenotype y must have the same number of samples n!");
}

void FaSTLMM::setCovariates(MatrixXd const& covs) throw (CGWASException) {
	__covs = covs;
	if (__covs.rows()!=__X.rows()) throw CGWASException("Covariates and genotype must have the same number of samples n!");
	__covs_set = true;
}

void FaSTLMM::setK(MatrixXd const& K) throw (CGWASException) {
	__K = K;
	//if (__K.rows()!=__X.rows()) throw CGWASException("Kernel and genotype must have the same number of samples n!");
	//if (__K.rows()!=__K.cols()) throw CGWASException("Kernel K must be a samples n x samples n matrix!");
}

void FaSTLMM::setIntercept(bool const& intercept) {
	__intercept = intercept;
}

void FaSTLMM::setREML(bool const& reml) {
	__REML = reml;
}

void FaSTLMM::setBrent(bool const& brent) {
	__use_brent = brent;
}

float64 FaSTLMM::getLogLikelihoodNullModel() {
	return __null_model.getLogLikelihood();
}

VectorXd FaSTLMM::getLogLikelihoodAlternativeModels() {
	return __logLikelihoods;
}

VectorXd FaSTLMM::getPValues() {
	return __p_values;
}

MatrixXd FaSTLMM::getBetas() {
	return __betas;
}

MatrixXd FaSTLMM::getSEBetas() {
	return __se_betas;
}

float64 FaSTLMM::getLogDelta0() {
	return __null_model.getLogDelta();
}

VectorXd FaSTLMM::getTestStatistics() {
	return __test_statistics;
}

GWASResults FaSTLMM::getResults() {
	GWASResults results;
	results.p_values = __p_values;
	results.betas = __betas;
	results.se_betas = __se_betas;
	results.test_statistics = __test_statistics;
	results.alternative_loglikelihoods = __logLikelihoods;
	results.null_loglikelihood = __null_model.getLogLikelihood();
	return results;
}

float64 FaSTLMM::getAIC() {
	return __null_model.getAIC();
}

float64 FaSTLMM::getBIC() {
	return __null_model.getBIC();
}

float64 FaSTLMM::getAICc() {
	return __null_model.getAICc();
}

VectorXd FaSTLMM::getPermutationPValues() {
	return __p_values_permutations;
}

void FaSTLMM::permutations() {
	permutations(1000);
}

void FaSTLMM::permutations(uint const& n_perm) {
	__n_permutations = n_perm;
	VectorXd original_y = __y;
	__p_values_permutations = VectorXd::Zero(__X.cols());
	__permutation_pool = VectorXd::Zero(__X.cols()*__n_permutations);
	//permute phenotypes
	MatrixXd Y = permuteVector(__y,__n_permutations-1);
	for(uint i=0;i<__n_permutations-1;i++) {
		__y = Y.col(i);
		test_associations();
		__permutation_pool.block(__X.cols()*i,0,__X.cols(),1) = __test_statistics;
    }
	__y = original_y;
	test_associations();
    __permutation_pool.block(__X.cols()*(__n_permutations-1),0,__X.cols(),1) = __test_statistics;
	
    //sort permutated p_values
	std::sort(__permutation_pool.data(),__permutation_pool.data()+__permutation_pool.rows(),std::greater<float64>());
	
    //Compute Permutation P-Values
    for (int64 i=0;i<__test_statistics.rows();i++) {
        uint64 n_elements = 0;
        for(int64 j=0; j<__permutation_pool.rows();j++) {
            if(__test_statistics(i)<=__permutation_pool(j)) {
                n_elements++;
            } else {
                break;
            }
        }
        float64 p_val=(float64)n_elements/((float64)__permutation_pool.rows());
        __p_values_permutations(i) = p_val;
    }
}

float64 FaSTLMM::getHeritabilityEstimate() {
    float64 sigma_g = exp(__null_model.getLogSigma());
    float64 sigma_e = exp(__null_model.getLogDelta())*sigma_g;
    return sigma_g/(sigma_e+sigma_g);
}

float64 FaSTLMM::getGeneticVariance() {
    return exp(__null_model.getLogSigma());
}

float64 FaSTLMM::getNoiseVariance() {
    float64 sigma_g = exp(__null_model.getLogSigma());
    return exp(__null_model.getLogDelta())*sigma_g;
}

VectorXd FaSTLMM::computeVarianceExplainedNullModel(uint const& folds) {
	MatrixXd covariates = __covs;
    if (__covs_set) {
        covariates = MatrixXd(__y.rows(),1+__covs.cols());
        covariates << VectorXd::Ones(__y.rows()), __covs;
    } else {
		covariates = VectorXd::Ones(__y.rows());
    }

    CCrossValidation cv(rand()%1000);
    //cv.kFold(folds,__X.rows());
    cv.ShuffleSplit(__X.rows(),folds,0.1);

    VectorXd y_estimated = VectorXd::Zero(__y.rows());
    
    VectorXd vtest = VectorXd::Zero(folds);

    for(uint k=0; k<folds;k++) {
        VectorXd tr_indices = cv.getTrainingIndices(k);
        VectorXd te_indices = cv.getTestingIndices(k);
        VectorXd y_train = sliceRowsMatrix(__y,tr_indices);
        VectorXd y_test = sliceRowsMatrix(__y,te_indices);
        MatrixXd cov_train = sliceRowsMatrix(covariates,tr_indices);
        MatrixXd cov_test = sliceRowsMatrix(covariates,te_indices);
        MatrixXd K_train = sliceRowsMatrix(__K,tr_indices);
        K_train = sliceColsMatrix(K_train,tr_indices);
        MatrixXd K_test = sliceRowsMatrix(__K,te_indices);
        K_test = sliceColsMatrix(K_test,tr_indices);
        CLinearMixedRegression null_model(false);
	    null_model.setInterval(100);
	    null_model.setLogDeltaMin(-5.0);
	    null_model.setLogDeltaMax(5.0);
	    null_model.setREML(__REML);
    	null_model.setBrent(__use_brent);
        null_model.fit(y_train,cov_train,K_train);

        VectorXd y_tmp;
        null_model.predict(&y_tmp,cov_test,K_test);
        insertColumnVectorAtIndices(&y_estimated,y_tmp,te_indices);       
        vtest(k) = 1.0 - CStats::varf(y_test.array() - y_tmp.array())/CStats::varf(y_test.array());
        //vtest(k) = pow(CStats::pearson_corr(y_test,y_tmp),2);
    }
    //float64 variance_explained = 1.0 - CStats::varf(__y.array() - y_estimated.array())/CStats::varf(__y);
    //return pow(CStats::pearson_corr(__y,y_estimated),2);
    //return variance_explained;
    return vtest;
}

void FaSTLMM::test_associations() {
	//Init data strucutres
    __p_values = VectorXd::Zero(__X.cols());
	__test_statistics = VectorXd::Zero(__X.cols());
	__logLikelihoods = VectorXd::Zero(__X.cols());
	if (__intercept)
		__n_features = 2;
	else
		__n_features = 1;
	Eigen::SelfAdjointEigenSolver<MatrixXd> solver(__K);
	MatrixXd S = solver.eigenvalues();
	MatrixXd U = solver.eigenvectors();
	MatrixXd Uy = U.transpose()*__y;
	MatrixXd Ucovs;
	MatrixXd Ux = U.transpose()*__X;
	if (__covs_set) {
		__n_features += __covs.cols();
		MatrixXd tmp(__y.rows(),1+__covs.cols());
        tmp << VectorXd::Ones(__y.rows()), __covs;
		Ucovs = U.transpose()*tmp;
	} else {
		Ucovs = U.transpose()*VectorXd::Ones(__y.rows());
	}
	__betas = MatrixXd::Zero(__X.cols(),__n_features);
	__se_betas = MatrixXd::Zero(__X.cols(),__n_features);
	//Compute NUll Model
	
	__null_model = CLinearMixedRegression(Ucovs,Uy,S,false);
	__null_model.setInterval(100);
	__null_model.setLogDeltaMin(-5.0);
	__null_model.setLogDeltaMax(5.0);
	__null_model.setREML(__REML);
	__null_model.setBrent(__use_brent);
	if (__covs_set) {
		MatrixXd tmp(__y.rows(),1+__covs.cols());
	        tmp << VectorXd::Ones(__y.rows()), __covs;
		__null_model.fit(__y,tmp,__K);
	} else {
		__null_model.fit(__y,VectorXd::Ones(__y.rows()),__K);
	}

	for(uint i=0; i < __X.cols(); i++) {
		//Compute Alternative Models
		CLinearMixedRegression a_model;
		if (__covs_set) {
			MatrixXd Utmp(__y.rows(),2+__covs.cols());
			Utmp << VectorXd::Ones(__y.rows()), __X.col(i), __covs;
			Utmp = U.transpose()*Utmp;
			a_model = CLinearMixedRegression(Utmp,Uy,S,false);
			a_model.setInterval(100);
			a_model.setLogDeltaMin(-5.0);
			a_model.setLogDeltaMax(5.0);
			a_model.setREML(__REML);
			a_model.setBrent(__use_brent);
			a_model.setLogDelta(__null_model.getLogDelta());
			MatrixXd tmp(__y.rows(),1+__covs.cols());
			tmp << __X.col(i), __covs;
			a_model.fit(__y,tmp,__K);
		} else {
			MatrixXd Utmp(__y.rows(),2);
			//Utmp << VectorXd::Ones(__y.rows()), Ux.col(i);
			Utmp << VectorXd::Ones(__y.rows()), __X.col(i);
			Utmp = U.transpose()*Utmp;
			a_model = CLinearMixedRegression(Utmp,Uy,S,false);
			a_model.setInterval(100);
			a_model.setLogDeltaMin(-5.0);
			a_model.setLogDeltaMax(5.0);
			a_model.setREML(__REML);
			a_model.setBrent(__use_brent);
			a_model.setLogDelta(__null_model.getLogDelta());
			a_model.fit(__y,__X.col(i),__K);
		}
		//store estimated data
		__betas.block(i,0,1,__n_features) = a_model.getBetas().transpose();
		__se_betas.block(i,0,1,__n_features) = a_model.getStdBetas().transpose();
		__logLikelihoods(i) = a_model.getLogLikelihood();
		
		//perform a statistical_test
		__test_statistics(i) = (a_model.getLogLikelihood()-__null_model.getLogLikelihood());
		__p_values(i) = CChi2::sf(2.0*__test_statistics(i),1);
	}
}

}; //namespace CSingleTraitGWAS
