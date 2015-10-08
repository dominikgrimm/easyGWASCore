#include "CRegression.h"
#include "../stats/CChi2.h"
#include "../stats/CFisherF.h"
#include "../utils/StringHelper.h"
#include "../utils/CMathHelper.h"
#include <unsupported/Eigen/IterativeSolvers>
#include <ctime>

/*
*CRegression Member methods
*/

CRegression::CRegression() {
	_x = MatrixXd::Zero(0,0);
	_y = VectorXd::Zero(0);
	_intercept = true;
}

CRegression::CRegression(bool const& intercept) {
	_x = MatrixXd::Zero(0,0);
	_y = VectorXd::Zero(0);
	_intercept = intercept;
}

void CRegression::_checkDimensions(uint64* n_samples,
				  uint64* n_features) throw (CRegressionException) {
	(*n_samples) = _y.rows();
	(*n_features) = _x.cols();
	if(_y.rows()!=_x.rows()) 
		throw CRegressionException("Y and X have different dimensions: Y=[" 
			+ StringHelper::to_string<uint64>(_y.rows()) + "," 
			+ StringHelper::to_string<uint64>(_y.cols()) + "], X=["
		       	+ StringHelper::to_string<uint64>(_x.rows()) + "," 
			+ StringHelper::to_string<uint64>(_x.cols()) + "]");
}

void CRegression::setX(MatrixXd const& X) {
	_x = X;
}

void CRegression::setY(VectorXd const& y) {
	_y = y;
}

float64 CRegression::getRSquared() const {
	return (_yhat.array()-_y.mean()).array().pow(2).sum()/
		(_y.array()-_y.mean()).array().pow(2).sum();
}

float64 CRegression::getAdjustedRSquared() const {
	return 1.0f - (_n_samples-1.0f)/(_n_samples-_rank)*(1.0f-getRSquared());
}

float64 CRegression::getLogLikelihood() const {
	return _loglikelihood;
}	

float64 CRegression::getAIC() const {
	return 2.0*_n_features - 2.0 * _loglikelihood;
}

float64 CRegression::getAICc() const {
	return getAIC() + (2.0*_n_features*(_n_features+1.0))/(_n_samples-_n_features-1.0);
}

float64 CRegression::getBIC() const {
	return _n_features*log(_n_samples) - 2.0 * _loglikelihood;
}

VectorXd CRegression::getBetas() const {
	return _betas;
}

VectorXd CRegression::getResiduals() const {
	return _residuals;
}

VectorXd CRegression::getYHat() {
	return _yhat;
} 

MatrixXd CRegression::getCovarianceBetas() {
	MatrixXd nM;
	return nM;
}

VectorXd CRegression::getStdBetas() {
	VectorXd nV;
	return nV;
}

uint CRegression::getDF() const {
	return _n_features;
}

void CRegression::print() {
	std::string formula, parameters;
	formula = "\n" + _model_info + ": y ~";
	parameters = "\t\tBetas\t\tSTD Betas\n";
	for(uint i=0; i<_n_features;i++) {
		if(_intercept && i==0){
		       	formula += " 1";
			parameters += "(Intercept)\t";
			parameters += StringHelper::to_string<float64>(_betas[i]) + "\t";
			parameters += StringHelper::to_string<float64>(getStdBetas()[i]) + "\n";
		}
		else if(_intercept==false && i==0) {
		       	formula += " x" + StringHelper::to_string<uint>(i+1);
			parameters += "x" + StringHelper::to_string<uint>(i+1) + "\t\t";
			parameters += StringHelper::to_string<float64>(_betas[i]) + "\t";
			parameters += StringHelper::to_string<float64>(getStdBetas()[i]) + "\n";
		}
		else {
			formula += " + x" + StringHelper::to_string<uint>(i);
			parameters += "x" + StringHelper::to_string<uint>(i) + "\t\t";
			parameters += StringHelper::to_string<float64>(_betas[i]) + "\t";
			parameters += StringHelper::to_string<float64>(getStdBetas()[i]) + "\n";
		}
	}
	
	std::cout << formula << "\n\n";
	std::cout << "Estimated Parameters: \n";
	std::cout << parameters << "\n\n";
	std::cout << "LogLikelihood:\t" << getLogLikelihood() << "\n";
	std::cout << "AIC:\t" << getAIC() << "\n";
	std::cout << "AICc:\t" << getAICc() << "\n";
	std::cout << "BIC:\t" << getBIC() << "\n";
	
}


/*
*CLinearRegression Member methods
*/

//Standard Constructor
CLinearRegression::CLinearRegression() {
	_x = MatrixXd::Zero(0,0);
	_y = VectorXd::Zero(0);
	_intercept = true;
	_model_info = "Linear Regression";
}

//Overloaded Constructor
CLinearRegression::CLinearRegression(bool const& intercept) : CRegression(intercept) {
	_x = MatrixXd::Zero(0,0);
	_y = VectorXd::Zero(0);
	_intercept = intercept;
	_model_info = "Linear Regression";
}

MatrixXd CLinearRegression::getCovarianceBetas() {
	return __variance*(_x.transpose()*_x).inverse();
}

VectorXd CLinearRegression::getStdBetas() {
	return getCovarianceBetas().diagonal().array().sqrt();
}

float64 CLinearRegression::getMSE() {
	return __variance;
}

float64 CLinearRegression::getRMSE() {
	return sqrt(__variance);
}

void CLinearRegression::_estimateBetas() {
	//_betas = (_x.transpose()*_x).inverse()*_x.transpose()*_y;
	MatrixXd pinv_x;
	pinv<MatrixXd>(_x,&pinv_x);
	_betas = pinv_x*_y;
}

void CLinearRegression::_estimateLogLikelihood() {
	_yhat = _x*_betas;
	_residuals = _y - _yhat;
    Eigen::FullPivLU<MatrixXd> luX(_x);
	_rank = luX.rank();
	__variance = (_residuals.transpose()*_residuals);
	__variance = __variance/(_n_samples-_rank);
	_loglikelihood = -(_n_samples/2.0f) * log(2.0f*PI*__variance) - 1.0f/(2.0f*__variance) * (_residuals.array().pow(2)).sum();
	
	/* This is only if matrix is not invertible
	MatrixXd pinv_x;
	pinv<MatrixXd>(x,&pinv_x);
	_betas = pinv_x*y;
	*/
	
}

void CLinearRegression::predict(VectorXd* results,MatrixXd const& x) throw (CRegressionException) {
	MatrixXd xtmp(_n_samples,_n_features);
	if(_intercept) {
		if(x.cols()!=(int64)_n_features-1) 
			throw CRegressionException("Prediction Dimension Error: #Features is wrong");
		xtmp << VectorXd::Ones(_n_samples,1), x;
		(*results) = xtmp*_betas;
	} else {
		if(x.cols()!=(int64)_n_features) 
			throw CRegressionException("Prediction Dimension Error: #Features is wrong");
		(*results) = x*_betas;
	}	
}

//Fit Linear Regression Model using maximum likelihood assuming Gaussian Noise
void CLinearRegression::fit(VectorXd const& y, MatrixXd const& x) throw (CRegressionException) {
	_x = x;
	_y = y;
	fit(true);
}

void CLinearRegression::fit(VectorXd const& y, MatrixXd const& x, bool const& estimateLL) throw (CRegressionException) {
	_x = x;
	_y = y;
	fit(estimateLL);
}
void CLinearRegression::fit() throw(CRegressionException) {
	fit(true);
}

void CLinearRegression::fit(bool const& estimateLL) throw(CRegressionException) {
	//checkDimensions
	_checkDimensions(&_n_samples,&_n_features);
	//Add intercept to x matrix if necessary
	if (_intercept) {
		_n_features += 1;
		MatrixXd tmp = _x;
		_x.resize(_n_samples,_n_features);
		MatrixXd intercept_v = VectorXd::Ones(_n_samples,1);
		_x << intercept_v, tmp;
	}

	//estimate betas
	_estimateBetas();
	//estimate log likelihood 
	if(estimateLL) _estimateLogLikelihood();
}


/*
*CLogisticRegression Member methods
*/

//Standard Constructor
CLogisticRegression::CLogisticRegression() {
	_x = MatrixXd::Zero(0,0);
	_y = VectorXd::Zero(0);
	_intercept = true;
	_model_info = "Logistic Regression";
	__epsilon = 1e-10;
}

//Overloaded Constructor
CLogisticRegression::CLogisticRegression(bool const& intercept) : CRegression(intercept) {
	_x = MatrixXd::Zero(0,0);
	_y = VectorXd::Zero(0);
	_intercept = intercept;
	_model_info = "Logistic Regression";
	__epsilon = 1e-10;
}

//Overloaded Constructor
CLogisticRegression::CLogisticRegression(bool const& intercept, float64 const& epsilon) : CRegression(intercept) {
	_x = MatrixXd::Zero(0,0);
	_y = VectorXd::Zero(0);
	_intercept = intercept;
	_model_info = "Logistic Regression";
	__epsilon = epsilon;
}

uint CLogisticRegression::getIterations() {
	return __iterations;
}

void CLogisticRegression::_estimateBetas() {
	MatrixXd beta0 = MatrixXd::Zero(_n_features,1);
	MatrixXd h = (_x*beta0).array().exp()/(1+(_x*beta0).array().exp());
	MatrixXd h_1 = 1.0-h.array();
	DiagXd Eye(_n_samples);
	Eye.diagonal() = VectorXd::Ones(_n_samples).cwiseProduct(h).cwiseProduct(h_1);
	MatrixXd W = Eye;
	__F = _x.transpose()*(W*_x);
	MatrixXd S = _x.transpose()*(_y-h);
	_betas = __F.colPivHouseholderQr().solve(S);
	_betas = beta0 + _betas;
	//_betas = beta0 + __F.inverse()*S;
	//Performe Fisher Scoring
	__iterations=0;
	while((_betas-beta0).norm()>__epsilon) {
		beta0 = _betas;
		h = (_x*beta0).array().exp()/(1.0+(_x*beta0).array().exp());
		h_1 = 1.0-h.array();
		Eye.diagonal() = VectorXd::Ones(_n_samples).cwiseProduct(h).cwiseProduct(h_1);
		W = Eye;
		__F = _x.transpose()*(W*_x);
		__F.array() += 1e-32;
		S = _x.transpose()*(_y-h);
		_betas = __F.colPivHouseholderQr().solve(S);
		_betas = beta0 + _betas;
		//_betas = beta0 + __F.inverse()*S;
		__iterations++;
		if(__iterations>=500) break;
	}
}

VectorXd CLogisticRegression::getYHat() {
	_yhat = _x*_betas;
	return _yhat;
} 

MatrixXd CLogisticRegression::getCovarianceBetas() {
	return __F.inverse();
}

VectorXd CLogisticRegression::getStdBetas() {
	return getCovarianceBetas().diagonal().array().sqrt();
}

void CLogisticRegression::_estimateLogLikelihood() {
	_yhat = _x*_betas;
	_residuals = _y - _yhat;
	_loglikelihood =  ((_y.cwiseProduct(_yhat)).array() - (_yhat.array().exp()+1.0).array().log().array()).array().sum();
}

void CLogisticRegression::predict(VectorXd* results,MatrixXd const& x) throw (CRegressionException) {
	MatrixXd xtmp(_n_samples,_n_features);
	if(_intercept) {
		if(x.cols()!=(int64)_n_features-1) 
			throw CRegressionException("Prediction Dimension Error: #Features is wrong");
		xtmp << VectorXd::Ones(_n_samples,1), x;
		(*results) = xtmp*_betas;
	} else {
		if(x.cols()!=(int64)_n_features) 
			throw CRegressionException("Prediction Dimension Error: #Features is wrong");
		(*results) = x*_betas;
	}	
}

void CLogisticRegression::fit(VectorXd const& y, MatrixXd const& x) throw (CRegressionException) {
	_x = x;
	_y = y;
	fit(true);
}

void CLogisticRegression::fit(VectorXd const& y, MatrixXd const& x, bool const& estimateLL) throw (CRegressionException) {
	_x = x;
	_y = y;
	fit(estimateLL);
}
void CLogisticRegression::fit() throw(CRegressionException) {
	fit(true);
}

//Fit Linear Regression Model using maximum likelihood assuming Gaussian Noise
void CLogisticRegression::fit(bool const& estimateLL) throw (CRegressionException) {
	//checkDimensions
	_checkDimensions(&_n_samples,&_n_features);

	//Add intercept to x matrix if necessary
	if (_intercept) {
		_n_features += 1;
		MatrixXd tmp = _x;
		_x.resize(_n_samples,_n_features);
		MatrixXd intercept_v = VectorXd::Ones(_n_samples,1);
		_x << intercept_v, tmp;
	}

	//estimate betas
	_estimateBetas();
	//estimate log likelihood 
	if(estimateLL) _estimateLogLikelihood();
}

/*
*CLinearMixedRegression Member methods
*/

//Standard Constructor
CLinearMixedRegression::CLinearMixedRegression() {
	_loglikelihood = 0;
	_intercept = true;
	__transform = true;
	_model_info = "LinearMixed Regression";
	__REML = false;
	__nInterval = 500;
	__logDeltaMin = -5;
	__logDeltaMax = 5;
	__logDeltaSet = false;
	__use_brent = true;
}

//Overloaded Constructor
CLinearMixedRegression::CLinearMixedRegression(bool const& intercept) : CRegression(intercept) {
	_loglikelihood = 0;
	_intercept = intercept;
	__transform = true;
	__REML = false;
	_model_info = "LinearMixed Regression";
	__nInterval = 100;
	__logDeltaMin = -5;
	__logDeltaMax = 5;
	__logDeltaSet = false;
	__use_brent = true;
}

CLinearMixedRegression::CLinearMixedRegression(MatrixXd const& Ux, MatrixXd const& Uy, MatrixXd const& S) {
	_loglikelihood = 0;
	_intercept = true;
	__transform = false;
	__REML = false;
	_model_info = "LinearMixed Regression";
	__Ux = Ux;
	__Uy = Uy;
	__S = S;
	__nInterval = 100;
	__logDeltaMin = -5;
	__logDeltaMax = 5;
	__logDeltaSet = false;
	__use_brent = true;
}

CLinearMixedRegression::CLinearMixedRegression(MatrixXd const& Ux, MatrixXd const& Uy, MatrixXd const& S, bool const& intercept) : CRegression(intercept) {
	_loglikelihood = 0;
	_intercept = intercept;
	__transform = false;
	__REML = false;
	_model_info = "LinearMixed Regression";
	__Ux = Ux;
	__Uy = Uy;
	__S = S;
	__nInterval = 100;
	__logDeltaMin = -5;
	__logDeltaMax = 5;
	__logDeltaSet = false;
	__use_brent = true;
}

float64 CLinearMixedRegression::getMSE() {
	return NAN;
}

float64 CLinearMixedRegression::getRMSE() {
	return NAN;
}

float64 CLinearMixedRegression::getLogDelta() {
	return __logDelta;
}

float64 CLinearMixedRegression::getLogSigma() {
	return __sigma;
}

void CLinearMixedRegression::setK(MatrixXd const& K) {
	__K = K;
}
void CLinearMixedRegression::setInterval(uint const& interval) {
	__nInterval = interval;
}
		
void CLinearMixedRegression::setLogDeltaMin(float64 const& d) {
	__logDeltaMin = d;
}
		
void CLinearMixedRegression::setLogDeltaMax(float64 const& d) {
	__logDeltaMax = d;
}

void CLinearMixedRegression::setLogDelta(float64 const& d) {
	__logDelta = d;
	__logDeltaSet = true;
}


void CLinearMixedRegression::setREML(bool const& REML) {
	__REML = REML;
}

void CLinearMixedRegression::setBrent(bool const& brent) {
	__use_brent = brent;
}

VectorXd CLinearMixedRegression::getYHat() {
	VectorXd fixed_component = _x*_betas;
	MatrixXd random_component = __K.array()+__logDelta;
	random_component = random_component.inverse();
	random_component = __K*random_component*(_y-_x*_betas);
	return fixed_component + random_component; 
} 

MatrixXd CLinearMixedRegression::getCovarianceBetas() {
	return __covarianceBetas;
}

VectorXd CLinearMixedRegression::getStdBetas() {
	return getCovarianceBetas().diagonal().array().sqrt();
}


void CLinearMixedRegression::__evaluateNLL(float64 const& ldelta, VectorXd* betas, float64* logLikelihood, float64* sigma) {
	//Transform data
	float64 delta = exp(ldelta);
	VectorXd Sd = (__S.array() + delta);
	Sd = Sd.transpose();
	float64 ldet = (log(Sd.array())).sum();
	Sd = 1.0/Sd.array();
	MatrixXd xSdi = __Ux.transpose().array().rowwise() * Sd.transpose().array();
	
	MatrixXd xSx = xSdi * __Ux;
	MatrixXd xSy = xSdi * __Uy;
	
	//Ordinary least squares
	//(*betas) = xSx.inverse()*xSy;
	(*betas) = xSx.colPivHouseholderQr().solve(xSy);
	//Evaluate residuals
	_residuals = (__Uy - (__Ux * (*betas)));
	_residuals = _residuals.array() * _residuals.array();
	_residuals = _residuals.array() * Sd.array();
	//compute sigma square 
	//use REML estimate if set otherwise use ML
	if(__REML) {
		MatrixXd XX = __Ux.transpose() * __Ux;
		Eigen::SelfAdjointEigenSolver<MatrixXd> reml_solver(XX);
		MatrixXd sXX = reml_solver.eigenvalues();
		float64 tmp = (log(sXX.array())).sum();
		Eigen::SelfAdjointEigenSolver<MatrixXd> solver(xSx);
		MatrixXd tmp_S = solver.eigenvalues();
		float64 ldetXSX = (log(tmp_S.array())).sum();
		//compute LogLikelihood
		(*sigma) = _residuals.sum()/(_n_samples-_n_features);
		(*logLikelihood) = 0.5*(_n_samples-_n_features)*log((*sigma)) + 0.5*((_n_samples-_n_features)*log(2.0*PI)+ldet+tmp+ldetXSX+(_n_samples-_n_features));
	} else {
		(*sigma) = _residuals.sum()/_n_samples;
		(*logLikelihood) = 0.5*(_n_samples*log(2.0*PI*(*sigma))+ldet+_n_samples);
	}

	__covarianceBetas = xSx.inverse();
}

float64 CLinearMixedRegression::__optimizeDelta() {
	_loglikelihood = std::numeric_limits<float64>::infinity();
	__logDeltas = VectorXd::LinSpaced(__nInterval,__logDeltaMin,__logDeltaMax);
	float64 LogLikelihood=0;
	float64 sigma=0;
	VectorXd betas;
    float64 logDelta=__logDeltaMin;
	VectorXd nll_grid(__logDeltas.rows());
	for(int i=0; i<__logDeltas.rows();i++) {
		//evaluate log likelihood
		__evaluateNLL(__logDeltas[i],&betas,&LogLikelihood,&sigma);
		nll_grid(i) = LogLikelihood;	
		if (LogLikelihood<_loglikelihood) {
			_loglikelihood = LogLikelihood;
			logDelta = __logDeltas[i];	
		}
	}
	if(__use_brent) {
		BrentFunction func(__Uy,__Ux,__S,__REML);
		for(int64 i=1; i < __logDeltas.rows()-1; i++) {
			if(nll_grid(i)<nll_grid(i+1) && nll_grid(i) < nll_grid(i-1)) {
				logDelta = CBrentOptimizer::solve(&func,__logDeltas(i-1),__logDeltas(i+1),1e-16,100);
			}
		}
	}
	return logDelta;
}

VectorXd CLinearMixedRegression::predict(MatrixXd const& x,MatrixXd const& K) throw (CRegressionException) {
    VectorXd results;
    predict(&results,x,K);
    return results;
}

void CLinearMixedRegression::predict(VectorXd* results,MatrixXd const& x,MatrixXd const& K) throw (CRegressionException) {
	MatrixXd xtmp(x.rows(),_n_features);
	if(_intercept) {
		if(x.cols()!=(int64)_n_features-1) 
			throw CRegressionException("Prediction Dimension Error: #Features is wrong");
        xtmp << VectorXd::Ones(x.rows(),1), x;
		VectorXd fixed_component = xtmp*_betas;
		//MatrixXd random_component = __K.array()+exp(__logDelta);
        MatrixXd identity = MatrixXd::Identity(__K.rows(),__K.cols());
		MatrixXd random_component = __K.array()+identity.array()*exp(__logDelta);
		random_component = random_component.inverse();
		random_component = K*random_component*(_y-_x*_betas);
		(*results) = fixed_component + random_component; 
	} else {
		if(x.cols()!=(int64)_n_features) 
			throw CRegressionException("Prediction Dimension Error: #Features is wrong");
        xtmp << x;
        VectorXd fixed_component = xtmp*_betas;
		//MatrixXd random_component = __K.array()+exp(__logDelta);
        //logging(STATUS,random_component.sum());
        MatrixXd identity = MatrixXd::Identity(__K.rows(),__K.cols());
		MatrixXd random_component = __K.array()+identity.array()*exp(__logDelta);
		random_component = random_component.inverse();
        random_component = K*random_component*(_y-_x*_betas);
		(*results) = fixed_component + random_component; 
	}	
}

void CLinearMixedRegression::fit(VectorXd const& y, MatrixXd const& X, MatrixXd const& K) throw (CRegressionException) {
	_x = X;
	_y = y;
	__K = K;
	fit();
}

//Fit LinearMixed Regression Model using maximum likelihood assuming Gaussian Noise
void CLinearMixedRegression::fit() throw (CRegressionException) {
	//checkDimensions
	_checkDimensions(&_n_samples,&_n_features);
	//Add intercept to x matrix if necessary
	if (_intercept) {
		_n_features += 1;
		MatrixXd tmp = _x;
		_x.resize(_n_samples,_n_features);
		MatrixXd intercept_v = VectorXd::Ones(_n_samples,1);
		_x << intercept_v, tmp;
	}
	if(__transform) {
		//Rotate data
		Eigen::SelfAdjointEigenSolver<MatrixXd> solver(__K);
		__S = solver.eigenvalues();
		__U = solver.eigenvectors();
		__Ux = __U.transpose() * _x;
		__Uy = __U.transpose() * _y;
	}
	//Optimize delta
	if (__logDeltaSet==false) {
		__logDelta = __optimizeDelta();
	}
	//estimate negative LogLikelihood
	__evaluateNLL(__logDelta,&_betas,&_loglikelihood,&__sigma);
	_loglikelihood = -_loglikelihood;
	
}

BrentFunction::BrentFunction(MatrixXd const& Uy, MatrixXd const& Ux, MatrixXd const& S, bool const& REML) {
	__Uy = Uy;
	__Ux = Ux;
	__S = S;
	__REML = REML;
}

float64 BrentFunction::evaluate(float64 const& ldelta) {
	//Transform data
	float64 n_samples = __Ux.rows();
	float64 n_features = __Ux.cols();
	float64 delta = exp(ldelta);
	VectorXd Sd = (__S.array() + delta);
	Sd = Sd.transpose();
	float64 ldet = (log(Sd.array())).sum();
	Sd = 1.0/Sd.array();
	MatrixXd xSdi = __Ux.transpose().array().rowwise() * Sd.transpose().array();
	
	MatrixXd xSx = xSdi * __Ux;
	MatrixXd xSy = xSdi * __Uy;
	
	//Ordinary least squares
	VectorXd betas = xSx.colPivHouseholderQr().solve(xSy);
	//Evaluate residuals
	VectorXd residuals = (__Uy - (__Ux * betas));
	//compute sigma square 
	residuals = residuals.array() * residuals.array();
	residuals = residuals.array() * Sd.array();
	//use REML estimate if set otherwise use MaximumLikelihood
	if(__REML) {
		MatrixXd XX = __Ux.transpose() * __Ux;
		Eigen::SelfAdjointEigenSolver<MatrixXd> reml_solver(XX);
		MatrixXd sXX = reml_solver.eigenvalues();
		float64 tmp = (log(sXX.array())).sum();
		Eigen::SelfAdjointEigenSolver<MatrixXd> solver(xSx);
		MatrixXd tmp_S = solver.eigenvalues();
		float64 ldetXSX = (log(tmp_S.array())).sum();
		//compute LogLikelihood
		float64 sigma = residuals.sum()/(n_samples-n_features);
		return 0.5*(n_samples-n_features)*log(sigma) + 0.5*((n_samples-n_features)*log(2.0*PI)+ldet+tmp+ldetXSX+(n_samples-n_features));
	} else {
		float64 sigma = residuals.sum()/n_samples;
		return 0.5*(n_samples*log(2.0*PI*sigma)+ldet+n_samples);
	}
}
