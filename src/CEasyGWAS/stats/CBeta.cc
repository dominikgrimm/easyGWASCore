#include <cmath>

#include "CBeta.h"

#include "CEasyGWAS/utils/CMathHelper.h"

#include "Cephes/cephes.h"

void CBeta::__checkParameters(float64 const& x,float64 const& alpha, float64 const& beta) throw (CBetaException) {
	if (x<0.0 || x>1.0) throw CBetaException("X has to be in the interval of [0,1]");
	if (alpha==0 || beta==0) throw CBetaException("Alpha and beta cannont be zero!");
	else if (alpha<0 || beta<0) throw CBetaException("Alpha and beta cannont be negative!");
}

float64 CBeta::cdf(float64 const& x, float64 const& alpha, float64 const& beta) throw (CBetaException) {
	__checkParameters(x,alpha,beta);
    return Cephes::incbet(alpha,beta,x);
    //return Special::regularizedIncompleteBetaFunction(x,alpha,beta);
}

float64 CBeta::logcdf(float64 const& x, float64 const& alpha, float64 const& beta) throw (CBetaException) {
	__checkParameters(x,alpha,beta);
	return log(cdf(x,alpha,beta));
}

float64 CBeta::sf(float64 const& x, float64 const& alpha, float64 const& beta) throw (CBetaException) {
	__checkParameters(x,alpha,beta);
	return 1.0-cdf(x,alpha,beta);
}

float64 CBeta::logsf(float64 const& x, float64 const& alpha, float64 const& beta) throw (CBetaException) {
	__checkParameters(x,alpha,beta);
	return log(sf(x,alpha,beta));
}

float64 CBeta::isf(float64 const& x, float64 const& alpha, float64 const& beta) throw (CBetaException) {
	__checkParameters(x,alpha,beta);
	return Cephes::incbi(alpha,beta,1.0-x);
}

float64 CBeta::ppf(float64 const& x, float64 const& alpha, float64 const& beta) throw (CBetaException) {
	__checkParameters(x,alpha,beta);
	return Cephes::incbi(alpha,beta,x);
}

float64 CBeta::pdf(float64 const& x, float64 const& alpha, float64 const& beta) throw (CBetaException) {
	__checkParameters(x,alpha,beta);
	return pow(x,alpha-1.0)*pow(1.0-x,beta-1.0)/tbeta(alpha,beta);
}

float64 CBeta::logpdf(float64 const& x, float64 const& alpha, float64 const& beta) throw (CBetaException) {
	__checkParameters(x,alpha,beta);
	return log(pdf(x,alpha,beta));
}

/*
BetaCF::BetaCF(float64 const& x, float64 const& alpha, float64 const& beta) {
	__x = x;
	__alpha = alpha;
	__beta = beta;
	_no_odd = false;
}	

cf_param BetaCF::cf_fraction(uint const& n) {
	cf_param param;
	if(n==0) {
		param.a_even = 1.0-(__alpha+__beta)*__x/(__alpha+1.0);
		param.b_even = 1.0;
	} else {
		param.a_odd = - ((__alpha+n)*(__alpha+__beta+n)*__x)/((__alpha+2*n)*(__alpha+1.0+2*n));
		param.b_odd = 1.0;
		param.a_even = (n*(__beta-n)*__x)/((__alpha-1.0+2*n)*(__alpha+2*n));
		param.b_even = 1.0;
	}
	return param;
}

float64 CBeta::Special::regularizedIncompleteBetaFunction(float64 const& x, float64 const& alpha, float64 const& beta) {
	__checkParameters(x,alpha,beta);
	float64 beta_f = exp(lgamma(beta+alpha) - lgamma(alpha) - lgamma(beta) + alpha * log(x)+beta*log(1-x));
	//Regular C-fraction for Ix(a,b) REF: "Handbook of Continued Fractions for Special Functions", p385, forumla: 18.5.17a
	if (x > (alpha+1.0)/(alpha+beta+2.0)) { //Speedup trick by using the symetrie of the function Ix(a,b) = I(1-x)(b,a)
		//float64 beta_f = pow(1.0-x,beta)*pow(x,alpha)/(beta*tbeta(beta,alpha));
		BetaCF cf(1.0-x,beta,alpha);
		//return 1.0-evaluateContinuedFraction(&cf)*beta_f;
		return 1.0-evaluateContinuedFraction(&cf)*beta_f/alpha;
	} else {
		//float64 beta_f = pow(x,alpha)*pow(1.0-x,beta)/(alpha*tbeta(alpha,beta));
		BetaCF cf(x,alpha,beta);
		//return evaluateContinuedFraction(&cf)*beta_f;
		return evaluateContinuedFraction(&cf)*beta_f/alpha;
	}
}
*/
