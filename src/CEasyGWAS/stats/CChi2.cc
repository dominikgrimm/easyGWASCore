#include <cmath>
#include <iostream>

#include "Cephes/cephes.h"

#include "CChi2.h"
//#include "CGamma.h"

void CChi2::__checkParameters(float64 const& k) throw (CChi2Exception) {
	if (k==0) throw CChi2Exception("Degress of freedom (k) cannont be zero!");
	else if (k<0) throw CChi2Exception("Degress of freedom (k) cannont be negative!");
}

float64 CChi2::cdf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
    return Cephes::chdtr(k,x);
    //Own Implementation
    /*if (k==2.0) { //Special case if k==2.0
		return 1.0 - exp(-0.5*x);
	} else { //for any k
		return (CGamma::Special::regularizedLowerIncompleteGamma(0.5*x,0.5*k));
	}*/	
}

float64 CChi2::logcdf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
	return log(cdf(x,k));
}

float64 CChi2::sf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
    if(x < 0) return 1;
    else return Cephes::chdtrc(k,x);
    //return 1-cdf(x,k); //is too less accurate! Use complemented approximation
    //OWN IMPlEMENTATION
	//return CGamma::Special::complementedIncompleteGamma(0.5*x,0.5*k);
}

float64 CChi2::isf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
	if(x>1.0 || x<0.0) throw CChi2Exception("x has to be in the range [0.0,1.0]!");
    return Cephes::chdtri(k,x);
}

float64 CChi2::logsf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
	return log(sf(x,k));
}

float64 CChi2::pdf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
	if(x<0.0) return 0.0;
	return pow(x,0.5*k-1)*exp(-0.5*x)/(pow(2.0,0.5*k)*tgamma(0.5*k));	
}

float64 CChi2::logpdf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
	return log(pdf(x,k));
}
