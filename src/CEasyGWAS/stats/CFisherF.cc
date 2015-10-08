#include <cmath>

#include "CFisherF.h"
#include "Cephes/cephes.h"

void CFisherF::__checkParameters(float64 const& x,int const& d1, int const& d2) throw (CFisherFException) {
	if (x<0.0) throw CFisherFException("X cannot be negative!");
	if (d1==0 || d2==0) throw CFisherFException("d1 and d2 cannot be zero!");
	else if (d1<0 || d2<0) throw CFisherFException("d1 and d2 cannot be negative!");
}

float64 CFisherF::cdf(float64 const& x, int const& d1, int const& d2) throw (CFisherFException) {
	__checkParameters(x,d1,d2);
	return Cephes::fdtr(d1,d2,x);
    //return CBeta::Special::regularizedIncompleteBetaFunction(d1*x/(d1*x+d2),0.5*d1,0.5*d2);
}

float64 CFisherF::logcdf(float64 const& x, int const& d1, int const& d2) throw (CFisherFException) {
	__checkParameters(x,d1,d2);
	return log(cdf(x,d1,d2));
}

float64 CFisherF::sf(float64 const& x, int const& d1, int const& d2) throw (CFisherFException) {
	__checkParameters(x,d1,d2);
	return Cephes::fdtrc(d1,d2,x); 
}

float64 CFisherF::isf(float64 const& x, int const& d1, int const& d2) throw (CFisherFException) {
	__checkParameters(x,d1,d2);
	return Cephes::fdtri(d1,d2,x); 
}

float64 CFisherF::ppf(float64 const& x, int const& d1, int const& d2) throw (CFisherFException) {
	__checkParameters(x,d1,d2);
	return Cephes::fdtri((int)d1,(int)d2,1.0-x); 
}

float64 CFisherF::logsf(float64 const& x, int const& d1, int const& d2) throw (CFisherFException) {
	__checkParameters(x,d1,d2);
	return log(sf(x,d1,d2));
}

float64 CFisherF::pdf(float64 const& x, int const& d1, int const& d2) throw (CFisherFException) {
	__checkParameters(x,d1,d2);
	float64 dh1 = 0.5*(float64)d1;
	float64 dh2 = 0.5*(float64)d2;
	return 1.0/tbeta(dh1,dh2)*pow((float64)d1/(float64)d2,dh1)*pow(x,dh1-1)*pow(1.0+(float64)d1/(float64)d2*x,-(float64)(d1+d2)*0.5);
}

float64 CFisherF::logpdf(float64 const& x, int const& d1, int const& d2) throw (CFisherFException) {
	__checkParameters(x,d1,d2);
	return log(pdf(x,d1,d2));
}
