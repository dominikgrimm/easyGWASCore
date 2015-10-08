#include <cmath>

#include "CStudentT.h"
#include "Cephes/cephes.h"

void CStudentT::__checkParameters(int const& d1) throw (CStudentTException) {
	if (d1<=0) throw CStudentTException("d1 (degrees of freedom) has to be a positive real number!");
}

float64 CStudentT::cdf(float64 const& x, int const& d1) throw (CStudentTException) {
	__checkParameters(d1);
	return Cephes::stdtr(d1,x);
}

float64 CStudentT::logcdf(float64 const& x, int const& d1) throw (CStudentTException) {
	__checkParameters(d1);
	return log(cdf(x,d1));
}

float64 CStudentT::sf(float64 const& x, int const& d1) throw (CStudentTException) {
	__checkParameters(d1);
	return Cephes::stdtr(d1,-x);
}

float64 CStudentT::logsf(float64 const& x, int const& d1) throw (CStudentTException) {
	__checkParameters(d1);
	return log(sf(x,d1));
}

float64 CStudentT::pdf(float64 const& x, int const& d1) throw (CStudentTException) {
	__checkParameters(d1);
	return 1.0/(sqrt((float64)d1)*tbeta(0.5,0.5*(float64)d1))*pow(1+pow(x,2)/(float64)d1,-0.5*(float64)(d1+1.0));
}

float64 CStudentT::logpdf(float64 const& x, int const& d1) throw (CStudentTException) {
	__checkParameters(d1);
	return log(pdf(x,d1));
}

float64 CStudentT::ppf(float64 const& x, int const& d1) throw (CStudentTException) {
	__checkParameters(d1);
	return Cephes::stdtri(d1,x);
}

float64 CStudentT::isf(float64 const& x, int const& d1) throw (CStudentTException) {
	__checkParameters(d1);
	return -Cephes::stdtri(d1,x);
}
