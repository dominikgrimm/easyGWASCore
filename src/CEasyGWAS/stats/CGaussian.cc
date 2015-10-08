#include <cmath>
#include "math.h"
#include <iostream>

#include "CEasyGWAS/utils/CMathHelper.h"
#include "CGaussian.h"
#include "Cephes/cephes.h"

void CGaussian::__checkParameters(float64 const& std) throw (CGaussianException) {
	if (std==0) throw CGaussianException("STD cannont be zero!");
	else if (std<0) throw CGaussianException("STD cannont be negative!");
}

float64 CGaussian::cdf(float64 const& x, float64 const& mean, float64 const& std) throw (CGaussianException) {
	__checkParameters(std);
	return Cephes::ndtr((x-mean)*1.0/std);
    //return 0.5 * (1.0 + erf((x-mean)/(sqrt(2.0*pow(std,2)))));
}

float64 CGaussian::logcdf(float64 const& x, float64 const& mean, float64 const& std) throw (CGaussianException) {
	__checkParameters(std);
	return log(cdf(x,mean,std));
}

float64 CGaussian::sf(float64 const& x, float64 const& mean, float64 const& std) throw (CGaussianException) {
	__checkParameters(std);
	return Cephes::ndtr(-((x-mean)*1.0/std));
}

float64 CGaussian::isf(float64 const& x, float64 const& mean, float64 const& std) throw (CGaussianException) {
	__checkParameters(std);
	return -Cephes::ndtri(x)*std + mean;
}

float64 CGaussian::logsf(float64 const& x, float64 const& mean, float64 const& std) throw (CGaussianException) {
	__checkParameters(std);
	return log(sf(x,mean,std));
}

float64 CGaussian::pdf(float64 const& x, float64 const& mean, float64 const& std) throw (CGaussianException) {
	__checkParameters(std);
    return 1.0f/(std*sqrt(2.0f*PI)) * exp(-(pow(x-mean,2))/(2.0f*pow(std,2)));
}

float64 CGaussian::logpdf(float64 const& x, float64 const& mean, float64 const& std) throw (CGaussianException) {
	__checkParameters(std);
	return log(pdf(x,mean,std));
}

float64 CGaussian::ppf(float64 const& p, float64 const& mean, float64 const& std) throw (CGaussianException) {
	__checkParameters(std);
	if (p>1 || p<0) return NAN;
	if (p==0) return -INFINITY;
	if (p==1) return INFINITY;
    return Cephes::ndtri(p)*std+mean;
}
