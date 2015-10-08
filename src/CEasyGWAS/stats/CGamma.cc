#include <cmath>

#include "CGamma.h"

#include "Cephes/cephes.h"

/*void CGamma::__checkParameters(float64 const& x, float64 const& k, float64 const& theta) throw (CGammaException) {
	if (theta==0 || k==0) throw CGammaException("Theta and k cannont be zero!");
	else if (theta<0 || k<0 || x<0) throw CGammaException("Theta, k and x cannont be negative!");
}*/

void CGamma::__checkParameters(float64 const& x, float64 const& k) throw (CGammaException) {
	if (k==0) throw CGammaException("Theta and k cannont be zero!");
	else if (k<0 || x<0) throw CGammaException("Theta, k and x cannont be negative!");
}

float64 CGamma::cdf(float64 const& x, float64 const& k) throw (CGammaException) {
    return Cephes::igam(k,x);
}

float64 CGamma::logcdf(float64 const& x, float64 const& k) throw (CGammaException) {
    return log(cdf(x,k));
}

float64 CGamma::sf(float64 const& x, float64 const& k) throw (CGammaException) {
	return Cephes::igamc(k,x);
}

float64 CGamma::logsf(float64 const& x, float64 const& k) throw (CGammaException) {
    return log(sf(x,k));
}

float64 CGamma::pdf(float64 const& x, float64 const& k) throw (CGammaException) {
	return (pow(x,k-1)*exp(-x))/tgamma(k);
}

float64 CGamma::logpdf(float64 const& x, float64 const& k) throw (CGammaException) {
    return log(pdf(x,k));
}

/*
float64 CGamma::cdf(float64 const& x, float64 const& k, float64 const& theta) throw (CGammaException) {
	__checkParameters(x,k,theta);
    return Cephes::igam(k,x);
    //OWN IMPLEMENTATION
    //return Special::regularizedLowerIncompleteGamma(x/theta,k);
}

float64 CGamma::logcdf(float64 const& x, float64 const& k, float64 const& theta) throw (CGammaException) {
	__checkParameters(x,k,theta);
	return log(cdf(x,k,theta));
}

float64 CGamma::sf(float64 const& x, float64 const& k, float64 const& theta) throw (CGammaException) {
	__checkParameters(x,k,theta);
	return 1.0-cdf(x,k,theta);
}

float64 CGamma::logsf(float64 const& x, float64 const& k, float64 const& theta) throw (CGammaException) {
	__checkParameters(x,k,theta);
	return log(sf(x,k,theta));
}
*/


/*
float64 CGamma::pdf(float64 const& x, float64 const& k, float64 const& theta=1.0) throw (CGammaException) {
	__checkParameters(x,k,theta);
	return (pow(x,k-1)*exp(-x/theta))/(pow(theta,k)*tgamma(k));
}
*/


/*
float64 CGamma::logpdf(float64 const& x, float64 const& k, float64 const& theta) throw (CGammaException) {
	__checkParameters(x,k,theta);
	return log(pdf(x,k,theta));
}
*/

/*
//OWN IMPLEMENTATIONS REPLACED BY CEPHES
float64 CGamma::Special::__regularizedIncompleteGamma(float64 const& x, float64 const& alpha){
	
	float64 gamma_f = exp(alpha*log(x) - lgamma(alpha) - x);
	if(x < alpha + 1.0) {
		float64 i = alpha;
		float64 tmp_sum = 1.0;
		float64 sum = tmp_sum;
		while(tmp_sum/sum > 1e-10) {
			i++;
			tmp_sum *= x/i;
			sum += tmp_sum;
		}
		return gamma_f*sum/alpha;
	} else { //Solve via evaluting continued fractions

		//BUG in this type of evaluation use implementation below
		//GammaCF cf(x,alpha);
		//float64 old = 1.0-evaluateContinuedFraction(&cf)*gamma_f;

		//Bugfree evaluation of continued fraction
		float64 a=1.0-alpha;
		float64 b=1+x+a;
		float64 pa1 = 1.0;
		float64 pb1 = x;
		float64 pa2 = x + 1.0;
		float64 pb2 = b*x;
		float64 func = pa2/pb2;
		float64 pa,pb,ratio,tmp;
		float64 i = 0;

		while(true) {
			i++;
			a++;
			b += 2.0;
			pa = b*pa2-a*i*pa1;
			pb = b*pb2-a*i*pb1;
			if(pb) {
				ratio = pa/pb;
				tmp = fabs((func-ratio));
				if(tmp<=1e-10*ratio) break;
				func = ratio;
			} else tmp=1.0;
			pa1=pa2;
			pb1=pb2;
			pa2=pa;
			pb2=pb;
			if(i>100) break;//Maximum number if iterations
		}
		return 1.0-func*gamma_f;
	}

}

float64 CGamma::Special::regularizedLowerIncompleteGamma(float64 const& x, float64 const& alpha) {
	if(x <= 0.0 || alpha <= 0.0) return 0.0;
	return __regularizedIncompleteGamma(x,alpha);
}

float64 CGamma::Special::regularizedUpperIncompleteGamma(float64 const& x, float64 const& alpha) {
	if(x <= 0.0 || alpha <= 0.0) return 0.0;
	return 1.0-__regularizedIncompleteGamma(x,alpha);
}

float64 CGamma::Special::complementedIncompleteGamma(float64 const& x, float64 const& alpha) {
	float64 gamma_f = exp(alpha*log(x) - lgamma(alpha) - x);
	if((x <= 0) || ( alpha <= 0))
		return 1.0;
	
	if((x < 1.0) || (x < alpha))
		return 1.0 - regularizedLowerIncompleteGamma(x,alpha);
	
	//continued fraction
	float64 y = 1.0 - alpha;
	float64 z = 1.0 + x + y;
	float64 c = 0.0;
	float64 pkm2 = 1.0;
	float64 qkm2 = x;
	float64 pkm1 = x + 1.0;
	float64 qkm1 = z * x;
	float64 func = pkm1/qkm1;
	float64 i = 0;
	float64 ratio,tmp,pk,qk;

	while(true) {
		i++;
		c += 1.0;
		y += 1.0;
		z += 2.0;
		pk = pkm1 * z  -  pkm2 * y*c;
		qk = qkm1 * z  -  qkm2 * y*c;
		if( qk != 0 ) {
			ratio = pk/qk;
			tmp = fabs( (func - ratio)/ratio );
			if(tmp<=1e-10*ratio) break;
			func = ratio;
		} else {
			tmp = 1.0;
		}
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;
		if( fabs(pk) > 1e32) {
			pkm2 *= 1e-32;
			pkm1 *= 1e-32;
			qkm2 *= 1e-32;
			qkm1 *= 1e-32;
		}
		if(i>100) break; //Max number of iterations
	}
	return func * gamma_f;
}
*/
/*
*deprecated: There is a small bug in the evaluation of the continued fraction!! 
GammaCF::GammaCF(float64 const& x, float64 const& alpha) {
	__x = x;
	__alpha = alpha;
	_no_odd = true;
}	

cf_param GammaCF::cf_fraction(uint const& n) {
	cf_param param;
	if(n==0) {
		param.a_even = __x+1.0-__alpha;
		param.b_even = 0.0;
	} else {
		param.a_even = -n*(n-__alpha);
		param.b_even = 2.0*n + (__x+1.0-__alpha);
	}
	return param;
}
*/ 
