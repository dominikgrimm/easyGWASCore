#ifndef CMATH_HELPER
#define CMATH_HELPER

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>

#include <cmath>
#include <iostream>

#include "CEasyGWAS/globals.h"
/*
*Compute the pseudo inverse of an matrix
*/
template<typename MatrixT> void pinv(MatrixT const& mat, MatrixT* result, float64 const& epsilon=1E-10) {
	Eigen::JacobiSVD<MatrixT> svd = mat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
	typename MatrixT::Scalar tolerance = epsilon * std::max(mat.cols(),mat.rows()) * svd.singularValues().array().abs().maxCoeff();
	(*result) = svd.matrixV() * MatrixT( (svd.singularValues().array().abs() > tolerance).select(
				svd.singularValues().array().inverse(),0)).asDiagonal()  * svd.matrixU().adjoint();
}	

/*
*Inline method to compute the factorial: n!
*/
inline uint64 factorial(uint64 const& n) {
	return (n==0 || n==1) ? 1 : factorial(n-1) * n;
}

/*
*tbeta returns the value of the Beta Function at x and y
*/
inline float64 tbeta(float64 const& x, float64 const& y) {
	if(x<0 || y<0) {
		std::cout << RED << "ERROR in tbeta: x and y cannot be negative" << BLACK << std::endl;
		throw(1);
	}
	return tgamma(x)*tgamma(y)/tgamma(x+y); 
}

/*
*lbeta returns the log of the Beta Function at x and y
*/
inline float64 lbeta(float64 const& x, float64 const& y) {
	return log(tbeta(x,y));
}

/*
*inverse erf function
*/
inline float64 erfinv(float64 const& p) {
	if(fabs(p)>=1) return NAN;
	if(p==1) return INFINITY;
	if(p==-1) return -INFINITY;
	float64 sign=1.0;
	if(p<0) sign=-1.0;
	float64 s = sqrt(PI)/2.0;
	float64 z = sqrt(-log(1.0-fabs(p)))*sign;
	uint iter=0;
	while(fabs(erf(z)-p)>1e-16*fabs(p)) {
		z = z - (erf(z)-p)*exp(pow(z,2))*s;
		if(iter>1000) {
			logging(WARNING,"Maxiterations reached in erfinv");
			break;
		}
		iter++;
	}
	return z;
}

/*
*isOdd checks if number is odd
*/
inline bool isOdd(uint const& n) {
	if (n%2==0) return false;
	else return true;
}

/*
*Compute a continued fraction a function of type func
*/
//Continued Fraction Class
class CF {
	protected:
		bool _no_odd;
	public:
		bool noOddStep() {return _no_odd;};

		virtual cf_param cf_fraction(uint const&) {
			cf_param param;
			return param;	
		};
};

//inline float64 evaluateContinuedFraction(cf_func func, uint const& max_iterations, float64 const& epsilon=1e-10) {
inline float64 evaluateContinuedFraction(CF* func, uint const& max_iterations, float64 const& epsilon=1e-10) {
	float64 underflow = 1e-30;
	float64 overflow = 1e30;

	cf_param params = func->cf_fraction(0);
	float64 delta = overflow;
	float64 d = params.a_even;
	if(fabs(d)<underflow) d = underflow;
	d = 1.0/d;
	float64 f = d;
	float64 c = 1.0;
	uint n = 1;
	while(fabs(delta-1.0)>=epsilon) {
		params = func->cf_fraction(n);
		//even stepp
		d = params.b_even + params.a_even*d;
		c = params.b_even + params.a_even/c;
		if(fabs(d)<underflow) d = underflow;
		if(fabs(c)<underflow) c = underflow;
		d = 1.0/d;
		delta = c*d;
		f *=delta;
		//Odd step
		if(!func->noOddStep()) {
			d = params.b_odd + params.a_odd*d;
			c = params.b_odd + params.a_odd/c;
			if(fabs(d)<underflow) d = underflow;
			if(fabs(c)<underflow) c = underflow;
			d = 1.0/d;
			delta = c*d;
			f *=delta;
		}
		n++;
		if(n==max_iterations) break;
	}
	return f;
}

inline float64 evaluateContinuedFraction(CF* func) {
	return evaluateContinuedFraction(func,100,1e-10);
}

#endif
