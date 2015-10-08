#include <math.h>

#include "CRootFinding.h"

/*
 *Brent's root finding method from http://en.wikipedia.org/wiki/Brent's_method
 */
float64 CBrentOptimizer::solve(CBrentFunction* func, float64 const& lower, float64 const& upper, float64 const& epsilon, uint const& max_iterations) {
	float64 a = lower;
	float64 b = upper;

	if(b <= a) {
		logging(ERROR,"ERROR: a must be smaller than b");
		throw(-1);
	}
	
	float64 fa = func->evaluate(a);
	float64 fb = func->evaluate(b);

	if(fa*fb >= 0) {
		if(fa<fb) return a;
		else return b;
	}

	float64 c = a;
	float64 fc = fa;
	bool mflag = true;
	uint i = 0;
	float64 s = 0;
	float64 fs = 0;
	float64 d = std::numeric_limits<float64>::infinity();

	while(!(fb==0) && (fabs(a-b)>epsilon)) {
		//Inverse quadratic interpolation
		if((fa!=fc) && (fb!=fc)) 
			s = a * fb * fc / (fa-fb) / (fa-fc) + b * fa * fc / (fb-fa) / (fb-fc) + c * fa * fb / (fc - fa) / (fc - fb);
		else //Secant Rule
			s = b - fb * (b - a) / (fb - fa);
		float64 tmp = (3*a+b)/4;
		if((!(((s>tmp) && (s<b)) || ((s<tmp) && (s>b)))) || (mflag && (fabs(s-b) >= (fabs(b-c)/2))) || (!mflag && (fabs(s-b) >= (fabs(c-d)/2)))) {
			s = (a+b)/2;
			mflag = true;
		} else {
			if ((mflag && (fabs(b-c)<epsilon)) || (!mflag && (fabs(c-d)<epsilon))) {
				s = (a+b)/2;
				mflag = true;
			} else {
				mflag = false;
			}
		}
		fs = func->evaluate(s);
		d = c;
		c = b;
		fc = fb;
		if(fa*fs < 0) {
			b = s;
			fb = fs;
		} else {
			a = s;
			fa = fs;
		}
		if(fabs(fa) < fabs(fb)) {
			float64 t = a;
			a = b;
			b = t;
			t = fa;
			fa = fb;
			fb = t;
		}
		i++;
		if(i>max_iterations) {
			logging(ERROR,"ERROR: Maximum number of iterations exeeded!");
			throw(-1);
		}
	}
	return b;
}
