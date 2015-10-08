#ifndef CROOTFINDING_CLASS
#define CROOTFINDING_CLASS

#include "CEasyGWAS/globals.h"

/*
*Brent's methods for root finding
*/

class CBrentFunction {
	public:
		virtual float64 evaluate(float64 const& x) {
			return x*0.0;
		}
};

class CBrentOptimizer {
	public:
		static float64 solve(CBrentFunction*, float64 const&, float64 const&, float64 const&, uint const&);
};

#endif //CROOTFINDING_CLASS
