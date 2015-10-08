#ifndef CGAMMA_CLASS
#define CGAMMA_CLASS

#include "CEasyGWAS/globals.h"
//#include "CEasyGWAS/utils/CMathHelper.h"

/*
*CGamma Exception Class
*/
class CGammaException {
	private:
		std::string __error_msg;
	public:
		CGammaException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CGamma Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

/*
* CGamma Class: Gamma Distribution class
*/
class CGamma {
	
	private:
		//static void __checkParameters(float64 const&, float64 const&, float64 const&) throw (CGammaException);
		static void __checkParameters(float64 const&, float64 const&) throw (CGammaException);

	public:
		static float64 cdf(float64 const&, float64 const&) throw (CGammaException);
		static float64 logcdf(float64 const&, float64 const&) throw (CGammaException);
		static float64 sf(float64 const&, float64 const&) throw (CGammaException);
		static float64 logsf(float64 const&, float64 const&) throw (CGammaException);
		static float64 pdf(float64 const&, float64 const&) throw (CGammaException);
		static float64 logpdf(float64 const&, float64 const&) throw (CGammaException);
	
};
		
        //static float64 cdf(float64 const&, float64 const&, float64 const&) throw (CGammaException);
		//static float64 logcdf(float64 const&, float64 const&, float64 const&) throw (CGammaException);
		//static float64 sf(float64 const&, float64 const&, float64 const&) throw (CGammaException);
		//static float64 isf(float64 const&, float64 const&, float64 const&) throw (CGammaException);
		//static float64 logsf(float64 const&, float64 const&, float64 const&) throw (CGammaException);
		//static float64 pdf(float64 const&, float64 const&, float64 const&) throw (CGammaException);
		//static float64 logpdf(float64 const&, float64 const&, float64 const&) throw (CGammaException);

    /* OWN IMPLEMENTATIONS
	class Special {

		private:
			//Regularized (Normalized) incomplete gamma function 
			static float64 __regularizedIncompleteGamma(float64 const&, float64 const&);
		public:
			//Regularized (Normalized) lower incomplete gamma function 
			static float64 regularizedLowerIncompleteGamma(float64 const&, float64 const&);
			//Regularized (Normalized) lower incomplete gamma function 
			static float64 regularizedUpperIncompleteGamma(float64 const&, float64 const&);
			static float64 complementedIncompleteGamma(float64 const&, float64 const&);
	};
    */

//Continued Fraction Class
/*
class GammaCF: public CF {
	private:
		float64 __x;
		float64 __alpha;
	public:
		GammaCF(float64 const&, float64 const&);

		cf_param cf_fraction(uint const&);
};
*/

#endif //CGAMMA_CLASS
