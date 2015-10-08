#ifndef CBETA_CLASS
#define CBETA_CLASS

#include "CEasyGWAS/globals.h"
#include "CEasyGWAS/utils/CMathHelper.h"

/*
*CBeta Exception Class
*/
class CBetaException {
	private:
		std::string __error_msg;
	public:
		CBetaException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CBeta Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

/*
* CBeta Class: Beta Distribution class
*/
class CBeta {
	
	private:
		static void __checkParameters(float64 const&, float64 const&, float64 const&) throw (CBetaException);

	public:
		static float64 cdf(float64 const&, float64 const&, float64 const&) throw (CBetaException);
		static float64 logcdf(float64 const&, float64 const&, float64 const&) throw (CBetaException);
		static float64 pdf(float64 const&, float64 const&, float64 const&) throw (CBetaException);
		static float64 logpdf(float64 const&, float64 const&, float64 const&) throw (CBetaException);
		static float64 sf(float64 const&, float64 const&, float64 const&) throw (CBetaException);
		static float64 isf(float64 const&, float64 const&, float64 const&) throw (CBetaException);
		static float64 logsf(float64 const&, float64 const&, float64 const&) throw (CBetaException);
		static float64 ppf(float64 const&, float64 const&, float64 const&) throw (CBetaException);

	/*
    class Special {
		public:
			static float64 regularizedIncompleteBetaFunction(float64 const&, float64 const&, float64 const&);
	};
    */
};

/*
//Continued Fraction Class
class BetaCF: public CF {
	private:
		float64 __x;
		float64 __alpha;
		float64 __beta;
	public:
		BetaCF(float64 const&, float64 const&, float64 const&);

		cf_param cf_fraction(uint const&);
};
*/


#endif //CBETA_CLASS
