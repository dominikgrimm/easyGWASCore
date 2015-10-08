#ifndef CGAUSSIAN_CLASS
#define CGAUSSIAN_CLASS

#include "CEasyGWAS/globals.h"

/*
*CGaussian Exception Class
*/
class CGaussianException {
	private:
		std::string __error_msg;
	public:
		CGaussianException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CGaussian Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

/*
* CGaussian Class: Gaussian Distribution class
*/
class CGaussian {
	
	private:
		static void __checkParameters(float64 const&) throw (CGaussianException);

	public:
		static float64 cdf(float64 const&, float64 const&, float64 const&) throw (CGaussianException);
		static float64 logcdf(float64 const&, float64 const&, float64 const&) throw (CGaussianException);
		static float64 pdf(float64 const&, float64 const&, float64 const&) throw (CGaussianException);
		static float64 logpdf(float64 const&, float64 const&, float64 const&) throw (CGaussianException);
		static float64 sf(float64 const&, float64 const&, float64 const&) throw (CGaussianException);
		static float64 isf(float64 const&, float64 const&, float64 const&) throw (CGaussianException);
		static float64 logsf(float64 const&, float64 const&, float64 const&) throw (CGaussianException);
		static float64 ppf(float64 const&, float64 const&, float64 const&) throw (CGaussianException);
};


#endif //CGAUSSIAN_CLASS
