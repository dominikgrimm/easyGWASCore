#ifndef CFISHERF_CLASS
#define CFISHERF_CLASS

#include "CEasyGWAS/globals.h"
#include "CEasyGWAS/utils/CMathHelper.h"

/*
*CFisherF Exception Class
*/
class CFisherFException {
	private:
		std::string __error_msg;
	public:
		CFisherFException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CFisherF Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

/*
* CFisherF Class: FisherF Distribution class
*/
class CFisherF {
	
	private:
		static void __checkParameters(float64 const&, int const&, int const&) throw (CFisherFException);

	public:
		static float64 cdf(float64 const&, int const&, int const&) throw (CFisherFException);
		static float64 logcdf(float64 const&, int const&, int const&) throw (CFisherFException);
		static float64 pdf(float64 const&, int const&, int const&) throw (CFisherFException);
		static float64 ppf(float64 const&, int const&, int const&) throw (CFisherFException);
		static float64 logpdf(float64 const&, int const&, int const&) throw (CFisherFException);
		static float64 sf(float64 const&, int const&, int const&) throw (CFisherFException);
		static float64 isf(float64 const&, int const&, int const&) throw (CFisherFException);
		static float64 logsf(float64 const&, int const&, int const&) throw (CFisherFException);

};
	

#endif //CFISHERF_CLASS
