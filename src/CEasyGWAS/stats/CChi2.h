#ifndef CCHI2_CLASS
#define CCHI2_CLASS

#include "CEasyGWAS/globals.h"


/*
*CChi2 Exception Class
*/
class CChi2Exception {
	private:
		std::string __error_msg;
	public:
		CChi2Exception(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CChi2 Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

/*
* CChi2 Class: Chi2 Distribution class
*/
class CChi2 {
	
	private:
		static void __checkParameters(float64 const&) throw (CChi2Exception);

	public:
		static float64 cdf(float64 const&, float64 const&) throw (CChi2Exception);
		static float64 logcdf(float64 const&, float64 const&) throw (CChi2Exception);
		static float64 pdf(float64 const&, float64 const&) throw (CChi2Exception);
		static float64 logpdf(float64 const&, float64 const&) throw (CChi2Exception);
		static float64 sf(float64 const&, float64 const&) throw (CChi2Exception);
		static float64 isf(float64 const&, float64 const&) throw (CChi2Exception);
		static float64 logsf(float64 const&, float64 const&) throw (CChi2Exception);
};


#endif //CCHI2_CLASS
