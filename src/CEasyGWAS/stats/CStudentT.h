#ifndef CSTUDENTT_CLASS
#define CSTUDENTT_CLASS

#include "CEasyGWAS/globals.h"
#include "CEasyGWAS/utils/CMathHelper.h"

/*
*CStudentT Exception Class
*/
class CStudentTException {
	private:
		std::string __error_msg;
	public:
		CStudentTException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CStudentT Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

/*
* CStudentT Class: StudentT Distribution class
*/
class CStudentT {
	
	private:
		static void __checkParameters(int const&) throw (CStudentTException);

	public:
		static float64 cdf(float64 const&, int const&) throw (CStudentTException);
		static float64 logcdf(float64 const&, int const&) throw (CStudentTException);
		static float64 pdf(float64 const&, int const&) throw (CStudentTException);
		static float64 logpdf(float64 const&, int const&) throw (CStudentTException);
		static float64 sf(float64 const&, int const&) throw (CStudentTException);
		static float64 isf(float64 const&, int const&) throw (CStudentTException);
		static float64 ppf(float64 const&, int const&) throw (CStudentTException);
		static float64 logsf(float64 const&, int const&) throw (CStudentTException);

};
	

#endif //CSTUDENTT_CLASS
