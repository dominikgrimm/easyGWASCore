#ifndef CSCONESIO_CLASS
#define CSCONESIO_CLASS

#include "CEasyGWAS/globals.h"
#include "CEasyGWAS/gwas/CGWASData.h"
#include "CEasyGWAS/gwas/CScones.h"

#include <fstream>
#include <string>

/*
*CSconesIO Exception Class
*/
class CSconesIOException {
	private:
		std::string __error_msg;
	public:
		CSconesIOException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CSconesIO Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

class CSconesIO {
	public:
		static void readSparseNetworkFile(std::string const&,
						 	    GWASData*) throw (CSconesIOException);
        static void writeOutput(std::string const&, GWASData const&, VectorXd const&, float64 const&, float64 const&);
        static void writeCMatrix(std::string const&, MatrixXd const&, CSconesSettings const&);
};

#endif //CSCONESIO_CLASS
