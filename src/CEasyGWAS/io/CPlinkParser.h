#ifndef CPLINKPARSER_CLASS
#define CPLINKPARSER_CLASS

#include "CEasyGWAS/globals.h"
#include "CEasyGWAS/utils/StringHelper.h"
#include "CEasyGWAS/gwas/CGWASData.h"

#include <fstream>
#include <vector>
#include <string>
#include <map>

/*
*CPlinkParser Exception Class
*/
class CPlinkParserException {
	private:
		std::string __error_msg;
	public:
		CPlinkParserException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CPlinkParser Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};


class CPlinkParser {
	private:

		static std::map<std::string,char> __init_map() {
			std::map<std::string,char> iupac_map;
			iupac_map["AA"] = 'A';
			iupac_map["GG"] = 'G';
			iupac_map["TT"] = 'T';
			iupac_map["CC"] = 'C';
			iupac_map["AG"] = 'R';
			iupac_map["GA"] = 'R';
			iupac_map["CT"] = 'Y';
			iupac_map["TC"] = 'Y';
			iupac_map["GC"] = 'S';
			iupac_map["CG"] = 'S';
			iupac_map["AT"] = 'W';
			iupac_map["TA"] = 'W';
			iupac_map["GT"] = 'K';
			iupac_map["TG"] = 'K';
			iupac_map["AC"] = 'M';
			iupac_map["CA"] = 'M';
            //iupac_map["--"] = 'D';
			return iupac_map;
		}

		const static std::map<std::string,char> __iupac_map;
	public:
		static void readPEDFile(std::string const&,
					GWASData*) 
					throw (CPlinkParserException);
		static void readMAPFile(std::string const&,
					GWASData*) 
					throw (CPlinkParserException);
		static void readPhenotypeFile(std::string const&,
					      GWASData*) 
					      throw (CPlinkParserException);
	
};

#endif //CPLINKPARSER_CLASS
