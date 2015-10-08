#ifndef CLOGGING_CLASS
#define CLOGGING_CLASS

#include "CEasyGWAS/globals.h"

#include <string>
#include <fstream>

class CLogging {

	private:
		std::string __filename;
		std::ofstream __ofs;
		bool file_logging;

	public:
		CLogging();
		CLogging(std::string const&);
		~CLogging();

		void log(std::string const&, std::string const&);

};

#endif //CLOGGING_CLASS
