#ifndef CIOPROGRESS_CLASS
#define CIOPROGRESS_CLASS

#include "CEasyGWAS/globals.h"

#include <fstream>

class CIOProgress {
	private:	
		float64 __progress;
		float64 __progress_per;
		float64 __progress_step;
		uint64 __file_size;

		uint64 __getFileSize(std::ifstream&);
	public:
		CIOProgress(std::ifstream&, float64 const&);
		void printProgress(std::ifstream&);
		uint64 getFileSize();
};

#endif //CIOPROGRESS_CLASS
