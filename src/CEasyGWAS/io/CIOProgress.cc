#include "CIOProgress.h"

//#include <iostream>
#include <stdio.h>

CIOProgress::CIOProgress(std::ifstream& ifs, float64 const& step) {
	__progress_step = step;
	__progress = step;
	__progress_per = 0.0;
	__file_size = __getFileSize(ifs);
}

uint64 CIOProgress::__getFileSize(std::ifstream& ifs) {
	ifs.seekg(0,std::ios::end);
	uint64 r = ifs.tellg();
	ifs.seekg(0,std::ios::beg);
	__progress_per = 100.0/((float64)r);
	return r;
}

void CIOProgress::printProgress(std::ifstream& ifs) {
	float64 current_pos = ifs.tellg();
	float64 current_per = __progress_per*current_pos;
    if(current_per >= __progress) {
	    __progress += __progress_step;
		printf("\rProcessing file: %.2f %%",current_per);
        fflush(stdout);
		//std::cout.precision(4);
        //std::cout << "\rProcessing file: " << current_per << " %";
	    //std::cout.flush();
	}
}

uint64 CIOProgress::getFileSize() {
	return __file_size;
}
