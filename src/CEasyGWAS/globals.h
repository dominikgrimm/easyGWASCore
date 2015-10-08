#ifndef GLOBALS
#define GLOBALS

#include <iostream>

#include "types.h"

#define RED	"\e[0;31m"
#define BLACK	"\e[0m"
#define GREEN	"\e[0;32m"
#define BLUE	"\e[0;34m"
#define YELLOW	"\e[0;33m"

#define LOG 	""
#define FERROR	"FATAL ERROR"
#define ERROR	"ERROR"
#define INFO	"INFORMATION"
#define WARNING "WARNING"
#define DEBUG	"DEBUG"
#define STATUS	"STATUS"
#define ATTENTION	"ATTENTION"

#define PI 3.14159265359f

//Writing log files
#define logging(B,C) {\
	time_t rt; struct tm* ct;\
	time(&rt);\
	ct = localtime(&rt);\
	if((std::string)B==ERROR || (std::string)B==FERROR){\
		std::cerr << RED << "[" << ct->tm_mday << "." << ct->tm_mon+1 << "." << ct->tm_year + 1900\
		     << "," << ct->tm_hour << ":" << ct->tm_min << ":" << ct->tm_sec << "] " << B << " in "\
		     << __FILE__ << " at line " << __LINE__  << ": "\
		     << C << BLACK << "\n";\
	} else if((std::string)B==WARNING) {\
		std::cerr << YELLOW <<  C << BLACK << "\n";\
	} else if((std::string)B==INFO) {\
		std::cout << GREEN <<  C << BLACK << "\n";\
	} else if((std::string)B==STATUS) {\
		std::cout << BLUE <<  C << BLACK << "\n";\
	} else if((std::string)B==ATTENTION) {\
		std::cout << RED <<  C << BLACK << "\n";\
	} else {\
		std::cout << C << "\n";\
	}\
} 


//Continued fraction struct
typedef struct cf_param {
	float64 a_odd;
	float64 b_odd;
	float64 a_even;
	float64 b_even;
} cf_param;

#endif //Globals
