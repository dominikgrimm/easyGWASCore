#ifndef CCROSSVALIDATION_HELPER
#define CCROSSVALIDATION_HELPER

#include "CEasyGWAS/globals.h"

#include <vector>

/*
*CCrossValidation Exception Class
*/
class CCrossValidationException {
	private:
		std::string __error_msg;
	public:
		CCrossValidationException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CCrossValidation Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

class CCrossValidation {
	private:
		float64 __seed;
		float64 __ratio;
		uint __k;
		uint64 __n;

		std::vector<VectorXd> __trainingData;
		std::vector<VectorXd> __testingData;

	public:
		CCrossValidation();
		CCrossValidation(float64 const&);

		//Split data into train and test set
		void train_test_split(uint64 const&, float64 const&);
		//Split data into k folds
		void kFold(uint const&, uint64 const&);
		void ShuffleSplit(uint64 const&, uint64 const&, float64 const&);
		//Split data into k stratified folds
		void stratifiedKFold(uint const&, VectorXd const&);

		VectorXd getTrainingIndices(uint const&) const throw (CCrossValidationException);
		VectorXd getTestingIndices(uint const&) const throw (CCrossValidationException);

		uint size();
};

#endif //CCROSSVALIDATION_HELPER
