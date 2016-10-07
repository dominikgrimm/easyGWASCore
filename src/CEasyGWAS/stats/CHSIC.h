#ifndef CHSIC_CLASS
#define CHSIC_CLASS

#include <string>
#include "CEasyGWAS/globals.h"

/*
*CHSIC Exception Class
*/
class CHSICException {
    private:
        std::string __error_msg;

    public:
        CHSICException(std::string const& error_msg) : __error_msg(error_msg) {
            std::cout << RED << "CHSIC Esception: " << error_msg << BLACK << "\n";
        }

        std::string what() {
            return __error_msg;
        }
};

/*
* CHSIC Class
* To compute HSIC and p-values via permutations
*/
class CHSIC {
    
    private:
        float64 __hsic_score;
        float64 __p_value;
        int64 __min_permutations;
        int64 __max_permutations;
        int64 __seed;

    public:
        CHSIC();
        CHSIC(int64 const&, int64 const&);
        CHSIC(int64 const&, int64 const&, int64 const&);

        float64 computeHSIC(MatrixXd const&, MatrixXd const&) throw (CHSICException);
        float64 computePValue(VectorXd const&, MatrixXd const&) throw (CHSICException);
        float64 getTestStatistic();
        float64 getPValue();
};

#endif //CHSIC_CLASS
