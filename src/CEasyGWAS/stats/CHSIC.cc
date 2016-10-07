#include <cmath>
#include "CHSIC.h"
#include "CEasyGWAS/kernel/CKernels.h"
#include "CEasyGWAS/utils/CMatrixHelper.h"

/*
*CHSIC Constructors
*/
CHSIC::CHSIC() {
    __hsic_score = std::numeric_limits<int>::quiet_NaN();
    __p_value = std::numeric_limits<int>::quiet_NaN();
    __min_permutations = 10000;
    __max_permutations = 100000;
    __seed = 0;
}

CHSIC::CHSIC(int64 const& min, int64 const& max) {
    if(max<min) 
        throw CHSICException("max cannot be smaller than min");
    __hsic_score = std::numeric_limits<int>::quiet_NaN();
    __p_value = std::numeric_limits<int>::quiet_NaN();
    __min_permutations = min;
    __max_permutations = max;
    __seed = 0;
}

CHSIC::CHSIC(int64 const& min, int64 const& max, int64 const& seed) {
    if(max<min) 
        throw CHSICException("max cannot be smaller than min");
    if(seed<0) 
        throw CHSICException("seed cannot be negative");
    __hsic_score = std::numeric_limits<int>::quiet_NaN();
    __p_value = std::numeric_limits<int>::quiet_NaN();
    __min_permutations = min;
    __max_permutations = max;
    __seed = seed;
}

/*
*Compute HSIC
*/
float64 CHSIC::computeHSIC(MatrixXd const& X, MatrixXd const& Y) throw (CHSICException) {
    MatrixXd centeredY = CKernels::centerKernel(Y);
    __hsic_score = pow(X.rows()-1,-2) * (X*centeredY.transpose()).trace();
    return __hsic_score;
}

float64 CHSIC::computePValue(VectorXd const& x, MatrixXd const& Y) throw (CHSICException) {
    MatrixXd centeredY = CKernels::centerKernel(Y);
    MatrixXd X = CKernels::linearKernel(x);
    __hsic_score = pow(X.rows()-1,-2) * (X*centeredY.transpose()).trace();
    VectorXd scores = VectorXd::Zero(__min_permutations);
    MatrixXd permX = permuteVector(x,__min_permutations-1);
    for(int64 i=0; i<__min_permutations-1; i++) {
        VectorXd xp = permX.col(i);
        MatrixXd Xp = CKernels::linearKernel(xp);
        scores(i) = pow(Xp.rows()-1,-2) * (Xp*centeredY.transpose()).trace();
    }
    scores(__min_permutations-1) = __hsic_score;
    //compute p_value
    __p_value = (scores.array() >= __hsic_score).count()/((float64)__min_permutations);

    //check if max_permutations is set
    if(__max_permutations > __min_permutations) {
        if(__p_value<=1.0/__min_permutations) {
            permX = permuteVector(x,__max_permutations);
            VectorXd tscores = VectorXd::Zero(__max_permutations);
            tscores.block(0,0,__min_permutations,1) = scores;
            for(int64 i=__min_permutations; i<__max_permutations-__min_permutations;i++) {
                VectorXd xp = permX.col(i);
                MatrixXd Xp = CKernels::linearKernel(xp);
                tscores(i) = pow(Xp.rows()-1,-2) * (Xp*centeredY.transpose()).trace();
            }
            __p_value = (tscores.array() >= __hsic_score).count()/((float64)__max_permutations);
        }
    }
    return __p_value;
}

/*
*get test statistic
*/
float64 CHSIC::getTestStatistic() {
    return __hsic_score;
}

/*
*get p-value
*/
float64 CHSIC::getPValue() {
    return __p_value;
}
