#include "CCrossValidation.h"
#include "math.h"
#include "StringHelper.h"

#include <iostream>

CCrossValidation::CCrossValidation() {
	__seed = 0;
	__k = 0;
	srand(__seed);
}

CCrossValidation::CCrossValidation(float64 const& seed) {
	__seed = seed;
	__k = 0;
	srand(__seed);
}

VectorXd CCrossValidation::getTrainingIndices(uint const& k) const throw (CCrossValidationException) {
	if(k>__trainingData.size()) throw CCrossValidationException("No TrainingIndices for fold k" + StringHelper::to_string<uint>(k));
	return __trainingData[k];
}

VectorXd CCrossValidation::getTestingIndices(uint const& k) const throw (CCrossValidationException) {
	if(k>__testingData.size()) throw CCrossValidationException("No TestingIndices for fold k" + StringHelper::to_string<uint>(k));
	return __testingData[k];
}

uint CCrossValidation::size() {
	return __k;
}	

void CCrossValidation::train_test_split(uint64 const& n, float64 const& ratio) {
	__ratio = ratio;
	__n = n;
	__k = 1;
	VectorXd tmp = VectorXd::LinSpaced(__n,0,__n-1);
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(__n);
	perm.setIdentity();
	std::random_shuffle(perm.indices().data(),perm.indices().data()+perm.indices().size());
	tmp = perm*tmp;
	uint sratio = ceil(__ratio*__n);
	__trainingData.push_back(tmp.segment(0,__n-sratio));
	__testingData.push_back(tmp.segment(__n-sratio,sratio));
}

void CCrossValidation::kFold(uint const& k, uint64 const& n) {
	__ratio = 0.0;
	__n = n;
	__k = k;
	VectorXd indices = VectorXd::LinSpaced(__n,0,__n-1);
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(__n);
	perm.setIdentity();
	std::random_shuffle(perm.indices().data(),perm.indices().data()+perm.indices().size());
	indices = perm*indices;
	float64 split = floor((float)__n/((float)__k));
	for(uint i=0; i<__k;i++) {
		if(! (i==__k-1)) {
			__testingData.push_back(indices.segment(i*split,split));
			VectorXd tmp = VectorXd::Zero(__n - __testingData[i].rows());
			if(i==0){
				tmp << indices.segment(i*split+split,__n-(i*split+split));
			}else{
				tmp << indices.segment(0,i*split), indices.segment(i*split+split,__n-(i*split+split));
			}
			__trainingData.push_back(tmp);
		} else {
			__testingData.push_back(indices.segment(i*split,__n-i*split));
			VectorXd tmp = VectorXd::Zero(__n-__testingData[i].rows());
			tmp << indices.segment(0,i*split);
			__trainingData.push_back(tmp);
		
		}
	}
}

void CCrossValidation::ShuffleSplit(uint64 const& n, uint64 const& k, float64 const& ratio) {
	__ratio = ratio;
	__n = n;
	__k = k;
	uint sratio = ceil(__ratio*__n);
	for(uint64 i=0; i<__k;i++) {
        VectorXd tmp = VectorXd::LinSpaced(__n,0,__n-1);
	    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(__n);
	    perm.setIdentity();
	    std::random_shuffle(perm.indices().data(),perm.indices().data()+perm.indices().size());
	    tmp = perm*tmp;
	    __trainingData.push_back(tmp.segment(0,__n-sratio));
	    __testingData.push_back(tmp.segment(__n-sratio,sratio));
    }
}
