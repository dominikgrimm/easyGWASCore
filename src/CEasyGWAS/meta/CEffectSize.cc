#include "CEffectSize.h"

#include "CEasyGWAS/stats/CStats.h"

namespace CEffectSize {

const uint MeanEffectSize::equalPopulationStandardDeviation;
const uint MeanEffectSize::nonEqualPopulationStandardDeviation;
const uint MeanEffectSize::independetGroups;
const uint MeanEffectSize::matchedGroups;
const uint MeanEffectSize::standardizedCorrected;
	
MeanEffectSize::MeanEffectSize(VectorXd const& x1, VectorXd const& x2,
			       VectorXd const& se1, VectorXd const& se2,
			       VectorXd const& n1, VectorXd const& n2) {
	__x1 = x1;
	__x2 = x2;
	__se1 = se1;
	__se2 = se2;
	__n1 = n1;
	__n2 = n2;
}

VectorXd MeanEffectSize::getD() {
	__D = __x1-__x2;
	return __D;
}

VectorXd MeanEffectSize::getVariance() {
	if(__standardized==MeanEffectSize::standardizedCorrected) 
		return getVariance(__standardized);
	else
		return getVariance(MeanEffectSize::nonEqualPopulationStandardDeviation);
}

VectorXd MeanEffectSize::getVariance(uint const& mode) {
	if(mode==MeanEffectSize::equalPopulationStandardDeviation) {
		VectorXd S1_2 = __se1.array().pow(2);
		VectorXd S2_2 = __se2.array().pow(2);
		__Swithin = (((__n1.array()-1.0).array()*S1_2.array()+(__n2.array()-1.0)*S2_2.array()).array()/((__n1+__n2).array()-2.0)).array().sqrt();
		__variance = __Swithin.array().pow(2).array()*(__n1+__n2).array()/(__n1+__n2).array();
	} else if (mode==MeanEffectSize::nonEqualPopulationStandardDeviation) {
		VectorXd S1_2 = __se1.array().pow(2);
		VectorXd S2_2 = __se2.array().pow(2);
		__Swithin = (((__n1.array()-1.0).array()*S1_2.array()+(__n2.array()-1.0)*S2_2.array()).array()/((__n1+__n2).array()-2.0)).array().sqrt();
		__variance = S1_2.array()/__n1.array() + S2_2.array()/__n2.array();
	} else if (mode==MeanEffectSize::standardizedCorrected) {
		__variance = __J.array().pow(2).array()*__variance.array();
	} else {
		//TODO implement matching groups
	}
	return __variance;
}

VectorXd MeanEffectSize::getSE() {
	return __variance.array().sqrt();
}

VectorXd MeanEffectSize::getJ() {
	return __J;
}

VectorXd MeanEffectSize::getHedgesG() {
	VectorXd S1_2 = __se1.array().pow(2);
	VectorXd S2_2 = __se2.array().pow(2);
	__Swithin = (((__n1.array()-1.0).array()*S1_2.array()+(__n2.array()-1.0)*S2_2.array()).array()/((__n1+__n2).array()-2.0)).array().sqrt();
	__D = __x1-__x2;
	__D.array() /= __Swithin.array();
	__variance = (__n1+__n2).array()/(__n1.array()*__n2.array()).array() + 
		     __D.array().pow(2).array()/(2*(__n1+__n2).array()).array();
	//This is an approximation of J
	__J = 1.0 - (3.0/(4.0*((__n1+__n2).array()-2.0)-1.0).array());
	//This is the exact J
	//VectorXd a = (__n1+__n2).array()-2.0;
	//__J.resize(a.rows());
	//for(uint i=0; i<a.rows(); i++) {
	//	__J(i) = tgamma(a(i)/2.0) / (sqrt(a(i)/2.0)*tgamma((a(i)-1.0)/2.0));
	//}
	__g = __J.array() * __D.array();
	__standardized = MeanEffectSize::standardizedCorrected;
	return __g;
}

};
