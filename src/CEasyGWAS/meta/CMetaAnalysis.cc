#include "math.h"

#include "CEasyGWAS/meta/CMetaAnalysis.h"
#include "CEasyGWAS/stats/CGaussian.h"
#include "CEasyGWAS/stats/CChi2.h"

namespace CMetaAnalysis {

/*
*Method for combining p-values
*/
float64 CombinedPvalues::FisherMethod(VectorXd const& pvalues) {
	float64 chi2 = -2.0 * pvalues.array().log().sum();
	return CChi2::sf(chi2,2*pvalues.rows());
}

float64 CombinedPvalues::StoufferZ(VectorXd const& pvalues) {
	return StoufferZ(pvalues,true);
}

float64 CombinedPvalues::StoufferZ(VectorXd const& pvalues, bool const& return_pval) {
	VectorXd z_scores = VectorXd::Ones(pvalues.rows());
	for(uint i=0; i<pvalues.rows();i++) {
		z_scores(i) = CGaussian::ppf(1.0-pvalues(i),0,1);
	}
	float64 zscore = z_scores.sum()/sqrt(pvalues.rows());
	if(return_pval) {
		return StoufferPval(zscore);
	} else {
		return zscore;
	}
}

float64 CombinedPvalues::StoufferZWeighted(VectorXd const& pvalues,
					   VectorXd const& weights) {
	VectorXd effects = VectorXd::Ones(pvalues.rows());
	return StoufferZWeighted(pvalues,weights,effects);
}

float64 CombinedPvalues::StoufferZWeighted(VectorXd const& pvalues,
					   VectorXd const& weights,//for example sample_size
					   VectorXd const& effect_direction) {
	VectorXd z_scores = VectorXd::Ones(pvalues.rows());
	for(uint i=0; i<pvalues.rows();i++) {
		float64 sign = 1.0;
		if(effect_direction(i)<0) sign = -1.0;
		z_scores(i) = CGaussian::ppf(1.0-pvalues(i)/2.0,0,1) * sign * weights(i);
	}
	float64 zscore = z_scores.sum()/sqrt(weights.array().pow(2).sum());
	return zscore;
}

/*
*Retrive pvalue for stouffer z-score
*/
float64 CombinedPvalues::StoufferPval(float64 const& zscore) {
	return StoufferPval(zscore,false);
}

float64 CombinedPvalues::StoufferPval(float64 const& zscore, bool const& two_tailed) {
	if(two_tailed) {
		return 2.0 * (CGaussian::sf(fabs(zscore),0,1));
	} else {
		return (CGaussian::sf(fabs(zscore),0,1));
	}
}

/*
*Effect Model Class: Base class for all Effect Models
*/	
EffectModel::EffectModel(VectorXd const& effects, VectorXd const& v) {
	_effects = effects;
	_variance_effects = v;
}

void EffectModel::setEffects(VectorXd const& effects) {
	_effects = effects;
}

void EffectModel::setVarianceEffects(VectorXd const& v) {
	_variance_effects = v;
}

VectorXd EffectModel::getWeights() {
	return _weights;
}

float64 EffectModel::getZvalue() {
	return _z;
}

float64 EffectModel::getPvalue() {
	return getPvalue(false);
}

float64 EffectModel::getPvalue(bool const& twosided) {
	if(twosided)
		return 2.0*(1-CGaussian::cdf(fabs(_z),0,1));
	else
		return 1-CGaussian::cdf(fabs(_z),0,1);
}

float64 EffectModel::getLowerLimit() {
	return _lower_limit;
}

float64 EffectModel::getUpperLimit() {
	return _upper_limit;
}

float64 EffectModel::getWeightedMean() {
	return _weighted_mean;
}

float64 EffectModel::getVarianceSummary(){
	return _variance_summary;
}

float64 EffectModel::getSESummary() {
	return _se_summary;
}

/*
 *FixedEffectModel Class
 */
FixedEffectModel::FixedEffectModel(VectorXd const& effects, VectorXd const& v) {
	_effects = effects;
	_variance_effects = v;
}


void FixedEffectModel::process() {
	_weights = 1.0/_variance_effects.array();
	_weighted_mean = (_weights.array()*_effects.array()).sum()/_weights.sum();
	_variance_summary = 1.0/_weights.sum();
	_se_summary = sqrt(_variance_summary);
	_lower_limit = _weighted_mean-1.96*_se_summary;
	_upper_limit = _weighted_mean+1.96*_se_summary;
	_z = _weighted_mean/_se_summary;
}

/*
 *RandomEffectModel Class
 */
RandomEffectModel::RandomEffectModel(VectorXd const& effects, VectorXd const& v) {
	_effects = effects;
	_variance_effects = v;
}


void RandomEffectModel::process() {
	_weights = 1.0/_variance_effects.array();
	
	__Q = (_weights.array()*_effects.array().pow(2).array()).sum() -
	      pow((_weights.array()*_effects.array()).sum(),2)/(_weights.sum());
	__C = _weights.sum()-_weights.array().pow(2).sum()/_weights.sum();
	__T_2 = (__Q-(_effects.rows()-1.0))/__C;
	if(__T_2<0) __T_2 = 0.0;
	VectorXd variance = _variance_effects.array() + __T_2; 
	_weights = 1.0/variance.array();

	_weighted_mean = (_weights.array()*_effects.array()).sum()/_weights.sum();
	_variance_summary = 1.0/_weights.sum();
	_se_summary = sqrt(_variance_summary);
	_lower_limit = _weighted_mean-1.96*_se_summary;
	_upper_limit = _weighted_mean+1.96*_se_summary;
	_z = _weighted_mean/_se_summary;
}

float64 RandomEffectModel::getQ() {
	return __Q;
}

float64 RandomEffectModel::getC() {
	return __C;
}

float64 RandomEffectModel::getT2() {
	return __T_2;
}

};//END Namespace
