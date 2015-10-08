#include "CEasyGWAS/gwas/CMetaGWAS.h"
#include "CEasyGWAS/utils/StringHelper.h"
#include "CEasyGWAS/meta/CMetaAnalysis.h"

#include "math.h"

CMetaGWAS::CMetaGWAS() {
	__study_counts = 0;
}

void CMetaGWAS::addPValuesStudy(VectorXd const& p_values,std::vector<std::string> const& chromosomes, VectorXd const& positions) {
	std::string id = ""; 
	for(int64 i=0; i<p_values.rows();i++) {
		id = chromosomes[i] + "_" + StringHelper::to_string<uint64>(positions[i]);
		__study_it = __study_map.find(id);
		CMetaData data;
		data.p_value = p_values[i];
		data.n_samples = NAN;
		data.beta = NAN;
		data.beta_se = NAN;
		data.position = positions[i];
		data.chromosome = chromosomes[i];
		if(__study_it != __study_map.end()) {
			std::vector<CMetaData>* tmp = &__study_it->second;
			tmp->push_back(data);
		} else {
			std::vector<CMetaData> new_v;
			new_v.push_back(data);
			__study_map[id] = new_v;
		}
	}
	__study_counts++;
}

void CMetaGWAS::addWeightedPValuesStudy(VectorXd const& p_values, 
	                            		std::vector<std::string> const& chromosomes,
			                            VectorXd const& positions,
	    		                        uint64 const& samples) {
	std::string id = ""; 
	for(int64 i=0; i<p_values.rows();i++) {
		id = chromosomes[i] + "_" + StringHelper::to_string<uint64>(positions[i]);
		__study_it = __study_map.find(id);
		CMetaData data;
		data.p_value = p_values[i];
		data.n_samples = sqrt(samples);
		data.beta = NAN;
		data.beta_se = NAN;
		data.position = positions[i];
		data.chromosome = chromosomes[i];
		if(__study_it != __study_map.end()) {
			std::vector<CMetaData>* tmp = &__study_it->second;
			tmp->push_back(data);
		} else {
			std::vector<CMetaData> new_v;
			new_v.push_back(data);
			__study_map[id] = new_v;
		}
	}
	__study_counts++;
}

void CMetaGWAS::addEffectSizeStudy(std::vector<std::string> const& chromosomes, 
			                       VectorXd const& positions,
	                        	   VectorXd const& betas,
			                       VectorXd const& betas_se) {
	std::string id = ""; 
	for(int64 i=0; i<betas.rows();i++) {
		id = chromosomes[i] + "_" + StringHelper::to_string<uint64>(positions[i]);
		__study_it = __study_map.find(id);
		CMetaData data;
		data.p_value = NAN;
		data.n_samples = NAN;
        data.beta = betas[i];
		data.beta_se = betas_se[i];
		data.position = positions[i];
		data.chromosome = chromosomes[i];
		if(__study_it != __study_map.end()) {
			std::vector<CMetaData>* tmp = &__study_it->second;
			tmp->push_back(data);
		} else {
			std::vector<CMetaData> new_v;
			new_v.push_back(data);
			__study_map[id] = new_v;
		}
	}
	__study_counts++;
}

void CMetaGWAS::addStudy(VectorXd const& p_values, 
			 std::vector<std::string> const& chromosomes, 
			 std::vector<uint64> const& positions) {
	addStudy(p_values,chromosomes,positions,VectorXd::Ones(p_values.rows()));
}

void CMetaGWAS::addStudy(VectorXd const& p_values, 
			 std::vector<std::string> const& chromosomes,
			 std::vector<uint64> const& positions,
			 VectorXd const& samples) {
	addStudy(p_values,chromosomes,positions,samples,VectorXd::Ones(p_values.rows()),VectorXd::Ones(p_values.rows()));
}

void CMetaGWAS::addStudy(VectorXd const& p_values, 
			 std::vector<std::string> const& chromosomes, 
			 std::vector<uint64> const& positions,
			 VectorXd const& samples,
			 VectorXd const& betas,
			 VectorXd const& betas_se) {
	std::string id = "";
	for(int64 i=0; i<p_values.rows();i++) {
		id = chromosomes[i] + "_" + StringHelper::to_string<uint64>(positions[i]);
		__study_it = __study_map.find(id);
		CMetaData data;
		data.p_value = p_values[i];
		data.n_samples = sqrt(samples[i]);
		data.beta = betas[i];
		data.beta_se = betas_se[i];
		data.position = positions[i];
		data.chromosome = chromosomes[i];
		if(__study_it != __study_map.end()) {
			std::vector<CMetaData>* tmp = &__study_it->second;
			tmp->push_back(data);
		} else {
			std::vector<CMetaData> new_v;
			new_v.push_back(data);
			__study_map[id] = new_v;
		}
	}
	if(__study_counts==0) {
		__study_counts = 2;
	} else {
		__study_counts++;
	}
}

void CMetaGWAS::performFisherMethod() {
	performFisherMethod(true);
}

void CMetaGWAS::performFisherMethod(bool const& ignore_single_studies) {
	uint64 size = 0;
	if(ignore_single_studies==true) {
		for(__study_it = __study_map.begin(); __study_it!=__study_map.end(); __study_it++)
			if(__study_it->second.size()>1) size++;
	} else { 
		size = __study_map.size();
	}
	__p_values = VectorXd::Zero(size);
	//Clear vectors and reinitialize
	std::vector<std::string>().swap(__chromosomes);
	std::vector<uint64>().swap(__positions);
	__chromosomes.resize(size);
	__positions.resize(size);
	uint64 counter=0;
	for(__study_it = __study_map.begin(); __study_it!=__study_map.end(); __study_it++) {
		//std::vector<float64> pvals_v = __study_it->second;
		std::vector<CMetaData> data_v = __study_it->second;
		if(ignore_single_studies==true) {
			if(data_v.size()==1) continue;
		}
		VectorXd pvals = VectorXd::Zero(data_v.size());
		for(uint i=0; i<data_v.size();i++) {
			pvals[i] = data_v[i].p_value;
		}
		__p_values[counter] = CMetaAnalysis::CombinedPvalues::FisherMethod(pvals);
		__chromosomes[counter] = data_v[0].chromosome;
		__positions[counter] = data_v[0].position;
		counter++;
	}
}	

void CMetaGWAS::performStoufferZ() {
	performStoufferZ(true,false);
}

void CMetaGWAS::performStoufferZ(bool const& two_tailed) {
	performStoufferZ(true,two_tailed);
}

void CMetaGWAS::performStoufferZ(bool const& ignore_single_studies, bool const& two_tailed) {
	uint64 size = 0;
	if(ignore_single_studies==true) {
		for(__study_it = __study_map.begin(); __study_it!=__study_map.end(); __study_it++)
			if(__study_it->second.size()>1) size++;
	} else { 
		size = __study_map.size();
	}
	__p_values = VectorXd::Zero(size);
	__z_values = VectorXd::Zero(size);
	//Clear vectors and reinitialize
	std::vector<std::string>().swap(__chromosomes);
	std::vector<uint64>().swap(__positions);
	__chromosomes.resize(size);
	__positions.resize(size);
	uint64 counter=0;
	for(__study_it = __study_map.begin(); __study_it!=__study_map.end(); __study_it++) {
		std::vector<CMetaData> data_v = __study_it->second;
		if(ignore_single_studies==true) {
			if(data_v.size()==1) continue;
		}
		VectorXd pvals = VectorXd::Zero(data_v.size());
		for(uint i=0; i<data_v.size();i++) {
			pvals[i] = data_v[i].p_value;
		}
		__z_values[counter] = CMetaAnalysis::CombinedPvalues::StoufferZ(pvals,false);
		__p_values[counter] = CMetaAnalysis::CombinedPvalues::StoufferPval(__z_values[counter],two_tailed);
		__chromosomes[counter] = data_v[0].chromosome;
		__positions[counter] = data_v[0].position;
		counter++;
	}
}	

void CMetaGWAS::performStoufferZWeighted() {
	performStoufferZWeighted(true,false);
}

void CMetaGWAS::performStoufferZWeighted(bool const& two_tailed) {
	performStoufferZWeighted(true,two_tailed);
}

void CMetaGWAS::performStoufferZWeighted(bool const& ignore_single_studies, bool const& two_tailed) {
	uint64 size = 0;
	if(ignore_single_studies==true) {
		for(__study_it = __study_map.begin(); __study_it!=__study_map.end(); __study_it++)
			if(__study_it->second.size()>1) size++;
	} else { 
		size = __study_map.size();
	}
	__p_values = VectorXd::Zero(size);
	__z_values = VectorXd::Zero(size);
	//Clear vectors and reinitialize
	std::vector<std::string>().swap(__chromosomes);
	std::vector<uint64>().swap(__positions);
	__chromosomes.resize(size);
	__positions.resize(size);
	uint64 counter=0;
	for(__study_it = __study_map.begin(); __study_it!=__study_map.end(); __study_it++) {
		std::vector<CMetaData> data_v = __study_it->second;
		if(ignore_single_studies==true) {
			if(data_v.size()==1) continue;
		}
		VectorXd pvals = VectorXd::Zero(data_v.size());
		VectorXd samples = VectorXd::Zero(data_v.size());
		for(uint i=0; i<data_v.size();i++) {
			pvals[i] = data_v[i].p_value;
			samples[i] = data_v[i].n_samples;
		}
		__z_values[counter] = CMetaAnalysis::CombinedPvalues::StoufferZWeighted(pvals,samples);
		__p_values[counter] = CMetaAnalysis::CombinedPvalues::StoufferPval(__z_values[counter],two_tailed);
		__chromosomes[counter] = data_v[0].chromosome;
		__positions[counter] = data_v[0].position;
		counter++;
	}
}	

void CMetaGWAS::performFixedEffectModel() {
	performFixedEffectModel(true,false);
}

void CMetaGWAS::performFixedEffectModel(bool const& two_tailed) {
	performFixedEffectModel(true,two_tailed);
}

void CMetaGWAS::performFixedEffectModel(bool const& ignore_single_studies, bool const& two_tailed) {
	uint64 size = 0;
	if(ignore_single_studies==true) {
		for(__study_it = __study_map.begin(); __study_it!=__study_map.end(); __study_it++)
			if(__study_it->second.size()>1) size++;
	} else { 
		size = __study_map.size();
	}
	__p_values = VectorXd::Zero(size);
	__z_values = VectorXd::Zero(size);
	__weighted_means = VectorXd::Zero(size);
	__variances = VectorXd::Zero(size);
	__lower_limits = VectorXd::Zero(size);
	__upper_limits = VectorXd::Zero(size);
	//Clear vectors and reinitialize
	std::vector<std::string>().swap(__chromosomes);
	std::vector<uint64>().swap(__positions);
	__chromosomes.resize(size);
	__positions.resize(size);
	uint64 counter=0;
	for(__study_it = __study_map.begin(); __study_it!=__study_map.end(); __study_it++) {
		std::vector<CMetaData> data_v = __study_it->second;
		if(ignore_single_studies==true) {
			if(data_v.size()==1) continue;
		}
		VectorXd betas = VectorXd::Zero(data_v.size());
		VectorXd betas_se = VectorXd::Zero(data_v.size());
		for(uint i=0; i<data_v.size();i++) {
			betas[i] = data_v[i].beta;
			betas_se[i] = data_v[i].beta_se;
		}
		CMetaAnalysis::FixedEffectModel model(betas,betas_se);
		model.process();
		__p_values[counter] = model.getPvalue(two_tailed);
		__z_values[counter] = model.getZvalue();
		__weighted_means[counter] = model.getWeightedMean();
	        __variances[counter] = model.getVarianceSummary();
		__lower_limits[counter] = model.getLowerLimit();
		__upper_limits[counter] = model.getUpperLimit();
		__chromosomes[counter] = data_v[0].chromosome;
		__positions[counter] = data_v[0].position;
		counter++;
	}
}	

void CMetaGWAS::performRandomEffectModel() {
	performRandomEffectModel(true,false);
}

void CMetaGWAS::performRandomEffectModel(bool const& two_tailed) {
	performRandomEffectModel(true,two_tailed);
}

void CMetaGWAS::performRandomEffectModel(bool const& ignore_single_studies, bool const& two_tailed) {
	uint64 size = 0;
	if(ignore_single_studies==true) {
		for(__study_it = __study_map.begin(); __study_it!=__study_map.end(); __study_it++)
			if(__study_it->second.size()>1) size++;
	} else { 
		size = __study_map.size();
	}
	__p_values = VectorXd::Zero(size);
	__z_values = VectorXd::Zero(size);
	__weighted_means = VectorXd::Zero(size);
	__variances = VectorXd::Zero(size);
	__lower_limits = VectorXd::Zero(size);
	__upper_limits = VectorXd::Zero(size);
	__q = VectorXd::Zero(size);
	__c = VectorXd::Zero(size);
	__t2 = VectorXd::Zero(size);
	//Clear vectors and reinitialize
	std::vector<std::string>().swap(__chromosomes);
	std::vector<uint64>().swap(__positions);
	__chromosomes.resize(size);
	__positions.resize(size);
	uint64 counter=0;
	for(__study_it = __study_map.begin(); __study_it!=__study_map.end(); __study_it++) {
		std::vector<CMetaData> data_v = __study_it->second;
		if(ignore_single_studies==true) {
			if(data_v.size()==1) continue;
		}
		VectorXd betas = VectorXd::Zero(data_v.size());
		VectorXd betas_se = VectorXd::Zero(data_v.size());
		for(uint i=0; i<data_v.size();i++) {
			betas[i] = data_v[i].beta;
			betas_se[i] = data_v[i].beta_se;
		}
		CMetaAnalysis::RandomEffectModel model(betas,betas_se);
		model.process();
		__p_values[counter] = model.getPvalue(two_tailed);
		__z_values[counter] = model.getZvalue();
		__weighted_means[counter] = model.getWeightedMean();
	        __variances[counter] = model.getVarianceSummary();
		__lower_limits[counter] = model.getLowerLimit();
		__upper_limits[counter] = model.getUpperLimit();
		__q[counter] = model.getQ();
		__c[counter] = model.getC();
		__t2[counter] = model.getT2();
		__chromosomes[counter] = data_v[0].chromosome;
		__positions[counter] = data_v[0].position;
		counter++;
	}
}

/*
*Getter methods
*/
VectorXd CMetaGWAS::getPValues() {
	return __p_values;
}

VectorXd CMetaGWAS::getZValues() {
	return __z_values;
}

VectorXd CMetaGWAS::getWeightedMeans() {
	return __weighted_means;
}

VectorXd CMetaGWAS::getVariances() {
	return __variances;
}

VectorXd CMetaGWAS::getLowerLimits() {
	return __lower_limits;
}

VectorXd CMetaGWAS::getUpperLimits() {
	return __upper_limits;
}

VectorXd CMetaGWAS::getQ() {
	return __q;
}

VectorXd CMetaGWAS::getC() {
	return __c;
}

VectorXd CMetaGWAS::getT2() {
	return __t2;
}

std::vector<std::string> CMetaGWAS::getChromosomes() {
	return __chromosomes;
}

std::vector<uint64> CMetaGWAS::getPositions() {
	return __positions;
}

CMetaResults CMetaGWAS::getResults() {
	CMetaResults results;
	results.chromosomes = __chromosomes;
	results.positions = __positions;
	results.p_values = __p_values;
	results.z_values = __z_values;
	results.weighted_means = __weighted_means;
	results.variances = __variances;
	results.lower_limits = __lower_limits;
	results.upper_limits = __upper_limits;
	results.q = __q;
	results.c = __c;
	results.t2 = __t2;
	results.study_counts = __study_counts;
	return results;
}
