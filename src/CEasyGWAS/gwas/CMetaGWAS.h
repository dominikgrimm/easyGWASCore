#ifndef CMETAGWAS_CLASS
#define CMETAGWAS_CLASS

#include "CEasyGWAS/globals.h"

#include <map>
#include <vector>

class CMetaData {
	public:
		std::string chromosome;
		uint64 position;
		float64 p_value;
		uint64 n_samples;
		float64 beta;
		float64 beta_se;
};

class CMetaResults {
	public:
		std::vector<std::string> chromosomes;
		std::vector<uint64> positions;
		VectorXd p_values;
		VectorXd z_values;
		VectorXd weighted_means;
		VectorXd variances;
		VectorXd lower_limits;
		VectorXd upper_limits;
		VectorXd q;
		VectorXd c;
		VectorXd t2;
		float64 study_counts;
};

class CMetaGWAS {
	
	private:
		std::vector<std::string> __chromosomes;
		std::vector<uint64> __positions;
		VectorXd __p_values;
		VectorXd __z_values;
		VectorXd __weighted_means;
		VectorXd __variances;
		VectorXd __lower_limits;
		VectorXd __upper_limits;
		VectorXd __q;
		VectorXd __c;
		VectorXd __t2;
		float64 __study_counts;

		std::map<std::string,std::vector<CMetaData> > __study_map;
		std::map<std::string,std::vector<CMetaData> >::iterator __study_it;
	public:
		CMetaGWAS();
	
		void addStudy(VectorXd const&, std::vector<std::string> const&, std::vector<uint64> const&);
		void addStudy(VectorXd const&, std::vector<std::string> const&, std::vector<uint64> const&, VectorXd const&);
		void addStudy(VectorXd const&, std::vector<std::string> const&, std::vector<uint64> const&, VectorXd const&, VectorXd const&, VectorXd const&);
		
        //Additional methods espacially for Python wrapping
        void addPValuesStudy(VectorXd const&,std::vector<std::string> const&,VectorXd const&);
        void addWeightedPValuesStudy(VectorXd const&,std::vector<std::string> const&,VectorXd const&,uint64 const&);
        void addEffectSizeStudy(std::vector<std::string> const&,VectorXd const&,VectorXd const&,VectorXd const&);

		//Meta analysis methods
		void performFisherMethod();
		void performFisherMethod(bool const&);
		void performStoufferZ();
		void performStoufferZ(bool const&);
		void performStoufferZ(bool const&, bool const&);
		void performStoufferZWeighted();
		void performStoufferZWeighted(bool const&);
		void performStoufferZWeighted(bool const&, bool const&);
		void performFixedEffectModel();
		void performFixedEffectModel(bool const&);
		void performFixedEffectModel(bool const&, bool const&);
		void performRandomEffectModel();
		void performRandomEffectModel(bool const&);
		void performRandomEffectModel(bool const&, bool const&);
		
		//Getter methods
		VectorXd getPValues();
		VectorXd getZValues();
		VectorXd getWeightedMeans();
		VectorXd getVariances();
		VectorXd getLowerLimits();
		VectorXd getUpperLimits();
		VectorXd getQ();
		VectorXd getC();
		VectorXd getT2();
		std::vector<std::string> getChromosomes();
		std::vector<uint64> getPositions();		

		//Get all results
		CMetaResults getResults();
};

#endif //CMETAGWAS_CLASS
