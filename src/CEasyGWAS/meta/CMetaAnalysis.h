#ifndef CMETAANALYSIS_CLASS
#define CMETAANALYSIS_CLASS

#include "CEasyGWAS/globals.h"

namespace CMetaAnalysis {
/*
*CMetaAnalysis Exception Class
*/
class CMetaAnalysisException {
	private:
		std::string __error_msg;
	public:
		CMetaAnalysisException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CMetaAnalysis Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

class CombinedPvalues {
	public:
		static float64 FisherMethod(VectorXd const&);
		static float64 StoufferZ(VectorXd const&);
		static float64 StoufferZ(VectorXd const&, bool const&);
		static float64 StoufferPval(float64 const&);
		static float64 StoufferPval(float64 const&, bool const&);
		static float64 StoufferZWeighted(VectorXd const&, VectorXd const&);
		static float64 StoufferZWeighted(VectorXd const&, VectorXd const&, VectorXd const&);
};

class EffectModel {
	protected:
		VectorXd _effects;
		VectorXd _variance_effects;
		VectorXd _weights;
		float64 _weighted_mean;
		float64 _variance_summary;
		float64 _se_summary;
		float64 _lower_limit;
		float64 _upper_limit;
		float64 _z;
	public:
		EffectModel() {};
		EffectModel(VectorXd const&, VectorXd const&);
		~EffectModel() {};

		void setEffects(VectorXd const&);
		void setVarianceEffects(VectorXd const&);
		
		virtual void process() {};
		
		float64 getZvalue();
		float64 getPvalue();
		float64 getPvalue(const bool&);
		float64 getLowerLimit();
		float64 getUpperLimit();
		float64 getVarianceSummary();
		float64 getSESummary();
		float64 getWeightedMean();
		VectorXd getWeights();

};

class FixedEffectModel: public EffectModel{

	public:
		FixedEffectModel(VectorXd const&, VectorXd const&);

		void process();

};

class RandomEffectModel: public EffectModel{
	private:
		float64 __Q;
		float64 __C;
		float64 __T_2;

	public:
		RandomEffectModel(VectorXd const&, VectorXd const&);

		void process();
		
		float64 getQ();
		float64 getC();
		float64 getT2();
};

}; //Namespace CMetaAnalysis


#endif //CMETAANALYSIS_CLASS
