#ifndef CEFFECTSIZE_CLASS
#define CEFFECTSIZE_CLASS

#include "CEasyGWAS/globals.h"

namespace CEffectSize {

class MeanEffectSize {

	private:
		VectorXd __x1;
		VectorXd __x2;
		VectorXd __se1;
		VectorXd __se2;
		VectorXd __n1;
		VectorXd __n2;
		VectorXd __variance;
		VectorXd __D;
		VectorXd __Swithin;
		VectorXd __J;
		VectorXd __g;
		uint __standardized;

	public:
		static const uint equalPopulationStandardDeviation = 0;
		static const uint nonEqualPopulationStandardDeviation = 1;
		static const uint independetGroups = 3;
		static const uint matchedGroups = 4;
		static const uint standardizedCorrected = 6;

		MeanEffectSize(VectorXd const&, VectorXd const&,
			       VectorXd const&, VectorXd const&,
			       VectorXd const&, VectorXd const&);

		VectorXd getD();
		VectorXd getVariance();
		VectorXd getVariance(uint const&);
		VectorXd getSE();
		VectorXd getJ();
		VectorXd getHedgesG();
};

};

#endif //CEFFECTSIZE_CLASS
