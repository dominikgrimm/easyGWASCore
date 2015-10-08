#ifndef C_FASTANOVA_CLASS
#define C_FASTANOVA_CLASS

#include "CEasyGWAS/globals.h"

#include <map>
#include <list>
#include <vector>

/*
*FastANOVA Algorithm 
*Xiang Zhang, Fei Zou and Wei Wang: FastANOVA: an Efficient Algorithm for Genome-Wide Association Study, KDD 08
*/

namespace CEpistasis {
	
	/*
	 *CEpistasisException Class
	 */
	class CEpistasisException {
		private:
			std::string __error_msg;
		public:
			CEpistasisException(std::string const& error_msg) : __error_msg(error_msg) {
				std::cout << RED << "CEpistasis Exception: " << error_msg << BLACK << "\n";
			}

			std::string what() {
				return __error_msg;
			}
	};

	typedef struct episnp {
		uint64 i;
		uint64 j;
		uint64 phenotype_id;
		float64 statistic;
		float64 pvalue;
		bool operator<(episnp const& esnp) const {
			if(statistic>esnp.statistic)
				return true;
			else
				return false;
		}
	} episnp;

	/*
	 *FastANOVA Class
	 */
	class CFastANOVA {
		private:
			uint64 __n_samples;
			uint64 __n_snps;
			uint __permutations;
			uint __alphaK;
			float64 __seed;
			float64 __fwer;

			float64 __tm;
			float64 __SST;
			float64 __theta;
			float64 __theta_pval;
			float64 __df1;
			float64 __df2;
			uint64 __visited;		
			
			VectorXd __pvalues;
			VectorXd __test_statistics;
			MatrixXd __snp_pairs;

			std::vector<VectorXd> __contingencyTable2;

			//index array
			std::map<float64,std::map<float64,std::list<uint64> > > __indexArray;

			VectorXd __y;
			MatrixXd __Y; //Permuted phenotypes
			MatrixXd __X;
		
			float64 __computeFStatistic(VectorXd const&, VectorXd const&, VectorXd const&);
			void __permutePhenotypes();
			void __compute_na1_nb1(uint64 const&, uint64 const&, float64*, float64*);
			void __updateIndexArray(float64 const&, float64 const&, uint64 const&);

			void __computeContingencyTable2(VectorXd const&, uint64 const&);
			void __computeCandidateList(float64 const&, std::list<uint>*);
			float64 __computeUpperBound(VectorXd const&, float64 const&);
			float64 __computeSSB();
		
			void __updateTopList(std::list<episnp>&, episnp const&);

			void __computeCriticalFalpha();
			void __computeSignificantSNPs();

			void __checkdata() throw (CEpistasisException);
			
		public:
			CFastANOVA();
			CFastANOVA(VectorXd const&, MatrixXd const&) throw (CEpistasisException);
			CFastANOVA(VectorXd const&, MatrixXd const&, uint const&) throw (CEpistasisException);
			CFastANOVA(VectorXd const&, MatrixXd const&, uint const&, float64 const&) throw (CEpistasisException);
			CFastANOVA(VectorXd const&, MatrixXd const&, uint const&, float64 const&, float64 const&) throw (CEpistasisException);

			void test_associations();
			void test_associations(float64 const&);
			
			//getter and setter methods
			void setPhenotype(VectorXd const&) throw (CEpistasisException);
			void setGenotype(MatrixXd const&);
			void setFWER(float64 const&) throw (CEpistasisException);
			void setPermutations(uint const&);
			void setSeed(float64 const&);

			VectorXd getPValues();
			VectorXd getTestStatistics();
			MatrixXd getSNPPairs();
			float64 getFalpha();
			float64 getAlphaPvalue();
			uint64 getNumVisitedPairs();
	};

};

#endif //C_FASTANOVA_CLASS
