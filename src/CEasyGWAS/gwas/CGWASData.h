#ifndef CGWASDATA_CLASS
#define CGWASDATA_CLASS

#include "CEasyGWAS/globals.h"
#include "CEasyGWAS/utils/StringHelper.h"

#include <fstream>
#include <vector>
#include <string>
#include <map>

#include <Eigen/Dense>

/*
*CGWASData Exception Class
*/
class CGWASDataException {
	private:
		std::string __error_msg;
	public:
		CGWASDataException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CGWASData Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};

/*
*Data object
*/
class GWASData {
	public:
		std::vector<std::vector<char> > raw_snps;
		std::vector<std::string> sample_ids;
		std::vector<std::string> chromosomes;
		std::vector<uint64> positions;
		std::vector<uint64> removed_snp_indices;
		std::vector<float64> snp_distance;
		std::vector<std::string> snp_identifiers;
		std::vector<std::string> phenotype_names;
		std::vector<std::string> family_ids;
		std::vector<std::string> paternal_ids;
		std::vector<std::string> maternal_ids;
		std::vector<uint> sex;
		VectorXd MAF;
		MatrixXd Y;
		MatrixXd X;
		MatrixXd K;
		SparseMatrixXd network;
		VectorXd snp_hash;
		uint64 n_samples;
		uint64 n_snps;
		uint64 n_unique_snps;
		std::string genotype_data_type;
		std::string genotype_encoding;
};

/*
*Result data class
*/ 
class GWASResults {
	public:
        std::vector<std::string> chromosomes;
		std::vector<uint64> positions;
		VectorXd samples;
		VectorXd p_values;
		MatrixXd betas;
		MatrixXd se_betas;
		VectorXd test_statistics;
		VectorXd alternative_loglikelihoods;
		float64 null_loglikelihood;
};

class CGWASDataHelper {
	
	private:
		static VectorXd __updateMAF(MatrixXd const&);
		GWASData __data;

	public:
		static const uint additive = 0;
		static const uint recessive = 1;
		static const uint dominant = 2;
		static const uint codominant = 3;
		
		//Class methods for more flexible interactions
		void encodeHomozygousData(std::vector< std::vector<char> > const&,
					  uint64 const&, uint64 const&) throw (CGWASDataException);
		void encodeHeterozygousData(std::vector< std::vector<char> > const&,
					  uint64 const&, uint64 const&, uint const&) throw (CGWASDataException);

		void releaseMemory();
		MatrixXd getEncodedData();
		VectorXd getMAF();

		//static member functions for C++ only
		static void encodeHomozygousData(GWASData*) throw (CGWASDataException);
		static void encodeHeterozygousData(GWASData*) throw (CGWASDataException);
		static void encodeHeterozygousData(GWASData*,uint const&) throw (CGWASDataException);
		static void filterSNPsByMAF(GWASData*,float64 const&) throw (CGWASDataException);
		static void filterNonInformativeSNPs(GWASData*) throw (CGWASDataException);
		//static void filterSNPsBySmallIndel(GWASData*,int const&) throw (CGWASDataException);
		static void filterUniqueSNPs(GWASData*) throw (CGWASDataException);
		static void createSNPHash(GWASData*) throw (CGWASDataException);
		static GWASData removeSamples4MissingData(GWASData const&, uint const&) throw (CGWASDataException);
		static GWASData removeSamples4MissingData(GWASData const&, uint const&, bool const) throw (CGWASDataException);
};

#endif //CGWASDATA_CLASS
