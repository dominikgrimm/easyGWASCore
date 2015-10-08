#include <iostream>
#include <string>
#include <fstream>
#include <Eigen/Dense>

#include "CEasyGWAS/stats/CStats.h"
#include "CEasyGWAS/stats/CGaussian.h"
#include "CEasyGWAS/stats/CGamma.h"
#include "CEasyGWAS/stats/CBeta.h"
#include "CEasyGWAS/stats/CChi2.h"
#include "CEasyGWAS/stats/CFisherF.h"
#include "CEasyGWAS/stats/CStudentT.h"
#include "CEasyGWAS/regression/CRegression.h"
#include "CEasyGWAS/utils/CMathHelper.h"
#include "CEasyGWAS/utils/CMatrixHelper.h"
#include "CEasyGWAS/gwas/CSingleTraitGWAS.h"
#include "CEasyGWAS/gwas/CFastANOVA.h"
#include "CEasyGWAS/gwas/CScones.h"
#include "CEasyGWAS/utils/StringHelper.h"

using namespace std;

int main() {

	VectorXd y = VectorXd::Zero(500);
	MatrixXd X = MatrixXd::Zero(500,1000);
	SparseMatrixXd L(1000,1000);
	
	string pheno = "data/testing/scones/y.txt";
	string geno = "data/testing/scones/X.txt";
	string network = "data/testing/scones/L.txt";
	
	//read phenotype
	ifstream ifs;
	string line;
	logging(INFO,"Reading Phenotype...");
	ifs.open(pheno.c_str(),ifstream::in);
	if(!ifs.is_open()) {
		logging(ERROR,"Opening file " + pheno);
		return 0;
	}
	int64 j=0;
	float64 num;
	while(ifs >> num) {
		y(j) = num;
		j++;
	}
	ifs.close();
	
	logging(INFO,"Reading Genotype...");
	ifs.open(geno.c_str(),ifstream::in);
	if(!ifs.is_open()) {
		logging(ERROR,"Opening file " + pheno);
		return 0;
	}
	j=0;
	while(ifs >> line) {
		std::vector<string> sv = StringHelper::split(line,",");
		for(uint i=0; i<sv.size();i++) {
			string t = sv[i];
			X(j,i) = StringHelper::string_to<uint>(t);
		}
		j++;
	}
	ifs.close();
	
	logging(INFO,"Reading Sparse Network File...");
	ifs.open(network.c_str(),ifstream::in);
	if(!ifs.is_open()) {
		logging(ERROR,"Opening file " + network);
		return 0;
	}
	typedef Eigen::Triplet<float64> T;
	std::vector<T> triplet;
	triplet.reserve(1000*1000);
	//triplet.reserve(100000*100000);
	while(ifs >> line) {
		std::vector<string> sv = StringHelper::split(line,",");
		triplet.push_back(T(StringHelper::string_to<uint64>(sv[0])-1,
				    StringHelper::string_to<uint64>(sv[1])-1,
				    StringHelper::string_to<uint64>(sv[2])));
		//uint64 t1 = StringHelper::string_to<uint64>(sv[0]);
		//uint64 t2 = StringHelper::string_to<uint64>(sv[1]);
		//cout << t1 << " " << t2 << endl;
		//L.insert(StringHelper::string_to<uint64>(sv[0]), 
		//	 StringHelper::string_to<uint64>(sv[1])) = 1;
	}
	L.setFromTriplets(triplet.begin(),triplet.end());
	logging(INFO,"Running SConES");
	CScones scones(y,X,L);
	scones.test_associations();
	logging(STATUS,"Indicator Vector:");
	logging(INFO,scones.getIndicatorVector().transpose());
	logging(INFO,"Best Eta: " + StringHelper::to_string<float64>(scones.getBestEta()));
	logging(INFO,"Best Lambda: " + StringHelper::to_string<float64>(scones.getBestLambda()));

	return 0;
}
