#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "CEasyGWAS/stats/CStats.h"
#include "CEasyGWAS/stats/CGaussian.h"
#include "CEasyGWAS/stats/CGamma.h"
#include "CEasyGWAS/stats/CBeta.h"
#include "CEasyGWAS/stats/CChi2.h"
#include "CEasyGWAS/stats/CFisherF.h"
#include "CEasyGWAS/stats/CStudentT.h"
#include "CEasyGWAS/regression/CRegression.h"
#include "CEasyGWAS/utils/CMathHelper.h"
#include "CEasyGWAS/utils/CCrossValidation.h"
#include "CEasyGWAS/gwas/CSingleTraitGWAS.h"
#include "CEasyGWAS/gwas/CFastANOVA.h"
#include "CEasyGWAS/gwas/CScones.h"
#include "CEasyGWAS/utils/CMatrixHelper.h"
#include "CEasyGWAS/io/CPlinkParser.h"
#include "CEasyGWAS/kernel/CKernels.h"
#include <vector>
#include <string>

using namespace std;

/*
*Print some summary of the regression model to demonstrate that CRegression is the parent
*class of the CLinearRegression and CLogisticRegression class. 
*/
void printSummary(CRegression* regression) {
    std::cout << "\n*******************************************" << std::endl;
    regression->print(); //Print results from regression model
    std::cout << "R2:\t" << regression->getRSquared() << std::endl; //Return R2 measure
    std::cout << "DF:\t" << regression->getDF() << std::endl; //Return degrees of freedom
}

int main() {
	MatrixXd X = MatrixXd::Random(12,2); //Initialise a Random Matrix with 12 samples and 2 features
    VectorXd y(12); //Initialise target vector with 12 samples
    y << 1,2,3,4.5,5,6,7.5,8,9,10.5,11,12;
    CLinearRegression lm;
    lm.fit(y,X);
    printSummary(&lm); //Print something to illustrate inheritence
                    
    VectorXd y_binary(12); //Initialise binary target vector with 12 samples
    y << 0,0,0,0,1,0,1,1,1,0,1,1;
    CLogisticRegression lg;
    lg.fit(y,X);
    printSummary(&lg); //Print something to illustrate inheritence
    return 0;
    
    /*
	VectorXd y(12);
	VectorXd x(12);
	y << 0,0,0,0,1,0,1,1,1,0,1,1;
	x << 2100,2300,2500,2700,2900,3100,3300,3500,3700,3900,4100,4300;
	CSingleTraitGWAS::LinearRegression gwas1(y,x);
	//gwas1.test_associations();
	gwas1.permutations(100000);
	logging(INFO,gwas1.getPValues());
	logging(INFO,gwas1.getPermutationPValues());
	

	y << 1,2,3,4.5,5,6,7.5,8,9,10.5,11,12;
	x << 1,2,3,4,5,6,7,8,9,10,11,12;
	VectorXd y1 = VectorXd::Random(12);
	MatrixXd tmp(x.rows(),2);
        tmp << x,x;
	MatrixXd K = CKernels::realizedRelationshipKernel(tmp);
	VectorXd test = CStats::std(x);
	//float64 t1 = CStats::std(x.col(0));
	logging(STATUS,test);
	//logging(INFO,t1);
	
	DiagXd K1 = DiagXd(12);
	K1.diagonal() = VectorXd::Ones(12);
	
	CSingleTraitGWAS::FaSTLMM gwas(y1,tmp,K);
	gwas.permutations(10000);
	logging(INFO,gwas.getPValues().transpose());
	logging(STATUS,gwas.getPermutationPValues().transpose());
	
    logging(INFO,"COMPUTING VARIANCE EXPLAINED");
    float64 var = gwas.computeVarianceExplainedNullModel(3);
    logging(STATUS,"Running Finished");
	logging(ERROR,var);

	exit(0);
	GWASData data;
	logging(STATUS,"Reading Genotype file...");
	CPlinkParser::readPEDFile("/agbs/agkbshare/data/hybrid_danelle/hybrid_genotypes_30x30_seg_0pm.ped",&data);
	logging(STATUS,"Reading Mapping file...");
	CPlinkParser::readMAPFile("/agbs/agkbshare/data/hybrid_danelle/hybrid_genotypes_30x30_seg_0pm.map",&data);
	logging(STATUS,"Reading Phenotype file...");
	CPlinkParser::readPhenotypeFile("/agbs/agkbshare/data/hybrid_danelle/hybrid_pheno_allad.txt",&data);
	logging(STATUS,"Encoding SNP data...");
	CGWASDataHelper::encodeHeterozygousData(&data);
	MatrixXd Xn = data.X;
	for(int64 i=0; i<data.X.cols(); i++) {
		VectorXd stdv = data.X.col(i).array() - data.X.col(i).mean();
		float64 std = sqrt(stdv.array().pow(2).sum());
		Xn.col(i) = (data.X.col(i).array()-data.X.col(i).mean()).array()/std;
	}
	MatrixXd K = Xn*Xn.transpose();
	logging(STATUS,"Running EMMAX");
	CSingleTraitGWAS::EMMAX gwas(data.Y,data.X,K);
	//CSingleTraitGWAS::LinearRegression gwas(data.Y,data.X);
	//gwas.setTestStatistic(1);
	gwas.test_associations();
	logging(STATUS,"Running Finished");
	
	
	srand(rand());
	MatrixXd x = MatrixXd::Random(200,10000);
	VectorXd y = VectorXd::Random(200);
	for(int i=0; i<10000; i++) {
		for(int j=0; j<10000; j++) {
			double res = CStats::pearson_corr(x.col(i),y);
    			CStats::pearson_pval(res,x.rows());
		}
	}
	//MatrixXd X = MatrixXd::Random(200,2);
	//MatrixXd X = MatrixXd::DiagonalMatrix(200);
	//Eigen::DiagonalMatrix<float64,Eigen::Dynamic> X(200);
	*/
	
	//GWASData data;
	//logging(STATUS,"Reading Genotype file...");
	//CPlinkParser::readPEDFile("/agbs/agkbshare/data/hybrid_danelle/hybrid_genotypes_30x30_seg_0pm.ped",&data);
	//logging(STATUS,"Reading Mapping file...");
	//CPlinkParser::readMAPFile("/agbs/agkbshare/data/hybrid_danelle/hybrid_genotypes_30x30_seg_0pm.map",&data);
	//logging(STATUS,"Reading Phenotype file...");
	//CPlinkParser::readPhenotypeFile("/agbs/agkbshare/data/hybrid_danelle/hybrid_pheno_allad.txt",&data);
	//logging(STATUS,"Encoding SNP data...");
	//CGWASDataHelper::encodeHeterozygousData(&data);
	
	/*
	//read phenotype
	ifstream ifs;
	string line;
	VectorXd y = VectorXd::Zero(172);
	MatrixXd X = MatrixXd::Zero(172,1);
	string pheno = "p1.txt";	
	ifs.open(pheno.c_str(),ifstream::in);
	if(!ifs.is_open()) {
		logging(ERROR,"Opening file " + pheno);
		return 0;
	}
	float64 num;
	float64 j=0;
	while(ifs >> num) {
		y(j) = num;
		j++;
	}
	ifs.close();
	pheno = "test.txt";	
	ifs.open(pheno.c_str(),ifstream::in);
	if(!ifs.is_open()) {
		logging(ERROR,"Opening file " + pheno);
		return 0;
	}
	j=0;
	while(ifs >> num) {
		X(j,0) = num;
		j++;
	}
	ifs.close();
	
	CLogisticRegression lm;
	lm.fit(y,X);
	lm.print();
	
	exit(0);
	*/

	/*
	MatrixXd Xn = data.X;
	for(int64 i=0; i<data.X.cols(); i++) {
		VectorXd stdv = data.X.col(i).array() - data.X.col(i).mean();
		float64 std = sqrt(stdv.array().pow(2).sum());
		Xn.col(i) = (data.X.col(i).array()-data.X.col(i).mean()).array()/std;
	}
	MatrixXd K = Xn*Xn.transpose();
	logging(STATUS,"Running EMMAX");
	CSingleTraitGWAS::EMMAX gwas(data.Y,data.X,K);
	//CSingleTraitGWAS::LinearRegression gwas(data.Y,data.X);
	//gwas.setTestStatistic(1);
	gwas.test_associations();
	logging(STATUS,"Running Finished");
	
	VectorXd pvals = gwas.getPValues();
	MatrixXd betas = gwas.getBetas();
	for(int i=0; i<data.n_snps; i++) {
		std::cout << data.chromosomes[i] << " " << data.positions[i] << " " << pvals(i) << " " << betas(i,1) << " " << std::endl;
	}
	logging(INFO,CFisherF::sf(75.8111,1,3));
	logging(INFO,1-CFisherF::cdf(75.8111,1,3));
	logging(INFO,CFisherF::cdf(5,2,1));
	logging(INFO,CFisherF::cdf(0.5,4,8));
	logging(INFO,CFisherF::cdf(100,4,8));
	exit(-1);
	VectorXd y(12);
	VectorXd x(12);
	y << 0,0,0,0,1,0,1,1,1,0,1,1;
	x << 2100,2300,2500,2700,2900,3100,3300,3500,3700,3900,4100,4300;
	
	
	double r = CStats::pearson_corr(x,y);
	logging(INFO,r);
	logging(INFO,CStats::pearson_pval(r,x.rows()));
	logging(INFO,x.rows());

	CLogisticRegression lm;
	lm.fit(y,x);
	lm.print();
	
	CLinearRegression lm1;
	lm1.fit(y,x);
	lm1.print();
	//logging(LOG,lm1.getResiduals());
	logging(LOG,lm1.getLogLikelihood());
	
	ifstream ifs;
	string line;
	VectorXd Y = VectorXd::Zero(172);
	string pheno = "p1.txt";	
	ifs.open(pheno.c_str(),ifstream::in);
	if(!ifs.is_open()) {
		logging(ERROR,"Opening file " + pheno);
		return 0;
	}
	float64 num;
	float64 j=0;
	while(ifs >> num) {
		Y(j) = num;
		j++;
	}
	ifs.close();
	
	MatrixXd X = MatrixXd::Random(172,10000);
	VectorXd y(12);
	VectorXd x(12);
	y << 0,0,0,0,1,0,1,1,1,0,1,1;
	x << 2100,2300,2500,2700,2900,3100,3300,3500,3700,3900,4100,4300;
	CSingleTraitGWAS::LogisticRegression gwas1(Y,X);
	gwas1.test_associations();
	
	*/
	/*
	DiagXd K = DiagXd(12);
	K.diagonal() = VectorXd::Ones(12);

	MatrixXd K1 = x*x.transpose();

	CLinearMixedRegression lmm;
	lmm.setREML(false);
	lmm.setInterval(10);
	lmm.fit(y,x,K);
	lmm.print();

	MatrixXd x1 = MatrixXd::Random(200,25000000);
	VectorXd y1 = VectorXd::Random(200);
	CSingleTraitGWAS::LinearRegression gwas(y1,x1);
	logging(INFO,gwas.getPValues());
	logging(INFO,gwas.getBetas());
	logging(LOG,"--------------");
	CSingleTraitGWAS::LinearRegression gwas1(y,x);
	gwas1.test_associations();
	logging(INFO,gwas1.getPValues());
	logging(INFO,gwas1.getBetas());
	
	logging(WARNING,1-CChi2::cdf(10*2,1));
	logging(WARNING,1-CChi2::cdf(0.2*2,1));
	logging(WARNING,1-CChi2::cdf(0.1*2,1));
	logging(WARNING,1-CChi2::cdf(15*2,1));
	logging(WARNING,1-CChi2::cdf(3.00465*2,1));
	logging(INFO,1-CGamma::Special::regularizedLowerIncompleteGamma(3.00465,0.5));
	
	logging(WARNING,CStudentT::cdf(0.00160324,500));
	logging(WARNING,CStudentT::cdf(0.00160324,100));
	logging(WARNING,CStudentT::cdf(0.1,2));

	MatrixXd X(x.rows(),2);
	X << x,x;
	SparseMatrixXd L(X.cols(),X.cols());
	L.insert(0,1) = 1;
	L.insert(1,0) = 1;

	CScones scones(y,X,L);
	scones.test_associations();
	
	*/
	//cout << test2 << endl;
	//logging(INFO,gwas.getSEBetas());
	/*
	//logging(WARNING,CBeta::cdf(0.5,0.3,0.5));
	logging(WARNING,CStudentT::cdf(0.9,3));
	logging(WARNING,CStudentT::cdf(10,2));
	logging(WARNING,CStudentT::cdf(-2.4,3));
	logging(WARNING,CStudentT::pdf(-10,4));
	logging(WARNING,CStudentT::pdf(2,5));
	*/

	/*
	CLinearRegression lm_null(false);
	VectorXd ones = VectorXd::Ones(12);
	lm_null.fit(y,ones);
	logging(WARNING,lm_null.getLogLikelihood());
	logging(INFO,lm_null.getAIC());
	
	MatrixXd x1 = MatrixXd::Random(12,1000);
	VectorXd y1 = VectorXd::Random(12);
	for(int i=0; i<1; i++) {
		CLinearRegression l1(true);
		//lm.fit(y,x1.col(i));
		l1.fit(y,x1);
		l1.print();
    	}
	return 0;
	*/
}
