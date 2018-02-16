#include <iostream>
#include <string>
#include <time.h>

#include "CEasyGWAS/gwas/CScones.h"
#include "CEasyGWAS/io/CPlinkParser.h"
#include "CEasyGWAS/kernel/CKernels.h"
#include "CEasyGWAS/io/CGWASDataIO.h"
#include "CEasyGWAS/io/CSconesIO.h"
#include "CEasyGWAS/stats/CStats.h"
#include "CEasyGWAS/utils/CMatrixHelper.h"

using namespace std;

int main(int argc, char* argv[]) {
	//get command line arguments
	if(argc<5) {
		logging(WARNING,"Wrong number of argurments:");
		logging(WARNING,"");
		logging(WARNING,"scones <plink genotype file> <plink phenotype file> <sparse network file> <float minor allele frequency filter> <outdir> <snp_encoding> <PC>");
		logging(WARNING,"");
		exit(-1);
	}

	float64 total_data = clock();
    float64 total = clock();
	string genotype_str = argv[1];
	string phenotype_str = argv[2];
	string network_str = argv[3];
	float64 maf = StringHelper::string_to<float64>(argv[4]);
	string outfolder_str = argv[5];
    string snp_encoding = argv[6];
	int64 pcs = StringHelper::string_to<int64>(argv[7]);
    
    uint encoding = 0;
    if(snp_encoding=="additive") encoding = 0;
    else if(snp_encoding=="recessive") encoding = 1;
    else if(snp_encoding=="dominant") encoding = 2;
    else if(snp_encoding=="overdominant") encoding = 3;
    else {
        logging(ERROR,"Encoding does not exist!");
        exit(-1);
    }
	GWASData data;
	
	float64 begin;
	begin = clock();
	logging(STATUS,"Reading Genotype file...");
	CPlinkParser::readPEDFile(genotype_str + ".ped",&data);
	logging(INFO,"Number of SNPs: " + StringHelper::to_string<uint64>(data.n_snps));
	logging(INFO,"Number of Samples: " + StringHelper::to_string<uint64>(data.n_samples));
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(STATUS,"Reading Mapping file...");
	CPlinkParser::readMAPFile(genotype_str + ".map",&data);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	
	begin = clock();
	logging(STATUS,"Reading Phenotype file...");
	CPlinkParser::readPhenotypeFile(phenotype_str,&data);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	
	begin = clock();
	logging(STATUS,"Encoding SNP data...");
	CGWASDataHelper::encodeHeterozygousData(&data,encoding);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	
	/*
	begin = clock();
	logging(STATUS,"Filter unique SNPs...");
	CGWASDataHelper::filterUniqueSNPs(&data);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	*/
	
	begin = clock();
	logging(STATUS,"Filter SNPs by MAF...");
	CGWASDataHelper::filterSNPsByMAF(&data,maf);
	logging(INFO,"Number of SNPs (after MAF): " + StringHelper::to_string<uint64>(data.n_snps));
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	
	begin = clock();
	logging(STATUS,"Loading and filtering network file...");
	CSconesIO::readSparseNetworkFile(network_str,&data);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	

	for(uint i=0; i<data.phenotype_names.size(); i++) {

		begin = clock();
		logging(STATUS,"Remove samples with missing values...");
		GWASData tmpData = CGWASDataHelper::removeSamples4MissingData(data,i);
		if(tmpData.n_samples!=data.n_samples) {
			logging(INFO,"#Samples removed: " + StringHelper::to_string<uint64>(data.n_samples-tmpData.n_samples));
		}
		logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
		logging(WARNING,"Total Data: " + StringHelper::to_string<float64>(float64(clock()-total_data)/CLOCKS_PER_SEC) + " sec\n");
	
		//OPTIONAL SET DIFFERENT SConES Settings
		begin = clock();
        //COMPUTE PRINCIPLE COMPONENTS
        if (pcs>0) {
            begin = clock();
            logging(STATUS,"Computing Realized Relationship Kernel and Principle Components...");
            tmpData.K = CKernels::realizedRelationshipKernel(tmpData.X);
            MatrixXd PCs = CStats::principle_components(tmpData.K);
            PCs = sliceColsMatrix(PCs,VectorXd::LinSpaced(pcs,0,pcs-1));
            logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

            logging(STATUS,"Running Scones for phenotype: " + data.phenotype_names[i]);
		    CSconesSettings settings;
            CScones scones(tmpData.Y.col(i),tmpData.X,tmpData.network,PCs,settings);
            scones.test_associations();
            logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
        
            begin = clock();
            string output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.out.txt";
            logging(STATUS,"Writing output to " + output_str);
            CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),scones.getBestLambda(),scones.getBestEta());
            output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.pmatrix.txt";
            logging(STATUS,"Writing pmatrix to " + output_str);
            CSconesIO::writeCMatrix(output_str, scones.getCMatrix(),scones.getSettings());
            logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
        } else {
		    logging(STATUS,"Running Scones for phenotype: " + data.phenotype_names[i]);
            CScones scones(tmpData.Y.col(i),tmpData.X,tmpData.network);
            scones.test_associations();
            logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
        
            begin = clock();
            string output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.out.txt";
            logging(STATUS,"Writing output to " + output_str);
            CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),scones.getBestLambda(),scones.getBestEta());
            output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.pmatrix.txt";
            logging(STATUS,"Writing pmatrix to " + output_str);
            CSconesIO::writeCMatrix(output_str, scones.getCMatrix(),scones.getSettings());
            logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

        }


        //logging(STATUS,"Indicator Vector:");
		//logging(INFO,scones.getIndicatorVector().transpose());
		//logging(INFO,"Best Eta: " + StringHelper::to_string<float64>(scones.getBestEta()));
		//logging(INFO,"Best Lambda: " + StringHelper::to_string<float64>(scones.getBestLambda()));
	}

	logging(STATUS,"Finished all coomputations in " + StringHelper::to_string<float64>((clock()-total)/CLOCKS_PER_SEC) + " sec\n");

	return 0;
}
