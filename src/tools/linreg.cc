#include <iostream>
#include <string>
#include <time.h>

#include "CEasyGWAS/io/CLogging.h"
#include "CEasyGWAS/gwas/CSingleTraitGWAS.h"
#include "CEasyGWAS/io/CPlinkParser.h"
#include "CEasyGWAS/kernel/CKernels.h"
#include "CEasyGWAS/io/CGWASDataIO.h"

int main(int argc, char* argv[]) {
	//get command line arguments
	if(argc<5) {
		logging(WARNING,"Wrong number of argurments:");
		logging(WARNING,"");
		logging(WARNING,"linreg <plink genotype file> <plink phenotype file> <float minor allele frequency filter> <output_folder> [permutations - optional]");
		logging(WARNING,"");
		exit(-1);
	}

	clock_t total = clock();
	string genotype_str = argv[1];
	string phenotype_str = argv[2];
	float64 maf = StringHelper::string_to<float64>(argv[3]);
	string outfolder_str = argv[4];
	
	int permutations = 0;
	if(argc==6) {
		permutations = StringHelper::string_to<int>(argv[5]);
	}

	GWASData data;
	CLogging logger(outfolder_str + "/run.log");
	
	clock_t begin;
	begin = clock();
	logger.log(STATUS,"Reading Genotype file...");
	CPlinkParser::readPEDFile(genotype_str + ".ped",&data);
	logger.log(INFO,"Number of SNPs:\t\t" + StringHelper::to_string<uint64>(data.n_snps));
	logger.log(INFO,"Number of Samples:\t" + StringHelper::to_string<uint64>(data.n_samples));
	logger.log(INFO,"MAF Filter:\t\t" + StringHelper::to_string<float64>(maf));
	logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
    
	begin = clock();
	logger.log(STATUS,"Reading Mapping file...");
	CPlinkParser::readMAPFile(genotype_str + ".map",&data);
	logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	
	begin = clock();
	logger.log(STATUS,"Reading Phenotype file...");
	CPlinkParser::readPhenotypeFile(phenotype_str,&data);
	logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	
	begin = clock();
	logger.log(STATUS,"Encoding SNP data...");
	CGWASDataHelper::encodeHeterozygousData(&data);
	logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	
	/*
    begin = clock();
	logger.log(STATUS,"Filter unique SNPs...");
	CGWASDataHelper::filterUniqueSNPs(&data);
	logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	*/

	
	for(uint i=0; i<data.phenotype_names.size(); i++) {
		begin = clock();
		logger.log(STATUS,"Restrict Genotype to Samples found in phenotype file: " + data.phenotype_names[i]);
		GWASData tmpData = CGWASDataHelper::removeSamples4MissingData(data,i);
		if(tmpData.n_samples!=data.n_samples) {
			logger.log(INFO,"#Samples removed: " + StringHelper::to_string<uint64>(data.n_samples-tmpData.n_samples));
			logger.log(INFO,"#Final Set of Samples: " + StringHelper::to_string<uint64>(tmpData.n_samples));
		}
		logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
        
        /*begin = clock();
	    logger.log(STATUS,"Filter Non Informative SNPs for phenotype: " + data.phenotype_names[i]);
	    CGWASDataHelper::filterNonInformativeSNPs(&tmpData);
	    logger.log(INFO,"Number of SNPs left: " + StringHelper::to_string<uint64>(tmpData.n_snps));
	    logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
        

        begin = clock();
	    logger.log(STATUS,"Filter SNPs by MAF for phenotype: " + data.phenotype_names[i]);
	    CGWASDataHelper::filterSNPsByMAF(&tmpData,maf);
	    logger.log(INFO,"Number of SNPs (after MAF): " + StringHelper::to_string<uint64>(tmpData.n_snps));
	    logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	    */

		begin = clock();
		logger.log(STATUS,"Create unique SNP hash...");
		CGWASDataHelper::createSNPHash(&tmpData);
		if(tmpData.n_unique_snps>0) {
			logger.log(ATTENTION, "This dataset contains " + StringHelper::to_string<uint64>(tmpData.n_unique_snps) + " unique SNPs (after encoding)!");
			logger.log(ATTENTION, "\t--> Thus " + StringHelper::to_string<uint64>(tmpData.n_snps-tmpData.n_unique_snps) + " SNP(s) are non-unique (same encoding as at least one of the other SNPs)!");
		}
		logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
		

		begin = clock();
		logger.log(STATUS,"Running Linear Regression for phenotype: " + data.phenotype_names[i]);
		CSingleTraitGWAS::LinearRegression gwas(tmpData.Y.col(i),tmpData.X);
		if(permutations==0) {
			gwas.test_associations();
		} else {
			logger.log(INFO,"\tPerforming Permutation Test with " + StringHelper::to_string<int>(permutations) + " permutations...");
			gwas.permutations(permutations);
		}

		GWASResults results = gwas.getResults();
		logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
		
		//write output 
		begin = clock();
		string output_str = outfolder_str + "/" + data.phenotype_names[i] + ".linreg.out.txt";
		logger.log(STATUS,"Writing output to " + output_str);
		CGWASDataIO::writeSummaryOutput(output_str, tmpData, results);
		logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	
	}

	logger.log(STATUS,"Finished all computations in " + StringHelper::to_string<clock_t>((clock()-total)/CLOCKS_PER_SEC) + " sec\n");

	return 0;
}
