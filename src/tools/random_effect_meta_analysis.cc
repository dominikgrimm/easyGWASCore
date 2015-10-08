#include <iostream>
#include <string>
#include <time.h>

#include "CEasyGWAS/io/CLogging.h"
#include "CEasyGWAS/gwas/CSingleTraitGWAS.h"
#include "CEasyGWAS/io/CPlinkParser.h"
#include "CEasyGWAS/kernel/CKernels.h"
#include "CEasyGWAS/io/CGWASDataIO.h"
#include "CEasyGWAS/gwas/CMetaGWAS.h"

using namespace std;

int main(int argc, char* argv[]) {
	//get command line arguments
	if(argc<5) {
		logging(WARNING,"Wrong number of argurments:");
		logging(WARNING,"");
		logging(WARNING,"random_effect_meta_analysis <number studies (n)> <outfile_prefix> <study file1> <study file2> ... <study file n>");
		logging(WARNING,"");
		exit(-1);
	}

	clock_t total = clock();
	clock_t begin;
	uint n_studies = StringHelper::string_to<uint>(argv[1]);
	string outfolder_str = argv[2];

	CLogging logger(outfolder_str + "_run.log");
	
	CMetaGWAS meta_gwas;
		
	for(uint i=0; i<n_studies; i++) {
		begin = clock();
		string study_file = argv[3+i];
		logger.log(STATUS,"Reading data from study " + study_file);
		GWASResults results = CGWASDataIO::readGWASResults(study_file);
		logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
		
		meta_gwas.addStudy(results.p_values,results.chromosomes,results.positions,results.samples,results.betas,results.se_betas);
			
	}
	
	begin = clock();
	logger.log(STATUS,"Performing Random Effect Model meta analysis...");
	//meta_gwas.performRandomEffectModel();
	//meta_gwas.performFixedEffectModel();
	meta_gwas.performStoufferZWeighted();
	logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	
	string output_str = outfolder_str + ".meta.out.txt";
	logger.log(STATUS,"Writing output to " + output_str);
	CGWASDataIO::writeMetaResultsFile(output_str, meta_gwas.getResults());
	logger.log(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	
	logger.log(STATUS,"Finished all coomputations in " + StringHelper::to_string<clock_t>((clock()-total)/CLOCKS_PER_SEC) + " sec\n");

	return 0;
}
