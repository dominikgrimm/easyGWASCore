#include <iostream>
#include <string>
#include <time.h>

#include "CEasyGWAS/gwas/CSingleTraitGWAS.h"
#include "CEasyGWAS/io/CPlinkParser.h"
#include "CEasyGWAS/kernel/CKernels.h"
#include "CEasyGWAS/io/CGWASDataIO.h"

using namespace std;

int main(int argc, char* argv[]) {
	//get command line arguments
	if(argc<3) {
		logging(WARNING,"Wrong number of argurments:");
		logging(WARNING,"");
		logging(WARNING,"filter_snps <plink genotype file> <float minor allele frequency filter> <out_prefix>");
		logging(WARNING,"");
		exit(-1);
	}

	clock_t total = clock();
	string genotype_str = argv[1];
	float64 maf = StringHelper::string_to<float64>(argv[2]);
	string outfolder_str = argv[3];

	GWASData data;
	
	clock_t begin;
	begin = clock();
	logging(STATUS,"Reading Genotype file...");
	CPlinkParser::readPEDFile(genotype_str + ".ped",&data);
	logging(INFO,"Number of SNPs: " + StringHelper::to_string<uint64>(data.n_snps));
	logging(INFO,"Number of Samples: " + StringHelper::to_string<uint64>(data.n_samples));
	logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(STATUS,"Reading Mapping file...");
	CPlinkParser::readMAPFile(genotype_str + ".map",&data);
	logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(STATUS,"Encoding SNP data...");
	CGWASDataHelper::encodeHeterozygousData(&data);
	logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	
	begin = clock();
	logging(STATUS,"Filter SNPs by MAF...");
	CGWASDataHelper::filterSNPsByMAF(&data,maf);
	logging(INFO,"Number of SNPs (after MAF): " + StringHelper::to_string<uint64>(data.n_snps));
	
	/*begin = clock();
	logging(STATUS,"Filter SNPs by SMALL INDELs...");
	CGWASDataHelper::filterSNPsBySmallIndel(&data,0);
	logging(INFO,"Number of SNPs (after indel removal): " + StringHelper::to_string<uint64>(data.n_snps));
    */

	CGWASDataIO::writeFilteredPlinkFile(outfolder_str, data);
	
	logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	return 0;
}
