#include "CGWASDataIO.h"

#include "CEasyGWAS/utils/StringHelper.h"

#include <fstream>
#include <iterator>
#include <sstream>

void CGWASDataIO::writeSummaryOutput(std::string const& outfile, GWASData const& data, GWASResults const& results) {
	std::ofstream ofs;
	ofs.open(outfile.c_str());
	if(!ofs.is_open()) {
		logging(ERROR,"Writing output failed!");
		exit(-1);
	}
	ofs << "SNP ID\tCHR\tPositions\t#Samples\tPValue\tTestStatistic\tBeta\tSEBeta\tLL0\tLLAlt\tSNPHash" << endl;
	for(uint i=0; i<results.p_values.rows();i++) {
		ofs << data.snp_identifiers[i] << "\t" 
		    << data.chromosomes[i] << "\t" 
		    << data.positions[i] << "\t" 
		    << data.n_samples << "\t" 
		    << results.p_values(i) << "\t"
		    << results.test_statistics(i) << "\t"
		    << results.betas(i,1) << "\t" 
		    << results.se_betas(i,1) << "\t" 
		    << results.null_loglikelihood << "\t" 
		    << results.alternative_loglikelihoods(i) << "\t" 
		    << data.snp_hash(i) << "\t"
		    << endl;
	}
	ofs.close();
	
}

void CGWASDataIO::writeFilteredPlinkFile(std::string const& outfile, GWASData const& data) {
	std::ofstream ofs;
	std::string ped_file = outfile + ".ped";
	ofs.open(ped_file.c_str());
	if(!ofs.is_open()) {
		logging(ERROR,"Writing output failed!");
		exit(-1);
	}
	uint64 snp_id = -1;
	for(uint64 i=0; i<data.n_samples;i++) {
		ofs << data.family_ids[i] << " " << data.sample_ids[i] << " " << data.paternal_ids[i] << " " << data.maternal_ids[i]
		    << " " << data.sex[i] << " 0 ";	
		for(uint64 j=0; j<data.n_snps;j++) {
			if(j==data.n_snps-1) {
				snp_id = data.removed_snp_indices[j]; 
				ofs << data.raw_snps[i][snp_id] << " " << data.raw_snps[i][snp_id] << "\n";
			} else {
				snp_id = data.removed_snp_indices[j]; 
				ofs << data.raw_snps[i][snp_id] << " " << data.raw_snps[i][snp_id] << " ";
			}
		}
		//snp_id = data.removed_snp_indices[data.n_snps]; 
		//ofs << data.raw_snps[i][snp_id] << " " << data.raw_snps[i][snp_id] << "\n";
	}
	ofs.close();
	
	//Write MAP file
	std::string map_file = outfile + ".map";
	ofs.open(map_file.c_str());
	if(!ofs.is_open()) {
		logging(ERROR,"Writing output failed!");
		exit(-1);
	}
	for(uint64 j=0; j<data.n_snps;j++) {
		ofs << data.chromosomes[j] << " " << data.snp_identifiers[j] << " " << data.snp_distance[j] << " " << data.positions[j] << "\n";
	}
	ofs.close();
}

GWASResults CGWASDataIO::readGWASResults(std::string const& filename) {
	GWASResults results;
	std::ifstream ifs;
	ifs.open(filename.c_str(), std::ifstream::in);
	
	if(!ifs.is_open()) {
		logging("ERROR", "Opening Result file " + filename);
		exit(0);
	}
	std::vector<float64> samples;
	std::vector<float64> p_values;
	std::vector<float64> betas;
	std::vector<float64> betas_se;
	std::string line;
	uint64 counter=0;
	while(ifs.good()) {
		getline(ifs,line);
		if (counter==0) {
			counter++;
			continue;
		}
		std::vector<std::string> sv = StringHelper::split(line,"\t");
		if(sv.size() < 7) continue;
		line = StringHelper::trim(line);
		results.chromosomes.push_back(sv[1]);
		results.positions.push_back(StringHelper::string_to<uint64>(sv[2]));
		samples.push_back(StringHelper::string_to<float64>(sv[3]));
		p_values.push_back(StringHelper::string_to<float64>(sv[4]));	
		betas.push_back(StringHelper::string_to<float64>(sv[6]));	
		betas_se.push_back(StringHelper::string_to<float64>(sv[7]));	
		counter++;
	}
	results.samples = VectorXd::Map(&samples[0],samples.size());
	results.p_values = VectorXd::Map(&p_values[0],p_values.size());
	results.betas = VectorXd::Map(&betas[0],betas.size());
	results.se_betas = VectorXd::Map(&betas_se[0],betas_se.size());
	
	ifs.close();
	return results;
}

void CGWASDataIO::writeMetaResultsFile(std::string const& outfile, CMetaResults const& results) {
	std::ofstream ofs;
	ofs.open(outfile.c_str());
	if(!ofs.is_open()) {
		logging(ERROR,"Writing output failed!");
		exit(-1);
	}
	ofs << "SNP ID\tCHR\tPositions\tPValue" << endl;
	for(uint i=0; i<results.p_values.rows();i++) {
		ofs << results.chromosomes[i] + "_" + StringHelper::to_string<uint64>(results.positions[i]) << "\t"
		    << results.chromosomes[i] << "\t"
	    	    << results.positions[i] << "\t"
		    << results.p_values(i)	    
		    << endl;
	}
	ofs.close();
}
