#include "CPlinkParser.h"

#include "CEasyGWAS/io/CIOProgress.h"

#include <math.h>
#include <algorithm>
#include <iterator>
#include <sstream>

const std::map<std::string,char> CPlinkParser::__iupac_map = CPlinkParser::__init_map();

void CPlinkParser::readPEDFile(std::string const& file,
			       GWASData* data)
				throw (CPlinkParserException){
	std::map<std::string,char>::const_iterator iupac_iterator;
	std::ifstream ifs;
	data->family_ids.clear();
	data->sample_ids.clear();
	data->paternal_ids.clear();
	data->maternal_ids.clear();
	data->sex.clear();
	ifs.open(file.c_str(), std::ifstream::in);
	
	if(!ifs.is_open()) {
		throw CPlinkParserException("ERROR opening PED file " + file);
	}

	CIOProgress progress(ifs,1);
	uint64 fsize = progress.getFileSize();
	logging(INFO,"File Size: " + StringHelper::to_string<float64>(((float64)fsize)/1024/1024) + " MB");

	std::string line;
	uint64 i=0;
	uint64 size=0;
	uint64 lcounter=0;
	std::vector<std::string> sv;
	while(ifs.good()) {
		getline(ifs,line);
		line = StringHelper::trim(line);
		std::istringstream iss(line);
		sv.clear();
		std::copy(std::istream_iterator<std::string>(iss),
		       	  std::istream_iterator<std::string>(),
			  std::back_inserter<std::vector<std::string> >(sv));
		if(sv.size() < 8) continue;
		data->family_ids.push_back(sv[0]);
		data->sample_ids.push_back(sv[1]);
		data->paternal_ids.push_back(sv[2]);
		data->maternal_ids.push_back(sv[3]);
		data->sex.push_back(StringHelper::string_to<int>(sv[4]));
		if((sv.size()-6)%2!=0) 
			throw CPlinkParserException("ERROR: PED file has wrong file format.");
		std::vector<char> snps;
		size = (sv.size()-6.0)/2.0;
		snps.resize(size);
		lcounter=0;
		for(i=6; i<sv.size(); i++) {
			std::string iupac = sv[i] + sv[i+1];
			iupac_iterator = __iupac_map.find(iupac);
			if(iupac_iterator!=__iupac_map.end()) {
			       	snps.at(lcounter) = (char)(iupac_iterator->second);
			} else {
				throw CPlinkParserException("ERROR in PED Parser: Wrong Nucleotide encoding");
			}
			lcounter++;
			i++;
		}
		data->raw_snps.push_back(snps);
		progress.printProgress(ifs);
	}
	progress.printProgress(ifs);
	logging(LOG,"");
	data->n_snps = size;
	data->n_samples = data->raw_snps.size();
	ifs.close();
}

void CPlinkParser::readMAPFile(std::string const& file,
			       GWASData* data) 
				throw (CPlinkParserException) {
	std::ifstream ifs;
	ifs.open(file.c_str(), std::ifstream::in);
	
	if(!ifs.is_open()) {
		throw CPlinkParserException("ERROR opening MAP file " + file);
	}
	std::string line;
	uint64 i=0;
	bool ped=false;
	if(data->raw_snps.size()>0) {
		data->chromosomes.resize(data->raw_snps[0].size());
		data->positions.resize(data->raw_snps[0].size());
		data->snp_identifiers.resize(data->raw_snps[0].size());
		data->snp_distance.resize(data->raw_snps[0].size());
		ped = true;
	}
	std::vector<std::string> sv;
	while(ifs.good()) {
		getline(ifs,line);
		line = StringHelper::trim(line);
		std::istringstream iss(line);
		sv.clear();
		std::copy(std::istream_iterator<std::string>(iss),
		       	  std::istream_iterator<std::string>(),
			  std::back_inserter<std::vector<std::string> >(sv));
		if(sv.size() > 4 || sv.size() == 0) continue;
		if(ped) {
			data->chromosomes.at(i) = sv[0];
			data->snp_distance.at(i)  = StringHelper::string_to<float64>(sv[2]);
			data->positions.at(i)  = StringHelper::string_to<uint64>(sv[3]);
			data->snp_identifiers.at(i) = sv[1];
		} else {
			data->chromosomes.push_back(sv[0]);
			data->positions.push_back(StringHelper::string_to<uint64>(sv[3]));
			data->snp_distance.push_back(StringHelper::string_to<float64>(sv[2]));
			data->snp_identifiers.push_back(sv[1]);
		}
		i++;
	}
	ifs.close();
}

void CPlinkParser::readPhenotypeFile(std::string const& file,
  				     GWASData* data) 
				     throw (CPlinkParserException) {
	std::ifstream ifs;
	ifs.open(file.c_str(), std::ifstream::in);
	
	if(!ifs.is_open()) {
		throw CPlinkParserException("ERROR opening PED file " + file);
	}
	std::string line;
	uint64 i=0;
	std::vector<std::string> sv;
	std::vector<std::string>::iterator it;
	if(data->sample_ids.size()==0)
		throw CPlinkParserException("WARNING: Load SNP data befor the phenotype file!");
	while(ifs.good()) {
		getline(ifs,line);
		line = StringHelper::trim(line);
		std::istringstream iss(line);
		sv.clear();
		std::copy(std::istream_iterator<std::string>(iss),
		       	  std::istream_iterator<std::string>(),
			  std::back_inserter<std::vector<std::string> >(sv));
		if(sv.size() < 3) continue;
		if(i==0) {
			for(uint j=2; j<sv.size(); j++) {
				data->phenotype_names.push_back(sv[j]);
			}
			data->Y = MatrixXd::Zero(data->sample_ids.size(),sv.size()-2);
			data->Y.array() *= NAN;
		} else {
			for(uint j=2; j<sv.size(); j++) {
				it = std::find(data->sample_ids.begin(),data->sample_ids.end(),sv[1]);
				if(it!=data->sample_ids.end()) {
					uint64 index = it - data->sample_ids.begin();
					if(sv[j]=="nan") continue;
					data->Y(index,j-2) = StringHelper::string_to<float64>(sv[j]);
				}
			}
		}
		i++;
	}
	ifs.close();

}
