#include "CSconesIO.h"

#include "CEasyGWAS/io/CIOProgress.h"
#include "CEasyGWAS/utils/StringHelper.h"

#include <sstream>
#include <vector>
#include <iterator>
#include <algorithm>

#include <Eigen/Dense>

void CSconesIO::readSparseNetworkFile(std::string const& file, GWASData* data) throw (CSconesIOException){
	std::ifstream ifs;
	ifs.open(file.c_str(),std::ifstream::in);

	if(!ifs.is_open()) 
		throw CSconesIOException("Error opening sparse network file " + file);
	CIOProgress progress(ifs,1);
	uint fsize = progress.getFileSize();
	logging(INFO,"File Size: " + StringHelper::to_string<float64>(((float64)fsize)/1024.0/1024.0) + " MB");
	//Create temprory positions map
	std::map<std::string,uint64> position_map;
	std::map<std::string,uint64>::iterator piter;
	for(uint64 i=0; i<data->positions.size(); i++) {
		std::string idp = data->chromosomes[i] + "_" +StringHelper::to_string<uint64>(data->positions[i]);
		position_map[idp] = i;
	}

	std::string line;
	std::vector<std::string> sv;
	SparseMatrixXd L(data->n_snps,data->n_snps);
	sparse_triplet_vector triplet;
	uint64 pos1;
	uint64 pos2;
	std::string idp1;
	std::string idp2;
	uint64 index1;

	while(ifs.good()) {
		getline(ifs,line);
		line = StringHelper::trim(line);
		std::istringstream iss(line);
		sv.clear();
		std::copy(std::istream_iterator<std::string>(iss),
			  std::istream_iterator<std::string>(),
			  std::back_inserter<std::vector<std::string> >(sv));
		if(sv.size()!=4) continue;
		pos1 = StringHelper::string_to<uint64>(sv[1]);
		idp1 = sv[0] + "_" + StringHelper::to_string<uint64>(pos1);
		pos2 = StringHelper::string_to<uint64>(sv[3]);
		idp2 = sv[2] + "_" + StringHelper::to_string<uint64>(pos2);
		//check if positions in GWASData
		piter = position_map.find(idp1);
		if(piter!=position_map.end()) {
			index1 = piter->second;
			piter = position_map.find(idp2);
			if(piter!=position_map.end()) {
				//if pos positions exist store interaction in sparse Matrix
				triplet.push_back(eigen_triplet(index1,
								piter->second,
								1));
			}
		}
		progress.printProgress(ifs);
	}
	//Storing data
	L.setFromTriplets(triplet.begin(),triplet.end());
	data->network = L;
	logging(LOG,"");
	ifs.close();
}

void CSconesIO::writeOutput(std::string const& outfile, GWASData const& data, VectorXd const& indicator, float64 const& best_lambda, float64 const& best_eta) {
	std::ofstream ofs;
	ofs.open(outfile.c_str());
	if(!ofs.is_open()) {
		logging(ERROR,"Writing output failed!");
		exit(-1);
	}
    ofs << "#Best Lambda:\t" << best_lambda << "\n";
    ofs << "#Best Eta:\t" << best_eta << "\n";
	ofs << "#Selected SNP ID\tCHR\tPositions" << "\n";
    for(uint i=0; i<indicator.rows();i++) {
        if(indicator(i)>=1) {
            ofs << data.snp_identifiers[i] << "\t"
                << data.chromosomes[i] << "\t"
                << data.positions[i] << "\n";
        }
    }
	ofs.close();
}

void CSconesIO::writeCMatrix(std::string const& outfile, MatrixXd const& cmat, CSconesSettings const& settings) {
	std::ofstream ofs;
	ofs.open(outfile.c_str());
	if(!ofs.is_open()) {
		logging(ERROR,"Writing output failed!");
		exit(-1);
	}
    ofs << "\t";
    for(int j=0;j<settings.lambdas.rows();j++) {
        if(j==settings.lambdas.rows()-1) ofs << settings.lambdas(j);
        else ofs << settings.lambdas(j) << "\t";
    }
    ofs << "\n";
    for(int i=0; i<cmat.rows(); i++) {
        ofs << settings.etas(i) << "\t";
        for(int j=0; j<cmat.cols();j++) {
            if(j==cmat.cols()-1) ofs << cmat(i,j);
            else ofs << cmat(i,j) << "\t";
        }
        ofs << "\n";
    }
    ofs.close();
}
