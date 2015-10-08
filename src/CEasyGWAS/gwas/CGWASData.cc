#include "CGWASData.h"
#include "CEasyGWAS/utils/CMatrixHelper.h"
#include <locale>
//#include <tr1/random>
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <cctype>

const uint CGWASDataHelper::dominant;
const uint CGWASDataHelper::codominant;
const uint CGWASDataHelper::recessive;
const uint CGWASDataHelper::additive;

void CGWASDataHelper::releaseMemory() {
    __data.X.resize(0,0);
    __data.Y.resize(0,1);
    __data.K.resize(0,0);
    __data.MAF.resize(0,1);
    __data.snp_hash.resize(0,1);
    std::vector<std::vector<char> >().swap(__data.raw_snps);
    std::vector<std::string>().swap(__data.sample_ids);
    std::vector<std::string>().swap(__data.chromosomes);
    std::vector<uint64>().swap(__data.positions);
    std::vector<float64>().swap(__data.snp_distance);
    std::vector<std::string>().swap(__data.snp_identifiers);
    std::vector<std::string>().swap(__data.family_ids);
    std::vector<std::string>().swap(__data.maternal_ids);
    std::vector<std::string>().swap(__data.paternal_ids);
    std::vector<uint>().swap(__data.sex);
}

MatrixXd CGWASDataHelper::getEncodedData() {
    return __data.X;
}

VectorXd CGWASDataHelper::getMAF() {
    return __data.MAF;
}

void CGWASDataHelper::encodeHomozygousData(std::vector< std::vector<char> > const& raw_snps, 
                       uint64 const& n_snps, 
                       uint64 const& n_samples) throw (CGWASDataException) {
    if(raw_snps.size()==0)
        throw CGWASDataException("Homozygous encoding error: GWASData object not initialized. raw_snps is empty.");
    std::locale loc;
    std::vector<uint> index_T;
    std::vector<uint> index_C;
    std::vector<uint> index_G;
    std::vector<uint> index_A;
    //Init X and MAF
    __data.n_samples = n_samples;
    __data.n_snps = n_snps;
    __data.genotype_data_type = "Homozygous";
    __data.genotype_encoding = "additive (0: Major, 2: Minor)";
    __data.X = MatrixXd::Ones(__data.n_samples,__data.n_snps);
    __data.X.array() *= 2;
    __data.MAF.resize(__data.n_snps);
    for(uint64 i=0; i<__data.n_snps; i++) {
        index_T.clear();
        index_C.clear();
        index_G.clear();
        index_A.clear();
        for(uint j=0; j<__data.n_samples; j++) {
            char nuc = std::toupper(raw_snps[j][i]);
            if(nuc == 'T') {
                index_T.push_back(j);
            } else if(nuc == 'C') {
                index_C.push_back(j);
            } else if(nuc == 'G') {
                index_G.push_back(j);
            } else if(nuc == 'A') {
                index_A.push_back(j);
            } else {
                std::string tmp(1,nuc);
                throw CGWASDataException("[Encode Homozygous Data]: Unknown nucleotide: " + tmp + " (Line: " + StringHelper::to_string(j+1) + " , SNP-Position: " + StringHelper::to_string(i));
            }
        }
        vector<uint64> max_elements;
        max_elements.push_back(index_T.size());
        max_elements.push_back(index_C.size());
        max_elements.push_back(index_G.size());
        max_elements.push_back(index_A.size());
        int index=0;
        uint64 max = 0;
        int n_elements = 0;
        for(uint j=0; j<max_elements.size(); j++) {
            if(max_elements[j]>0) {
                n_elements += 1;
            }
            if(max_elements[j]>max) {
                max = max_elements[j];
                index = j;
            }
        }
        if(n_elements>2)
            throw CGWASDataException("[Encode Homozygous Data]: Non-biallelic SNP found at SNP-Position: " + StringHelper::to_string(i));
        __data.MAF(i) = 1.0 - ((float64)max)/((float64)__data.n_samples);
        if(index==0) {
            for(uint64 j=0; j < index_T.size(); j++) {
                __data.X(index_T[j],i) = 0;
            }
        } else if(index==1) {
            for(uint64 j=0; j < index_C.size(); j++) {
                __data.X(index_C[j],i) = 0;
            }
        } else if(index==2) {
            for(uint64 j=0; j < index_G.size(); j++) {
                __data.X(index_G[j],i) = 0;
            }
        } else {
            for(uint64 j=0; j < index_A.size(); j++) {
                __data.X(index_A[j],i) = 0;
            }
        }
    }
}

void CGWASDataHelper::encodeHeterozygousData(std::vector< std::vector<char> > const& raw_snps, 
                         uint64 const& n_snps, 
                         uint64 const& n_samples,
                         uint const& encoding) throw (CGWASDataException) {
    if(raw_snps.size()==0)
        throw CGWASDataException("Heterozygous encoding error: GWASData object not initialized. raw_snps is empty.");
    std::locale loc;
    std::vector<uint> index_T;
    std::vector<uint> index_C;
    std::vector<uint> index_G;
    std::vector<uint> index_A;
    std::vector<uint> index_R;
    std::vector<uint> index_Y;
    std::vector<uint> index_S;
    std::vector<uint> index_W;
    std::vector<uint> index_K;
    std::vector<uint> index_M;
    //Init X and MAF
    __data.n_samples = n_samples;
    __data.n_snps = n_snps;
    __data.X = MatrixXd::Ones(__data.n_samples,__data.n_snps);
    uint het_encoding = 1;
    __data.genotype_data_type = "Heterozygous";
    if(encoding==CGWASDataHelper::additive) {
               __data.genotype_encoding = "additive (0: Major, 1: Heterozygous, 2: Minor)";
        __data.X.array() *= 2;
        het_encoding = 1;
    } else if(encoding==CGWASDataHelper::dominant) {
               __data.genotype_encoding = "dominant (0: Major, 1: Heterozygous, 1: Minor)";
        __data.X.array() *= 1;
        het_encoding = 1;
    } else if(encoding==CGWASDataHelper::codominant) {
        __data.genotype_encoding = "co-dominant (0: Major, 1: Heterozygous, 0: Minor)";
        __data.X.array() *= 0;
        het_encoding = 1;
    } else if(encoding==CGWASDataHelper::recessive) {
               __data.genotype_encoding = "recessive (0: Major, 0: Heterozygous, 1: Minor)";
        __data.X.array() *= 1;
        het_encoding = 0;
    } else {
        throw CGWASDataException("Encoding is wrong. Has to be [0: additive, 1: recessive, 2: dominant, 3: codominant]");
    }
    
    __data.MAF.resize(__data.n_snps);
    for(uint64 i=0; i<__data.n_snps; i++) {
        index_T.clear();
        index_C.clear();
        index_G.clear();
        index_A.clear();
        index_R.clear();
        index_Y.clear();
        index_S.clear();
        index_W.clear();
        index_K.clear();
        index_M.clear();
        for(uint j=0; j<__data.n_samples; j++) {
            char nuc = std::toupper(raw_snps[j][i]);
            if(nuc == 'T') {
                index_T.push_back(j);
            } else if(nuc == 'C') {
                index_C.push_back(j);
            } else if(nuc == 'G') {
                index_G.push_back(j);
            } else if(nuc == 'A') {
                index_A.push_back(j);
            } else if(nuc == 'R') {
                index_R.push_back(j);
            } else if(nuc == 'Y') {
                index_Y.push_back(j);
            } else if(nuc == 'S') {
                index_S.push_back(j);
            } else if(nuc == 'W') {
                index_W.push_back(j);
            } else if(nuc == 'K') {
                index_K.push_back(j);
            } else if(nuc == 'M') {
                index_M.push_back(j);
            } else {
                std::string tmp(1,nuc);
                throw CGWASDataException("[Encode Heterozygous Data]: Unknown nucleotide: " + tmp + " (Line: " + StringHelper::to_string(j+1) + " , SNP-Position: " + StringHelper::to_string(i));
            }
        }
        vector<uint64> max_elements;
        max_elements.push_back(index_T.size());
        max_elements.push_back(index_C.size());
        max_elements.push_back(index_G.size());
        max_elements.push_back(index_A.size());
        int index=0;
        uint64 max = 0;
        uint64 min = ULONG_MAX;
        int zero_counter=0;
        int n_elements = 0;
        for(uint j=0; j<max_elements.size(); j++) {
            if(max_elements[j]>0) {
                n_elements += 1;
            }
            if(max_elements[j]>max) {
                max = max_elements[j];
                index = j;
            }
            if(max_elements[j]==0) zero_counter++;
            else if(max_elements[j]<min) {
                min = max_elements[j];
            }
        }
        if(zero_counter==3) min = 0;
        uint64 het_count = __data.n_samples - max - min;
        if(n_elements>2)
            throw CGWASDataException("[Encode Heterozygous Data]: Non-biallelic SNP found at SNP-Position: " + StringHelper::to_string(i));
        //__data.MAF(i) = 1.0 - ((float64)max)/((float64)__data.n_samples);
        __data.MAF(i) = ((float64)min)/((float64)__data.n_samples) + 0.5 * ((float64)het_count)/((float64)__data.n_samples);
        if(index==0) {
            for(uint64 j=0; j < index_T.size(); j++) {
                __data.X(index_T[j],i) = 0;
            }
        } else if(index==1) {
            for(uint64 j=0; j < index_C.size(); j++) {
                __data.X(index_C[j],i) = 0;
            }
        } else if(index==2) {
            for(uint64 j=0; j < index_G.size(); j++) {
                __data.X(index_G[j],i) = 0;
            }
        } else if(index==3) {
            for(uint64 j=0; j < index_A.size(); j++) {
                __data.X(index_A[j],i) = 0;
            }
        }
        for(uint64 j=0; j < index_R.size(); j++) {
            __data.X(index_R[j],i) = het_encoding;
        }
        for(uint64 j=0; j < index_Y.size(); j++) {
            __data.X(index_Y[j],i) = het_encoding;
        }
        for(uint64 j=0; j < index_S.size(); j++) {
            __data.X(index_S[j],i) = het_encoding;
        }
        for(uint64 j=0; j < index_W.size(); j++) {
            __data.X(index_W[j],i) = het_encoding;
        }
        for(uint64 j=0; j < index_K.size(); j++) {
            __data.X(index_K[j],i) = het_encoding;
        }
        for(uint64 j=0; j < index_M.size(); j++) {
            __data.X(index_M[j],i) = het_encoding;
        }
    }

}

void CGWASDataHelper::encodeHomozygousData(GWASData* data) throw (CGWASDataException) {
    if(data->raw_snps.size()==0)
        throw CGWASDataException("Homozygous encoding error: GWASData object not initialized. raw_snps is empty.");
    std::locale loc;
    std::vector<uint> index_T;
    std::vector<uint> index_C;
    std::vector<uint> index_G;
    std::vector<uint> index_A;
    //Init X and MAF
    data->genotype_data_type = "Homozygous";
    data->genotype_encoding = "additive (0: Major, 2: Minor)";
    data->X = MatrixXd::Ones(data->n_samples,data->n_snps);
    data->X.array() *= 2;
    data->MAF.resize(data->n_snps);
    for(uint64 i=0; i<data->n_snps; i++) {
        index_T.clear();
        index_C.clear();
        index_G.clear();
        index_A.clear();
        for(uint j=0; j<data->n_samples; j++) {
            char nuc = std::toupper(data->raw_snps[j][i]);
            if(nuc == 'T') {
                index_T.push_back(j);
            } else if(nuc == 'C') {
                index_C.push_back(j);
            } else if(nuc == 'G') {
                index_G.push_back(j);
            } else if(nuc == 'A') {
                index_A.push_back(j);
            } else {
                std::string tmp(1,nuc);
                throw CGWASDataException("[Encode Homozygous Data]: Unknown nucleotide: " + tmp + " (Line: " + StringHelper::to_string(j+1) + " , SNP-Position: " + StringHelper::to_string(i));
            }
        }
        vector<uint64> max_elements;
        max_elements.push_back(index_T.size());
        max_elements.push_back(index_C.size());
        max_elements.push_back(index_G.size());
        max_elements.push_back(index_A.size());
        int index=0;
        uint64 max = 0;
        int n_elements = 0;
        for(uint j=0; j<max_elements.size(); j++) {
            if(max_elements[j]>0) {
                n_elements += 1;
            }
            if(max_elements[j]>max) {
                max = max_elements[j];
                index = j;
            }
        }
        if(n_elements>2)
            throw CGWASDataException("[Encode Homozygous Data]: Non-biallelic SNP found at SNP-Position: " + StringHelper::to_string(i));
        data->MAF(i) = 1.0 - ((float64)max)/((float64)data->n_samples);
        if(index==0) {
            for(uint64 j=0; j < index_T.size(); j++) {
                data->X(index_T[j],i) = 0;
            }
        } else if(index==1) {
            for(uint64 j=0; j < index_C.size(); j++) {
                data->X(index_C[j],i) = 0;
            }
        } else if(index==2) {
            for(uint64 j=0; j < index_G.size(); j++) {
                data->X(index_G[j],i) = 0;
            }
        } else {
            for(uint64 j=0; j < index_A.size(); j++) {
                data->X(index_A[j],i) = 0;
            }
        }
    }
}

void CGWASDataHelper::encodeHeterozygousData(GWASData* data) throw (CGWASDataException) {
    encodeHeterozygousData(data,CGWASDataHelper::additive);
}

void CGWASDataHelper::encodeHeterozygousData(GWASData* data, uint const& encoding) throw (CGWASDataException) {
    if(data->raw_snps.size()==0)
        throw CGWASDataException("Heterozygous encoding error: GWASData object not initialized. raw_snps is empty.");
    std::locale loc;
    std::vector<uint> index_T;
    std::vector<uint> index_C;
    std::vector<uint> index_G;
    std::vector<uint> index_A;
    std::vector<uint> index_R;
    std::vector<uint> index_Y;
    std::vector<uint> index_S;
    std::vector<uint> index_W;
    std::vector<uint> index_K;
    std::vector<uint> index_M;
    //Init X and MAF
    data->X = MatrixXd::Ones(data->n_samples,data->n_snps);
    uint het_encoding = 1;
    data->genotype_data_type = "Heterozygous";
    if(encoding==CGWASDataHelper::additive) {
               data->genotype_encoding = "additive (0: Major, 1: Heterozygous, 2: Minor)";
        data->X.array() *= 2;
        het_encoding = 1;
    } else if(encoding==CGWASDataHelper::dominant) {
               data->genotype_encoding = "dominant (0: Major, 1: Heterozygous, 1: Minor)";
        data->X.array() *= 1;
        het_encoding = 1;
    } else if(encoding==CGWASDataHelper::codominant) {
        data->genotype_encoding = "co-dominant (0: Major, 1: Heterozygous, 0: Minor)";
        data->X.array() *= 0;
        het_encoding = 1;
    } else if(encoding==CGWASDataHelper::recessive) {
               data->genotype_encoding = "recessive (0: Major, 0: Heterozygous, 1: Minor)";
        data->X.array() *= 1;
        het_encoding = 0;
    }
    
    data->MAF.resize(data->n_snps);
    for(uint64 i=0; i<data->n_snps; i++) {
        index_T.clear();
        index_C.clear();
        index_G.clear();
        index_A.clear();
        index_R.clear();
        index_Y.clear();
        index_S.clear();
        index_W.clear();
        index_K.clear();
        index_M.clear();
        for(uint j=0; j<data->n_samples; j++) {
            char nuc = std::toupper(data->raw_snps[j][i]);
            if(nuc == 'T') {
                index_T.push_back(j);
            } else if(nuc == 'C') {
                index_C.push_back(j);
            } else if(nuc == 'G') {
                index_G.push_back(j);
            } else if(nuc == 'A') {
                index_A.push_back(j);
            } else if(nuc == 'R') {
                index_R.push_back(j);
            } else if(nuc == 'Y') {
                index_Y.push_back(j);
            } else if(nuc == 'S') {
                index_S.push_back(j);
            } else if(nuc == 'W') {
                index_W.push_back(j);
            } else if(nuc == 'K') {
                index_K.push_back(j);
            } else if(nuc == 'M') {
                index_M.push_back(j);
            } else {
                std::string tmp(1,nuc);
                throw CGWASDataException("[Encode Heterozygous Data]: Unknown nucleotide: " + tmp + " (Line: " + StringHelper::to_string(j+1) + " , SNP-Position: " + StringHelper::to_string(i));
            }
        }
        vector<uint64> max_elements;
        max_elements.push_back(index_T.size());
        max_elements.push_back(index_C.size());
        max_elements.push_back(index_G.size());
        max_elements.push_back(index_A.size());
        int index=0;
        uint64 max = 0;
        uint64 min = ULONG_MAX;
        int zero_counter=0;
        int n_elements = 0;
        for(uint j=0; j<max_elements.size(); j++) {
            if(max_elements[j]>0) {
                n_elements += 1;
            }
            if(max_elements[j]>max) {
                max = max_elements[j];
                index = j;
            }
            if(max_elements[j]==0) zero_counter++;
            else if(max_elements[j]<min) {
                max = max_elements[j];
                index = j;
            }
        }
        if(zero_counter==3) min = 0;
        uint64 het_count = data->n_samples - max - min;
        if(n_elements>2)
            throw CGWASDataException("[Encode Heterozygous Data]: Non-biallelic SNP found at SNP-Position: " + StringHelper::to_string(i));
        data->MAF(i) = ((float64)min)/((float64)data->n_samples) + 0.5 * ((float64)het_count)/((float64)data->n_samples);
        //data->MAF(i) = 1.0 - ((float64)max)/((float64)data->n_samples);
        if(index==0) {
            for(uint64 j=0; j < index_T.size(); j++) {
                data->X(index_T[j],i) = 0;
            }
        } else if(index==1) {
            for(uint64 j=0; j < index_C.size(); j++) {
                data->X(index_C[j],i) = 0;
            }
        } else if(index==2) {
            for(uint64 j=0; j < index_G.size(); j++) {
                data->X(index_G[j],i) = 0;
            }
        } else if(index==3) {
            for(uint64 j=0; j < index_A.size(); j++) {
                data->X(index_A[j],i) = 0;
            }
        }
        for(uint64 j=0; j < index_R.size(); j++) {
            data->X(index_R[j],i) = het_encoding;
        }
        for(uint64 j=0; j < index_Y.size(); j++) {
            data->X(index_Y[j],i) = het_encoding;
        }
        for(uint64 j=0; j < index_S.size(); j++) {
            data->X(index_S[j],i) = het_encoding;
        }
        for(uint64 j=0; j < index_W.size(); j++) {
            data->X(index_W[j],i) = het_encoding;
        }
        for(uint64 j=0; j < index_K.size(); j++) {
            data->X(index_K[j],i) = het_encoding;
        }
        for(uint64 j=0; j < index_M.size(); j++) {
            data->X(index_M[j],i) = het_encoding;
        }
    }
}

void CGWASDataHelper::filterNonInformativeSNPs(GWASData* data) throw (CGWASDataException) {
    if(data->X.cols()==0)
        throw CGWASDataException("GWAS Data object not initialized or SNP data not encoded");
    std::vector<uint64> indices_v;
    std::vector<std::string> chromosomes;
    std::vector<std::string> snp_identifiers;
    std::vector<uint64> positions;
    std::vector<float64> snp_distance;
    VectorXd snp_sum = data->X.colwise().sum().array()/2;
    for(uint i=0; i<snp_sum.rows(); i++) {
        if(!(snp_sum(i)==0 || snp_sum(i)==data->n_samples)) {
            indices_v.push_back(i);
            chromosomes.push_back(data->chromosomes[i]);
            snp_identifiers.push_back(data->snp_identifiers[i]);
            positions.push_back(data->positions[i]);
            snp_distance.push_back(data->snp_distance[i]);
        }
    }    
    
    MatrixXd X(data->n_samples,indices_v.size());
    VectorXd MAF(indices_v.size());
    for(uint64 i=0; i<indices_v.size();i++) {
        X.col(i) = data->X.col(indices_v[i]);
        MAF(i) = data->MAF(indices_v[i]);
    }
    std::vector<std::vector<char> > raw_snps;
    for(uint64 j=0; j<data->n_samples; j++) {
        std::vector<char> snps;
        snps.resize(indices_v.size());
        for(uint64 i=0; i<indices_v.size();i++) {
            snps.at(i) = data->raw_snps[j][indices_v[i]];
        }
        raw_snps.push_back(snps);
    }
    data->raw_snps = raw_snps;
    data->removed_snp_indices = indices_v;
    data->X = X;
    data->MAF = MAF;
    data->chromosomes = chromosomes;
    data->positions = positions;
    data->snp_distance = snp_distance;
    data->snp_identifiers = snp_identifiers;
    data->n_snps = indices_v.size();
}

void CGWASDataHelper::filterSNPsByMAF(GWASData* data, float64 const& maf) throw (CGWASDataException) {
    if(data->MAF.size()==0)
        throw CGWASDataException("Data object not initialized or SNP data not encoded!");
    std::vector<uint64> maf_indices_v;
    std::vector<std::string> chromosomes;
    std::vector<std::string> snp_identifiers;
    std::vector<uint64> positions;
    std::vector<float64> snp_distance;
    for(uint i=0; i<data->MAF.rows(); i++) {
        if(data->MAF(i)>=maf) {
            maf_indices_v.push_back(i);
            chromosomes.push_back(data->chromosomes[i]);
            snp_identifiers.push_back(data->snp_identifiers[i]);
            positions.push_back(data->positions[i]);
            snp_distance.push_back(data->snp_distance[i]);
        }
    }    
    MatrixXd X(data->n_samples,maf_indices_v.size());
    VectorXd MAF(maf_indices_v.size());
    for(uint64 i=0; i<maf_indices_v.size();i++) {
        X.col(i) = data->X.col(maf_indices_v[i]);
        MAF(i) = data->MAF(maf_indices_v[i]);
    }
    std::vector<std::vector<char> > raw_snps;
    for(uint64 j=0; j<data->n_samples; j++) {
        std::vector<char> snps;
        snps.resize(maf_indices_v.size());
        for(uint64 i=0; i<maf_indices_v.size();i++) {
            snps.at(i) = data->raw_snps[j][maf_indices_v[i]];
        }
        raw_snps.push_back(snps);
    }
    data->raw_snps = raw_snps;
    data->removed_snp_indices = maf_indices_v;
    data->X = X;
    data->MAF = MAF;
    data->chromosomes = chromosomes;
    data->positions = positions;
    data->snp_distance = snp_distance;
    data->snp_identifiers = snp_identifiers;
    data->n_snps = maf_indices_v.size();
}

/*
void CGWASDataHelper::filterSNPsBySmallIndel(GWASData* data, int const& indel) throw (CGWASDataException) {
    if(data->small_indel.size()==0)
        throw CGWASDataException("Data object not initialized or SNP data not encoded!");
    if(indel!=0 && indel!=1) 
        throw CGWASDataException("Indel flag can be only 0: Remove all Indels, or 1: Remove all SNPs");

    std::vector<uint64> indel_indices_v;
    std::vector<std::string> chromosomes;
    std::vector<std::string> snp_identifiers;
    std::vector<uint64> positions;
    std::vector<float64> snp_distance;
    for(uint i=0; i<data->small_indel.rows(); i++) {
        if(data->small_indel(i)==indel) {
            indel_indices_v.push_back(i);
            chromosomes.push_back(data->chromosomes[i]);
            snp_identifiers.push_back(data->snp_identifiers[i]);
            positions.push_back(data->positions[i]);
            snp_distance.push_back(data->snp_distance[i]);
        }
    }    
    MatrixXd X(data->n_samples,indel_indices_v.size());
    VectorXd MAF(indel_indices_v.size());
    VectorXd small_indel(indel_indices_v.size());
    for(uint64 i=0; i<indel_indices_v.size();i++) {
        X.col(i) = data->X.col(indel_indices_v[i]);
        MAF(i) = data->MAF(indel_indices_v[i]);
        small_indel(i) = data->small_indel(indel_indices_v[i]);
    }
    std::vector<std::vector<char> > raw_snps;
    for(uint64 j=0; j<data->n_samples; j++) {
        std::vector<char> snps;
        snps.resize(indel_indices_v.size());
        for(uint64 i=0; i<indel_indices_v.size();i++) {
            snps.at(i) = data->raw_snps[j][indel_indices_v[i]];
        }
        raw_snps.push_back(snps);
    }
    data->raw_snps = raw_snps;
    data->removed_snp_indices = indel_indices_v;
    data->X = X;
    data->MAF = MAF;
    data->small_indel = small_indel;
    data->chromosomes = chromosomes;
    data->positions = positions;
    data->snp_distance = snp_distance;
    data->snp_identifiers = snp_identifiers;
    data->n_snps = indel_indices_v.size();
}
*/
        
//TODO: Implement filtering method
void CGWASDataHelper::filterUniqueSNPs(GWASData* data) throw (CGWASDataException) {
    VectorXd randv = VectorXd::Random(data->X.rows());
    logging(INFO,randv.rows());
    logging(INFO,data->X.rows());
    VectorXd hash = (data->X.array()*randv.replicate(1,data->X.cols()).array()).colwise().sum();
    logging(STATUS,hash.rows());
}

void CGWASDataHelper::createSNPHash(GWASData* data) throw (CGWASDataException) {
    /*VectorXd snp_hash = VectorXd::Zero(data->X.cols());
    data->n_unique_snps = 0;
    snp_hash(0) = 1;
    for(int64 i=0; i<data->X.cols();i++) {
        for(int64 j=i+1; j<data->X.cols();j++) {
            if((data->X.col(i)-data->X.col(j)).sum()==0) {
                snp_hash(j)=snp_hash(i);
            } else {
                data->n_unique_snps++;    
                snp_hash(j)=data->n_unique_snps;
            }
        }
    }
    logging(INFO,snp_hash.transpose());
    logging(STATUS,data->n_unique_snps);
    */
    
    //Accurate but extremely slowely!!
    /*VectorXd snp_hash = VectorXd::Ones(data->X.cols()).array()*-1;
    for(int64 i=0; i<data->X.cols(); i++) {
        VectorXd tmpSNP = data->X.col(i);
        if (snp_hash(i)==-1) {
            snp_hash(i) = i;
            for(int64 j=i+1;j<data->X.cols();j++) {
                if(i==j) continue;
                if((tmpSNP.array()-data->X.col(j).array()).array().abs().array().sum()==0) {
                    snp_hash(j) = snp_hash(i);
                }
            }
        }
    }
    data->snp_hash = snp_hash;
    std::vector<float64> tmp(snp_hash.rows());
    Eigen::Map<VectorXd>(tmp.data(),snp_hash.rows(),1) = snp_hash;
    std::stable_sort(tmp.begin(),tmp.end());
    std::vector<float64>::iterator it;
    it = std::unique(tmp.begin(),tmp.end());
    data->n_unique_snps = std::distance(tmp.begin(),it);
    */
    
    //FAST: Please check if accurate
    VectorXd hash = VectorXd::Ones(data->X.rows());
    for(int64 i=0; i<data->X.rows(); i++) {
        //float64 r = randn();
        float64 r = fabs(rand()/(float64)(RAND_MAX) + 1);
        uint found = (hash.array()==r).any();
        while(found>0) {
            r = fabs(rand()/(float64)(RAND_MAX) + 1);
            found = (hash.array()==r).any();
        }
        hash(i) = r;
    }
    VectorXd snp_hash = (data->X.array()*hash.replicate(1,data->X.cols()).array()).colwise().sum();
    data->snp_hash = snp_hash;
    std::vector<float64> tmp(snp_hash.rows());
    Eigen::Map<VectorXd>(tmp.data(),snp_hash.rows(),1) = snp_hash;
    std::stable_sort(tmp.begin(),tmp.end());
    std::vector<float64>::iterator it;
    it = std::unique(tmp.begin(),tmp.end());
    data->n_unique_snps = std::distance(tmp.begin(),it);
}

GWASData CGWASDataHelper::removeSamples4MissingData(GWASData const& data, uint const& phenotype_id) throw (CGWASDataException) {
    return removeSamples4MissingData(data,phenotype_id,false);
}

GWASData CGWASDataHelper::removeSamples4MissingData(GWASData const& data, uint const& phenotype_id, bool const copy_raw) throw (CGWASDataException) {
    GWASData newData;
    std::vector<float64> splice_indices;
    std::vector<std::string> sample_ids;
    std::vector<std::string> family_ids;
    std::vector<std::string> paternal_ids;
    std::vector<std::string> maternal_ids;
    std::vector<uint> sex;
    std::vector<std::string> phenotype_names;
    for(uint i=0; i<data.X.rows(); i++) {
        if(isnan(data.Y(i,phenotype_id))==false) {
            splice_indices.push_back(i);
            sample_ids.push_back(data.sample_ids[i]);
            family_ids.push_back(data.family_ids[i]);
            paternal_ids.push_back(data.paternal_ids[i]);
            maternal_ids.push_back(data.maternal_ids[i]);
            sex.push_back(data.sex[i]);
        }
    }
    phenotype_names.push_back(data.phenotype_names[phenotype_id]);
    Eigen::Map<VectorXd> sindices(&splice_indices[0],splice_indices.size());
    newData.snp_identifiers = data.snp_identifiers;
    newData.chromosomes = data.chromosomes;
    newData.positions = data.positions;
    newData.snp_distance = data.snp_distance;
    newData.n_samples = splice_indices.size();
    newData.n_snps = data.n_snps;

    newData.X = sliceRowsMatrix(data.X,sindices);
    newData.Y = sliceRowsMatrix(data.Y,sindices);
    
    if (copy_raw==true) {
        std::vector<std::vector<char> > raw_snps;
        for(uint64 j=0; j<sindices.size(); j++) {
            std::vector<char> snps;
            snps.resize(data.n_snps);
            for(uint64 i=0; i<data.n_snps;i++) {
                snps.at(i) = data.raw_snps[sindices[j]][i];
            }
            raw_snps.push_back(snps);
        }
        newData.raw_snps = raw_snps;
    } else {
        newData.raw_snps = std::vector<std::vector<char> >();
    }

    if(newData.K.rows()>0 && newData.K.cols()>0) {
        newData.K = MatrixXd::Zero(newData.n_samples,newData.n_samples);
        newData.K = sliceRowsMatrix(data.K,sindices);
        newData.K = sliceColsMatrix(newData.K,sindices);
    }
    newData.network = data.network;
    newData.sample_ids = sample_ids;
    newData.family_ids = family_ids;
    newData.paternal_ids = paternal_ids;
    newData.maternal_ids = maternal_ids;
    newData.sex = sex;
    newData.phenotype_names = phenotype_names;
    newData.MAF = __updateMAF(newData.X);;
    return newData;
}

VectorXd CGWASDataHelper::__updateMAF(MatrixXd const& m) {
    uint n = m.cols();
    VectorXd maf(n);
    for(uint i=0; i<n;i++) {
        maf(i) = (m.col(i).array()==2).count()/100.0;
    }
    return maf;
}
