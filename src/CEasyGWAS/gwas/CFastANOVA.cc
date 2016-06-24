/*
Reimplementation of FastANOVA
Paper Link: http://www.ncbi.nlm.nih.gov/pubmed/20945829
*/

#include "CFastANOVA.h"

#include <list>
#include <vector>
#include <algorithm>
#include <time.h>

#include <fstream>

#include "CEasyGWAS/utils/StringHelper.h"
#include "CEasyGWAS/stats/CFisherF.h"

namespace CEpistasis {

CFastANOVA::CFastANOVA() {
	__tm = 0.0;
	__permutations = 100;
	__fwer = 0.05;
	__theta = 0;
	__alphaK = floor(__permutations*__fwer);
}

CFastANOVA::CFastANOVA(VectorXd const& y, MatrixXd const& X) throw (CEpistasisException) {
	__y = y;
	__X = X;
	__permutations = 100;
	__fwer = 0.05;
	__theta = 0;
	__alphaK = floor(__permutations*__fwer);
	__checkdata();
}

CFastANOVA::CFastANOVA(VectorXd const& y, MatrixXd const& X, uint const& perm) throw (CEpistasisException) {
	__y = y;
	__X = X;
	__permutations = perm;
	__fwer = 0.05;
	__theta = 0;
	__alphaK = floor(__permutations*__fwer);
	__checkdata();
}

CFastANOVA::CFastANOVA(VectorXd const& y, MatrixXd const& X, uint const& perm, float64 const& fwer) throw (CEpistasisException) {
	__y = y;
	__X = X;
	__permutations = perm;
	__fwer = fwer;
	__theta = 0;
	__alphaK = floor(__permutations*__fwer);
	__checkdata();
}

CFastANOVA::CFastANOVA(VectorXd const& y, MatrixXd const& X, uint const& perm, float64 const& fwer, float64 const& seed) throw (CEpistasisException) {
	__y = y;
	__X = X;
	__permutations = perm;
	__fwer = fwer;
	__theta = 0;
	__alphaK = floor(__permutations*__fwer);
	__seed = seed;
	srand(__seed);
	__checkdata();
}

void CFastANOVA::__checkdata() throw (CEpistasisException) {
	if(__y.cols()>1) throw CEpistasisException("Phenotype y has wrong dimensions! (n x 1)");
	if(__X.rows() != __y.rows()) throw CEpistasisException("Genotype X and Phenotype y must have the sam enumer of samples n!");
	if(__fwer <= 0 || __fwer >=1) throw CEpistasisException("Family-wise error rate fwer has to be between 0 and 1 and cannot be 0!");
	std::string msg = "\t<Too few iterations to determine a threshold F_alpha! For a FWER-Threshold of " + 
				StringHelper::to_string<float>(__fwer) + " one needs at least " + StringHelper::to_string<float>(1.0/__fwer) + " permutations!" +
				" Bruteforce method is performed!>";
	if(__alphaK==0) logging(WARNING, msg); 
	//Check if why is continuous 
	
	//Check if X is binary
}

//Getter and Setter
void CFastANOVA::setPhenotype(VectorXd const& y) throw (CEpistasisException) {
	__y = y;
	if(__y.cols()>1) throw CEpistasisException("Phenotype y has wrong dimensions! (n x 1)");
}

void CFastANOVA::setGenotype(MatrixXd const& X) {
	__X = X;
}

void CFastANOVA::setFWER(float64 const& fwer) throw (CEpistasisException) {
	__fwer = fwer;
	if(__fwer <= 0 || __fwer >=1) throw CEpistasisException("Family-wise error rate fwer has to be between 0 and 1 and cannot be 0!");
}

void CFastANOVA::setPermutations(uint const& permutations) {
	__permutations = permutations;
	__alphaK = floor(__permutations*__fwer);
	std::string msg = "\t<Too few iterations to determine a threshold F_alpha! For a FWER-Threshold of " + 
				StringHelper::to_string<float>(__fwer) + " one needs at least " + StringHelper::to_string<float>(1.0/__fwer) + " permutations!" +
				" Bruteforce method is performed!>";
	if(__alphaK==0) logging(WARNING, msg); 
}

void CFastANOVA::setSeed(float64 const& seed) {
	__seed = seed;
	srand(__seed);
}

VectorXd CFastANOVA::getPValues() {
	return __pvalues;
}

VectorXd CFastANOVA::getTestStatistics() {
	return __test_statistics;
}

MatrixXd CFastANOVA::getSNPPairs() {
	return __snp_pairs;
}

float64 CFastANOVA::getFalpha() {
	return __theta;
}

float64 CFastANOVA::getAlphaPvalue() {
	return __theta_pval;
}

uint64 CFastANOVA::getNumVisitedPairs() {
	return __visited;
}

/*
*Compute F-Statistic for two interacting SNPs
*/
float64 CFastANOVA::__computeFStatistic(VectorXd const& y, VectorXd const& x1, VectorXd const& x2) {
	float64 T[4] = {0};
	float64 n[4] = {0};
	for(uint i=0; i<__n_samples; i++) {
		if(x1(i)==1) {
			if(x2(i)==1) {
				T[0] += y(i);
				n[0]++;
			} else {
				T[1] += y(i);
				n[1]++;
			}
		} else {
			if(x2(i)==1) {
				T[2] += y(i);
				n[2]++;
			} else {
				T[3] += y(i);
				n[3]++;
			}
		}
	}

	//compute SSB, between group variance
	float64 SSB = 0.0;
	__df1 = 4;
	for(uint i=0; i<4; i++) {
		if(n[i]>0) SSB += (T[i]*T[i])/n[i];
		else __df1 = 3;
	}
	__df2 = __n_samples-__df1;
	__df1--;
	SSB = SSB - __tm;
	return (__SST-SSB!=0)?(__df2)/(__df1) * SSB/(__SST-SSB):0;
}

/*
*Compute Contingency Table for 2 groups 
*/
void CFastANOVA::__computeContingencyTable2(VectorXd const& phenoP, uint64 const& i) {
	__contingencyTable2.clear();
	VectorXd x1 = __X.col(i);
	std::vector<float64> g1;
	std::vector<float64> g2;
	for(uint i=0; i<__n_samples;i++) {
		if(x1(i)==1) g1.push_back(phenoP(i));
		if(x1(i)==0) g2.push_back(phenoP(i));
	}
	//sort in ascending order	
	std::sort(g1.begin(),g1.end());
	std::sort(g2.begin(),g2.end());
	//Map stl vector into a eigen matrix
	Eigen::Map<VectorXd> g1v(&g1[0],g1.size());
	Eigen::Map<VectorXd> g2v(&g2[0],g2.size());
	//Store Eigen vectors in contingency table
	__contingencyTable2.push_back(g1v);
	__contingencyTable2.push_back(g2v);
}

float64 CFastANOVA::__computeSSB() {
	float64 SSB = 0.0;
	for(uint i=0; i<2;i++) {
		if(__contingencyTable2[i].rows()!=0) {
			VectorXd tmp = __contingencyTable2[i];
			float64 tmp_sum = tmp.sum();
			SSB += (tmp_sum*tmp_sum)/tmp.rows();
		}
	}
	return SSB - __tm;
}

void CFastANOVA::__permutePhenotypes() {
	__Y = MatrixXd::Zero(__n_samples,__permutations);
	Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(__n_samples);
	perm.setIdentity();
	for(uint i=0; i<__permutations;i++) {
		std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
		__Y.block(0,i,__n_samples,1) = perm*__y;
	}
	
	/*
	ifstream ifs;
	__Y = MatrixXd::Zero(__n_samples,__permutations);
	logging(INFO,"Reading Permutations Phenotype...");
	ifs.open("perm.txt",ifstream::in);
	if(!ifs.is_open()) {
		logging(ERROR,"Opening file perm.txt");
	}
	uint64 p=0;
	std::string line;
	while(!ifs.eof()) {
		if(p==100) break;
		getline(ifs,line);
		vector<string> sp = StringHelper::split(line," ");
		for(uint i=0; i<sp.size();i++) {
			__Y(i,p) = StringHelper::string_to<float64>(sp[i]);
		}
		p++;
	}
	ifs.close();
	*/
}

//Compute indices n_a1 and n_b1
void CFastANOVA::__compute_na1_nb1(uint64 const& i, uint64 const& j, float64* n_a1, float64* n_b1) {
	VectorXd a1 = (__X.col(i).array().colwise()*__X.col(j).array());
       	(*n_a1) = a1.sum();	
	(*n_b1) = (__X.col(j)-a1).sum();
}

//init and update index Array
void CFastANOVA::__updateIndexArray(float64 const& n_a1, float64 const& n_b1, uint64 const& k) {
	std::map<float64,std::map<float64, std::list<uint64> > >::iterator it1;
	it1 = __indexArray.find(n_a1);
	if(it1 == __indexArray.end()) {
		std::map<float64, std::list<uint64> > tmp_map;
		std::list<uint64> tmp_list;
		tmp_list.push_back(k);
		tmp_map.insert(std::pair<float64, std::list<uint64> >(n_b1,tmp_list));
		__indexArray.insert(std::pair<float64, std::map<float64, std::list<uint64 > > >(n_a1,tmp_map));
	} else {
		std::map<float64, std::list<uint64> >::iterator it2;
		std::map<float64, std::list<uint64> > *tmp = &it1->second;
		it2 = tmp->find(n_b1);
		if(it2 != tmp->end()) {
			it2->second.push_back(k);
		} else {
			std::list<uint64> tmp_list;
			tmp_list.push_back(k);
			tmp->insert(std::pair<float64,std::list<uint64> >(n_b1,tmp_list));
		}
	}
}

float64 CFastANOVA::__computeUpperBound(VectorXd const& plist, float64 const& r) {
	float64 nA = plist.rows();
	if(nA == 0) return 0;
	if(((nA-r)==0) || (r==0) ) return 0;
	float64 tsum = plist.sum();
	//compute l for the first r elements
	float64 l = plist.segment(0,r).sum();
	//compute u for the last r elements
	float64 u = plist.segment(nA-r,r).sum();
	return (std::max(pow(nA*l-r*tsum,2),pow(nA*u-r*tsum,2)))/(r*nA*(nA-r));
}

void CFastANOVA::__computeCandidateList(float64 const& SSB_XiY, std::list<uint>* candidateList) {
	float64 r1,r2;

	float64 t4 = __SST / ((__n_samples-4)/(3*__theta) + 1.0);
	float64 t3 = __SST / ((__n_samples-3)/(2*__theta) + 1.0);

	std::map<float64,std::map<float64, std::list<uint64> > >::iterator it1;
	std::map<float64, std::list<uint64> >::iterator it2;
	std::list<uint64>::iterator it3;

	for(it1=__indexArray.begin();it1!=__indexArray.end();it1++) {
		for(it2=it1->second.begin();it2!=it1->second.end();it2++) {
			//compute upper bounds R1 and r2 
			r1 = __computeUpperBound(__contingencyTable2[0],it1->first);
			r2 = __computeUpperBound(__contingencyTable2[1],it2->first);
			if(r1==0 || r2==0) {
				//store candidate SNPs in list
				if(SSB_XiY + r1 + r2 >= t3) {
					for(it3=it2->second.begin();it3!=it2->second.end();it3++) {
						candidateList->push_back((*it3));
					}
				}
			} else {
				//store candidate SNPs in list
				if(SSB_XiY + r1 + r2 >= t4) {
					for(it3=it2->second.begin();it3!=it2->second.end();it3++) {
						candidateList->push_back((*it3));
					}
				}
			}
		}
	}
}

void CFastANOVA::__updateTopList(std::list<episnp>& topList,episnp const& esnp) {
	//Update top k list of significant values
	if(topList.back().statistic>=esnp.statistic && topList.size()>=__alphaK) return;
	list<episnp>::iterator eit;
	list<episnp>::iterator itp=topList.end();
	list<episnp>::iterator iti=topList.end();
	for(eit=topList.begin();eit!=topList.end();eit++) {
		if(eit->phenotype_id == esnp.phenotype_id) itp = eit;
		if(iti==topList.end() && eit->statistic<=esnp.statistic) iti = eit;
	}
	if(topList.size()<__alphaK) {
		if(itp!=topList.end()) {
			if(itp->statistic<esnp.statistic) {
				topList.insert(iti,esnp);
				topList.erase(itp);
			}
		} else topList.insert(iti,esnp);
	} else {
		if(itp!=topList.end()) {
			if(itp->statistic<esnp.statistic) {
				topList.insert(iti,esnp);
				topList.erase(itp);
			}
		} else {
			topList.insert(iti,esnp);
			topList.pop_back();
		}
	}
}

void CFastANOVA::__computeCriticalFalpha() {
	//Init top F alpha list with zero values 
	std::list<episnp> topList;
	__theta = 0.0;
	__visited = 0;
	float64 n_a1, n_b1;

	for(uint64 i=0; i<__n_snps; i++) {
		//clear index array 
		__indexArray.clear();
		//compute Index array for all SNPs
		//clock_t time1 = std::clock();
		for(uint64 j=i+1;j<__n_snps;j++) {
			//compute n_a1 and n_b1
			__compute_na1_nb1(i,j,(&n_a1),(&n_b1));
			//update index array
			__updateIndexArray(n_a1,n_b1,j);
		}
		//clock_t time2 = std::clock() - time1;
		//std::cout << "Runtime: " << (float64)time2/CLOCKS_PER_SEC << " s" << std::endl;

		//for all phenotype permutations compute candidate SNPs
		for(uint64 p=0;p<__permutations;p++) {
			__computeContingencyTable2(__Y.col(p),i);
			float64 SSB_XiY = __computeSSB();
			//compute candidate list
			std::list<uint> candidateList;
			__computeCandidateList(SSB_XiY,&candidateList);
			//Count how many SNP have been visited
			__visited += candidateList.size();
			
			//Compute F-Statistic for all candidate SNPs 
			std::list<uint>::iterator it1;
			for(it1=candidateList.begin();it1!=candidateList.end();it1++) {
				episnp esnp;
				esnp.statistic = __computeFStatistic(__Y.col(p),__X.col(i),__X.col(*it1));
				esnp.phenotype_id = p;
				//update theta for candidate screening
				__updateTopList(topList,esnp);
				__theta = (topList.back()).statistic;
				__theta_pval = CFisherF::sf((topList.back()).statistic,__df1,__df2);
			}
		}
	}
}

void CFastANOVA::__computeSignificantSNPs() {
	//Init top F alpha list with zero values 
	std::vector<episnp> snplist;
	float64 n_a1, n_b1, SSB_XiY;
	for(uint64 i=0; i<__n_snps; i++) {
		//clear index array 
		__indexArray.clear();
		//compute Index array for all SNPs
		for(uint64 j=i+1;j<__n_snps;j++) {
			//compute n_a1 and n_b1
			__compute_na1_nb1(i,j,(&n_a1),(&n_b1));
			//update index array
			__updateIndexArray(n_a1,n_b1,j);
		}
		__computeContingencyTable2(__y,i);
		SSB_XiY = __computeSSB();
		//compute candidate list
		std::list<uint> candidateList;
		__computeCandidateList(SSB_XiY,&candidateList);
		//Compute F-Statistic for all candidate SNPs 
		std::list<uint>::iterator it1;
		for(it1=candidateList.begin();it1!=candidateList.end();it1++) {	
			float64 fstat = __computeFStatistic(__y,__X.col(i),__X.col(*it1));
			if(fstat>=__theta) {
				episnp esnp;
				esnp.i = i;
				esnp.j = *it1;
				esnp.statistic = fstat;
				esnp.pvalue = CFisherF::sf(fstat,__df1,__df2);
				snplist.push_back(esnp);
			}
		}
	}
	std::sort(snplist.begin(),snplist.end());
	__pvalues = VectorXd::Zero(snplist.size());
	__test_statistics = VectorXd::Zero(snplist.size());
	__snp_pairs = MatrixXd::Zero(snplist.size(),2);
	for(uint i=0; i<snplist.size();i++) {
		__pvalues(i) = snplist[i].pvalue;
		__test_statistics(i) = snplist[i].statistic;
		__snp_pairs(i,0) = snplist[i].i; 
		__snp_pairs(i,1) = snplist[i].j; 
		//logging(INFO, "SNP i " + StringHelper::to_string<uint64>(snplist[i].i) + " SNP j " + StringHelper::to_string<uint64>(snplist[i].j) + " : " + StringHelper::to_string<float64>(snplist[i].statistic) + "\t" + StringHelper::to_string<float64>(snplist[i].pvalue));
	}
}

void CFastANOVA::test_associations() {
	__n_samples = __y.rows();
	__n_snps = __X.cols();
	//compute __tm and global __SST
	__tm  = pow(__y.sum(),2)/__n_samples;
	__SST = __y.array().pow(2).sum() - __tm;
	//Permute Phenotype
	__permutePhenotypes();
	//compute critical value Falpha (upper bound)
	__computeCriticalFalpha();
	//compute significant SNPs for a given threshold theta
	__computeSignificantSNPs();
}

void CFastANOVA::test_associations(float64 const& theta) {
	__n_samples = __y.rows();
	__n_snps = __X.cols();
	//compute __tm and global __SST
	__tm  = pow(__y.sum(),2)/__n_samples;
	__SST = __y.array().pow(2).sum() - __tm;
	//__theta = CFisherF::sf(1-fwer_pval,__n_samples-4,3);
	__theta = theta;
	//compute significant SNPs for a given threshold theta
	__computeSignificantSNPs();
}	

};
