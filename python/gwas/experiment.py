import sys,os
sys.path.append("bin/" + sys.platform + "/interfaces/python/")

from transformations import zeroMean,zeroMeanUnitVarianz,logTransformation,sqrtTransformation,boxcoxTransformation,createDummyVariables
import fdr_helper as fdr
import fisher_exact as fe

#TODO IMPLEMENT IN C++
#from pymap.pca.kernelPCA import kernelPCA

#from dataio.tasks import *

#TODO IMPLEMENT RANKSUM IN C++
#import pymap.asso as asso
#import pymap.linear.ranksum as ranksum

import CEasyGWAS as gwas_core
import scipy as sp
import scipy.linalg as linalg
import h5py

class GWACovariates():
    def __init__(self):
        self.covariate_id = None
        self.covariate_transformation = "none"

class GWASettings():
    def __init__(self):
        self.algorithm = "linear" #others: EMMAX, FaSTLMM, EMMAXperm,logit,linearperm,logitperm
        self.output_path = ""
        self.output_file = "out.hdf5"
        self.hdf5_file = None
        self.snp_encoding = "additive" #other options: recessive, dominant, codominant
        self.phenotype_id = None
        self.phenotype_transformation = "boxcox"
        self.covariates = []
        self.homozygous = True
        self.maf = 0.1
        self.principle_components = 0
        self.unique_snps_only = True

class GWAExperiment():
    
    def __init__(self):
        self.__dbfile = None
        self.__y = None
        self.__sample_indices = None
        self.__x = None
        self.__x_additive = None
        self.__cov = None
        self.__ass = None
        self.__selected_indices = None
        self.__sample_ids = None
        self.__opts = None
        self.__number_chrs = None
        self.__algorithm = None
        self.__boxcox_lambda = None
        self.__phenotype_indices = None
        self.__limix_lmm = None
        self.__algo_model = None
        self.__permutation = False
        self.__perms = 1000
        self.__phenotype_name = None
        self.__transformed = False
        self.__raw = None
        #Init logger for periodic logger


    def openHDF5DB(self, dbfilename=None, mode='r'):
        self.__dbfile = h5py.File(dbfilename,mode)

    def closeHDF5DB(self):
        if not self.__dbfile == None:
            self.__dbfile.close()

    """
    Transform Nucleotides into binary matrix
    return: encoded matrix
    """
    def encodeHomozygousData(self,raw_data=None):
        gwas_data = gwas_core.CGWASDataHelper()
        gwas_data.encodeHomozygousData(raw_data,raw_data.shape[1],raw_data.shape[0])
        encoded = gwas_data.getEncodedData()
        maf_data = gwas_data.getMAF()
        gwas_data.releaseMemory()
        return [encoded,maf_data]
    
    """
    Transform Nucleotides into 0,1,2 matrix
    return: encoded matrix
    """
    def encodeHeterozygousData(self,raw_data=None,snp_encoding="additive"):
        gwas_data = gwas_core.CGWASDataHelper()
        if snp_encoding=="recessive":
            encoding = gwas_data.recessive
        elif snp_encoding=="dominant":
            encoding = gwas_data.dominant
        elif snp_encoding=="overdominant":
            encoding = gwas_data.codominant
        else:
            encoding = gwas_data.additive
        gwas_data.encodeHeterozygousData(raw_data,raw_data.shape[1],raw_data.shape[0],encoding)
        encoded = gwas_data.getEncodedData()
        maf_data = gwas_data.getMAF()
        gwas_data.releaseMemory()
        return [encoded,maf_data]

    #compute realized relationship matrix
    def computeRealizedRelationshipKernel(self,genotype=None):
        X = zeroMeanUnitVarianz(data=genotype)
        K = 1.0/X.shape[1]*sp.dot(X,X.T)
        #if K contains negative values set these values to 0
        #ind = sp.where(K<0)
        #K[ind] = 0
        return K
    
    def computePCA(self,X=None,number_pcs=5):
        K = self.computeRealizedRelationshipKernel(X)
        [evalue,PC] = linalg.eig(K)
        indices = sp.argsort(-evalue)
        PC = PC[:,indices]
        return PC[:,0:number_pcs]

    def getSampleIDs(self):
        return self.__sample_ids
    
    #load data for computation
    def loadData(self,settings=None):
        '''
        Load Phenotype Data
        '''
        try:
            if settings.phenotype_file==None:
                self.__y = self.__dbfile['Phenotypes/' + str(settings.phenotype_id) + '/y'][:]
                self.__sample_ids = self.__dbfile['Phenotypes/' + str(settings.phenotype_id) + '/sample_ids'][:]
                self.__phenotype_name = self.__dbfile['Phenotypes/' + str(settings.phenotype_id) + "/name"].value 
            else:
                f = h5py.File(settings.phenotype_file,'r')
                self.__y = f['Phenotypes/' + str(settings.phenotype_id) + '/y'][:]
                self.__sample_ids = f['Phenotypes/' + str(settings.phenotype_id) + '/sample_ids'][:]
                self.__phenotype_name = f['Phenotypes/' + str(settings.phenotype_id) + "/name"].value 
                f.close()
        except:
            print "[ERROR] Loading Phenotype went wrong"
            quit()
        #remove missing values
        ind = sp.where(~sp.isnan(self.__y))[0]
        self.__y = self.__y[ind]
        self.__sample_ids = self.__sample_ids[ind]
        #transform phenotypes
        self.__y = self.transformData(self.__y,settings.phenotype_transformation)
        
        '''
        Load Covariate Data and restrict samples
        '''
        self.__cov = None
        covariates = settings.covariates
        for covariate in covariates:
            try:
                cov = self.__dbfile['Covariates/' + str(covariate.covaraite_id) + '/y'][:]
                sample_ids = self.__dbfile['Covariates/' + str(covariate.covariate_id) + '/sample_ids'][:]
            except:
                print "[ERROR] Loading Covariate went wrong"
                quit()
            #transform covariates
            cov = self.transformData(cov,covariate.covariate_transformation)
            #match samples
            sample_indices = (sp.reshape(sample_ids,(sample_ids.shape[0],1))==self.__sample_ids).nonzero()
            sample_ids = sample_ids[sample_indices[0]]
            if self.__cov==None:
                self.__cov = cov[sample_indices[0]]
            else:
                self.__cov = sp.column_stack([self.__cov,cov])
            self.__y = self.__y[sample_indices[1]]
            self.__sample_ids = self.__sample_ids[sample_indices[1]]
        
        '''
        Load Genotype Data and restrict samples
        '''
        sample_ids_file = self.__dbfile['Genotype/sample_ids'][:]
        raw_data = self.__dbfile['Genotype/raw'][:]
        self.__chr_index = self.__dbfile['Genotype/chr_index'][:]
        self.__pos_index = self.__dbfile['Genotype/position_index'][:]
        sample_indices = (sp.reshape(self.__sample_ids,(self.__sample_ids.shape[0],1))==sample_ids_file).nonzero()
        self.__sample_ids = self.__sample_ids[sample_indices[0]]
        self.__y = self.__y[sample_indices[0]]
        if self.__cov != None:
            self.__cov = self.__cov[sample_indices[0]]
        raw_data = raw_data[sample_indices[1],:]
        self.__raw = raw_data
        if settings.homozygous==True:
            [self.__x, self.__maf_data] = self.encodeHomozygousData(raw_data)
        else:
            [self.__x, self.__maf_data] = self.encodeHeterozygousData(raw_data,settings.snp_encoding)
            if settings.snp_encoding!="additive":
                [self.__x_additive, self.__maf_data] = self.encodeHeterozygousData(raw_data)
                #This was experimental to use an additve Kinship in case a other encoding was selected, now we only use the MAF filtering based
                #on an additve model! Comment the next line out if you wish to use the additve kinship matrix again
                self.__x_additive = None
        if settings.maf > 0.0:
            self.filter_mAF(settings.maf)
        self.filterNonInformativeSNPs()
        
        if settings.principle_components > 0:
            if not self.__x_additive is None:
                cov = sp.real(self.computePCA(X=self.__x_additive,number_pcs=settings.principle_components))
            else:
                cov = sp.real(self.computePCA(X=self.__x,number_pcs=settings.principle_components))
            if self.__cov == None:
                self.__cov = cov
            else:
                self.__cov = sp.column_stack([self.__cov,cov])

        if not self.__x_additive is None:
            tmpX = self.__x_additive.T
        else:            
            tmpX = self.__x.T
        usnps,self.__snp_hash=sp.unique(sp.ascontiguousarray(tmpX).view(sp.dtype((sp.void,tmpX.dtype.itemsize*tmpX.shape[1]))),return_inverse=True)
        
        #compute kinship kernel if necessary
        if self.__algorithm == "FaSTLMM" or self.__algorithm=="EMMAX" or self.__algorithm=="EMMAXperm":
            if settings.unique_snps_only:
                tmp,uindex = sp.unique(self.__snp_hash,return_index=True)
                if not self.__x_additive is None:
                    self.__ass.setK(gwas_core.CKernels.realizedRelationshipKernel(self.__x_additive[:,uindex]))
                    #self.__ass.setK(self.computeRealizedRelationshipKernel(genotype=self.__x_additive[:,uindex]))
                else:
                    self.__ass.setK(gwas_core.CKernels.realizedRelationshipKernel(self.__x[:,uindex]))
                    #self.__ass.setK(self.computeRealizedRelationshipKernel(genotype=self.__x[:,uindex]))
            else:
                if not self.__x_additive is None:
                    self.__ass.setK(gwas_core.CKernels.realizedRelationshipKernel(self.__x_additive))
                    #self.__ass.setK(self.computeRealizedRelationshipKernel(genotype=self.__x_additive))
                else:
                    self.__ass.setK(gwas_core.CKernels.realizedRelationshipKernel(self.__x))
                    #self.__ass.setK(self.computeRealizedRelationshipKernel(genotype=self.__x))
        

    def transformData(self,data=None,transformation="none"):
        data = data.flatten()
        if transformation=="zeroMean":
            return zeroMean(data)
        elif transformation=="unitVariance":
            return zeroMeanUnitVarianz(data)
        elif transformation=="sqrt":
            return sqrtTransformation(data)
        elif transformation=="log10":
            return logTransformation(data)
        elif transformation=="boxcox":
            new_data = boxcoxTransformation(data)
            if new_data.sum()==data.sum():
                self.__transformed = False
            else:
                self.__transformed = True
            return new_data
        elif transformation=="dummyVariable":
            return createDummyVariables(data)
        else:
            return data
    
    def filterNonInformativeSNPs(self):
        tmp = sp.where((self.__x==2).sum(axis=0)!=self.__x.shape[0])[0]
        if not tmp.shape[0]==self.__x.shape[0]:
            self.__x = self.__x[:,tmp]
            self.__chr_index = self.__chr_index[tmp]
            self.__pos_index = self.__pos_index[tmp]
            self.__maf_data = self.__maf_data[tmp]
            self.__raw = self.__raw[:,tmp]
            if not self.__x_additive is None:
                self.__x_additive = self.__x_additive[:,tmp]
        tmp = sp.where((self.__x==1).sum(axis=0)!=self.__x.shape[0])[0]
        if not tmp.shape[0]==self.__x.shape[0]:
            self.__x = self.__x[:,tmp]
            self.__chr_index = self.__chr_index[tmp]
            self.__pos_index = self.__pos_index[tmp]
            self.__maf_data = self.__maf_data[tmp]
            self.__raw = self.__raw[:,tmp]
            if not self.__x_additive is None:
                self.__x_additive = self.__x_additive[:,tmp]
        tmp = sp.where((self.__x==0).sum(axis=0)!=self.__x.shape[0])[0]
        if not tmp.shape[0]==self.__x.shape[0]:
            self.__x = self.__x[:,tmp]
            self.__chr_index = self.__chr_index[tmp]
            self.__pos_index = self.__pos_index[tmp]
            self.__maf_data = self.__maf_data[tmp]
            self.__raw = self.__raw[:,tmp]
            if not self.__x_additive is None:
                self.__x_additive = self.__x_additive[:,tmp]

    #filter by minimum allele frequency
    def filter_mAF(self,maf_filter):
        selected_indices = sp.where(self.__maf_data>maf_filter)[0]
        self.__x = self.__x[:,selected_indices]
        self.__chr_index = self.__chr_index[selected_indices]
        self.__pos_index = self.__pos_index[selected_indices]
        self.__maf_data = self.__maf_data[selected_indices]
        self.__raw = self.__raw[:,selected_indices]
        if not self.__x_additive is None:
            self.__x_additive = self.__x_additive[:,selected_indices]


    def selectAlgorithm(self,algo_model=None):
        if algo_model=="linear":
            self.__ass = gwas_core.LinearRegression()
        elif algo_model=="logit":
            self.__ass = gwas_core.LogisticRegression()
        elif algo_model=="MWUrt":
            self.__ass = ranksum.RankSumAsso(test="MannWhitneyU")
        elif algo_model=="WCrt":
            self.__ass = ranksum.RankSumAsso(test="Wilcoxon")
        elif algo_model=="ttest":
            self.__ass = ranksum.RankSumAsso(test="ttest")
        elif algo_model=="FaSTLMM":
            self.__ass = gwas_core.FaSTLMM()
        elif algo_model=="EMMAX":
            self.__ass = gwas_core.EMMAX()
        elif algo_model=="EMMAXperm":
            self.__ass = gwas_core.EMMAX()
            self.__permutation = True
        elif algo_model=="linearperm":
            self.__ass = gwas_core.LinearRegression()
            self.__permutation = True
        elif algo_model=="logitperm":
            self.__ass = gwas_core.LogisticRegression()
            self.__permutation = True
        elif algo_model=="fisher":
            self.__ass = fe.FisherExact()
        self.__algorithm = algo_model
            
    def saveResults(self, filename=None):
        if self.__algorithm=="MWUrt" or self.__algorithm=="WCrt" or self.__algorithm=="fisher":
            pval = self.__ass.getPvalues()
        else:
            if self.__permutation==True:
                pval = self.__ass.getPermutationPValues()
            else:
                pval = self.__ass.getPValues()
        result_file = h5py.File(filename,'w')
        
        ind = sp.where(~sp.isnan(pval))[0]
        pval= pval[ind]
        self.__pos_index = self.__pos_index[ind]
        self.__chr_index = self.__chr_index[ind]
        
        if self.__algorithm=="MWUrt" or self.__algorithm=="WCrt" or self.__algorithm=="fisher":
            scores = self.__ass.getScores()[ind]
            betas = sp.array([])
            betas_se = sp.array([])
        else:
            scores = self.__ass.getTestStatistics()[ind]
            betas = self.__ass.getBetas()[ind]
            betas_se = self.__ass.getSEBetas()[ind]
        chromosome_names = sp.unique(self.__chr_index) 
        
        bf_thresholds = []
        try:
            bf_thresholds.append(0.1/float(pval.shape[0]))
            bf_thresholds.append(0.05/float(pval.shape[0]))
            bf_thresholds.append(0.01/float(pval.shape[0]))
        except:
            pass
        bf_thresholds = sp.array(bf_thresholds)
        
        #Benjamini and Hochberg FDR
        bh_thresholds = []
        [thr,pvals_adjusted,sort_idx] = fdr.benjamini_hochberg(p_values=pval, q_value=0.1, return_sort_idx=True)
        bh_thresholds.append(thr)
        bh_p_values = pvals_adjusted
        [thr,pvals_adjusted] = fdr.benjamini_hochberg(p_values=pval, q_value=0.05, sort_idx=sort_idx)
        bh_thresholds.append(thr)
        [thr,pvals_adjusted] = fdr.benjamini_hochberg(p_values=pval, q_value=0.01, sort_idx=sort_idx)
        bh_thresholds.append(thr)
        bh_thresholds = sp.array(bh_thresholds)

        #Benjamini, Hochberg and Yekutieli
        bhy_thresholds = []
        [thr,pvals_adjusted] = fdr.benjamini_hochberg_yekutieli(p_values=pval, q_value=0.1, sort_idx=sort_idx)
        bhy_thresholds.append(thr)
        bhy_p_values = pvals_adjusted
        [thr,pvals_adjusted] = fdr.benjamini_hochberg_yekutieli(p_values=pval, q_value=0.05, sort_idx=sort_idx)
        bhy_thresholds.append(thr)
        [thr,pvals_adjusted] = fdr.benjamini_hochberg_yekutieli(p_values=pval, q_value=0.01, sort_idx=sort_idx)
        bhy_thresholds.append(thr)
        bhy_thresholds = sp.array(bhy_thresholds)
        
        #Compute Q-values as in Storey and Robert Tibshirani
        q_values = fdr.storey_tibishirani(p_values=pval, sort_idx=sort_idx)
        
        #for chr in range(self.__number_chrs):
        #for chrom in chromosome_names:
        #    ind_chr = sp.where(self.__chr_index==chrom)[0]
        #    try:
        #        hd5_chr = result_file.create_group(chrom)
        #    except:
        #        hd5_chr = result_file['hd5_chr']
        result_file.create_dataset("chromosomes",data=self.__chr_index,compression="gzip",compression_opts=9,chunks=True)
        result_file.create_dataset("maf",data=self.__maf_data,compression="gzip",compression_opts=9,chunks=True)
        result_file.create_dataset("snp_hash",data=self.__snp_hash,compression="gzip",compression_opts=9,chunks=True)
        result_file.create_dataset("positions",data=self.__pos_index,compression="gzip",compression_opts=9,chunks=True)
        result_file.create_dataset("p_values",data=pval,compression="gzip",compression_opts=9,chunks=True)
        result_file.create_dataset("q_values",data=q_values,compression="gzip",compression_opts=9,chunks=True)
        result_file.create_dataset("bh_p_values",data=bh_p_values,compression="gzip",compression_opts=9,chunks=True)
        result_file.create_dataset("bhy_p_values",data=bhy_p_values,compression="gzip",compression_opts=9,chunks=True)
        if self.__algorithm=="MWUrt" or self.__algorithm=="WCrt" or self.__algorithm=="fisher":
            result_file.create_dataset('scores',data=scores,chunks=True,compression="gzip", compression_opts=9)
        else:
            result_file.create_dataset('scores',data=scores,chunks=True,compression="gzip", compression_opts=9)
            result_file.create_dataset('betas',data=betas,chunks=True,compression="gzip", compression_opts=9)
            result_file.create_dataset('betas_se',data=betas_se,chunks=True,compression="gzip", compression_opts=9)

        result_file.create_dataset('sample_ids',data=self.__sample_ids,compression="gzip",compression_opts=9,chunks=True)
        result_file.create_dataset('n_samples',data=self.__sample_ids.shape[0])
        result_file.create_dataset('snps_left',data=pval.shape[0])
        result_file.create_dataset('phenotype_name',data=self.__phenotype_name)
        if not (self.__algorithm=="MWUrt" or self.__algorithm=="WCrt" or self.__algorithm=="fisher"):
            result_file.create_dataset('AIC',data=self.__ass.getAIC())
            result_file.create_dataset('AICc',data=self.__ass.getAICc())
            result_file.create_dataset('BIC',data=self.__ass.getBIC())
            result_file.create_dataset('contains_betas',data="True")
            result_file.create_dataset('contains_betas_se',data="True")
        result_file.create_dataset('contains_scores',data="True")
        result_file.create_dataset('contains_maf',data="True")
        result_file.create_dataset('contains_snp_hash',data="True")
        result_file.create_dataset('transformed',data=self.__transformed)
        result_file.create_dataset('chromosome_names',data=chromosome_names,compression="gzip",compression_opts=9,chunks=True)
        result_file.create_dataset('bhy_thresholds',data=bhy_thresholds,chunks=True,compression="gzip", compression_opts=9)
        result_file.create_dataset('bf_thresholds',data=bf_thresholds,chunks=True,compression="gzip", compression_opts=9)
        result_file.create_dataset('bh_thresholds',data=bh_thresholds,chunks=True,compression="gzip", compression_opts=9)
        
        if self.__algorithm=="FaSTLMM" or self.__algorithm=="EMMAX":
            variance_null_model = self.__ass.computeVarianceExplainedNullModel(1000)
            result_file.create_dataset('variance_null_model',data=variance_null_model.mean())
            result_file.create_dataset('variance_null_model_sem',data=variance_null_model.std()/sp.sqrt(variance_null_model.shape[0]))
            result_file.create_dataset('estimated_heritability',data=self.__ass.getHeritabilityEstimate())
            result_file.create_dataset('estimated_genetic_variance',data=self.__ass.getGeneticVariance())
            result_file.create_dataset('estimated_noise_variance',data=self.__ass.getNoiseVariance())
            #narrow sense heritability TODO
        if self.__permutation==True:
            result_file.create_dataset('n_permutations_per_snp',data=self.__perms)        
        result_file.close()

    def saveFiles(self,path=None):
        f = open(os.path.join(path,"genotype_" + self.__phenotype_name.replace(" ","_").replace("<i>","").replace("</i>","") + ".txt"),'w')
        for i in xrange(self.__x.shape[1]):
            string = ""
            for j in xrange(self.__x.shape[0]):
                if self.__x[j,i]==2:
                    string += "1 "
                else:
                    string += "0 "
                #string += str(int(self.__x[j,i])) + " "
            f.write(string[:-1] + "\n")
        f.close()
        
        f = open(os.path.join(path,"plink_" + self.__phenotype_name.replace(" ","_").replace("<i>","").replace("</i>","") + ".ped"),'w')
        for i in xrange(self.__raw.shape[0]):
            string = self.__sample_ids[i] + " "
            for j in xrange(self.__raw.shape[1]):
                string += self.__raw[i,j] + " " + self.__raw[i,j] + " "
                #string += str(int(self.__x[j,i])) + " "
            f.write(string[:-1] + "\n")
        f.close()

        f = open(os.path.join(path,"phenotype_" + self.__phenotype_name.replace(" ","_").replace("<i>","").replace("</i>","") + ".txt"),'w')
        y = self.__y.flatten()
        for i in xrange(y.shape[0]):
            f.write(str(int(y[i])) + "\n")
        f.close()

        f = open(os.path.join(path,"plink_" + self.__phenotype_name.replace(" ","_").replace("<i>","").replace("</i>","") + ".map"),'w')
        for i in xrange(self.__chr_index.shape[0]):
            f.write(str(self.__chr_index[i]) + " " + str(self.__chr_index[i]) + "_" + str(self.__pos_index[i]) + " 0 " + str(self.__pos_index[i]) + "\n")
        #for i in xrange(self.__chr_index.shape[0]):
        #    f.write(str(self.__chr_index[i]) + "\t" + str(self.__pos_index[i]) + "\n")
        f.close()

    def train(self):
        if self.__algo_model=="MWUrt" or self.__algo_model=="WCrt":
            data = asso.MatrixData(x=self.__x,y=self.__y,covariates=self.__cov)
            self.__ass.setData(data)
            self.__ass.train()
        else:
            self.__ass.setPhenotype(self.__y)
            self.__ass.setGenotype(self.__x)

            if not self.__cov is None:
                self.__ass.setCovariates(sp.column_stack([sp.ones(self.__y.shape),self.__cov]))
            if self.__permutation==True:
                if self.__x.shape[1]<1000:
                    self.__perms = 1000000
                elif self.__x.shape[1]<10000:
                    self.__perms = 100000
                elif self.__x.shape[1]<100000:
                    self.__perms = 10000
                elif self.__x.shape[1]<1000000:
                    self.__perms = 1000
                elif self.__x.shape[1]<10000000:
                    self.__perms = 100
                else:
                    self.__perms = 10
                print "SNPS: ", self.__x.shape[1]
                print "Perms: ", self.__perms
                self.__ass.permutations(self.__perms)
            else:
                self.__ass.test_associations()
