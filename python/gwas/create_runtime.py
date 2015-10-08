import h5py,sys,os
sys.path.append("bin/" + sys.platform + "/interfaces/python/")
import scipy as sp

import CEasyGWAS as gwas_core

def encodeHomozygousData(raw_data=None):
    gwas_data = gwas_core.CGWASDataHelper()
    gwas_data.encodeHomozygousData(raw_data,raw_data.shape[1],raw_data.shape[0])
    encoded = gwas_data.getEncodedData()
    maf_data = gwas_data.getMAF()
    gwas_data.releaseMemory()
    return [encoded,maf_data]

if __name__ in "__main__":
    if len(sys.argv)<4:
        print "USAGE: HDF5_input_file output_dir nr_samples nr_snps" 
        print

    f = h5py.File(sys.argv[1],'r')
    output_dir = sys.argv[2]

    nr_samples = float(sys.argv[3])
    nr_snps = float(sys.argv[4])

    raw = f['Genotype/raw'][:]
    sample_ids = f['Genotype/sample_ids'][:]
    chr_index = f['Genotype/chr_index'][:]
    pos_index = f['Genotype/position_index'][:]
    f.close()
    
    #PRINT
    print "Samples Total: ", sample_ids.shape[0]
    print "SNPS Total: ", raw.shape[1]
    
    #restrict data
    if nr_snps>0:
        raw = raw[:,0:nr_snps]
        chr_index = chr_index[0:nr_snps]
        pos_index = pos_index[0:nr_snps]
    if nr_samples>0:
        raw = raw[0:nr_samples,:]
        sample_ids = sample_ids[0:nr_samples]
    
    #write PED file
    f = open(os.path.join(output_dir,"samples_" + str(nr_samples) + "_snps_" + str(nr_snps) + ".ped"),'w')
    for i in xrange(sample_ids.shape[0]):
        string = sample_ids[i] + " " + sample_ids[i] + " -9 -9 -9 -9 " # + str(int(y[i])) + " "
        for j in xrange(raw.shape[1]):
            string += raw[i,j] + " " + raw[i,j] + " "
        f.write(string[:-1] + "\n")
    f.close()

    #write MAP file
    f = open(os.path.join(output_dir,"samples_" + str(nr_samples) + "_snps_" + str(nr_snps) + ".map"),'w')
    for i in xrange(chr_index.shape[0]):
        chr_id = chr_index[i].replace("Chr","").replace("chr","")
        f.write(chr_id + " " + chr_id + "_" + str(int(pos_index[i])) + " 0 " + str(int(pos_index[i])) + "\n")
    f.close()

    #write phenotype
    f = open(os.path.join(output_dir,"samples_" + str(nr_samples) + "_snps_" + str(nr_snps) + ".pheno"),'w')
    f.write("FID IID random_phenotype\n")
    for i in xrange(sample_ids.shape[0]):
        string = sample_ids[i] + " " + sample_ids[i] + " " + str(sp.absolute(sp.random.randn()))#str(int(sp.random.binomial(1,0.5)))
        f.write(string + "\n")
    f.close()
