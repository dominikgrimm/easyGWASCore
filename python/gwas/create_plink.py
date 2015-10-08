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
        print "USAGE: HDF5_input_file output_dir phenotype_id maf" 
        print

    f = h5py.File(sys.argv[1],'r')
    output_dir = sys.argv[2]
    phenotype_id = sys.argv[3]
    maf = float(sys.argv[4])

    y = f['Phenotypes'][phenotype_id]['y'][:]
    y = y.flatten()
    sample_ids_y = f['Phenotypes'][phenotype_id]['sample_ids'][:]
    phenotype_name = f['Phenotypes'][phenotype_id]['name'].value
    phenotype_name = phenotype_name.replace(" ","_").replace("<i>","").replace("</i>","")
    raw = f['Genotype/raw'][:]
    sample_ids = f['Genotype/sample_ids'][:]
    chr_index = f['Genotype/chr_index'][:]
    pos_index = f['Genotype/position_index'][:]
    sample_indices = (sp.reshape(sample_ids,(sample_ids.shape[0],1))==sample_ids_y).nonzero()
    y = y[sample_indices[1]]
    sample_ids = sample_ids[sample_indices[0]]
    raw = raw[sample_indices[0],:]
    [x,maf_data] = encodeHomozygousData(raw)
    if maf > 0.0:
        selected_indices = sp.where(maf_data>maf)[0]
        x = x[:,selected_indices]
        raw = raw[:,selected_indices]
        chr_index = chr_index[selected_indices]
        pos_index = pos_index[selected_indices]
        maf_data = maf_data[selected_indices]
    f.close()

    #write PED file
    f = open(os.path.join(output_dir,phenotype_name + ".ped"),'w')
    for i in xrange(sample_ids.shape[0]):
        string = sample_ids[i] + " " + sample_ids[i] + " -9 -9 -9 " + str(int(y[i])) + " "
        for j in xrange(raw.shape[1]):
            string += raw[i,j] + " " + raw[i,j] + " "
        f.write(string[:-1] + "\n")
    f.close()

    #write MAP file
    f = open(os.path.join(output_dir,phenotype_name + ".map"),'w')
    for i in xrange(chr_index.shape[0]):
        chr_id = chr_index[i].replace("Chr","").replace("chr","")
        f.write(chr_id + " " + chr_id + "_" + str(int(pos_index[i])) + " 0 " + str(int(pos_index[i])) + "\n")
    f.close()

    #write phenotype file
    f = open(os.path.join(output_dir,phenotype_name + ".pheno"),'w')
    f.write("FID IID " + phenotype_name + "\n")
    for i in xrange(sample_ids.shape[0]):
        string = sample_ids[i] + " " + sample_ids[i] + " " + str(int(y[i]))
        f.write(string + "\n")
    f.close()
