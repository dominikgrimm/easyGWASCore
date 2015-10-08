import h5py,sys
import scipy as sp

def itemfreq(a):
    items,ind, inv = sp.unique(a, return_inverse=True,return_index=True)
    freq = sp.bincount(inv)
    return sp.array([ind, items, freq]).T

if __name__ in "__main__":
    f = h5py.File(sys.argv[1])

    raw = f['Genotype']['raw'][:]
    sample_ids = f['Genotype']['sample_ids'][:]
    chr_index = f['Genotype']['chr_index'][:]
    position_index = f['Genotype']['position_index'][:]

    rawT = raw.T
    snp_strings = sp.ascontiguousarray(rawT).view(sp.dtype((sp.void,rawT.dtype.itemsize * rawT.shape[1])))
    frequencies = itemfreq(snp_strings)

    ind = sp.where(frequencies[:,2]<=int(sys.argv[3]))[0]
    indices = sp.array(frequencies[ind,0],dtype="int")
    
    print "Number of Samples:\t\t\t", raw.shape[0]
    print "Number of SNPs before Filtering:\t", raw.shape[1]
    print "Number of truly unique SNPs:\t\t", sp.where(frequencies[ind,2]==1)[0].shape[0]
    print sp.histogram(frequencies[:,2],bins=[1,2,5,10,50,100,500,1000,2000])

    out = h5py.File(sys.argv[2],'w')
    g = out.create_group("Genotype")
    g.create_dataset("raw",data=raw[:,indices],chunks=True)
    g.create_dataset("sample_ids",data=sample_ids,chunks=True)
    g.create_dataset("chr_index",data=chr_index[indices],chunks=True)
    g.create_dataset("position_index",data=position_index[indices],chunks=True)

    #copy phenotype data
    p = out.create_group("Phenotypes")
    for key in f['Phenotypes'].keys():
        pp = p.create_group(str(key))
        for element in f['Phenotypes'][key].keys():
            try:
                pp.create_dataset(element,data=f['Phenotypes'][key][element][:],chunks=True)
            except:
                try:
                    pp.create_dataset(element,data=f['Phenotypes'][key][element].value)
                except:
                    pass
    out.close()

    f.close()
