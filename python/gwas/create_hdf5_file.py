import h5py,sys,os
import scipy as sp

iupac_map = {"AA":"A","GG":"G","TT":"T","CC":"C","AG":"R","GA":"R","CT":"Y","TC":"Y","GC":"S","CG":"S",
             "AT":"W","TA":"W","GT":"K","TG":"K","AC":"M","CA":"M"}

if __name__ in "__main__":
    if len(sys.argv)<4:
        print
        print "USAGE: ./create_hdf5_files.py <genotype prefix> <phenotype folder> <output file>"
        print
        quit()
    
    f = open(sys.argv[1] + '.map','r')
    chromosomes = []
    positions = []
    for line in f:
        sv = line.strip().split("\t")
        chromosomes.append(sv[0].strip())
        positions.append(int(sv[-1].strip()))
    f.close()
    chromosomes = sp.array(chromosomes)
    positions = sp.array(positions)

    f = open(sys.argv[1] + '.ped','r')
    sample_ids = []
    matrix = []
    for line in f:
        sv = line.strip().split("\t")
        sample_ids.append(sv[1].strip())
        snps = []
        j = 6
        while j < len(sv)-1:
            snps.append(iupac_map[sv[j] + sv[j+1]])
            j = j+2
        matrix.append(snps)
    f.close()
    sample_ids = sp.array(sample_ids)
    matrix = sp.array(matrix)

    #restrict input
    if len(sys.argv)==5:
        f = open(sys.argv[4],'r')
        chrids = []
        for line in f:
            sv = line.strip().split("\t")
            chrids.append(sv[-1].split(",")[0] + "_" + sv[-1].split(",")[1])
        f.close()
        indices = []
        for i in xrange(chromosomes.shape[0]):
            if not chromosomes[i] + "_" + str(positions[i]) in chrids:
                indices.append(i)
        indices = sp.array(indices)
        chromosomes = chromosomes[indices]
        positions = positions[indices]
        matrix = matrix[:,indices]

    hd5 = h5py.File(sys.argv[3])

    #save genotype
    genotype = hd5.create_group("Genotype")
    genotype.create_dataset("chr_index",data=chromosomes,chunks=True)
    genotype.create_dataset("position_index",data=positions,chunks=True)
    genotype.create_dataset("raw",data=matrix,chunks=True)
    genotype.create_dataset("sample_ids",data=sample_ids,chunks=True)

    #read phenotypes per file in folder and store data
    counter = 0
    phenotypes = hd5.create_group("Phenotypes")
    for subdir, dirs, files in os.walk(sys.argv[2]):
        for file in files:
                #processing file
                print os.path.join(subdir, file)
                f = open(os.path.join(subdir,file),'r')
                sample_ids = []
                phenotype_names = []
                Y = []
                for i,line in enumerate(f):
                    sv = line.strip().split("\t")
                    if i==0:
                        for j in xrange(2,len(sv)):
                            phenotype_names.append(sv[j].strip())
                        continue
                    sample_ids.append(sv[1].strip())
                    yline = []
                    for j in xrange(2,len(sv)):
                        yline.append(float(sv[j].strip()))
                    Y.append(yline)
                f.close()
                sample_ids = sp.array(sample_ids)
                Y = sp.array(Y)
                
                for i,phenotype in enumerate(phenotype_names):
                    ph = phenotypes.create_group(str(counter))
                    ph.create_dataset("name",data=phenotype)
                    ph.create_dataset("sample_ids",data=sample_ids,chunks=True)
                    ph.create_dataset("y",data=Y[:,i])
                    counter += 1
    hd5.close()
