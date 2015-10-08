from gwas import start_gwas
from experiment import GWASettings
import sys,h5py
import multiprocessing

if __name__ in "__main__":

    if len(sys.argv)<3:
        print "\nUsage: InputFile OutputDir <ALGORITHM: linear (default), EMMAX, FaSTLMM, logit, EMMAXperm, logitperm, linearperm, fisher> <SNPENCODING: additive (default), recessive, dominant, codominant>"
        quit()

    hdf5_file = sys.argv[1]
    output_dir = sys.argv[2]
    try:
        algorithm = sys.argv[3]
    except:
        algorithm = "linear"
    try:
        encoding = sys.argv[4]
    except:
        encoding = "additive"
    
    phenotype_ids = []
    #phenotype_ids.append(int(sys.argv[5]))
    
    try:
        phenotype_file = sys.argv[5]
    except:
        phenotype_file = None

    if phenotype_file==None:
        f = h5py.File(hdf5_file,'r')
        phenotype_ids = f['Phenotypes'].keys()
        f.close()
    else:
        f = open(sys.argv[6],'r')
        phenotype_ids = []
        for line in f:
            phenotype_ids.append(line.strip())
        f.close()
    

    jobs = []
    for pid in phenotype_ids:
        settings = GWASettings()
        settings.algorithm = algorithm
        settings.output_path = output_dir
        settings.output_file = str(pid) + ".hdf5"
        settings.hdf5_file = hdf5_file
        #settings.phenotype_file = phenotype_file
        settings.phenotype_file = hdf5_file
        settings.snp_encoding = encoding
        settings.phenotype_id = pid
        settings.phenotype_transformation = "none"
        settings.homozygous = False
        settings.maf = 0.0
        settings.unique_snps_only = False
        #start_gwas(settings)
        #quit()
        p = multiprocessing.Process(target=start_gwas, args=(settings,))
        jobs.append(p)
        p.start()
