import os,time
from experiment import GWAExperiment 

def start_gwas(settings=None):   
    if settings.principle_components>0:
        print "Started GWAS for Phenotype ID: " + str(settings.phenotype_id) + " (PCs: " + str(settings.principle_components) + ")"
    else:
        print "Started GWAS for Phenotype ID: ", settings.phenotype_id
    execution_time = time.time()
    #start experiment
    experiment = GWAExperiment()
    experiment.selectAlgorithm(algo_model=settings.algorithm)    
    experiment.openHDF5DB(dbfilename=settings.hdf5_file)
    experiment.loadData(settings=settings)
    experiment.train()
    
    if not os.path.exists(settings.output_path):
        os.makedirs(settings.output_path)
    status_file_path = os.path.join(settings.output_path,settings.output_file)
       
    experiment.saveResults(status_file_path)
    #experiment.saveFiles(settings.output_path)
        
    experiment.closeHDF5DB()
        
    print "Finished Experiment in: ", str(time.time()-execution_time)
