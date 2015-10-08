import h5py,sys,os
import scipy as sp
    
def create_csv(file_path=None,csv_filename=None):
    f = h5py.File(file_path,'r')
    csv_filename = os.path.join(csv_filename,f['phenotype_name'].value.replace(" ","_") + ".csv") 
    csv_header = None
    csv_matrix = None
    
    if "betas" in f.keys():
        csv_header = sp.array(["CHR","Positions","P-Value","TestStatistic","Q-Value","Benjamini-Hochberg-P-Value","Benjamini-Hochberg-Yekutieli-P-Value","Beta0","SEBeta0","Beta1","SEBeta1","MAF","SNP-Hash"])
    else:
        csv_header = sp.array(["CHR","Positions","P-Value","TestStatistic","Q-Value","Benjamini-Hochberg-P-Value","Benjamini-Hochberg-Yekutieli-P-Value","MAF","SNP-Hash"])
    tmp_matrix = None
    if "betas" in f.keys():
        tmp_matrix = sp.column_stack([sp.array(f["chromosomes"],dtype="S50"),
                                          f['positions'],
                                          f['p_values'],
                                          f['scores'],
                                          f['q_values'],
                                          f['bh_p_values'],
                                          f['bhy_p_values'],
                                          f['betas'][:,0],
                                          f['betas_se'][:,0],
                                          f['betas'][:,1],
                                          f['betas_se'][:,1],
                                          f['maf'],
                                          f['snp_hash']])
    else:
        tmp_matrix = sp.column_stack([sp.array(f["chromosomes"],dtype="S50"),
                              f['positions'],
                              f['p_values'],
                              f['scores'],
                              f['q_values'],
                              f['bh_p_values'],
                              f['bhy_p_values'],
                              f['maf'],
                              f['snp_hash']])
     
    csv_matrix = tmp_matrix
    mf = open(csv_filename,'w')
    string = ""
    for i in xrange(csv_header.shape[0]):
        string += str(csv_header[i]) + ","
    mf.write(string[:-1] + "\n")
    for i in xrange(csv_matrix.shape[0]):
        string = ""
        for j in xrange(csv_matrix.shape[1]):
            string += str(csv_matrix[i,j]) + ","
        mf.write(string[:-1] + "\n")
    mf.close()
    f.close()

if __name__ in "__main__":
    create_csv(sys.argv[1],sys.argv[2])
