import h5py,os,sys
import scipy as sp
import scipy.stats as stats
import sqlite3
import estimate_ld as ld
sys.path.append("bin/" + sys.platform + "/interfaces/python/")
import CEasyGWAS as gwas_core

#IUPAC Nucleotide MAP
iupac_map = {"AA":"A","GG":"G","TT":"T","CC":"C","AG":"R","GA":"R","CT":"Y","TC":"Y","GC":"S","CG":"S",
             "AT":"W","TA":"W","GT":"K","TG":"K","AC":"M","CA":"M"}
iupac_map_reverse = {'A':'A','C':'C','G':'G','T':'T','R':['A','G'],'Y':['C','T'],'S':['G','C'],'W':['A','T'],'K':['G','T'],'M':['A','C']}

'''
Helper function to create a item set with pattern frequencies
'''
def itemfreq(a):
    items,ind, inv = sp.unique(a, return_inverse=True,return_index=True)
    freq = sp.bincount(inv)
    return sp.array([ind, freq]).T, inv


'''
encode data
'''
def encodeHeterozygousData(raw_data=None,snp_encoding="additive"):
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

'''
get encoded data
'''
def getEncodedData(filename,encoding="additive",phenotype_id=None,maf=0.0):
    f = h5py.File(filename,'r')
    print phenotype_id
    if not phenotype_id==None:
        sample_ids = f['Genotype/sample_ids'][:]
        p_sample_ids = f['Phenotypes'][phenotype_id]['sample_ids'][:]
        y = f['Phenotypes'][phenotype_id]['y'][:]
        ind = sp.where(~sp.isnan(y))[0]
        y = y[ind]
        p_sample_ids = p_sample_ids[ind]
        ind = (sp.reshape(sample_ids,(sample_ids.shape[0],1))==p_sample_ids).nonzero()
        raw = f['Genotype/raw'][:]
        raw = raw[ind[0],:]
        [encoded,maf_v] = encodeHeterozygousData(raw)
        ind = sp.where(maf_v>=maf)[0]
        encoded = encoded[:,ind]
        identifiers = f['Genotype/identifiers'][:]
        identifiers = identifiers[ind]
        maf_v = maf_v[ind]
        f.close()
        return [encoded,maf_v,identifiers]
    if encoding=="additive":
        if 'encoded_additive' in f['Genotype'].keys():
            encoded = f['Genotype/encoded_additive'][:]
            maf_v = f['Genotype/global_maf'][:]
        else:
            [encoded,maf_v] = encodeHeterozygousData(f['Genotype/raw'][:])
    identifiers = f['Genotype/identifiers'][:]
    f.close()
    return [encoded,maf_v,identifiers]

'''
Read HDF5 Output File
'''
def read_HDF5_file(fileName=None):
    f = h5py.File(fileName,'r')
    pv = f['p_values'][:]
    positions = f['positions'][:]
    chromosomes = f['chromosomes'][:]
    hashs = f['snp_hash'][:]
    maf = f['maf'][:]
    name = f['phenotype_name'].value.replace(" ","_").replace("<i>","").replace("</i>","")
    f.close()
    
    tmp,u_ind = sp.unique(hashs,return_index=True)
    unique_pv = pv[u_ind]
    return [pv,positions,chromosomes,hashs,unique_pv,name]

'''
Read CSV Result File
'''
def read_CSV_file(fileName=None):
    f = open(fileName,'r')
    pv = []
    positions = []
    chromosomes = []
    hashs = []
    name = fileName.replace(".csv","").split("/")[-1]
    for i,line in enumerate(f):
        if line[0]=="#" or line[0]==">":
            continue
        sv = line.strip().split(",")
        chromosomes.append(sv[0].replace("Chr",""))
        positions.append(int(float(sv[1])))
        pv.append(float(sv[2]))
        hashs.append(sv[-1])
    f.close()
    pv = sp.array(pv)
    chromosomes = sp.array(chromosomes)
    positions = sp.array(positions)
    hashs = sp.array(hashs)

    tmp,u_ind = sp.unique(hashs,return_index=True)
    unique_pv = pv[u_ind]
    return [pv,positions,chromosomes,hashs,unique_pv,name]

'''
Parse Covariate Data and add to HDF5 file
'''
def addCovariates2HDF5(arguments):
    if not os.path.isfile(arguments.hdata):
        print "Argument --hdata " + arguments.hdata + "does not exist or is not a file\n"
        quit()
    
    covariate_list = []
    if arguments.addcovariates!=None:
        if os.path.isdir(arguments.addcovariates):
            for fn in os.listdir(arguments.addcovariates):
                filename = os.path.join(arguments.addcovariates,fn)
                if os.path.isfile(filename):
                    covariate_list.append(filename)
                else:
                    print "Argument --addcovariates " + filename + " is not a file\n"
                    quit()

        else:
            if os.path.isfile(arguments.addcovariates):
                covariate_list.append(arguments.addcovariates)
            else:
                print "Argument --addcovariates " + arguments.addcovariates + " does not exist or is not a file\n"
                quit()

    #open HDF5 file
    hd5 = h5py.File(arguments.hdata)

    #read covariates per file in folder and store data
    
    if "Covariates" in hd5.keys():
        covariates = hd5["Covariates"]
    else:
        covariates = hd5.create_group("Covariates")
    #find maximum covariate id in data file
    counter = 0
    for i in covariates.keys():
        if int(i)>counter:
            counter=int(i)
    counter += 1
    for filename in covariate_list:
        f = open(filename,'r')
        sample_ids = []
        covariate_names = []
        Y = []
        delimiter = "\t"
        for i,line in enumerate(f):
            sv = line.strip().split(delimiter)
            if len(sv)==1:
                delimiter = " "
                sv = line.strip().split(delimiter)

            if i==0:
                for j in xrange(2,len(sv)):
                    covariate_names.append(sv[j].strip())
                continue
            sample_ids.append(sv[1].strip())
            yline = []
            for j in xrange(2,len(sv)):
                yline.append(float(sv[j].strip()))
            Y.append(yline)
        f.close()
        sample_ids = sp.array(sample_ids)
        Y = sp.array(Y)
        
        for i,covariate in enumerate(covariate_names):
            ph = covariates.create_group(str(counter))
            ph.create_dataset("name",data=covariate)
            ph.create_dataset("sample_ids",data=sample_ids,chunks=True,compression="gzip",compression_opts=9)
            ph.create_dataset("y",data=Y[:,i],compression="gzip",compression_opts=9,chunks=True)
            counter += 1
    hd5.close()
    

'''
Plink 2 HDF5 file
'''
def plink2HDF5(arguments):
    if not os.path.isfile(arguments.plink_data + ".map"):
        print "Argument --plink_data " + arguments.plink_data + ".map does not exist or is not a file\n"
        quit()
    if not os.path.isfile(arguments.plink_data + ".ped"):
        print "Argument --plink_data " + arguments.plink_data + ".ped does not exist or is not a file\n"
        quit()
    if arguments.maf>1 or arguments.maf<0:
        print "Argument --maf " + str(arguments.maf) + " must be between 0 and 1\n"
        quit()
    phenotype_list = []
    if arguments.plink_phenotype!=None:
        if os.path.isdir(arguments.plink_phenotype):
            for fn in os.listdir(arguments.plink_phenotype):
                filename = os.path.join(arguments.plink_phenotype,fn)
                if os.path.isfile(filename):
                    phenotype_list.append(filename)
                else:
                    print "Argument --plink_phenotype " + filename + " is not a file\n"
                    quit()

        else:
            if os.path.isfile(arguments.plink_phenotype):
                phenotype_list.append(arguments.plink_phenotype)
            else:
                print "Argument --plink_phenotype " + arguments.plink_phenotype + " does not exist or is not a file\n"
                quit()
    exclude_snps = []
    if not arguments.exclude_snps==None:
        if not os.path.isfile(arguments.exclude_snps):
            print "Argument --exclude_snps " + arguments.exclude_snps + " does not exist or is not a file\n"
            quit()
        f = open(arguments.exclude_snps,'r')
        for line in f:
            exclude_snps.append(line.strip())
        f.close()
    exclude_snps = sp.array(exclude_snps)

    f = open(arguments.plink_data + '.map','r')
    chromosomes = []
    positions = []
    identifiers = []
    ref_allele = []
    delimiter = "\t"
    for line in f:
        sv = line.strip().split(delimiter)
        if len(sv)==1:
            delimiter = " "
            sv = line.strip().split(delimiter)
            
        chromosomes.append(sv[0].strip())
        identifiers.append(sv[1].strip())
        positions.append(int(sv[3].strip()))
        if len(sv)==5:
            ref_allele.append(sv[4].strip())
    f.close()
    chromosomes = sp.array(chromosomes)
    positions = sp.array(positions)
    identifiers = sp.array(identifiers)
    ref_allele = sp.array(ref_allele)
    if ref_allele.shape[0]>0:
        if ref_allele.shape[0] != identifiers.shape[0]:
            print "[ERROR] in *.map file: Column 5 (Reference allele) does not exist for all identifiers\n"
            quit()

    f = open(arguments.plink_data + '.ped','r')
    sample_ids = []
    family_ids = []
    paternal_ids = []
    maternal_ids = []
    sex = []
    matrix = []
    delimiter = "\t"
    for line in f:
        sv = line.strip().split(delimiter)
        if len(sv)==1:
            delimiter = " "
            sv = line.strip().split(delimiter)

        family_ids.append(sv[0].strip())
        sample_ids.append(sv[1].strip())
        paternal_ids.append(sv[2].strip())
        maternal_ids.append(sv[3].strip())
        sex.append(sv[4].strip())
        snps = []
        #TODO ADD EXCEPTION IF NUMBER OF SNPs are wrong
        j = 6
        while j < len(sv)-1:
            snps.append(iupac_map[sv[j] + sv[j+1]])
            j = j+2
        matrix.append(snps)
    f.close()
    family_ids = sp.array(family_ids)
    sample_ids = sp.array(sample_ids)
    paternal_ids = sp.array(paternal_ids)
    maternal_ids = sp.array(maternal_ids)
    matrix = sp.array(matrix)

    #Exclude SNPs
    print "SNPs before excluding SNPs:\t\t\t" + str(matrix.shape[1])
    ex_indices = []
    for i,ident in enumerate(identifiers):
        if not ident in exclude_snps:
            ex_indices.append(i)
    ex_indices = sp.array(ex_indices)
    if ex_indices.shape[0]!=matrix.shape[1]:
        print "SNPs after excluding SNPs:\t\t\t" + str(ex_indices.shape[0])
    matrix = matrix[:,ex_indices]
    identifiers = identifiers[ex_indices]
    chromosomes = chromosomes[ex_indices]
    positions = positions[ex_indices]
    if ref_allele.shape[0]>0:
        ref_allele = ref_allele[ex_indices]

    #FILTER MAF
    encoded = sp.array([])
    if arguments.maf>0:
        [encoded,maf] = encodeHeterozygousData(matrix)
        ind = sp.where(maf>=arguments.maf)[0]
        print "SNPs before MAF filtering:\t\t\t" + str(encoded.shape[1])
        encoded = encoded[:,ind]
        matrix = matrix[:,ind]
        maf = maf[ind]
        identifiers = identifiers[ind]
        chromosomes = chromosomes[ind]
        positions = positions[ind]
        if ref_allele.shape[0]>0:
            ref_allele = ref_allele[ind]
        print "SNPs after MAF filtering:\t\t\t" + str(encoded.shape[1])

    #distinct filtering
    if arguments.distinct_filter>0:
        if encoded.shape != matrix.shape:
            [encoded,maf] = encodeHeterozygousData(matrix)
        rawT = encoded.T
        snp_strings = sp.ascontiguousarray(rawT).view(sp.dtype((sp.void,rawT.dtype.itemsize * rawT.shape[1])))
        [frequencies, inverse] = itemfreq(snp_strings)
        if arguments.distinct_filter==1:
            ind = sp.where(frequencies[:,1]==int(arguments.distinct_filter))[0]
            indices = sp.array(frequencies[ind,0],dtype="int")
        else:
            ind = sp.where(frequencies[:,1]==1)[0]
            indices = sp.array(frequencies[ind,0],dtype="int")
            for i in sp.arange(2,int(arguments.distinct_filter)+1):
                tmp_indices = sp.where(frequencies[:,1]==i)[0]
                for tmp in tmp_indices:
                    ind = sp.where(inverse==tmp)[0]
                    chrom = chromosomes[ind]
                    un = sp.unique(chrom)
                    if un.shape[0]==1:
                        indices = sp.concatenate([indices,ind])
        print "Number of SNPs before distinct filtering:\t", matrix.shape[1]
        print "Number of truly unique SNPs:\t\t\t", sp.where(frequencies[:,1]==1)[0].shape[0]
        matrix = matrix[:,indices]
        chromosomes = chromosomes[indices]
        positions = positions[indices]
        identifiers = identifiers[indices]
        if ref_allele.shape[0]>0:
            ref_allele = ref_allele[indices]
        print "Number of SNPs after distinct filtering:\t", indices.shape[0]
        
    #sort data
    ind = sp.argsort(chromosomes)
    matrix = matrix[:,ind]
    identifiers = identifiers[ind]
    chromosomes = chromosomes[ind]
    positions = positions[ind]
    if ref_allele.shape[0]>0:
        ref_allele = ref_allele[ind]
    
    chrom_list = sp.unique(chromosomes)
    for chrom in chrom_list:
        ind = sp.where(chromosomes==chrom)[0]
        pos_tmp = positions[ind]
        chrom_tmp = chromosomes[ind]
        matrix_tmp = matrix[:,ind]
        ident_tmp = identifiers[ind]
        if ref_allele.shape[0]>0:
            ref_tmp = ref_allele[ind]
        #sort by position
        ind2 = sp.argsort(pos_tmp)
        positions[ind] = pos_tmp[ind2]
        chromosomes[ind] = chrom_tmp[ind2]
        matrix[:,ind] = matrix_tmp[:,ind2]
        identifiers[ind] = ident_tmp[ind2]
        if ref_allele.shape[0]>0:
            ref_allele[ind] = ref_tmp[ind2]

    #STORE DATA
    hd5 = h5py.File(arguments.hout)

    #save genotype
    genotype = hd5.create_group("Genotype")
    genotype.create_dataset("raw",data=matrix,chunks=True,compression="gzip",compression_opts=9)
    genotype.create_dataset("chr_index",data=chromosomes,chunks=True,compression="gzip",compression_opts=9)
    genotype.create_dataset("position_index",data=positions,chunks=True,compression="gzip",compression_opts=9)
    genotype.create_dataset("identifiers",data=identifiers,chunks=True,compression="gzip",compression_opts=9)
    genotype.create_dataset("sample_ids",data=sample_ids,chunks=True,compression="gzip",compression_opts=9)
    if ref_allele.shape[0]>0:
        genotype.create_dataset("ref_allele",data=ref_allele,chunks=True,compression="gzip",compression_opts=9)

    #read phenotypes per file in folder and store data
    counter = 0
    phenotypes = hd5.create_group("Phenotypes")
    for filename in phenotype_list:
        f = open(filename,'r')
        sample_ids = []
        phenotype_names = []
        Y = []
        delimiter = "\t"
        for i,line in enumerate(f):
            sv = line.strip().split(delimiter)
            if len(sv)==1:
                delimiter = " "
                sv = line.strip().split(delimiter)

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
            ph.create_dataset("sample_ids",data=sample_ids,chunks=True,compression="gzip",compression_opts=9)
            ph.create_dataset("y",data=Y[:,i],compression="gzip",compression_opts=9,chunks=True)
            counter += 1
    hd5.close()

def encodeData(arguments):
    if not os.path.isfile(arguments.hdata):
        print "Argument --hdata " + arguments.hdata + " is not a file\n"
        quit()
    f = h5py.File(arguments.hdata)
    if arguments.encode=="additive":
        if "encoded_additive" in f['Genotype'].keys(): 
            print "HDF5 File already contains an additive encoded genotype matrix\n"
            quit()
    if arguments.encode=="dominant":
        if "encoded_dominant" in f['Genotype'].keys(): 
            print "HDF5 File already contains a dominant encoded genotype matrix\n"
            quit()
    if arguments.encode=="recessive":
        if "encoded_recessive" in f['Genotype'].keys(): 
            print "HDF5 File already contains a recessive encoded genotype matrix\n"
            quit()
    if arguments.encode=="overdominant":
        if "encoded_overdominant" in f['Genotype'].keys(): 
            print "HDF5 File already contains an overdominant encoded genotype matrix\n"
            quit()
    raw = f['Genotype']['raw'][:]
    [encoding,maf] = encodeHeterozygousData(raw,arguments.encode)
    if arguments.encode=="additive":
        g = f['Genotype']
        g.create_dataset("encoded_additive",data=encoding,chunks=True,compression="gzip",compression_opts=9)
        g.create_dataset("global_maf",data=maf,chunks=True,compression="gzip",compression_opts=9)
    if arguments.encode=="dominant":
        g = f['Genotype']
        g.create_dataset("encoded_dominant",data=encoding,chunks=True,compression="gzip",compression_opts=9)
    if arguments.encode=="recessive":
        g = f['Genotype']
        g.create_dataset("encoded_recessive",data=encoding,chunks=True,compression="gzip",compression_opts=9)
    if arguments.encode=="overdominant":
        g = f['Genotype']
        g.create_dataset("encoded_overdominant",data=encoding,chunks=True,compression="gzip",compression_opts=9)
    f.close()

'''
Create VCF file from HDF5
'''
def convertHDF5_2_VCF(arguments):
    if not os.path.isfile(arguments.hdata):
        print "Argument --hdata " + arguments.hdata + " is not a file\n"
        quit()
    if os.path.isfile(arguments.vout):
        print "File in --vout " + arguments.vout + " already exists. Please specify a different file!\n"
        quit()
    f = h5py.File(arguments.hdata,'r')
    raw = f['Genotype/raw'][:]
    positions = f['Genotype']['position_index'][:]
    chromosomes = f['Genotype']['chr_index'][:]
    identifiers = f['Genotype']['identifiers'][:]
    major_ref = False
    if "ref_allele" in f['Genotype'].keys():
        ref_allele = f['Genotype']['ref_allele'][:]
    else:
        major_ref = True
        print "No reference allele list in HDF5 file! Reference allele is set to the major allele!"
    f.close()

    out = open(arguments.vout,'w')
    snps_one_allele = 0
    snps_more_alleles = 0
    out.write("##Hyrbid VCF\n")
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for i in xrange(raw.shape[1]):
        stmp = sp.unique(raw[:,i])
        snps = []
        for s in stmp:
            nucs = iupac_map_reverse[s]
            for nuc in nucs:
                snps.append(nuc)
        snps = sp.array(snps)
        snps = sp.unique(snps)
        if snps.shape[0]<2:
            snps_one_allele += 1
            continue
        if snps.shape[0]>=3:
            snps_more_alleles += 1
            continue
        if major_ref==False:
            ind = sp.where(snps!=ref_allele[i])[0]
            out.write(str(chromosomes[i]) + "\t" + str(positions[i]) + "\t" + identifiers[i] + "\t" + ref_allele[i] + "\t" + snps[ind][0] + "\t.\t.\t.\n")
        else: #Use major allele
            ind1 = sp.where(snps[0]==raw[:,i])[0][0]
            ind2 = sp.where(snps[1]==raw[:,i])[0][0]
            minor_index = sp.argmin(sp.array([ind1,ind2]))
            major_index = sp.argmax(sp.array([ind1,ind2]))
            out.write(str(chromosomes[i]) + "\t" + str(positions[i]) + "\t" + identifiers[i] + "\t" + snps[major_index][0] + "\t" + snps[minor_index][0] + "\t.\t.\t.\n")
    out.close()

'''
GFF Element INFO
'''
def getGFFInfoElements(info=None):
    elements = {}
    for element in info:
        sv = element.strip().split("=")
        if len(sv)==2:
            elements[sv[0].strip().upper()] = sv[1].strip()
    return elements

'''
GFF 2 SQLLite DB
'''
def createGFFDB(arguments):
    if not os.path.isfile(arguments.gfile):
        print "Argument --gfile " + arguments.gfile + " is not a file\n"
        quit()
    if os.path.isfile(arguments.sqlout):
        print "File in --sqlout " + arguments.sqlout + " already exists. Please specify a different file!\n"
        quit()
    gfff = open(arguments.gfile,'r')
    split_delimiter = " "
    tmp_chromosome_map = {}
    
    gene_set = []
    for i,line in enumerate(gfff):
        sv = line.strip().split(split_delimiter)
        if i==0:
            if len(sv)==1:
                sv = line.strip().split("\t")
                if len(sv)>1:
                    split_delimiter = "\t"
                else:
                    print "[GFF File]: Wrong file format in GFF file. Delimiter has to be either a whitespace or a tab."
                    quit()
        if len(sv)<8 or len(sv)==0:
            print "[GFF File]: Wrong file format. GFF file must contain 9 columns. See FAQ for details."
            quit()
        if sv[2].strip().upper() == "CHROMOSOME" or sv[2].strip().upper() == "CHROMOSOME_ARM":
            info = sv[8].strip().split(";")
            elements = getGFFInfoElements(info)
            if len(elements)==0:
                print "[GFF File]: Wrong file format. The attribute column has to less elements (Column 9). At least the ID field has to be specified. See FAQ for details."
                quit()
            if not sv[0].strip()==elements["ID"]:
                print "[GFF File]: Wrong file format. Chromosome ID in column 1 has to be the same identifier as the ID in column 9. See FAQ for details."
                print "[GFF File]: Chromosome IDs do not match. Unknown chromosome identifier found in GFF file: " + elements["ID"]
                quit()
        elif sv[2].strip().upper()=="GENE":
            info = sv[8].strip().split(";")
            elements = getGFFInfoElements(info)
            if len(elements)==0:
                print "[GFF File]: Wrong file format. The attribute column has to less elements (Column 9). At least the ID field has to be specified. See FAQ for details."
                quit()
            if not (sv[6].strip()=="+" or sv[6].strip()=="-" or sv[6].strip()=="."):
                print "[GFF File]: Strand has to be either + or -. To indicate missing information use . (dot). See FAQ for details."
                quit()
            if "NAME" in elements:
                name = elements["NAME"]
            else:
                name = ""
            if "NOTE" in elements:
                note = elements["NOTE"]
            else:
                note = ""
            gene = (str(elements["ID"]),str(sv[3].strip()),str(sv[4].strip()),str(sv[1].strip()),
                    str(sv[6].strip()),name,note,str(sv[0].strip()))
            gene_set.append(gene)
           
    gfff.close()
    #CREATE SQLITE3 DATABASE
    sqlite = sqlite3.connect(arguments.sqlout)
    sqlite_cursor = sqlite.cursor()
    #CREATE TABLE
    sqlite_cursor.execute('''
        CREATE TABLE "geneannotation" (
        "id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        "annotation_id" varchar(30) NOT NULL,
        "annotation_source" varchar(255) NOT NULL,
        "annotation_start" integer NOT NULL,
        "annotation_end" integer NOT NULL,
        "annotation_strand" varchar(1),
        "annotation_note" varchar(255),
        "annotation_name" varchar(255),
        "chromosome_id" varchar(255) NOT NULL
        )
    ''')
    sqlite_cursor.executemany('''INSERT INTO geneannotation (annotation_id,
                                                        annotation_start,
                                                        annotation_end,
                                                        annotation_source,
                                                        annotation_strand,
                                                        annotation_name,
                                                        annotation_note,
                                                        chromosome_id)
                                VALUES (?,?,?,?,?,?,?,?)''', gene_set)
    sqlite.commit()
    sqlite.close()

'''
READ SIFT4G Output
'''
class SIFT4GEffect():
    def __init__(self,chrom=None,pos=None,ref_allele=None,alt_allele=None,transcript_id=None,
                 gene_id=None,gene_name=None,region=None,variant_type=None,ref_amino=None,
                 alt_amino=None,sift_score=None,sift_prediction=None,amino_pos=None):
        self.pos = pos
        self.chrom = chrom
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.region = region
        self.variant_type = variant_type
        self.ref_amino = ref_amino
        self.alt_amino = alt_amino
        self.amino_pos = amino_pos
        self.sift_score = sift_score
        self.sift_prediction = sift_prediction

def readSIFT4GFile(arguments,read_all=True):
    if not os.path.isfile(arguments.sfile):
        print "Argument --sfile " + arguments.sfile + " is not a file\n"
        quit()
    f = open(arguments.sfile,'r')
    sift_map = {}
    for i,line in enumerate(f):
        if i==0:
            continue
        sv = line.strip().split("\t")
        if sv[-1] == "NA":
            continue
        if "DELETERIOUS" in sv[-1]:
            state = "DELETERIOUS"
        else:
            state = "TOLERANT"
        if read_all==False:
            if sv[8]=="NONSYNONYMOUS":
                effect =  SIFT4GEffect(chrom=sv[0],pos=int(sv[1]),ref_allele=sv[2],alt_allele=sv[3],transcript_id=sv[4],
                                   gene_id=sv[5],gene_name=sv[6],region=sv[7],variant_type=sv[8],ref_amino=sv[9],
                                   alt_amino=sv[10],amino_pos=int(sv[11]),sift_score=float(sv[12]),sift_prediction=state)
                sift_map[sv[0] + "_" + sv[1]] = effect
        else:
            effect =  SIFT4GEffect(chrom=sv[0],pos=int(sv[1]),ref_allele=sv[2],alt_allele=sv[3],transcript_id=sv[4],
                                   gene_id=sv[5],gene_name=sv[6],region=sv[7],variant_type=sv[8],ref_amino=sv[9],
                                   alt_amino=sv[10],amino_pos=int(sv[11]),sift_score=float(sv[12]),sift_prediction=state)
            sift_map[sv[0] + "_" + sv[1]] = effect
    f.close()
    return sift_map

class PathogenicityEffect():
    def __init__(self,chrom=None,pos=None,ref_allele=None,alt_allele=None,transcript_id=None,
                 gene_id=None,gene_name=None,region=None,variant_type=None,ref_amino=None,
                 alt_amino=None,prediction_score=None,prediction=None,amino_pos=None):
        self.pos = pos
        self.chrom = chrom
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.region = region
        self.variant_type = variant_type
        self.ref_amino = ref_amino
        self.alt_amino = alt_amino
        self.amino_pos = amino_pos
        self.prediction_score = prediction_score
        self.prediction = prediction

def readPathogenicityFile(arguments,read_all=True):
    if not os.path.isfile(arguments.sfile):
        print "Argument --sfile " + arguments.sfile + " is not a file\n"
        quit()
    f = open(arguments.sfile,'r')
    path_map = {}
    for i,line in enumerate(f):
        if i==0:
            continue
        sv = line.strip().split("\t")
        if sv[4] == "NA":
            continue
        if "DELETERIOUS" in sv[4]:
            state = "DELETERIOUS"
        else:
            state = "TOLERANT"
        if read_all==False:
            if sv[5]=="NONSYNONYMOUS":
                effect =  PathogenicityEffect(chrom=sv[1],pos=int(sv[2]),ref_allele=sv[6],alt_allele=sv[7],transcript_id=sv[8],
                                   gene_id=sv[9],gene_name=sv[10],region=sv[11],variant_type=sv[5],ref_amino=sv[12],
                                   alt_amino=sv[13],amino_pos=int(sv[14]),prediction_score=float(sv[3]),prediction=state)
                path_map[sv[0]] = effect
        else:
            effect =  PathogenicityEffect(chrom=sv[1],pos=int(sv[2]),ref_allele=sv[6],alt_allele=sv[7],transcript_id=sv[8],
                               gene_id=sv[9],gene_name=sv[10],region=sv[11],variant_type=sv[5],ref_amino=sv[12],
                               alt_amino=sv[13],amino_pos=int(sv[14]),prediction_score=float(sv[3]),prediction=state)
            path_map[sv[0]] = effect
    f.close()
    return path_map


'''
Create CSV file from HDF5
'''
def writeCSVOutput(hdf5File=None,cout=None):
    f = h5py.File(hdf5File,'r')

    csv_filename = os.path.join(cout,f['phenotype_name'].value.replace(" ","_") + ".csv") 
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

def convertHDF5_2_CSV(arguments):
    if os.path.isdir(arguments.hfile):
        for fn in os.listdir(arguments.hfile):
            filename = os.path.join(arguments.hfile,fn)
            if not os.path.isfile(filename):
                continue
            else:
                print "Creating CSV output for:\t" + fn,
                writeCSVOutput(filename,arguments.cout)
                print "\t\t\t[Done]"
    else:
        if not os.path.isfile(arguments.hfile):
            print "Argument --hfile " + arguments.hfile + " is not a file\n"
            quit()
        writeCSVOutput(arguments.hfile,arguments.cout)

'''
Write TOP genes to output file
'''
def writeTopXGenes2File(filename,sqlfile,outdir,top=1000):
    f = h5py.File(filename,'r')
    chromosomes = f['chromosomes'][:]
    positions = f['positions'][:]
    p_values = f['p_values'][:].flatten()
    name = f['phenotype_name'].value.replace(" ","_").replace("<i>","").replace("</i>","")
    ind = sp.argsort(p_values)[:-1]
    chromosomes = chromosomes[ind]
    positions = positions[ind]
    p_values = p_values[ind]
    chromosomes = chromosomes[0:top]
    positions = positions[0:top]
    p_values = p_values[0:top]
    f.close()

    sqlite = sqlite3.connect(sqlfile)
    sqlite_cursor = sqlite.cursor()

    out = open(os.path.join(outdir,name + ".csv"),"w")

    out.write("Chr,Pos,PVal,GeneID (closest),Distance (bp)\n")
    for i in xrange(chromosomes.shape[0]):
        sqlite_cursor.execute("SELECT * FROM geneannotation WHERE chromosome_id=? ORDER BY ABS(annotation_start - ?) LIMIT 1",(str(chromosomes[i]),int(positions[i])))
        annotation = sqlite_cursor.fetchall()
        #print annotation
        if len(annotation)==1:
            if positions[i] >= annotation[0][3] and positions[i] <= annotation[0][4]:
                distance = 0
            elif positions[i] > annotation[0][4]:
                distance = abs(positions[i]-annotation[0][4])
            else:
                distance = abs(positions[i]-annotation[0][3])
            out.write(chromosomes[i] + "," + str(int(positions[i])) + ",%.2e"%(p_values[i]) + "," + annotation[0][1] + "," + str(int(distance)) + "\n")
    sqlite.close()

def writeTopXGenes(arguments=None):
    if not os.path.isfile(arguments.sqlfile):
        print "File in --sqlfile " + arguments.sqlfile + " is not a file!\n"
        quit()
    if not os.path.isdir(arguments.gout):
        print "File in --gout " + arguments.gout + " is not a directory!\n"
        quit()
    if os.path.isdir(arguments.hfile):
        for fn in os.listdir(arguments.hfile):
            filename = os.path.join(arguments.hfile,fn)
            if not os.path.isfile(filename):
                continue
            else:
                print "Creating GENE output for top " + str(arguments.topx) +  " Genes:\t" + fn,
                writeTopXGenes2File(filename,arguments.sqlfile,arguments.gout,arguments.topx)
                print "\t\t\t[Done]"
    else:
        if not os.path.isfile(arguments.hfile):
            print "Argument --hdata " + arguments.hfile + " is not a file\n"
            quit()
        writeTopXGenes2File(arguments.hfile,arguments.sqlfile,arguments.gout,arguments.topx)

'''
Generate LD file
'''

def writeLDInfo2File(filename=None,nr_hypothesis=-1,
                     distinct=False,selected_snp='-1',
                     distance=10000,
                     r2_measure="excoffier_slatkin",
                     outdir=None,ignore=None,hdata=None):

    [pv,positions,chromosomes,hashs,unique_pv,fname] = read_HDF5_file(filename)
    [encoded,maf,identifiers] = getEncodedData(hdata)
    
    if ignore!=None:
        if ignore in fname:
            return

    if nr_hypothesis==-1:
        if distinct:
            bf_threshold = 0.05/unique_pv.shape[0]
        else:
            bf_threshold = 0.05/pv.shape[0]
    else:
        bf_threshold = 0.05/arguments.nr_hypothesis
    
    snp_distance = distance
    r2_measure = r2_measure
    
    #select different SNP
    if selected_snp!='-1':
        ind = sp.where(identifiers==selected_snp)[0]
        if ind.shape[0]==0:
            print "\nSNP " + selected_snp + " not found in dataset!"
            print "Please select a different SNP identifier!\n"
            quit()
        else:
            chrom_list = sp.array([selected_snp.split("_")[0]])
            __pos = sp.array([int(selected_snp.split("_")[1])])
    else:
        chrom_list = sp.unique(chromosomes)

    out = open(os.path.join(outdir,fname.replace(" ","_") + ".txt"),"w")
    out.write("Chr_Focal\tPosition_Focal\tSNP_ID_Focal\tMAF_Focal\tPValue_Focal\tChr_B\tPosition_B\tSNP_ID_B\tMAF_B\tPValue_B\tR2\tPearsonR2\n")
    
    for i,chrom in enumerate(chrom_list):
        idx = chromosomes==chrom
        _pos = positions[idx]
        _chrs = chromosomes[idx]
        _pv = pv[idx]

        if selected_snp=='-1':
            idx = sp.where(_pv<=bf_threshold)[0]
            __pos = _pos[idx]
            __chrs = _chrs[idx]
            __pv = _pv[idx]
        else:
            idx = sp.where(__pos[0]==_pos)[0]
            __pv = _pv[idx]
        
        ind = sp.argsort(_pos)
        _pos = _pos[ind]
        _pv = _pv[ind]
        for k,pp in enumerate(__pos):
            if pp-snp_distance > 0 and pp+snp_distance<=_pos.max():
                ranges = sp.where((_pos>=pp-snp_distance) & (_pos<=pp+snp_distance))[0]
            elif pp-snp_distance < 0 and pp+snp_distance<=_pos.max():
                ranges = sp.where(_pos<=pp+snp_distance)[0]
            elif pp-snp_distance > 0 and pp+snp_distance>_pos.max():
                ranges = sp.where(_pos>=pp-snp_distance)[0]
            
            idx = sp.where(str(chrom) + "_" + str(pp)==identifiers)[0][0]
            
            for l,sra in enumerate(ranges):
                string = str(chrom) + "\t" + str(pp) + "\t" + identifiers[idx] + "\t%.4f"%(maf[idx]) + "\t%.4e"%(__pv[k]) + "\t" + str(chrom) + "\t"
                if _pos[sra]==pp:
                    continue
                sind = sp.where(str(chrom) + "_" + str(_pos[sra])==identifiers)[0][0]
                if r2_measure=="excoffier_slatkin":
                    rr = ld.esem_r(sp.array(encoded[:,sind].flatten(),dtype="int"),sp.array(encoded[:,idx].flatten(),dtype="int"))**2
                elif r2_measure=="roger_huff":
                    rr = ld.get_r(sp.array(encoded[:,sind].flatten(),dtype="int"),sp.array(encoded[:,idx].flatten(),dtype="int"))**2
                elif r2_measure=="pearson_r2":
                    rr = stats.pearsonr(encoded[:,sind].flatten(),encoded[:,idx].flatten())[0]**2
                rrp = stats.pearsonr(encoded[:,sind].flatten(),encoded[:,idx].flatten())[0]**2
                string += str(_pos[sra]) + "\t" + identifiers[sind] + "\t%.4f"%(maf[sind]) + "\t%.4e"%(_pv[sra]) + "\t%.4f"%(rr) + "\t%.4f"%(rrp) + "\n"
                out.write(string)
    out.close()

def writeLDInfo(arguments=None):
    if not os.path.isdir(arguments.ldout):
        print "File in --ldout " + arguments.ldout + " is not a directory!\n"
        quit()
    if not os.path.isfile(arguments.hdata):
        print "Argument --hdata " + arguments.hdata + " is not a file\n"
        quit()
    if os.path.isdir(arguments.hfile):
        for fn in os.listdir(arguments.hfile):
            filename = os.path.join(arguments.hfile,fn)
            if not os.path.isfile(filename):
                continue
            else:
                print "Creating LD output for:\t" + fn,
                writeLDInfo2File(filename=filename,nr_hypothesis=arguments.nr_hypothesis,
                                 distinct=arguments.distinct,selected_snp=arguments.selected_snp,
                                 distance=arguments.distance, r2_measure=arguments.r2_measure,
                                 outdir=arguments.ldout,ignore=arguments.ignore,hdata=arguments.hdata)
                print "\t\t\t[Done]"
    else:
        if not os.path.isfile(arguments.hfile):
            print "Argument --hfile " + arguments.hfile + " is not a file\n"
            quit()
        writeLDInfo2File(filename=arguments.hfile,nr_hypothesis=arguments.nr_hypothesis,
                         distinct=arguments.distinct,selected_snp=arguments.selected_snp,
                         distance=arguments.distance, r2_measure=arguments.r2_measure,
                         outdir=arguments.ldout,ignore=arguments.ignore,hdata=arguments.hdata)

'''
CONVERT HDF5 to PLINK
'''
def convertHDF5_2_PLINK(arguments=None):
    if not os.path.isfile(arguments.hdata):
        print "Argument --hdata " + arguments.hdata + " is not a file\n"
        quit()
    hf = h5py.File(arguments.hdata,'r')
    output_dir = arguments.pout
    
    raw = hf['Genotype/raw'][:]
    sample_ids = hf['Genotype/sample_ids'][:]
    chr_index = hf['Genotype/chr_index'][:]
    pos_index = hf['Genotype/position_index'][:]
    identifiers = hf['Genotype/identifiers'][:]

    #write PED file
    f = open(output_dir + ".ped",'w')
    for i in xrange(sample_ids.shape[0]):
        string = sample_ids[i] + " " + sample_ids[i] + " 0 0 0 0 "
        for j in xrange(raw.shape[1]):
            snps = iupac_map_reverse[raw[i,j]]
            if len(snps)==2:
                string += snps[0] + " " + snps[1] + " "
            else:
                string += snps + " " + snps + " "
        f.write(string[:-1] + "\n")
    f.close()

    #write MAP file
    f = open(output_dir + ".map",'w')
    for i in xrange(chr_index.shape[0]):
        chr_id = chr_index[i].replace("Chr","").replace("chr","")
        f.write(chr_id + " " + identifiers[i] + " 0 " + str(int(pos_index[i])) + "\n")
    f.close()
    
    #write phenotype file
    phenotype_ids = hf['Phenotypes'].keys()
    for pid in phenotype_ids:
        phenotype = hf['Phenotypes'][pid]
        y = phenotype['y'][:]
        sample_ids = phenotype['sample_ids'][:]
        phenotype_name = phenotype['name'].value
        phenotype_name = phenotype_name.replace(" ","_").replace("<i>","").replace("</i>","")
        f = open(os.path.join(output_dir,phenotype_name + ".pheno"),'w')
        f.write("FID IID " + phenotype_name + "\n")
        for i in xrange(sample_ids.shape[0]):
            string = sample_ids[i] + " " + sample_ids[i] + " " + str(y[i])
            f.write(string + "\n")
        f.close()
    hf.close()

'''
CONVERT HDF5 to PLINK: TRAINING AND TESTING
'''
def convertHDF5_2_PLINK_Split(arguments=None):
    if not os.path.isfile(arguments.hdata):
        print "Argument --hdata " + arguments.hdata + " is not a file\n"
        quit()
    hf = h5py.File(arguments.hdata,'r')
    output_dir = arguments.spout
    
    o_raw = hf['Genotype/raw'][:]
    o_sample_ids = hf['Genotype/sample_ids'][:]
    chr_index = hf['Genotype/chr_index'][:]
    pos_index = hf['Genotype/position_index'][:]
    identifiers = hf['Genotype/identifiers'][:]
    #SPLIT DATA
    ind = sp.random.permutation(o_sample_ids.shape[0])
    ratio = sp.floor(float(o_sample_ids.shape[0])*arguments.ratio)
    train_indices = ind[0:(ind.shape[0]-ratio)]
    test_indices = ind[(ind.shape[0]-ratio):]
    

    #WRITE TRAINING DATA
    raw = o_raw[train_indices,:]
    sample_ids = o_sample_ids[train_indices]
    
    #write PED file
    f = open(output_dir + ".train.ped",'w')
    for i in xrange(sample_ids.shape[0]):
        string = sample_ids[i] + " " + sample_ids[i] + " 0 0 0 0 "
        for j in xrange(raw.shape[1]):
            snps = iupac_map_reverse[raw[i,j]]
            if len(snps)==2:
                string += snps[0] + " " + snps[1] + " "
            else:
                string += snps + " " + snps + " "
        f.write(string[:-1] + "\n")
    f.close()

    #write MAP file
    f = open(output_dir + ".train.map",'w')
    for i in xrange(chr_index.shape[0]):
        chr_id = chr_index[i].replace("Chr","").replace("chr","")
        f.write(chr_id + " " + identifiers[i] + " 0 " + str(int(pos_index[i])) + "\n")
    f.close()
    
    #WRITE TESTING DATA
    raw = o_raw[test_indices,:]
    sample_ids = o_sample_ids[test_indices]
    #write PED file
    f = open(output_dir + ".test.ped",'w')
    for i in xrange(sample_ids.shape[0]):
        string = sample_ids[i] + " " + sample_ids[i] + " 0 0 0 0 "
        for j in xrange(raw.shape[1]):
            snps = iupac_map_reverse[raw[i,j]]
            if len(snps)==2:
                string += snps[0] + " " + snps[1] + " "
            else:
                string += snps + " " + snps + " "
        f.write(string[:-1] + "\n")
    f.close()

    #write MAP file
    f = open(output_dir + ".test.map",'w')
    for i in xrange(chr_index.shape[0]):
        chr_id = chr_index[i].replace("Chr","").replace("chr","")
        f.write(chr_id + " " + identifiers[i] + " 0 " + str(int(pos_index[i])) + "\n")
    f.close()

    #write phenotype file
    phenotype_ids = hf['Phenotypes'].keys()
    for pid in phenotype_ids:
        phenotype = hf['Phenotypes'][pid]
        o_y = phenotype['y'][:]
        o_sample_ids = phenotype['sample_ids'][:]
        phenotype_name = phenotype['name'].value
        phenotype_name = phenotype_name.replace(" ","_").replace("<i>","").replace("</i>","")
        sample_ids = o_sample_ids[train_indices]
        y = o_y[train_indices]
        f = open(os.path.join(output_dir,phenotype_name + ".train.pheno"),'w')
        f.write("FID IID " + phenotype_name + "\n")
        for i in xrange(sample_ids.shape[0]):
            string = sample_ids[i] + " " + sample_ids[i] + " " + str(y[i])
            f.write(string + "\n")
        f.close()

        sample_ids = o_sample_ids[test_indices]
        y = o_y[test_indices]
        f = open(os.path.join(output_dir,phenotype_name + ".test.pheno"),'w')
        f.write("FID IID " + phenotype_name + "\n")
        for i in xrange(sample_ids.shape[0]):
            string = sample_ids[i] + " " + sample_ids[i] + " " + str(y[i])
            f.write(string + "\n")
        f.close()
    
    hf.close()
