'''
Copyright by Dominik G. Grimm
'''
import sys, argparse, os, h5py
import scipy as sp
import multiprocessing
from utils.dataio import read_HDF5_file,plink2HDF5,encodeData,convertHDF5_2_VCF,getEncodedData,createGFFDB
from utils.dataio import readSIFT4GFile,convertHDF5_2_CSV,writeTopXGenes,writeLDInfo,convertHDF5_2_PLINK
from utils.dataio import convertHDF5_2_PLINK_Split,readPathogenicityFile,read_CSV_file,addCovariates2HDF5
from utils.dataio import addPhenotypes2HDF5
from utils.plotting import ManhattanPlot,QQPlot,LDPlot
from gwas.experiment import GWASettings
from gwas.gwas import start_gwas

def main():
    #PARSE PLOTTING ARGUMENTS
    plot_parser = argparse.ArgumentParser(add_help=False)
    plot_parser.add_argument("--hfile",action="store",dest="hfile",help="HDF5 result input file or directory with several input files")
    plot_parser.add_argument("--csvfile",action="store",dest="csvfile",help="CSV result input file or directory with several input files")
    plot_parser.add_argument("--out",action="store",dest="out",help="Path to output directory (Required)",default="")
    plot_parser.add_argument("--iformat",action="store",dest="iformat",help="Image File Format (Default: png)",default="png",choices=('png', 'pdf','tiff'))
    plot_parser.add_argument("--nogc",action="store_false",default=True,help="Do not compute genomic control (GC) and display in plot",dest="gc")
    plot_parser.add_argument("--phenotype_id",action="store",dest="phenotype_id",type=int,help="Specifiy certain phenotype in HDF5 file. If not specified loop over all phenotypes (default=all)",default=-1)
    #mutgroup = plot_parser.add_mutually_exclusive_group()
    #mutgroup.add_argument("--cfile",action="store",dest="cfile",help="CSV delimited input file or directory")
    #mutgroup.add_argument("--tfile",action="store",dest="tfile",help="TAB delimited input file or directory")
    
    plotting_common_parser = plot_parser.add_argument_group("Common Plotting Parameters")
    plotting_common_parser.add_argument("--distinct","-d",action="store_true",dest="distinct",default=False,help="Use distinct SNPs only to compute multiple hypothesis threshold")
    plotting_common_parser.add_argument("--notitle",action="store_false",dest="title",default=True,help="Add title to plots including Phenotype name and lambda value")
    plotting_common_parser.add_argument("--nhypothesis","-n",action="store",default=-1,type=int,dest="nr_hypothesis",help="Number of hypothesis to correct for (default: all markers)")
    plotting_common_parser.add_argument("--nhypothesis2","-n2",action="store",default=-1,type=int,dest="nr_hypothesis2",help="Second Number of hypothesis to correct for (default: no markers)")
    plotting_common_parser.add_argument("--ignore",action="store",help="Ignore Phenotypes that contain a certain string (Optional)")

    manhattan_parser = plot_parser.add_argument_group("Manhattan-Plots")
    manhattan_parser.add_argument("--manhattan",action="store_true",default=False,help="Create a Manhattan plot")
    
    qq_parser = plot_parser.add_argument_group("QQ-Plots")
    qq_parser.add_argument("--qqplot",action="store_true",default=False,help="Create a Quantile-Quantile (QQ)-plot")
    qq_parser.add_argument("--estpv",action="store_true",default=False,help="Add the estimated theoretical destribution of p-values to the plot")
    
    ld_parser = plot_parser.add_argument_group("LD-Plots")
    ld_parser.add_argument("--ldplot",action="store_true",default=False,help="Create a Manhattan Linkage Plot")
    ld_parser.add_argument("--hdata",action="store",dest="hdata",help="HDF5 File containing the genotype data (Required)")
    ld_parser.add_argument("--snp",action="store",dest="selected_snp",help="SNP identifier that should be used for analysis (default: loop over all singificantly associated SNPs)",default="-1")
    ld_parser.add_argument("--distance",action="store",dest="distance",type=int,help="Distance in bp around the selected SNP (default=10000)",default=10000)
    ld_parser.add_argument("--r2-measure",action="store",dest="r2_measure",help="Choice of Linkage Disequilibrium measure (default: Excoffier-Slatkin)",default="excoffier_slatkin",choices=('excoffier_slatkin','pearson_r2','roger_huff'))
    ld_parser.add_argument("--sql_gene",action="store",dest="sql_gene",help="Add gene annotations to plot. Requires a SQL file generated from GFF file (see data manipulation commands) (Optional)")
    ld_parser.add_argument("--pathogenicity_scores",action="store",dest="sfile",help="Add pathogenicity scores to LD plot (Optional)")
    ld_parser.add_argument("--maf",action="store",dest="maf",type=float,help="Remove SNPs with a population based minor allele frequency (MAF) smaller than specified (default=0)",default=0)
    
    #GWAS PARSING OPTION
    gwas_parser = argparse.ArgumentParser(add_help=False)
    g_parser = gwas_parser.add_argument_group("General parameters shared by all algorithms")
    g_parser.add_argument("--out",action="store",dest="out",help="Path to output directory (Required)",default="")
    g_parser.add_argument("--hdata",action="store",dest="hdata",help="HDF5 File containing the genotype data (Required)")
    g_parser.add_argument("--maf",action="store",dest="maf",type=float,help="Remove SNPs with a minor allele frequency (MAF) smaller than specified (default=0)",default=0)
    g_parser.add_argument("--encoding",action="store",dest="encoding",default="additive",help="Encode Genotype with a different encoding (Default: additive)",choices=('additive','dominant','recessive','overdominant'))
    g_parser.add_argument("--transform",action="store",dest="transform",default=None,help="Transform Phenotype (default: No transformation)",choices=('sqrt','log10','boxcox','zeroMean','unitVariance'))
    g_parser.add_argument("--homozygous",action="store_true",dest="homozygous",help="Genotype is homozygous (default=False)",default=False)
    g_parser.add_argument("--phenotype_id",action="store",dest="phenotype_id",type=int,help="Specifiy certain phenotype in HDF5 file. If not specified loop over all phenotypes (default=all)",default=-1)
    g_parser.add_argument("--algorithm",action="store",dest="algorithm",default="linear",help="Select Algorithm",choices=('linear','logit','FaSTLMM','EMMAX','ttest','fisher','WCrt','MWUrt','linearperm','logitperm','EMMAXperm'))
    g_parser.add_argument("--threads",action="store",dest="threads",type=int,help="If mutliple phenotypes in file parallelize computations (default=Available CPUs - 1)",default=multiprocessing.cpu_count()-1)
    g_parser.add_argument("--pcs",action="store",dest="principle_components",type=int,help="Number of Principle Components (default=0)",default=0)
    g_parser.add_argument("--pc_iterative",action="store_true",dest="pc_iterative",help="Iterate through the number of 'pcs' (default=False)",default=False)

    algo_parser = gwas_parser.add_argument_group("Linear Mixed Model specific parameters")
    algo_parser.add_argument("--unique_snps",action="store_true",dest="unique_snps",help="Use unique SNPs only to compute kinship matrix (default=False)",default=False)
    

    #DATA PARSING OPTION
    data_parser = argparse.ArgumentParser(add_help=False)
    #mutgroup = resfile_parser.add_mutually_exclusive_group()
    #mutgroup.add_argument("--hdata",action="store",dest="hfile",help="easyGWASCore HDF5 input file")
    #mutgroup.add_argument("--plink_data",action="store",dest="cfile",help="Prefix of Plink input files (*.ped, *.map)")

    plink_parser = data_parser.add_argument_group("Convert Plink into easyGWASCore HDF5 file")
    plink_parser.add_argument("--plink2hdf5",action="store_true",default=False,help="Convert Plink into an easyGWASCore HDF5 file")
    plink_parser.add_argument("--plink_data",action="store",dest="plink_data",help="Prefix of Plink input files (*.ped, *.map)")
    plink_parser.add_argument("--hout",action="store",dest="hout",help="Filename for HDF5 output file")
    plink_parser.add_argument("--plink_phenotype",action="store",dest="plink_phenotype",help="Plink input file or directory with plink input files (Optinal)")
    plink_parser.add_argument("--maf",action="store",dest="maf",type=float,help="Remove SNPs with a population based minor allele frequency (MAF) smaller than specified (default=0)",default=0)
    plink_parser.add_argument("--exclude_snps",action="store",dest="exclude_snps",help="Remove a list of SNP identifiers (Optional)")
    plink_parser.add_argument("--distinct_filter",action="store",dest="distinct_filter",default=0,type=int,help="Exclude SNPs that are not distinct or share a certain pattern more often than x (default: no filtering)")
    

    resfile_parser = data_parser.add_argument_group("File Input Flags")
    resfile_parser.add_argument("--hdata",action="store",dest="hdata",help="Filename of HDF5 input file (needed for options {--encode, --vcf, --hdf5toplink, --addcovariates})")
    resfile_parser.add_argument("--hfile",action="store",dest="hfile",help="HDF5 result input file or directory with several input files (needed for options {--csv,--ld})")
    
    add_parser = data_parser.add_argument_group("Add additional data to HDF5 file")
    add_parser.add_argument("--adddata",action="store_true",default=False,help="Store additional data to HDF5 file")
    add_parser.add_argument("--addphenotypes",action="store",dest="addphenotypes",help="Add additional phenotypes to HDF5 file (either a folder with files or a single file)")
    add_parser.add_argument("--addcovariates",action="store",dest="addcovariates",help="Add additional covariates to HDF5 file (either a folder with files or a single file)")
    
    plink_parser = data_parser.add_argument_group("Convert HDF5 File into Plink files")
    plink_parser.add_argument("--hdf5toplink",action="store_true",default=False,help="Convert HDF5 file into PLINK files")
    plink_parser.add_argument("--pout",action="store",dest="pout",help="Path with file prefix to PLINK output folder")
    
    plink_parser = data_parser.add_argument_group("Convert HDF5 File into Plink files (SPLIT DATA INTO TRAIN/TEST)")
    plink_parser.add_argument("--hdf5toplink_split",action="store_true",default=False,help="Convert HDF5 file into PLINK files")
    plink_parser.add_argument("--spout",action="store",dest="spout",help="Path with file prefix to PLINK output folder")
    plink_parser.add_argument("--ratio",action="store",dest="ratio",type=float,help="Splitting Ratio Training:Testing Set (default=0.2, 80%% Training, 20%% Testing)",default=0.2)
    
    encode_parser = data_parser.add_argument_group("Encode data in HDF5 file")
    encode_parser.add_argument("--encode",action="store",dest="encode",default=None,help="Encode raw data matrix in HDF5 file and store encoded data in file",choices=('additive','dominant','recessive','overdominant'))
    
    vcf_parser = data_parser.add_argument_group("HDF5 to VCF File")
    vcf_parser.add_argument("--vcf",action="store_true",default=False,help="Convert HDF5 file to VCF output file")
    vcf_parser.add_argument("--vout",action="store",dest="vout",help="Filename for VCF output file")
    
    gff_parser = data_parser.add_argument_group("Gene GFF File to SQLLite Database")
    gff_parser.add_argument("--gff2sql",action="store_true",default=False,help="Store genes from GFF file in local SQLITE3 database")
    gff_parser.add_argument("--gfile",action="store",dest="gfile",help="GFF input file")
    gff_parser.add_argument("--sqlout",action="store",dest="sqlout",help="SQL output file")
    
    csv_parser = data_parser.add_argument_group("HDF5 Result file to CSV File")
    csv_parser.add_argument("--csv",action="store_true",default=False,help="Write HDF5 Results to CSV output file")
    csv_parser.add_argument("--cout",action="store",dest="cout",help="Path to CSV output folder")
    
    
    csv_parser = data_parser.add_argument_group("HDF5 Result file to Gene Annotation Output")
    csv_parser.add_argument("--agene",action="store_true",default=False,help="Write HDF5 Results to an Annotated Gene Output File")
    csv_parser.add_argument("--sqlfile",action="store",dest="sqlfile",help="SQL input file")
    csv_parser.add_argument("--gout",action="store",dest="gout",help="Path to output folder")
    csv_parser.add_argument("--topx",action="store",dest="topx",help="Write top x associations with genes to output file (default=1000)",default=1000,type=int)
    
    ld_parser = data_parser.add_argument_group("Create Linkage-Disequilibrium Output for HDF5 Result Files")
    ld_parser.add_argument("--ld",action="store_true",default=False,help="Create Linkage-Disequilibrium Files")
    ld_parser.add_argument("--ldout",action="store",dest="ldout",help="Path to output folder")
    ld_parser.add_argument("--snp",action="store",dest="selected_snp",help="SNP identifier that should be used for analysis (default: loop over all singificantly associated SNPs)",default="-1")
    ld_parser.add_argument("--distance",action="store",dest="distance",type=int,help="Distance in bp around the selected SNP (default=10000)",default=10000)
    ld_parser.add_argument("--r2-measure",action="store",dest="r2_measure",help="Choice of Linkage Disequilibrium measure (default: Excoffier-Slatkin)",default="excoffier_slatkin",choices=('excoffier_slatkin','pearson_r2','roger_huff'))
    ld_parser.add_argument("--nhypothesis","-n",action="store",default=-1,type=int,dest="nr_hypothesis",help="Number of hypothesis to correct for (default: all markers)")
    ld_parser.add_argument("--distinct","-d",action="store_true",dest="distinct",default=False,help="Use distinct SNPs only to compute multiple hypothesis threshold")
    ld_parser.add_argument("--ignore",action="store",help="Ignore Phenotypes that contain a certain string (Optional)")
    
    #GENERAL ARGUMENT PARSER
    parser = argparse.ArgumentParser(description="easGWASCore: Performing, visualizing and annotating Genome-wide Association Studies\n",
                                     version="1.0")
    subparsers = parser.add_subparsers(help="Subcommands: Please specify the command and use the flag -h to print the parameters for the different subcommands")
    plotting = subparsers.add_parser("plot",help="Create different plots (e.g. Manhattan Plot, QQ-Plot, LD-Plot etc.! To list all options please use 'plot -h')",parents=[plot_parser])
    plotting.set_defaults(which="plot")

    gwas = subparsers.add_parser("gwas",help="Perform a Genome-Wide Association Scan. To list all available options use 'gwas -h')",parents=[gwas_parser])
    gwas.set_defaults(which="gwas")
    
    data = subparsers.add_parser("data",help="Data processing, converting and manipulation methods. To list all available options use 'data -h')",parents=[data_parser])
    data.set_defaults(which="data")
    
    results = parser.parse_args()

    #Plotting
    if results.which=="data":
        enc = False
        if results.encode!=None:
            enc = True
        arg_list = sp.array([results.plink2hdf5,results.vcf,results.gff2sql,enc,results.csv,results.agene,results.ld,results.hdf5toplink,results.hdf5toplink_split,results.adddata])
        if sp.where(arg_list==True)[0].shape[0]>=2 or sp.where(arg_list==True)[0].shape[0]==0:
            data_parser.print_help()
            print "\n-------------------------------------------------------------"
            print "Please select ONE of the following options: {--plink2hdf5, --hdf5toplink, --encode, --vcf, --gff2sql, --csv, --agene, --ld,--hdf5toplink_split,--adddata}\n"
            quit()

        
        if results.plink2hdf5:
            if results.plink_data==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --plink_data has to be set!\n"
                quit()
            if results.hout==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --hout has to be set!\n"
                quit()
            plink2HDF5(results)
        elif enc!=False:
            if results.hdata==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --hdata has to be set!\n"
                quit()
            encodeData(results)
        elif results.hdf5toplink:
            if results.hdata==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --hdata has to be set!\n"
                quit()
            if results.pout==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --pout has to be set!\n"
                quit()
            convertHDF5_2_PLINK(results)
        elif results.hdf5toplink_split:
            if results.hdata==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --hdata has to be set!\n"
                quit()
            if results.spout==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --spout has to be set!\n"
                quit()
            convertHDF5_2_PLINK_Split(results)
        elif results.vcf:
            if results.hdata==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --hdata has to be set!\n"
                quit()
            if results.vout==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --vout has to be set!\n"
                quit()
            convertHDF5_2_VCF(results)
        elif results.gff2sql:
            if results.gfile==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --gfile has to be set!\n"
                quit()
            if results.sqlout==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --sqlout has to be set!\n"
                quit()
            createGFFDB(results)
        elif results.csv:
            if results.hfile==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --hfile has to be set!\n"
                quit()
            if results.cout==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --cout has to be set!\n"
                quit()
            convertHDF5_2_CSV(results)
        elif results.agene:
            if results.hfile==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --hfile has to be set!\n"
                quit()
            if results.gout==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --gout has to be set!\n"
                quit()
            if results.sqlfile==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --sqlfile has to be set!\n"
                quit()
            writeTopXGenes(results)
        elif results.ld:
            if results.hdata==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --hdata has to be set!\n"
                quit()
            if results.hfile==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --hfile has to be set!\n"
                quit()
            if results.ldout==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --ldout has to be set!\n"
                quit()
            writeLDInfo(results)
        elif results.adddata:
            if results.hdata==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Argument --hdata has to be set!\n"
                quit()
            if results.addphenotypes==None and results.addcovariates==None:
                data_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "Either --addphenotypes or --addcovariates has to be set!\n"
                quit()
            if results.addphenotypes!=None:
                addPhenotypes2HDF5(results)
            if results.addcovariates!=None:
                addCovariates2HDF5(results)
            quit()

    elif results.which=="plot":
        #CHECK Additional Paramters
        if results.manhattan==False and results.qqplot==False and results.ldplot==False:
            plot_parser.print_help()
            print "\n-------------------------------------------------------------"
            print "Either --manhattan, --qqplot or --ldplot has to be specified\n"
            quit()
        if results.hfile==None and results.csvfile==None:
            plot_parser.print_help()
            print "\n-------------------------------------------------------------"
            print "--hfile and --csvfile cannot be empty! Either one of both has to be set!\n"
            quit()

        if results.out==None:
            plot_parser.print_help()
            print "\n-------------------------------------------------------------"
            print "--out cannot be empty. Please specify a path for output files!\n"
            quit()
        elif not os.path.isdir(results.out):
            plot_parser.print_help()
            print "\n-------------------------------------------------------------"
            print "--out: " + results.out + " is not a directory!\n"
            quit()
        
        if results.ldplot:
            if not results.hdata==None:
                if not os.path.isfile(results.hdata):
                    plot_parser.print_help()
                    print "\n-------------------------------------------------------------"
                    print "--hdata " + results.hdata + " is not a file!\n" 
                    quit()
            else:
                plot_parser.print_help()
                print "\n-------------------------------------------------------------"
                print "--hdata cannot be empty!\n" 
                quit()
            pathogenicity_map = {}
            if results.sfile!=None:
                if os.path.isfile(results.sfile):
                    #sift_map = readSIFT4GFile(results,read_all=False)
                    pathogenicity_map = readPathogenicityFile(results,read_all=True)
                else:
                    plot_parser.print_help()
                    print "\n-------------------------------------------------------------"
                    print "--pathogenicity_scores " + results.sfile + " is not a file!\n" 
                    quit()

        #DATA READING
        if not (results.hfile==None and results.csvfile==None):
            if results.hfile!=None:
               hdf5_file_flag = True
               rfile = results.hfile
            else:
               hdf5_file_flag = False
               rfile = results.csvfile
            
            if os.path.isdir(rfile):
                for fn in os.listdir(rfile):
                    filename = os.path.join(rfile,fn)
                    if not os.path.isfile(filename):
                        continue
                    else:
                        print "Creating plots for " + str(filename),
                        if hdf5_file_flag:
                            [pv,positions,chromosomes,hashs,unique_pv,phenotype_name] = read_HDF5_file(filename)
                        else:
                            [pv,positions,chromosomes,hashs,unique_pv,phenotype_name] = read_CSV_file(filename)
                        if results.manhattan:
                            ManhattanPlot(results,pv,positions,chromosomes,hashs,unique_pv,phenotype_name)
                        if results.qqplot:
                            QQPlot(results,pv,unique_pv,phenotype_name)
                        if results.ldplot:
                            [encoded,maf,identifiers] = getEncodedData(results.hdata,phenotype_id=fn.split(".")[0],maf=results.maf)
                            LDPlot(results,identifiers,encoded,maf,pv,unique_pv,positions,chromosomes,phenotype_name,pathogenicity_map)
                        print  "\t[DONE]"
            else:
                if hdf5_file_flag:
                    [pv,positions,chromosomes,hashs,unique_pv,phenotype_name] = read_HDF5_file(rfile)
                else:
                    [pv,positions,chromosomes,hashs,unique_pv,phenotype_name] = read_CSV_file(rfile)
                if results.manhattan:
                    ManhattanPlot(results,pv,positions,chromosomes,hashs,unique_pv,phenotype_name)
                if results.qqplot:
                    QQPlot(results,pv,unique_pv,phenotype_name)
                if results.ldplot:
                    print results.hfile
                    [encoded,maf,identifiers] = getEncodedData(results.hdata,phenotype_id=results.hfile.split("/")[-1].split(".")[0],maf=results.maf)
                    LDPlot(results,identifiers,encoded,maf,pv,unique_pv,positions,chromosomes,phenotype_name,pathogenicity_map)
        else:
            plot_parser.print_help()
            print "\n-------------------------------------------------------------"
            print "--hfile and --csvfile cannot be empty! Either one of these paramters has to be set!\n"
            quit()
    
    #PERFORM a GWAS
    elif results.which=="gwas":
        #CHECK Additional Paramters
        if results.hdata==None:
            gwas_parser.print_help()
            print "\n-------------------------------------------------------------"
            print "--hdata cannot be empty!\n"
            quit()
        elif not os.path.isfile(results.hdata):
            gwas_parser.print_help()
            print "\n-------------------------------------------------------------"
            print "--hdata: " + results.hdata + " is not a file!\n"
            quit()

        if results.out==None:
            gwas_parser.print_help()
            print "\n-------------------------------------------------------------"
            print "--out cannot be empty. Please specify a path for output files!\n"
            quit()
        elif not os.path.isdir(results.out):
            gwas_parser.print_help()
            print "\n-------------------------------------------------------------"
            print "--out: " + results.out + " is not a directory!\n"
            quit()
        
        #if only one phenotype is selected
        if results.phenotype_id!=-1:
            if results.pc_iterative==True:
                pool = multiprocessing.Pool(processes=results.threads)
                for pc in xrange(results.principle_components):
                    settings = GWASettings()
                    settings.algorithm = results.algorithm
                    settings.output_path = results.out
                    settings.output_file = str(results.phenotype_id) + "_PCs_" + str(pc) + ".hdf5"
                    settings.hdf5_file = results.hdata
                    settings.phenotype_file = results.hdata
                    settings.snp_encoding = results.encoding
                    settings.phenotype_id = str(results.phenotype_id)
                    settings.phenotype_transformation = results.transform
                    settings.homozygous = results.homozygous
                    settings.principle_components = pc
                    settings.maf = results.maf
                    settings.unique_snps_only = results.unique_snps
                    pool.apply_async(start_gwas,args=(settings,))
                pool.close()
                pool.join()
            else:
                settings = GWASettings()
                settings.algorithm = results.algorithm
                settings.output_path = results.out
                settings.output_file = str(results.phenotype_id) + ".hdf5"
                settings.hdf5_file = results.hdata
                settings.phenotype_file = results.hdata
                settings.snp_encoding = results.encoding
                settings.phenotype_id = str(results.phenotype_id)
                settings.phenotype_transformation = results.transform
                settings.homozygous = results.homozygous
                settings.principle_components = results.principle_components
                settings.maf = results.maf
                settings.unique_snps_only = results.unique_snps
                start_gwas(settings)
        else:
            pool = multiprocessing.Pool(processes=results.threads)
            f = h5py.File(results.hdata,'r')
            try:
                phenotype_ids = f['Phenotypes'].keys()
                f.close()
            except:
                print "\nNo Phenotypes found in file --hdata: " + results.hdata
                f.close()
                quit()
            for pid in phenotype_ids:
                if results.pc_iterative==True:
                    for pc in xrange(results.principle_components):
                        settings = GWASettings()
                        settings.algorithm = results.algorithm
                        settings.output_path = results.out
                        settings.output_file = str(pid) + "_PCs_" + str(pc) + ".hdf5"
                        settings.hdf5_file = results.hdata
                        settings.phenotype_file = results.hdata
                        settings.snp_encoding = results.encoding
                        settings.phenotype_id = str(pid)
                        settings.phenotype_transformation = results.transform
                        settings.homozygous = results.homozygous
                        settings.principle_components = pc
                        settings.maf = results.maf
                        settings.unique_snps_only = results.unique_snps
                        pool.apply_async(start_gwas,args=(settings,))
                else:
                    settings = GWASettings()
                    settings.algorithm = results.algorithm
                    settings.output_path = results.out
                    settings.output_file = str(pid) + ".hdf5"
                    settings.hdf5_file = results.hdata
                    settings.phenotype_file = results.hdata
                    settings.snp_encoding = results.encoding
                    settings.phenotype_id = str(pid)
                    settings.phenotype_transformation = results.transform
                    settings.homozygous = results.homozygous
                    settings.principle_components = results.principle_components
                    settings.maf = results.maf
                    settings.unique_snps_only = results.unique_snps
                    pool.apply_async(start_gwas,args=(settings,))
            pool.close()
            pool.join()

if __name__ in "__main__":
    main()
