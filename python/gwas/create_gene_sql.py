import sqlite3,sys

def getGFFInfoElements(info=None):
    elements = {}
    for element in info:
        sv = element.strip().split("=")
        if len(sv)==2:
            elements[sv[0].strip().upper()] = sv[1].strip()
    return elements

def createGFFDB(gff_file=None,sqlite3_file=None):
    gfff = open(gff_file,'r')
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
                    return [None,"[GFF File]: Wrong file format in GFF file. Delimiter has to be either a whitespace or a tab."]
        if len(sv)<8 or len(sv)==0:
            return [None,"[GFF File]: Wrong file format. GFF file must contain 9 columns. See FAQ for details."]
        if sv[2].strip().upper() == "CHROMOSOME" or sv[2].strip().upper() == "CHROMOSOME_ARM":
            info = sv[8].strip().split(";")
            elements = getGFFInfoElements(info)
            if len(elements)==0:
                return [None,"[GFF File]: Wrong file format. The attribute column has to less elements (Column 9). At least the ID field has to be specified. See FAQ for details."]
            if not sv[0].strip()==elements["ID"]:
                return [None,"[GFF File]: Wrong file format. Chromosome ID in column 1 has to be the same identifier as the ID in column 9. See FAQ for details."]
                return [None,"[GFF File]: Chromosome IDs do not match. Unknown chromosome identifier found in GFF file: " + elements["ID"]]
        elif sv[2].strip().upper()=="GENE":
            info = sv[8].strip().split(";")
            elements = getGFFInfoElements(info)
            if len(elements)==0:
                return [None,"[GFF File]: Wrong file format. The attribute column has to less elements (Column 9). At least the ID field has to be specified. See FAQ for details."]
            if not (sv[6].strip()=="+" or sv[6].strip()=="-" or sv[6].strip()=="."):
                return [None,"[GFF File]: Strand has to be either + or -. To indicate missing information use . (dot). See FAQ for details."] 
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
    sqlite = sqlite3.connect(sqlite3_file)
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
    return [gene_set,""]

if __name__ in "__main__":
    if len(sys.argv)<2:
        print "USAGE: GFF_File OutputFile"
        print
        quit()

    [gene_set,message] = createGFFDB(sys.argv[1],sys.argv[2])
    if gene_set==None:
        print message
    else:
        print "\nDone\n"
