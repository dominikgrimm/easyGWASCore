import sys,sqlite3,h5py,os
import scipy as sp

if __name__ in "__main__":

    if len(sys.argv)<4:
        print "USAGE: ResultFile SQLITE3File TOPX output"
        quit()

    top = int(sys.argv[3])
    f = h5py.File(sys.argv[1],'r')
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

    sqlite = sqlite3.connect(sys.argv[2])
    sqlite_cursor = sqlite.cursor()

    out = open(os.path.join(sys.argv[4],name + ".csv"),"w")

    out.write("Chr,Pos,PVal,GeneID (closest),Distance (bp)\n")
    for i in xrange(chromosomes.shape[0]):
        sqlite_cursor.execute("SELECT * FROM geneannotation WHERE chromosome_id=? ORDER BY ABS(annotation_start - ?) LIMIT 1",(str("Chr" + chromosomes[i]),int(positions[i])))
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
