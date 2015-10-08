import sys,sqlite3,h5py,os
import scipy as sp

if __name__ in "__main__":

    if len(sys.argv)<3:
        print "USAGE: SQLITE3File DataFile output [window_around_gene, default=0bp]"
        quit()
    
    sqlite = sqlite3.connect(sys.argv[1])
    sqlite_cursor = sqlite.cursor()
    
    out_dir = sys.argv[3]
    if len(sys.argv)==3:
        window = float(sys.argv[4])
    else:
        window = 0
    
    #f = h5py.File(sys.argv[2],'r')
    #f.close()

    #out = open(os.path.join(sys.argv[4],name + ".csv"),"w")
    
    sqlite_cursor.execute("SELECT * FROM geneannotation")
    annotation = sqlite_cursor.fetchall()
    import pdb
    pdb.set_trace()

    #out.write("Chr,Pos,PVal,GeneID (closest),Distance (bp)\n")
    #for i in xrange(chromosomes.shape[0]):
    if len(annotation)==1:
            if positions[i] >= annotation[0][3] and positions[i] <= annotation[0][4]:
                distance = 0
            elif positions[i] > annotation[0][4]:
                distance = abs(positions[i]-annotation[0][4])
            else:
                distance = abs(positions[i]-annotation[0][3])
            out.write(chromosomes[i] + "," + str(int(positions[i])) + ",%.2e"%(p_values[i]) + "," + annotation[0][1] + "," + str(int(distance)) + "\n")
    sqlite.close()
