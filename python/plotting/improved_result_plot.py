import pylab as pl
import scipy as sp
import sys,os
import scipy.stats as stats
import h5py
import matplotlib.colors as mcol
import matplotlib.cm as cm
import matplotlib as mpl

sys.path.append("bin/" + sys.platform + "/interfaces/python/")
import CEasyGWAS as gwas_core

def encodeHomozygousData(raw_data=None):
    gwas_data = gwas_core.CGWASDataHelper()
    gwas_data.encodeHomozygousData(raw_data,raw_data.shape[1],raw_data.shape[0])
    encoded = gwas_data.getEncodedData()
    maf_data = gwas_data.getMAF()
    gwas_data.releaseMemory()
    return [encoded,maf_data]

def encodeHeterozygousData(raw_data=None,snp_encoding="additive"):
    gwas_data = gwas_core.CGWASDataHelper()
    if snp_encoding=="recessive":
        encoding = gwas_data.recessive
    elif snp_encoding=="dominant":
        encoding = gwas_data.dominant
    elif snp_encoding=="codominant":
        encoding = gwas_data.codominant
    else:
        encoding = gwas_data.additive
    gwas_data.encodeHeterozygousData(raw_data,raw_data.shape[1],raw_data.shape[0],encoding)
    encoded = gwas_data.getEncodedData()
    maf_data = gwas_data.getMAF()
    gwas_data.releaseMemory()
    return [encoded,maf_data]

def remove_border(axes=None, top=False, right=False, left=True, bottom=True):
    ax = axes or pl.gca()
    ax.spines['top'].set_visible(top)
    ax.spines['right'].set_visible(right)
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)
    
    #turn off all ticks
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')
    
    #now re-enable visibles
    if top:
        ax.xaxis.tick_top()
    if bottom:
        ax.xaxis.tick_bottom()
    if left:
        ax.yaxis.tick_left()

def read_result_files(fileName=None):
    f = h5py.File(sys.argv[1],'r')
    pv = f['p_values'][:]
    positions = f['positions'][:]
    chromosomes = f['chromosomes'][:]
    hashs = f['snp_hash'][:]
    name = f['phenotype_name'].value.replace(" ","_").replace("<i>","").replace("</i>","")
    #transformed = f['transformed'].value
    #if transformed==True:
    #    name = name + " (boxcox)"
    f.close()
    chromosomes = sp.array(chromosomes)
    positions = sp.array(positions)
    pv = sp.array(pv)
    hashs = sp.array(hashs)
    tmp,u_ind = sp.unique(hashs,return_index=True)
    unique_pv = pv[u_ind]
    #pv = unique_pv
    #positions = positions[u_ind]
    #chromosomes = chromosomes[u_ind]
    return [pv,positions,chromosomes,hashs,unique_pv,name]

if __name__ in '__main__':
    
    if len(sys.argv)<3:
        print
        print "Usage: result_plots.py <result_file_name> <output_folder> <snp_data> [use unique snp true|false] [encoding=default:additive]"
        print
        quit(0);
    
    font_size = 10
    mpl.rcParams['font.family']="sans-serif"
    mpl.rcParams['font.sans-serif']="Arial"
    mpl.rcParams['font.size']=font_size
    mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['font.weight']='medium'
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['patch.edgecolor'] = 'white'
    mpl.rcParams['grid.linestyle'] = '-'
    mpl.rcParams['grid.color'] = '#FFFFFF'

    [pv,positions,chromosomes,hashs,unique_pv,fname] = read_result_files(sys.argv[1])
    
    chrom_list = sp.unique(chromosomes)
    unsnps = False
    encoding = "additive"
    if len(sys.argv)==5:
        if sys.argv[4]=="true":
            unsnps = True
    if len(sys.argv)==6:
        if sys.argv[4]=="true":
            unsnps = True
        encoding = sys.argv[5]

    f = h5py.File(sys.argv[3],'r')
    raw = f['Genotype']['raw'][:]
    [encoded_tmp,maf] = encodeHeterozygousData(raw,encoding)
    real_pos_tmp = f['Genotype']['position_index'][:]
    real_chr_tmp = f['Genotype']['chr_index'][:]
    encoded = None
    real_pos = None
    real_chr = None
    
    for i,chrom in enumerate(chrom_list):
        idx = chromosomes==chrom
        _pos = positions[idx]

        idx = real_chr_tmp==chrom
        _real_chr = real_chr_tmp[idx]
        _real_pos = real_pos_tmp[idx]
        _encoded = encoded_tmp[:,idx]
        idx = (sp.reshape(_real_pos,(_real_pos.shape[0],1))==_pos).nonzero()[0]
        
        if i==0:
            real_chr = _real_chr[idx]
            real_pos = _real_pos[idx]
            encoded = _encoded[:,idx]
        else:
            real_chr = sp.concatenate([real_chr,_real_chr[idx]])
            real_pos = sp.concatenate([real_pos,_real_pos[idx]])
            encoded = sp.column_stack([encoded,_encoded[:,idx]])
    f.close()
    
    pl.ion()
    pl.figure(figsize=(15,4))

    #fname = sys.argv[1].split("90.")[1].split(".")[0] + "." + sys.argv[1].split("90.")[1].split(".")[1]

    #chrom_list = ['1','2','3','4','5']
    color_list = ['#F26C4F','#F68E55','#7CC576','#00BFF3','#605CA8','#F06EA9','#F26C4F','#F68E55']
    edge_list = ['#ED1C24','#F26522','#39B54A','#00AEEF','#2E3192','#EC008C','#ED1C24','#F26522']
    darker_list = ['#790000','#7B2E00','#005E20','#005B7F','#0D004C','#7B0046','#790000','#7B2E00']
    
    print pv.shape
    print unique_pv.shape
    if unsnps:
        bf_threshold = 0.05/unique_pv.shape[0]
    else:
        bf_threshold = 0.05/pv.shape[0]

    #bf_threshold = 0.05/(128558*2)
    
    ''' IGNORE IS NOT WORKING
    thr_snps = 0.05
    ind = sp.where(pv<0.05)[0]
    pv = pv[ind]
    positions = positions[ind]
    chromosomes = chromosomes[ind]
    hashs = hashs[ind]
    '''

    current_pos = 0

    max_y = (-sp.log10(pv[:])).max()
    if max_y<-sp.log10(bf_threshold):
        max_y = -sp.log10(bf_threshold)
    max_y += 1

    split_list = []
    xtick_list = []
    pl.axhline(-sp.log10(bf_threshold),color='#F68E55',linestyle='--')
    for i,chrom in enumerate(chrom_list):
        idx = chromosomes==chrom
        _pos = positions[idx]
        _chrs = chromosomes[idx]
        _encoded = encoded[:,idx]
        _pv = pv[idx]
        idx = sp.where(_pv>bf_threshold)[0]
        __pos = _pos[idx]
        __chrs = _chrs[idx]
        __pv = _pv[idx]
        __pv = -sp.log10(__pv)
        
        #pl.plot(sp.arange(0+current_pos,current_pos+__pv.shape[0]),__pv,'.',#color=color_list[i],alpha=0.2,
        pl.plot(__pos+current_pos,__pv,'.',#color=color_list[i],alpha=0.2,
                alpha=0.3,
                color="#A0A0A0",
                #markerfacecolor="#FFFFFF",
                #markeredgecolor="LightGray",
                #markeredgewidth=1.5,
                #markerfacecolor=color_list[i],
                #markerfacecoloralt=color_list[i],
                #markeredgecolor=edge_list[i],
                markersize=5,
                )
        
        idx = sp.where(_pv<=bf_threshold)[0]
        __pos = _pos[idx]
        __chrs = _chrs[idx]
        __pv = _pv[idx]
        __pv = -sp.log10(__pv)

        #pl.plot(current_pos+idx,__pv,'.',#color=darker_list[i],#alpha=0.6,
        pl.plot(__pos+current_pos,__pv,'.',#color=darker_list[i],#alpha=0.6,
                #markerfacecolor=edge_list[i],
                #markerfacecoloralt=edge_list[i],
                #markeredgecolor="#D7D7D7",
                #markeredgecolor="#0EBFE9",
                #markerfacecolor="#FFFFFF",
                #markeredgewidth=1.5,
                color="#b000ff",
                alpha=0.9,
                #markeredgewidth=0,
                markersize=8,
                )

        #xtick_list.append(float(current_pos) + float(_pv.shape[0])/2.0)
        xtick_list.append(float(current_pos) + float(_pos.max())/2.0)
        current_pos = _pos.max() + current_pos
        #current_pos += _pv.shape[0]+1
        #split_list.append(current_pos-1)
        split_list.append(current_pos)
        pl.xlim(0,current_pos)
        pl.ylabel("$-log10(p-value)$")
    
    pl.ylim(0,max_y)
    
    for i in xrange(len(split_list)):
        if i%2==0:
            try:
                pl.fill_between([split_list[i],split_list[i+1]],0,max_y,color="#D7D7D7",linewidth=0,alpha=0.5)
            except:
                pass
    s_list = split_list[:-1]

    for split in s_list:
        pl.axvline(split,color='k',linestyle='--')
    
    pl.xticks(xtick_list,chrom_list)
    pl.subplots_adjust(left=0.05,bottom=0.08,right=0.99,top=0.9,wspace=0.45)
    remove_border()
    if unsnps:
        pv = unique_pv
    gc = sp.median(stats.chi2.isf(pv,1))/0.456

    pl.title("Phenotype: %s, genomic control $\lambda=%0.2f$"%(fname,gc))
    pl.savefig(os.path.join(sys.argv[2],'manhattan_' + fname.replace(" (boxcox)","")+'.png') )
    

    for i,chrom in enumerate(chrom_list):
        idx = chromosomes==chrom
        _pos = positions[idx]
        _chrs = chromosomes[idx]
        _encoded = encoded[:,idx]
        _pv = pv[idx]
        
        idx = sp.where(_pv<=bf_threshold)[0]
        __pos = _pos[idx]
        __chrs = _chrs[idx]
        __pv = _pv[idx]
        __pv = -sp.log10(__pv)

        
        #plot LD
        cm1 = mcol.LinearSegmentedColormap.from_list("dgcolor",["#00A651","#39B54A","#8DC73F","#F7941D","#F26522","#ED1C24"])
        cnorm = mcol.Normalize(vmin=0.0,vmax=1)
        #cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
        cpick = cm.ScalarMappable(norm=cnorm,cmap=pl.get_cmap("jet"))
        cpick.set_array([])
        ind = sp.argsort(_pos)
        _pos = _pos[ind]
        _encoded = _encoded[:,ind]
        _pv = _pv[ind]
        nr_snps = 1000
        for pp in __pos:
            fig = pl.figure()
            pl.axhline(-sp.log10(bf_threshold),color='#F68E55',linestyle='--')
            #plot significant point
            pl.plot(__pos,__pv,'.',
                color="#b000ff",
                alpha=0.9,
                markersize=8,
                )
            idx = sp.where(pp==_pos)[0][0]
            if idx-nr_snps>0 and idx+nr_snps<=_encoded.shape[1]:
                ranges = sp.arange(idx-nr_snps,idx+nr_snps)
            elif idx-nr_snps<0 and idx+nr_snps<=_encoded.shape[1]:
                ranges = sp.arange(0,idx+nr_snps)
            elif idx-nr_snps>0 and idx+nr_snps>_encoded.shape[1]:
                ranges = sp.arange(idx-nr_snps,_encoded.shape[1])
            ind = sp.where(ranges==idx)[0]
            ranges = sp.delete(ranges,ind)
            for sind in ranges:
                rr = stats.pearsonr(_encoded[:,sind],_encoded[:,idx])[0]**2
                if rr>=0.7:
                    pl.plot(_pos[sind],-sp.log10(_pv[sind]),'.',color=cpick.to_rgba(rr),alpha=0.9,markersize=8)
                elif rr>=0.5:
                    pl.plot(_pos[sind],-sp.log10(_pv[sind]),'.',color=cpick.to_rgba(rr),alpha=0.9,markersize=8)
                elif rr>=0.3:
                    pl.plot(_pos[sind],-sp.log10(_pv[sind]),'.',color=cpick.to_rgba(rr),alpha=0.9,markersize=7)
                else:
                    pl.plot(_pos[sind],-sp.log10(_pv[sind]),'.',color=cpick.to_rgba(rr),alpha=0.4,markersize=5)
                #if _pv[sind]<1.0 and rr>=0.3:
                #if _pv[sind]<1.0 and rr>=0.3:
                    #pl.plot(_pos[sind]+current_pos,-sp.log10(_pv[sind]),'.',color=cpick.to_rgba(rr),alpha=1,markersize=6)
                #    pl.plot(_pos[sind],-sp.log10(_pv[sind]),'.',color=cpick.to_rgba(rr),alpha=1,markersize=6)
                #elif _pv[sind]<0.005:
                #    pl.plot(_pos[sind],-sp.log10(_pv[sind]),'.',color=cpick.to_rgba(rr),alpha=1,markersize=6)
            pl.subplots_adjust(left=0.05,bottom=0.08,right=0.99,top=0.9,wspace=0.45)
            remove_border()
            pl.ylabel("$-log10(p-value)$")
            pl.colorbar(cpick,label="SNP r^2",shrink=0.5,drawedges=False)
            pl.ylim(0,max_y)
            pl.title("Phenotype: %s"%(fname))
            pl.savefig(os.path.join(sys.argv[2],'ld_plot_' + str(pp) + "_" + fname.replace(" (boxcox)","")+'.png') )
            
    #pl.fill_between(x=[split_list[0],split_list[1]],y1=0,y2=max_y,color="#E5FFF0",linewidth=0,alpha=0.5)
    #pl.grid(True)
    #compute gc
    

    #print QQ-plot 
    pl.figure(figsize=(5,5))
    pv_uni = (sp.arange(1.0/float(pv.shape[0]),1,1.0/float(pv.shape[0]+1)))
    pl.plot(-sp.log10(pv_uni),-sp.log10(sp.sort(pv_uni)),'b--')
    pl.ylim(0,(-sp.log10(pv[:])).max()+1)
    pl.plot(-sp.log10(pv_uni),-sp.log10(sp.sort(pv[:],axis=0)),'.',color="#F68E55",markersize=4,markeredgewidth=0,alpha=0.8)
    #pl.grid(True)
    pl.title("Phenotype: %s"%(fname))
    pl.xlabel('Expected $-log10(p-value)$')
    pl.ylabel('Observed $-log10(p-value)$')
    pl.text(4,1,"$\lambda=%.2f$"%(gc))
    remove_border()
    pl.subplots_adjust(left=0.1,bottom=0.15,right=0.98,top=0.9,wspace=0.45)
    pl.savefig(os.path.join(sys.argv[2],'qqplot_' + fname.replace(" (boxcox)","")+'.png') )
    
    '''
    ####new qq plot######
    gc = sp.median(stats.chi2.isf(pv,1))/0.456
    gc_null = (0.5 + sp.arange(pv.shape[0]))/float(pv.shape[0])

    pl.figure(figsize=(5,5))
    qq_null = -sp.log10(gc_null)
    qq = -sp.log10(sp.sort(pv[:],axis=0))

    pl.plot(qq_null,qq,'.',color="#F68E55",markersize=4,markeredgewidth=0,alpha=0.8)
    pl.plot([0,qq_null.max()],[0,qq_null.max()],'b--')
    pl.title("Phenotype: %s"%(fname))
    pl.xlabel('Expected $-log10(p-value)$')
    pl.ylabel('Observed $-log10(p-value)$')
    pl.text(4,1,"$\lambda=%.2f$"%(gc))
    
    #plot theoretical expectations
    mRange=10**(sp.arange(sp.log10(0.5),sp.log10(pv.shape[0]-0.5)+0.1,0.1))
    numPts=len(mRange)
    betaalphaLevel=sp.zeros(numPts)
    betaOneMinusalphaLevel=sp.zeros(numPts)
    betaInvHalf=sp.zeros(numPts)
    for n in xrange(numPts):
        m=mRange[n]
        betaInvHalf[n]=stats.beta.ppf(0.5,m,pv.shape[0]-m)
        betaalphaLevel[n]=stats.beta.ppf(0.05,m,pv.shape[0]-m)
        betaOneMinusalphaLevel[n]=stats.beta.ppf(1-0.05,m,pv.shape[0]-m)
        betaDown=betaInvHalf-betaalphaLevel
        betaUp=betaOneMinusalphaLevel-betaInvHalf
        theoreticalPvals=mRange/pv.shape[0]

    lower = -sp.log10(theoreticalPvals-betaDown)
    upper = -sp.log10(theoreticalPvals+betaUp)
    pl.fill_between(-sp.log10(theoreticalPvals),lower,upper,color='#00BFF3',alpha=0.4,linewidth=0)
    
    
    remove_border()
    pl.subplots_adjust(left=0.1,bottom=0.15,right=0.98,top=0.9,wspace=0.45)
    pl.savefig(os.path.join(sys.argv[2],'qqplot_new_' + fname.replace(" (boxcox)","")+'.png') )
    '''
