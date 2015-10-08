import pylab as pl
import scipy as sp
import sys,os
import scipy.stats as stats
import h5py
import matplotlib.colors as mcol
import matplotlib.cm as cm
import matplotlib as mpl
import estimate_ld as ld
import sqlite3

from general import remove_border


def ManhattanPlot(arguments,pv,positions,chromosomes,hashs,unique_pv,fname): 
    font_size = 14
    mpl.rcParams['font.family']="sans-serif"
    mpl.rcParams['font.sans-serif']="Arial"
    mpl.rcParams['font.size']=font_size
    mpl.rcParams['figure.dpi'] = 300
    mpl.rcParams['font.weight']='medium'
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['patch.edgecolor'] = 'white'
    mpl.rcParams['grid.linestyle'] = '-'
    mpl.rcParams['grid.color'] = 'LightGray'
    
    if arguments.ignore!=None:
        if arguments.ignore in fname:
            return

    chrom_list = sp.unique(chromosomes)
    unsnps = arguments.distinct

    pl.ion()
    pl.figure(figsize=(12,3))

    if arguments.nr_hypothesis==-1:
        if unsnps:
            bf_threshold = 0.05/unique_pv.shape[0]
        else:
            bf_threshold = 0.05/pv.shape[0]
    else:
        bf_threshold = 0.05/arguments.nr_hypothesis

    current_pos = 0

    max_y = (-sp.log10(pv[:])).max()
    if max_y<-sp.log10(bf_threshold):
        max_y = -sp.log10(bf_threshold)
    max_y += 1

    split_list = []
    xtick_list = []
    pl.axhline(-sp.log10(bf_threshold),color='#F68E55',linestyle='--')
    if not arguments.nr_hypothesis2==-1: 
        pl.axhline(-sp.log10(0.05/arguments.nr_hypothesis2),color='#00BFF3',linestyle='--')
        factor = sp.around(float(arguments.nr_hypothesis2)/float(arguments.nr_hypothesis))
        leg = pl.legend(['Bonferroni','Bonferroni ' + str(factor) + 'x'],ncol=2,fancybox=True,prop={'size':12},bbox_to_anchor=(1,1.05))
        leg.get_frame().set_alpha(0.8)
        leg.get_frame().set_linewidth(0.2)

    for i,chrom in enumerate(chrom_list):
        idx = chromosomes==chrom
        _pos = positions[idx]
        _chrs = chromosomes[idx]
        _pv = pv[idx]
        idx = sp.where(_pv>bf_threshold)[0]
        __pos = _pos[idx]
        __chrs = _chrs[idx]
        __pv = _pv[idx]
        __pv = -sp.log10(__pv)
        
        pl.plot(__pos+current_pos,__pv,'.',
                alpha=0.3,
                color="#A0A0A0",
                markersize=5,
                )
        
        idx = sp.where(_pv<=bf_threshold)[0]
        __pos = _pos[idx]
        __chrs = _chrs[idx]
        __pv = _pv[idx]
        __pv = -sp.log10(__pv)

        pl.plot(__pos+current_pos,__pv,'.',
                color="#b000ff",
                alpha=0.9,
                markersize=8,
                )

        xtick_list.append(float(current_pos) + float(_pos.max())/2.0)
        current_pos = _pos.max() + current_pos
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
    pl.subplots_adjust(left=0.09,bottom=0.09,right=0.99,top=0.9,wspace=0.45)
    remove_border()
    if unsnps:
        pv = unique_pv

    if arguments.title:
        if arguments.gc:
            gc = sp.median(stats.chi2.isf(pv,1))/0.456
            pl.title("Phenotype: %s, genomic control $\hat \lambda=%0.2f$"%(fname,gc))
        else:
            pl.title("Phenotype: %s"%(fname))
    pl.savefig(os.path.join(arguments.out,'manhattan_' + fname + '.' + arguments.iformat) )
    pl.close()

def QQPlot(arguments,pv,unique_pv,fname):
    font_size = 18
    mpl.rcParams['font.family']="sans-serif"
    mpl.rcParams['font.sans-serif']="Arial"
    mpl.rcParams['font.size']=font_size
    mpl.rcParams['figure.dpi'] = 300
    mpl.rcParams['font.weight']='medium'
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['patch.edgecolor'] = 'white'
    mpl.rcParams['grid.linestyle'] = '-'
    mpl.rcParams['grid.color'] = 'LightGray'
    if arguments.ignore!=None:
        if arguments.ignore in fname:
            return 
    
    if arguments.distinct:
        pv = unique_pv

    pl.figure(figsize=(5,5))
    pv_uni = (sp.arange(1.0/float(pv.shape[0]),1,1.0/float(pv.shape[0]+1)))
    pl.plot(-sp.log10(pv_uni),-sp.log10(sp.sort(pv_uni)),'b--')
    pl.ylim(0,(-sp.log10(pv[:])).max()+1)
    pl.plot(-sp.log10(pv_uni),-sp.log10(sp.sort(pv[:],axis=0)),'.',color="#F68E55",markersize=12,markeredgewidth=0,alpha=1)
    #plot theoretical expectations
    if arguments.estpv:
        datapoints=10**(sp.arange(sp.log10(0.5),sp.log10(pv.shape[0]-0.5)+0.1,0.1))
        beta_alpha=sp.zeros(datapoints.shape[0])
        beta_nalpha=sp.zeros(datapoints.shape[0])
        beta_tmp=sp.zeros(datapoints.shape[0])
        for n in xrange(datapoints.shape[0]):
            m=datapoints[n]
            beta_tmp[n]=stats.beta.ppf(0.5,m,pv.shape[0]-m)
            beta_alpha[n]=stats.beta.ppf(0.05,m,pv.shape[0]-m)
            beta_nalpha[n]=stats.beta.ppf(1-0.05,m,pv.shape[0]-m)
        estimated_pvals=datapoints/pv.shape[0]
        lower_bound = -sp.log10(estimated_pvals-(beta_tmp-beta_alpha))
        upper_bound = -sp.log10(estimated_pvals+(beta_nalpha-beta_tmp))
        pl.fill_between(-sp.log10(estimated_pvals),lower_bound,upper_bound,color='#00BFF3',alpha=0.4,linewidth=0)
    if arguments.title:
        pl.title("Phenotype: %s"%(fname))
    pl.xlabel('Expected $-log10(p-value)$')
    pl.ylabel('Observed $-log10(p-value)$')
    if arguments.gc:
        gc = sp.median(stats.chi2.isf(pv,1))/0.456
        pl.text(4,1,"$\hat \lambda=%.2f$"%(gc))
    remove_border()
    pl.subplots_adjust(left=0.14,bottom=0.13,right=0.97,top=0.95,wspace=0.45)
    pl.savefig(os.path.join(arguments.out,'qqplot_' + fname + '.' + arguments.iformat) )
    pl.close()

def LDPlot(arguments,identifiers,encoded,maf,pv,unique_pv,positions,chromosomes,fname,pathogenicity_map):
    font_size = 16
    mpl.rcParams['font.family']="sans-serif"
    mpl.rcParams['font.sans-serif']="Arial"
    mpl.rcParams['font.size']=font_size
    #mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['font.weight']='medium'
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['patch.edgecolor'] = 'white'
    mpl.rcParams['grid.linestyle'] = '-'
    mpl.rcParams['grid.color'] = 'LightGray'

    if arguments.ignore!=None:
        if arguments.ignore in fname:
            return

    if not arguments.sql_gene==None:
        sqlite = sqlite3.connect(arguments.sql_gene)
        sqlite_cursor = sqlite.cursor()
    else:
        sqlite = None

    if arguments.nr_hypothesis==-1:
        if arguments.distinct:
            bf_threshold = 0.05/unique_pv.shape[0]
        else:
            bf_threshold = 0.05/pv.shape[0]
    else:
        bf_threshold = 0.05/arguments.nr_hypothesis
    
    snp_distance = arguments.distance
    r2_measure = arguments.r2_measure

    pl.ion()
    pl.figure(figsize=(12,4))

    color_list = ['#F26C4F','#F68E55','#7CC576','#00BFF3','#605CA8','#F06EA9','#F26C4F','#F68E55']
    
    #select different SNP
    if arguments.selected_snp!='-1':
        ind = sp.where(identifiers==arguments.selected_snp)[0]
        if ind.shape[0]==0:
            print "\nSNP " + arguments.selected_snp + " not found in dataset!"
            print "Please select a different SNP identifier!\n"
            quit()
        else:
            chrom_list = sp.array([arguments.selected_snp.split("_")[0]])
            __pos = sp.array([int(arguments.selected_snp.split("_")[1])])
    else:
        chrom_list = sp.unique(chromosomes)

    for i,chrom in enumerate(chrom_list):
        idx = chromosomes==chrom
        _pos = positions[idx]
        _chrs = chromosomes[idx]
        _pv = pv[idx]

        if arguments.selected_snp=='-1':
            idx = sp.where(_pv<=bf_threshold)[0]
            __pos = _pos[idx]
            __chrs = _chrs[idx]
            __pv = _pv[idx]
            __pv = -sp.log10(__pv)
        else:
            idx = sp.where(__pos[0]==_pos)[0]
            __pv = _pv[idx]
            __pv = -sp.log10(__pv)
        #plot LD
        cnorm = mcol.Normalize(vmin=0.0,vmax=1)
        cpick = cm.ScalarMappable(norm=cnorm,cmap=pl.get_cmap("jet"))
        cpick.set_array([])
        ind = sp.argsort(_pos)
        _pos = _pos[ind]
        _pv = _pv[ind]
        for k,pp in enumerate(__pos):
            fig = pl.figure(figsize=(12,4))
            if sqlite==None:
                ax1 = pl.subplot2grid((5,5),(0,0),colspan=5,rowspan=4)
            else:
                ax1 = pl.subplot2grid((6,5),(0,0),colspan=5,rowspan=4)
            ax1.axhline(-sp.log10(bf_threshold),color='#F68E55',linestyle='--')
            marker = '.'
            if str(chrom) + "_" + str(pp) in pathogenicity_map:
                if pathogenicity_map[str(chrom) + "_" + str(pp)].prediction == "DELETERIOUS":
                    marker = '^' 
                else:
                    marker = 'v'
            ax1.plot(pp,__pv[k],marker,
                color="#b000ff",
                alpha=0.9,
                markersize=10,
                )
            if pp-snp_distance > 0 and pp+snp_distance<=_pos.max():
                ranges = sp.where((_pos>=pp-snp_distance) & (_pos<=pp+snp_distance))[0]
            elif pp-snp_distance < 0 and pp+snp_distance<=_pos.max():
                ranges = sp.where(_pos<=pp+snp_distance)[0]
            elif pp-snp_distance > 0 and pp+snp_distance>_pos.max():
                ranges = sp.where(_pos>=pp-snp_distance)[0]
            
            pp_idx = sp.where(pp==_pos)[0][0]
            vline_map = {}
            rr_color_list = []
            maf_list = []
            idx = sp.where(str(chrom) + "_" + str(pp)==identifiers)[0]
            for l,sra in enumerate(ranges):
                sind = sp.where(str(chrom) + "_" + str(_pos[sra])==identifiers)[0]
                maf_list.append(maf[sind])

                if r2_measure=="excoffier_slatkin":
                    rr = ld.esem_r(sp.array(encoded[:,sind].flatten(),dtype="int"),sp.array(encoded[:,idx].flatten(),dtype="int"))**2
                elif r2_measure=="roger_huff":
                    rr = ld.get_r(sp.array(encoded[:,sind].flatten(),dtype="int"),sp.array(encoded[:,idx].flatten(),dtype="int"))**2
                elif r2_measure=="pearson_r2":
                    rr = stats.pearsonr(encoded[:,sind].flatten(),encoded[:,idx].flatten())[0]**2
                rr_color_list.append(cpick.to_rgba(rr))
                if _pos[sra]==pp:
                    ax1.vlines(pp, 0, __pv[k],color='#b000ff',linestyle='--',alpha=0.8,linewidth=0.5)
                    vline_map[pp] = "#b000ff"
                    continue
                marker = '.'
                if str(chrom) + "_" + str(_pos[sra]) in pathogenicity_map:
                    if pathogenicity_map[str(chrom) + "_" + str(_pos[sra])].prediction == "DELETERIOUS":
                        marker = '^' 
                    else:
                        marker = 'v'
                if rr>=0.7:
                    ax1.plot(_pos[sra],-sp.log10(_pv[sra]),marker,color=cpick.to_rgba(rr),alpha=0.9,markeredgecolor='#DDDDDD',markersize=10)#8
                elif rr>=0.5:
                    ax1.plot(_pos[sra],-sp.log10(_pv[sra]),marker,color=cpick.to_rgba(rr),alpha=0.9,markeredgecolor='#DDDDDD',markersize=9)#7
                elif rr>=0.3:
                    ax1.plot(_pos[sra],-sp.log10(_pv[sra]),marker,color=cpick.to_rgba(rr),alpha=0.9,markeredgecolor='#DDDDDD',markersize=8)#6
                else:
                    ax1.plot(_pos[sra],-sp.log10(_pv[sra]),marker,color=cpick.to_rgba(rr),alpha=0.4,markeredgecolor='#DDDDDD',markersize=7)#4
                if _pv[sra]<=bf_threshold:
                    ax1.vlines(_pos[sra], 0, -sp.log10(_pv[sra]),color=cpick.to_rgba(rr),linestyle='--',alpha=0.8,linewidth=0.5)
                    vline_map[_pos[sra]] = cpick.to_rgba(rr)
            remove_border()
            ax1.set_ylabel("$-log10(p-value)$")
            ax1.set_ylim(0,sp.maximum(-sp.log10(bf_threshold)+1,sp.maximum((-sp.log10(_pv[ranges])).max()+1,-sp.log10(_pv[pp_idx])+1)))
            if arguments.title:
                ax1.set_title("Phenotype: %s, SNP: %d"%(fname,pp))
            ax1.set_xticks([])
            ax1.set_yticks(sp.arange(1,sp.ceil(sp.maximum(-sp.log10(bf_threshold)+1,sp.maximum((-sp.log10(_pv[ranges])).max()+1,-sp.log10(_pv[pp_idx])+1))),3))
            ax1.set_xlim(_pos[ranges[0]]-snp_distance*0.01,_pos[ranges[-1]]+snp_distance*0.01)
            ax1.yaxis.grid()
            if not pathogenicity_map==None:
                deleterious = pl.Line2D(range(1), range(1), marker='^', color="white",linestyle="None")
                benign = pl.Line2D(range(1), range(1), marker='v', color="white",linestyle="None")
                leg = pl.legend([benign,deleterious],['Benign Missense Mutation','Deleterious Missense Mutation'],frameon=True,scatterpoints=1,prop={'size':12},
                                fancybox=True,bbox_to_anchor=(1,1.2),numpoints=1,ncol=2)
                leg.get_frame().set_alpha(0.5)
                leg.get_frame().set_linewidth(0.2)
                leg.get_frame().set_facecolor("#DDDDDD")

            if sqlite==None:
                ax3 = pl.subplot2grid((5,5),(4,0),colspan=5)
                ax3.set_xlabel("Genomic positions on chromosome: " + str(chrom)) 
            else:
                ax3 = pl.subplot2grid((6,5),(4,0),colspan=5)
            ax3.set_ylim(0,0.6)
            ax3.set_yticks([0.25,0.5])
            ax3.yaxis.grid()
            ax3.set_ylabel("MAF")
            ax3.set_xticks([])
            remove_border(top=True)
            ax3.set_xlim(_pos[ranges[0]]-snp_distance*0.01,_pos[ranges[-1]]+snp_distance*0.01)
            for key in vline_map:
                ax3.axvline(key,color=vline_map[key],linestyle='--',alpha=0.8,linewidth=0.5)
            for i in xrange(ranges.shape[0]):
                ax3.plot(_pos[ranges[i]],maf_list[i],'x',color=rr_color_list[i])

            if sqlite!=None:
                ax2 = pl.subplot2grid((6,5),(5,0),colspan=5)
                remove_border(top=True)
                ax2.set_xlabel("Genomic positions on chromosome: " + str(chrom)) 
                ax2.set_yticks([])
                ax2.set_xlim(_pos[ranges[0]]-snp_distance*0.01,_pos[ranges[-1]]+snp_distance*0.01)
                for key in vline_map:
                    ax2.axvline(key,color=vline_map[key],linestyle='--',alpha=0.8,linewidth=0.5)
        
                sqlite_cursor.execute("SELECT * FROM geneannotation WHERE chromosome_id=? AND annotation_end >=? AND annotation_start<=?",(str(chrom),int(_pos[ranges[0]]),int(_pos[ranges[-1]])))
                annotations = sqlite_cursor.fetchall()
                for g,annotation in enumerate(annotations):
                    alpha = 0.6
                    start = int(annotation[3]) 
                    length = int(annotation[4]) - int(annotation[3])
                    head_width=0.25
                    head_length=length*0.15
                    width=0.1
                    name = annotation[1]

                    if annotation[5]=="+":
                        shape = "right"
                        arrow_params={'length_includes_head':True, 'shape':shape,'head_starts_at_zero':True}
                        if g%2==0:
                            y_pos = 0.6
                        else:
                            y_pos = 0.5
                        y_text = y_pos+2*width
                        ax2.arrow(start,y_pos,length,0,head_width=head_width,head_length=head_length,fc=color_list[0],ec=color_list[0],alpha=alpha,width=width,**arrow_params)
                    else:
                        shape = "left"
                        arrow_params={'length_includes_head':True, 'shape':shape,'head_starts_at_zero':True}
                        if g%2==0:
                            y_pos = 0.15
                        else:
                            y_pos = 0.35
                        y_text = y_pos-1*width
                        ax2.arrow(start+length,y_pos,-length,0,head_width=head_width,head_length=head_length,fc=color_list[3],ec=color_list[3],alpha=alpha,width=width,**arrow_params)
                    ax2.text(start+length/2.0, y_text,name, size=12, ha='center', va='center', color="k")

            
            if not sqlite==None:
                cax = fig.add_axes([0.93, 0.4, 0.01, 0.5])
                pl.colorbar(cpick,label="SNP r^2",cax=cax)
                pl.subplots_adjust(left=0.07,bottom=0.15,right=0.92,top=0.9,wspace=0,hspace=0)
            else:
                pl.subplots_adjust(left=0.07,bottom=0.15,right=0.99,top=0.9,wspace=0,hspace=0)
            if not arguments.title:
                pl.subplots_adjust(top=0.9)

            pl.savefig(os.path.join(arguments.out,'ld_plot_chrom_' + str(chrom) + "_pos_" + str(pp) + "_" + fname +'.' + arguments.iformat) )
            pl.close()
