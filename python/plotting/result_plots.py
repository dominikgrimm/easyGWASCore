import pylab as pl
import scipy as sp
import sys,os
import scipy.stats as stats

def read_result_files(fileName=None):
	f = open(sys.argv[1])
	pv = []
	positions = []
	chromosomes = []
	hashs = []
	for i,line in enumerate(f):
		if i==0: continue
		sv = line.strip().split("\t")
		chromosomes.append(sv[1])
		positions.append(int(sv[2]))
		pv.append(float(sv[4]))
		hashs.append(float(sv[-1].strip()))
	f.close()
	chromosomes = sp.array(chromosomes)
	positions = sp.array(positions)
	pv = sp.array(pv)
	hashs = sp.array(hashs)
	tmp,u_ind = sp.unique(hashs,return_index=True)
	unique_pv = pv[u_ind]
	return [pv,positions,chromosomes,hashs,unique_pv]

if __name__ in '__main__':
	
	if len(sys.argv)<3:
		print
		print "Usage: result_plots.py <result_file_name> <output_folder> [use unique snp true|false]"
		print
		quit(0);

	[pv,positions,chromosomes,hashs,unique_pv] = read_result_files(sys.argv[1])

	unsnps = False
	if len(sys.argv)==4:
		if sys.argv[3]=="true":
			unsnps = True
	
	pl.ion()
	pl.figure(figsize=(15,4))

	fname = sys.argv[1].split("/")[-1].split(".")[0]
	#fname = sys.argv[1].split("90.")[1].split(".")[0] + "." + sys.argv[1].split("90.")[1].split(".")[1]

	#chrom_list = ['1','2','3','4','5']
	chrom_list = sp.unique(chromosomes)
	color_list = ['#F26C4F','#F68E55','#7CC576','#00BFF3','#605CA8','#F06EA9']
	edge_list = ['#ED1C24','#F26522','#39B54A','#00AEEF','#2E3192','#EC008C']
	darker_list = ['#790000','#7B2E00','#005E20','#005B7F','#0D004C','#7B0046']
	
	if unsnps:
		bf_threshold = 0.05/unique_pv.shape[0]
	else:
		bf_threshold = 0.05/pv.shape[0]

	current_pos = 0

	max_y = (-sp.log10(pv[:])).max()
	if max_y<-sp.log10(bf_threshold):
		max_y = -sp.log10(bf_threshold)
	max_y += 1

	split_list = []
	xtick_list = []
	pl.axhline(-sp.log10(bf_threshold),color='r',linestyle='--')
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

		pl.plot(sp.arange(0+current_pos,current_pos+__pv.shape[0]),__pv,'.',color=color_list[i],alpha=0.8,
				markerfacecolor=color_list[i],
				markerfacecoloralt=color_list[i],
				markeredgecolor=edge_list[i],
				markeredgewidth=1,
				)
		
		idx = sp.where(_pv<=bf_threshold)[0]
		__pos = _pos[idx]
		__chrs = _chrs[idx]
		__pv = _pv[idx]
		__pv = -sp.log10(__pv)
		pl.plot(current_pos+idx,__pv,'.',color=darker_list[i],#alpha=0.6,
				#markerfacecolor=edge_list[i],
				#markerfacecoloralt=edge_list[i],
				#markeredgecolor="#D7D7D7",
				markeredgewidth=0,
				markersize=10,
				)


		xtick_list.append(float(current_pos) + float(_pv.shape[0])/2.0)
		current_pos += _pv.shape[0]+1
		split_list.append(current_pos-1)
		pl.xlim(0,current_pos)
		pl.ylabel("$-log10(p-value)$")

	pl.ylim(0,max_y)
	s_list = split_list[:-1]
	
	for split in s_list:
		pl.axvline(split,color='k',linestyle='--')
	
	pl.xticks(xtick_list,chrom_list)
	pl.grid(True)
	#compute gc
	if unsnps:
		pv = unique_pv
	gc = sp.median(stats.chi2.isf(pv,1))/0.456

	pl.title("Phenotype: %s, genomic control $\lambda=%0.2f$"%(fname,gc))
	pl.savefig(os.path.join(sys.argv[2],'manhattan_' + fname+'.png') )

	#print QQ-plot 
	pl.figure(figsize=(5,5))
	pv_uni = (sp.arange(1.0/float(pv.shape[0]),1,1.0/float(pv.shape[0]+1)))
	pl.plot(-sp.log10(pv_uni),-sp.log10(sp.sort(pv_uni)),'b--')
	pl.ylim(0,(-sp.log10(pv[:])).max()+1)
	pl.plot(-sp.log10(pv_uni),-sp.log10(sp.sort(pv[:],axis=0)),'r.',markeredgewidth=0,alpha=0.8)
	pl.grid(True)
	pl.title("Phenotype: %s"%(fname))
	pl.xlabel('Expected $-log10(p-value)$')
	pl.ylabel('Observed $-log10(p-value)$')
	pl.text(4,1,"$\lambda=%.2f$"%(gc))
	pl.savefig(os.path.join(sys.argv[2],'qqplot_' + fname+'.png') )
