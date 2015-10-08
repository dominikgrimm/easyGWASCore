import scipy as sp
import sys, os
import pylab as pl
from scipy.stats import shapiro
from scipy.stats.morestats import boxcox

color_t = ['#F7977A','#FDC68A','#A2D39C','#6ECFF6','#8493CA','#BC8DBF','#F6989D','#FFF79A','#998675','#A4A4A4']

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
    if right:
        ax.yaxis.tick_right()

def read_data(filename=None):
	f = open(filename,'r')
	phenotype_names = []
	sample_ids = []
	fid = []
	y = []
	
	for i,line in enumerate(f):
		sv = line.strip().split(" ")
		if i==0:
			for j in range(2,len(sv)):
				phenotype_names.append(sv[j])
		else:	
			fid.append(sv[0])
			sample_ids.append(sv[1])
			tmp = []
			for j in range(2,len(sv)):
				tmp.append(float(sv[j]))
			y.append(tmp)
	y = sp.array(y)
	phenotype_names = sp.array(phenotype_names,dtype="S100")
	sample_ids = sp.array(sample_ids)
	fid = sp.array(fid)
	f.close()
	return [y,phenotype_names,sample_ids,fid]

def plotHistogram(y=None,phenotype_name=None,transform='sqrt',outdir=None):
	ind = sp.where(~sp.isnan(y))
	y = y[ind]
	pl.figure(figsize=(12,6))
	pl.subplot(121)
	[test_statistic, p_value] = shapiro(y)
	pl.hist(y,bins=30,color=color_t[1],label="Original")
	pl.title(phenotype_name.replace("_"," ") + ", Shapiro: %.2e"%(p_value))
	leg = pl.legend(fancybox=True)
	leg.get_frame().set_alpha(0.2)
	#leg.get_frame().set_edgecolor("none")
	remove_border()
	
	pl.subplot(122)

	if transform=="all":
		p_vals = []
		transformations = sp.array(['boxcox','sqrt','log','log10'])
		for t in transformations:
			zeros = sp.where(y==0)[0]
			if zeros.shape[0]==0:
				if t=='sqrt':
					tmpy = sp.sqrt(y)
				elif t=="boxcox":
					[tmpy,b_lambda] = boxcox(y)
				elif t=="log":
					tmpy = sp.log(y)
				elif t=="log10":
					tmpy = sp.log10(y)
				[test_statistic, pv] = shapiro(tmpy)
			else:
				pv = 0.0
			p_vals.append(pv)
		p_vals = sp.array(p_vals)
		ind = sp.argmax(p_vals)
		transform = transformations[ind]
			
	ind = sp.where(y==0)[0]
	if ind.shape[0]>0:
		print "IMPORTANT: y contains 0 -> transformation changed to SQRT"
		transform = "sqrt"
	
	if transform=='sqrt':
		y = sp.sqrt(y)
	elif transform=="boxcox":
		[y,b_lambda] = boxcox(y)
	elif transform=="log":
		y = sp.log(y)
	elif transform=="log10":
		y = sp.log10(y)
	[test_statistic, p_value_t] = shapiro(y)
	pl.hist(y,bins=30,color=color_t[4],label=transform)
	pl.title(phenotype_name.replace("_"," ") + ", Shapiro: %.2e"%(p_value_t))
	leg = pl.legend(fancybox=True)
	leg.get_frame().set_alpha(0.2)
	#leg.get_frame().set_edgecolor("none")
	remove_border()
	pl.subplots_adjust(left=0.03,bottom=0.05,right=0.99,top=0.94,wspace=0.07,hspace=0.34)
	pl.savefig(os.path.join(outdir,phenotype_name + ".pdf"))
	if(p_value>p_value_t):
		return "original"
	else:
		return transform

if __name__ in "__main__":
	if len(sys.argv)<3:
		print
		print "Usage: phenotype_plots.py <PhenotypeFile> <OutputDirectory> [transform: default sqrt]"
		print
		quit(0)

	transform = 'all'
	if len(sys.argv)==4:
		transform = sys.argv[3]
	
	[y,phenotype_names,sample_ids,fid] = read_data(sys.argv[1])
	output_dir = sys.argv[2]

	for i,phenotype in enumerate(phenotype_names):
		selected_transform = plotHistogram(y=y[:,i],phenotype_name=phenotype,transform=transform,outdir=output_dir)
		ind = sp.where(~sp.isnan(y[:,i]))[0]
		if selected_transform=='sqrt':
			phenotype_names[i] = "sqrt_" + phenotype
			y[ind,i] = sp.sqrt(y[ind,i])
		elif selected_transform=="boxcox":
			phenotype_names[i] = "boxcox_" + phenotype 
			tmp = y[ind,i]
			[y[ind,i],b_lambda] = boxcox(tmp)
		elif selected_transform=="log":
			phenotype_names[i] = "log_" + phenotype
			y[ind,i] = sp.log(y[ind,i])
		elif selected_transform=="log10":
			phenotype_names[i] = "log10_" + phenotype
			y[ind,i] = sp.log10(y[ind,i])
	
	f = open(os.path.join(output_dir, "transformed_phenotypes.txt"),'w')
	f.write("FID IID ")
	string = ""
	for phenotype in phenotype_names:
		string += phenotype + " "
	f.write(string[:-1] + "\n")
	for i in range(fid.shape[0]):
		f.write(fid[i] + " " + sample_ids[i] + " ")
		string = ""
		for j,phenotype in enumerate(phenotype_names):
			string += str(y[i,j]) + " "
		f.write(string[:-1] + "\n")
	f.close()
