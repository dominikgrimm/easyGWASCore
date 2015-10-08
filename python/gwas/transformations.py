import scipy as sp
from scipy.stats.morestats import boxcox

#zero mean data
def zeroMean(data=None,x=True):
    if x:
        return data-data.mean(axis=0)
    else:
        mean = data[~sp.isnan(data)].mean(axis=0)
        return data - mean

#zero mean and unit variance
def zeroMeanUnitVarianz(data=None,x=True):
    if x:
        return (data-data.mean(axis=0))/data.std(axis=0)
    else:
        mean = data[~sp.isnan(data)].mean(axis=0)
        std = data[~sp.isnan(data)].std(axis=0)
        return (data - mean)/std

#log10 transformation
def logTransformation(data=None):
    return sp.log10(data)
        
#SQRT transformation
def sqrtTransformation(data=None):
    return sp.sqrt(data)
        
def boxcoxTransformation(data=None):
    ind = sp.where(data.flatten()<=0.0)[0]
    if ind.shape[0]>0:
        return data
    if data.ndim==1:
        [data, boxcox_lambda] = boxcox(data)           
    elif data.shape[1]==1:
        [tmp_y, boxcox_lambda] = boxcox(data[:,0])
        data = sp.zeros((tmp_y.shape[0],1))
        data[:,0] = tmp_y
    elif data.shape[0]==1:
        [tmp_y, boxcox_lambda] = boxcox(data[0,:])
        data = sp.zeros((1,tmp_y.shape[0]))
        data[0,:] = tmp_y
    return data
    
def createDummyVariables(data=None):
    data = data.ravel()
    variables_list = sp.unique(data)
    n_variables = variables_list.shape[0]-1
    dummy_matrix = sp.zeros([data.shape[0],n_variables])
    for i in xrange(n_variables):
        ind = sp.where(variables_list[i+1]==data)[0]
        dummy_matrix[ind,i] = 1
    return dummy_matrix
