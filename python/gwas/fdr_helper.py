import scipy as sp
import scipy.interpolate as interpolate

def benjamini_hochberg(p_values=None,q_value=0.05,sort_idx=None,return_sort_idx=False):
    p_values = p_values.ravel()
    if sort_idx is None:
        sort_idx = sp.argsort(p_values)
        p_values = p_values[sort_idx]
    else:
        sort_idx = sort_idx.ravel()
        p_values = p_values[sort_idx]
    m = p_values.shape[0]
    idx_line = sp.arange(1,m+1)
    thr_line = (idx_line*q_value)/float(m);
    thr_ind = sp.where(p_values<=thr_line)[0]
    if thr_ind.shape[0]==0:
        thr = 0.0;
    else:
        thr = p_values[thr_ind.max()]
    #adjust p_values
    p_values_adjusted = sp.ones(m)
    prev = 1.0;
    for i in range(m,0,-1):
        p_values_adjusted[i-1] = sp.minimum(prev,float(m)/float(i)*p_values[i-1])
        if p_values_adjusted[i-1]>1:
            p_values_adjusted[i-1]=1
        prev = p_values_adjusted[i-1]
    #resort pvalues
    p_tmp = p_values_adjusted.copy()
    p_values_adjusted[sort_idx] = p_tmp
    if return_sort_idx==True:
        return [thr,p_values_adjusted,sort_idx]        
    else:
        return [thr,p_values_adjusted]
    
def benjamini_hochberg_yekutieli(p_values=None,q_value=0.05,sort_idx=None,return_sort_idx=False):
    p_values = p_values.ravel()
    if sort_idx is None:
        sort_idx = sp.argsort(p_values)
        p_values = p_values[sort_idx]
    else:
        sort_idx = sort_idx.ravel()
        p_values = p_values[sort_idx]
    m = p_values.shape[0]
    idx_line = sp.arange(1,m+1)
    cV = (1.0/idx_line).sum()
    thr_line = (idx_line*q_value*cV)/float(m);
    thr_ind = sp.where(p_values<=thr_line)[0]
    if thr_ind.shape[0]==0:
        thr = 0.0;
    else:
        thr = p_values[thr_ind.max()]
    #adjust p_values
    p_values_adjusted = sp.ones(m)
    prev = 1.0
    for i in range(m,0,-1):
        p_values_adjusted[i-1] = sp.minimum(prev,p_values[i-1]*float(m)*cV/float(i))
        if p_values_adjusted[i-1]>1:
            p_values_adjusted[i-1]=1
        prev = p_values_adjusted[i-1]
    #resort pvalues
    p_tmp = p_values_adjusted.copy()
    p_values_adjusted[sort_idx] = p_tmp
    if return_sort_idx==True:
        return [thr,p_values_adjusted,sort_idx]        
    else:
        return [thr,p_values_adjusted]
    
def storey_tibishirani(p_values=None,sort_idx=None,return_sort_idx=False):
    p_values = p_values.ravel()
    if sort_idx is None:
        sort_idx = sp.argsort(p_values)
        p_values = p_values[sort_idx]
    else:
        sort_idx = sort_idx.ravel()
        p_values = p_values[sort_idx]
    m = p_values.shape[0]
    if m<100: #if number if tests is too small use pi0=1
        pi0=1.0
    else: # otherwise estimate pi0 using a natural cubic spline
        #evaluate pi0 for a set of lambdas
        pi0 = []
        lambdas = sp.arange(0.01,0.96,0.01)
        counts = []
    
        for __lambda in lambdas:
            counts.append((p_values>__lambda).sum())
        counts = sp.array(counts)
        for i in xrange(lambdas.shape[0]):
            pi0.append(counts[i]/float(m*(1.0-lambdas[i])))
        pi0 = sp.array(pi0)
        
        splrep = interpolate.splrep(lambdas,pi0,k=3)
        pi0 = interpolate.splev(lambdas[-1],splrep)
    
        if pi0>1.0:
            pi0 = 1.0

    q_values = pi0*p_values;
    #q_values[-1] = sp.minimum(q_values[-1],1.0)
    for i in xrange(m-2,-1,-1):
        q_values[i] = sp.minimum(pi0*m*p_values[i]/float(i+1.0),q_values[i+1])
    #resort q_values
    q_tmp = q_values.copy()
    q_values[sort_idx] = q_tmp
    if return_sort_idx==True:
        return [q_values,sort_idx]        
    else:
        return q_values
