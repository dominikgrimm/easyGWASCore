import sys
sys.path.append("bin/" + sys.platform + "/interfaces/python/")
import CEasyGWAS as gwas
import scipy as sp

#create dummy genotype wiht 100 sampkes and 10 SNPs. 
#Just to demonstrate how to call some stuff in python
X = sp.random.randn(100,10)
Y = sp.random.randn(100)

K = gwas.CKernels.realizedRelationshipKernel(X)

lmm = gwas.EMMAX()
lmm.setPhenotype(Y)
lmm.setGenotype(X)
lmm.setK(K)

lmm.test_associations()

print lmm.getPValues()
