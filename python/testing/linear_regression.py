import sys
sys.path.append("bin/" + sys.platform + "/interfaces/python/")
import CEasyGWAS as gwas
import scipy as sp

#create dummy genotype wiht 100 sampkes and 10 SNPs. 
#Just to demonstrate how to call some stuff in python
X = sp.random.randn(100,10)
Y = sp.random.randn(100)

linreg = gwas.LinearRegression()
linreg.setPhenotype(Y)
linreg.setGenotype(X)

linreg.test_associations()

print linreg.getPValues()
