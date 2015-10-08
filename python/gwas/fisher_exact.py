import scipy.stats as stats
import scipy as sp

class FisherExact():
    
    def __init__(self):
        self.__y = None
        self.__X = None
        self.__p_values = None
        self.__scores = None

    def setPhenotype(self,y):
        self.__y = y

    def setGenotype(self,X):
        self.__X = X

    def getPvalues(self):
        return self.__p_values

    def getScores(self):
        return self.__scores

    def test_associations(self):
        self.__p_values = sp.ones(self.__X.shape[1])
        self.__scores = sp.ones(self.__X.shape[1])
        for i in xrange(self.__X.shape[1]):
            x = self.__X[:,i]
            ind = sp.where(x==2)[0]
            if ind.shape[0]>0:
                x[ind] = 1
            y = self.__y.flatten()
            table = sp.zeros((2,2))
            ind_cases = sp.where(y==1)[0]
            ind_controls = sp.where(y==0)[0]
            table[0,0] = ind_cases.shape[0]-x[ind_cases].sum()
            table[1,0] = x[ind_cases].sum()
            table[0,1] = ind_controls.shape[0]-x[ind_controls].sum()
            table[1,1] = x[ind_controls].sum()
            [odds_ratio,pval] = stats.fisher_exact(table)
            self.__p_values[i] = pval
            self.__scores[i] = odds_ratio
