import sys
sys.path.append("bin/" + sys.platform + "/interfaces/python/")
import CEasyGWAS as gwas
import scipy as sp

means1 = sp.array([94, 98, 98, 94, 98, 96])
means2 = sp.array([92, 92, 88, 82, 88, 92])

sd1 = sp.array([22, 21, 28, 19, 21, 21])
sd2 = sp.array([20, 22, 26, 17, 22, 22])

n1 = sp.array([ 60,  65,  40, 200,  50,  85])
n2 = sp.array([ 60,  65,  40, 200,  45,  85])

mEffect = gwas.MeanEffectSize(means1,means2,sd1,sd2,n1,n2)
effects = mEffect.getHedgesG()
variance = mEffect.getVariance()

print effects, variance
tmp = sp.array([0.129297,0.193972])
tmp1 = sp.array([ 0.265015,0.341087])
random_model = gwas.RandomEffectModel(tmp,tmp1)

random_model.process()

weights = random_model.getWeights()

print weights

mean = random_model.getWeightedMean()
Z = random_model.getZvalue()
T2 = random_model.getT2()
Q = random_model.getQ()
C = random_model.getC()
p_value = random_model.getPvalue()
p_value2 = random_model.getPvalue(True)


print
print "Mean\tZ\tT2\tQ\tC\t\tp_value\tp_value (two sided)"
print str(mean) + "\t" + str(Z) + "\t" + str(T2) + "\t" + str(Q) + "\t" + str(C) + "\t" + str(p_value) + "\t" + str(p_value2)
print

print (Q-5)/C
