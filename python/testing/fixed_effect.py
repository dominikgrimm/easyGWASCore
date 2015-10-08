import sys
print "bin/" + sys.platform + "/interfaces/python/"
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

fixed_model = gwas.FixedEffectModel(effects,variance)

fixed_model.process()

weights = fixed_model.getWeights()

print weights

mean = fixed_model.getWeightedMean()
Z = fixed_model.getZvalue()
p_value = fixed_model.getPvalue()
p_value2 = fixed_model.getPvalue(True)

print
print "Mean\tZ\tp_value\tp_value (two sided)"
print str(mean) + "\t" + str(Z) + "\t" + str(p_value) + "\t" + str(p_value2)
print

