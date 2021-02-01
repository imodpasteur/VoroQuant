from Voro import *


v=Voro()

print('simulation of a sphere with 1000 localization of radius R=10')
print('theoretical surface = ',(4.*3.1415*10.*10.),'nm^2')
print('theoretical volume = ',((4./3.)*3.1415*10.*10.*10.),'nm^3')

v.simulate3Dsphere(1000)

v.addUniformNoiseN(300)#add noise with twice less data


v.makeVoroDiagram()



rank=1
v.computeVoroStatistic(rank)

factor=1#threshold at twice the average density

v.threshold(rank,label='density',factor=factor)




v.splitClusters()



v.plotcluster(2)





res=v.getFullStat(alpha=100,rank=rank,show=False)

surface=res[0]
volume=res[1]



print('surface=',surface,'nm^2')
print('volume=',volume,'nm^3')


print('CONCAV VOLUMES WITH ALPHA=100 AND ALPHA=500 SHOULD BE THE SAME FOR A SPHERE')
print('volume concav alpha=100:',res[13])
print('volume concav alpha=500:',res[21])
print('SMOOTHNESS (%) :',100*res[13]/res[21])
