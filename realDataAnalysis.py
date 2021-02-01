
import sys
from Voro import *
import time


if (len(sys.argv)==2):
    print("python batch file: ",sys.argv[1])
    
    f=sys.argv[1]
    
    


    
    #copy from here to run on tars
    
    print(f)
    for i in range(1):



        v=Voro()
        v.loadTable(f)


        print('file loaded')

        print(len(v.data))

        if (len(v.data)>500):
            
            
            fcrop = open(f[:-3]+'crop', "r")
            scrop=fcrop.readline()
            vcrop=scrop.split(',')
            cropx=float(vcrop[0])
            cropy=float(vcrop[1])
            cropw=float(vcrop[2])
            croph=float(vcrop[3])
            print(cropx,cropy,cropw,croph)
            
            length=len(v.data)

            #v.showImage(100)

            v.makeVoroDiagram(100)

            v.removefromcells('adjacency')#to save memory: romove unuseful data





            rank=2
            v.computeVoroStatistic(rank)


            cropx+=v.miniBound[0]
            cropy+=v.miniBound[1]

            
            

            th=v.thresholdBackckgroundRatio(rank,cropx,cropy,cropw,croph,ratio=4)


            print("threshold ok")

            v.splitClusters()
            
            v.filterLargeClusters(.02)#remove small regions that may appear in the background (multiple localizations of the same emitter in a very small area)

            print("split cluster ok")

            #v.showCluster(100)

            #set 'show=True' to plot concave hull at alpha=100 and alpha=500
            res=v.getFullStat(alpha=100,rank=rank,z_weight=0.5,show=False)
            print(i,res)

            #fff.write("%d,%f\r\n" % (i,res))
            
            v.saveClusterTable(f+'_rank='+str(rank)+'_cluster1.csv')



            
            

            #save result:
            fff= open(f+'_rank='+str(rank)+'_cluster1.result',"w+")
            fff.write("%d,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f\r\n" % (i,res[0],res[1],res[2],res[3],res[4],res[5],res[6],res[7],length,th,res[8],res[9],res[10],res[11],res[12],res[13],res[14],res[15],res[16],res[17],res[18],res[19],res[20],res[21],res[22],res[23],res[24],res[25],res[26],res[27],res[28],res[29],res[30],res[31]))
            

            fff.close()
            
            
            
            
            #sum up:
            print('---------------------------------------------')
            print('gyration radius:',np.sqrt(res[3])/1000,'um')
            print('volume concav alpha=100:',res[13],'nm^3')
            print('volume concav alpha=500:',res[21],'nm^3')
            print('SMOOTHNESS (%) :',100*res[13]/res[21])
            print('---------------------------------------------')
            
            
            

            

