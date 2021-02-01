
import sys
from Voro import *
import os
import pickle


times=1


theLength=49800


increasedSampling=3#5#density

kilobasePerBead=5



#PLEASE, run makeRandomCuts.py first
print("PLEASE, run makeRandomCuts.py first")
            










if (len(sys.argv)==2):

    
    print("python batch file: ",sys.argv[1])
    
    f=sys.argv[1]
    
    
    thefilename=os.path.basename(f)
    p1=thefilename.find('_timepoint')
    if (p1<=6):
        print("ERROR, coord timepoint not found in filename ",f)
    thenumber=thefilename[6:p1]

    with open('randomCutting/Coord_'+thenumber+'_timepoint', 'rb') as rc:
        chrRandomCut = pickle.load(rc)
    
    
    start=chrRandomCut[0]
    end=chrRandomCut[1]
    startKb=chrRandomCut[2]
    endKb=chrRandomCut[3]


            
    print("start/end",start,end)



        
        

    split=len(start)#split number
    
    for i in range(split):
        
        
        #copy from here to run on tars
        
        print(f)
        
        beadSize=45#nm 45 --> 5.5um nucleus size
        #dataNumber=int(splitposit*1000)
        
        dataNumber=(end[i]-start[i])*increasedSampling
        dataNumberKb=endKb[i]-startKb[i]
        
        
        
        v=Voro()
        v.load3DMatrix(f,' ')
        
        

        
        print('file loaded')
        
        
        number=len(v.data)
        
        print('!!!!!!!!!!!!!!!!  nanometer conversion  bead=50nm !!!!!!!!!!!!!!!!!!!')
        for numb in range(number):
            v.data[numb][0]*=beadSize
            v.data[numb][1]*=beadSize
            v.data[numb][2]*=beadSize
            
        #print("SPLIT ",i*(number/split),min((i+1)*(number/split),number),dataNumber)
        #v.selectData(i*(number/split),min((i+1)*(number/split),number),dataNumber)
        '''if i==0:
            v.selectData(0,int(float(number)*splitposit/100),int(splitposit*1000))
        if i==1:
            v.selectData(int(float(number)*splitposit/100),len(v.data),int((100-splitposit)*1000))
        '''
        
        print('START',start[i],'END',end[i],dataNumber)
        v.selectData(start[i],end[i],dataNumber)
        print('startKb',startKb[i],'endKb',endKb[i],dataNumberKb)
        
            
        #if i==1:
        #    v.selectData(int(float(number)*splitposit/100),len(v.data),int((100-splitposit)*1000))
            
        
        
        finalsize=dataNumber
        v.saveTable(f+'_'+str(split)+'split='+str(i)+'_GT.csv')
        
        
        
        #number=v.selectData(c[0],c[1],length)
        number=len(v.data)
        v.addLocalizationError(65,65,130)
        
        
        ###############################added...
        mini=np.min(v.data,0)
        maxi=np.max(v.data,0)
        
        diff=np.mean(np.subtract(maxi,mini))
        
        gap=1500 ############################### gap=1.5 um
        density=0.00000015#0.00000015
        v.addUniformNoiseD(density,gap,50)#add non specific noise
        #v.addUniformNoiseN(number,gap)#add non specific noise
        
        cropx=mini[0]-gap+gap*.05
        cropy=mini[1]-gap+gap*.05
        cropw=gap*.5
        croph=maxi[1]-mini[1]+gap*2*.9#the full left part
        
        ###############################
        
        v.saveTable(f+'_'+str(split)+'times='+str(times)+'split='+str(i)+'_GT_wnoise.csv')
        
        print("newlength",len(v.data))

        if (len(v.data)<1000000):####################################
            #v.showImage(100)
            length=len(v.data)
            v.makeVoroDiagram(100)

            v.removefromcells('adjacency')#to save memory: remove unuseful data


            


            rank=2
            v.computeVoroStatistic(rank)

           


            threshold=0.0000001 #  = 1000 emitters in 1000.000.000nm^3
            


            #th=v.thresholdFixed(rank,label='density',threshold=threshold)
            th=v.thresholdBackckgroundRatio(rank,cropx,cropy,cropw,croph,ratio=4)

            print('WARNING:  THRESHOLD=',th)

            print("threshold ok")

            v.splitClusters()
            
            v.filterLargeClusters(.01)

            print("split cluster ok")

            #v.showCluster(100)

            res=v.getFullStat(alpha=100,rank=rank,z_weight=0.5)

            print(res)
            fff= open(f+'_'+str(split)+'times='+str(times)+'split='+str(i)+'_rank='+str(rank)+'_cluster1.result',"w+")
            #fff.write("%d,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\r\n" % (i,res[0],res[1],res[2],res[3],res[4],res[5],res[6],res[7],length,th,res[8],finalsize,dataNumberKb))
            fff.write("%d,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f,%.20f\r\n" % (i,res[0],res[1],res[2],res[3],res[4],res[5],res[6],res[7],length,th,res[8],res[9],res[10],res[11],res[12],res[13],res[14],res[15],res[16],res[17],res[18],res[19],res[20],res[21],res[22],res[23],res[24],res[25],res[26],res[27],res[28],res[29],res[30],res[31],finalsize,dataNumberKb))
            fff.close()
            fff.close()
            
            
            
            v.saveClusterTable(f+'_'+str(split)+'times='+str(times)+'split='+str(i)+'_rank='+str(rank)+'_cluster1.csv')
            
            
            
            
            #sum up:
            print('---------------------------------------------')
            print('SCE split number:',i)
            print('gyration radius:',np.sqrt(res[3])/1000,'um')
            print('volume concav alpha=100:',res[13],'nm^3')
            print('volume concav alpha=500:',res[21],'nm^3')
            print('SMOOTHNESS (%) :',100*res[13]/res[21])
            print('---------------------------------------------')
            
            
