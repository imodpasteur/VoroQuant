
# coding: utf-8

# In[1]:


import sys
import resource
resource.setrlimit(resource.RLIMIT_STACK, [0x10000000, resource.RLIM_INFINITY])
sys.setrecursionlimit(0x10000000)

import re
import numpy as np
import pandas as pd
import pyvoro
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from shapely.geometry import Polygon


from collections import defaultdict


import random
#random.seed(1)
import copy
from scipy.spatial import ConvexHull

from scipy.spatial import Delaunay

from skimage import measure

from scipy.interpolate import interp1d

from bisect import bisect_left#search closest higher value

from collections import namedtuple




#%matplotlib inline

class Voro():
    
    is3D=True
    data=None
    cells=None
    volume=0
    density=0
    miniBound=None
    maxiBound=None
    
    clusterLabels=[]
    
    def __init__(self,is3D=True,data=None):
        self.is3D=is3D
        self.data=None
        self.cells=None
        self.volume=0
        self.density=0
        self.miniBound=None
        self.maxiBound=None
        self.clusterLabels=[]
        if data!=None:
            self.data=np.array(data)
        
    def simulate3Dsphere(self,number):
        gap=250.;
        radius=500.;
        self.data=[]
        for i in range(0,number):
            rad=radius*(pow(random.random(), 1.0/3.0))
            a=random.random()*3.141592*2.;
            b=random.random()*3.141592*2.;
            x=rad*np.sin(a)*np.cos(b)
            y=rad*np.sin(a)*np.sin(b)
            z=rad*np.cos(a)
            self.data.append([x+gap/2.,y+gap/2.,z+gap/2.])
        
        self.data=np.array(self.data)
        self.is3D=True 
        return radius
        
    
    def loadTable(self,path= "/home/benoit/data/chromatin/tmp2.csv",x_id="x",y_id="y",z_id="z"):
        
        delimiter=','
        
        
        #read identifiers x, y and z
        x_posit=None
        y_posit=None
        z_posit=None
        x_tmp=None
        y_tmp=None
        z_tmp=None
        with open(path) as f:
            first_line = f.readline()
            first_line_list = first_line.split(delimiter)
            for i in range(len(first_line_list)):
                identifier=first_line_list[i]
                identifier=re.sub(r'\W+', '', identifier)
                if identifier.startswith(x_id):
                    x_tmp=i
                if identifier.startswith(y_id):
                    y_tmp=i
                if identifier.startswith(z_id):
                    z_tmp=i
                if identifier==x_id:
                    x_posit=i
                if identifier==y_id:
                    y_posit=i
                if identifier==z_id:
                    z_posit=i
        if x_posit==None:
            x_posit=x_tmp
        if y_posit==None:
            y_posit=y_tmp
        if z_posit==None:
            z_posit=z_tmp
        
        
        if x_posit==None:
            raise NameError('identifier x not found')
        if y_posit==None:
            raise NameError('identifier y not found')
        if self.is3D&z_posit==None:
            raise NameError('identifier z not found')
        
        
                
        
        
        #if (self.is3D):
        #    self.data=(np.loadtxt(path, delimiter=delimiter, skiprows=1,usecols = (x_posit,y_posit,z_posit)))
        #else :
        #    self.data=(np.loadtxt(path, delimiter=delimiter, skiprows=1,usecols = (x_posit,y_posit)))
            
        if (self.is3D):
            self.data=(pd.read_csv(path, delimiter, skiprows=1,usecols = (x_posit,y_posit,z_posit))).values
        else :
            self.data=(pd.read_csv(path, delimiter, skiprows=1,usecols = (x_posit,y_posit))).values
        
        
        
    def load3DMatrix(self,path,delimiter=','):
        
        
        x_posit=0
        y_posit=1
        z_posit=2
        
        
        
        if (self.is3D):
            self.data=(np.loadtxt(path, delimiter=delimiter, skiprows=0,usecols = (x_posit,y_posit,z_posit)))
        else :
            self.data=(np.loadtxt(path, delimiter=delimiter, skiprows=0,usecols = (x_posit,y_posit)))
    
        
    
    #add points on the surface of a cube around the structure to segment
    #each axis is separated in "step" points
    #Then, a fixed threshold with th=0 allows to remove them and keep the structure
    #This is useful since thresholding is mandatory using voronoi diagram
    #doing addContourMesh + thresholding is the same as doing NO THRESHOLD
    def addContourMesh(self,step=5,gap=-1):
        
        mini=np.min(self.data,0)
        maxi=np.max(self.data,0)
        if gap==-1:
            gap=np.max(np.subtract(maxi,mini))/float(step)#gap equals max size
            gap*=2
        mini = [x - gap for x in mini]
        maxi = [x + gap for x in maxi]

        
        for x in range(1,step):
            for y in range(1,step):
                zz=mini[2]
                xx=(float(x)/float(step))*(maxi[0]-mini[0])+mini[0]
                yy=(float(y)/float(step))*(maxi[1]-mini[1])+mini[1]
                self.data=np.append(self.data,[[xx,yy,zz]],axis=0)
                zz=maxi[2]
                xx=(float(x)/float(step))*(maxi[0]-mini[0])+mini[0]
                yy=(float(y)/float(step))*(maxi[1]-mini[1])+mini[1]
                self.data=np.append(self.data,[[xx,yy,zz]],axis=0)
                xx=mini[0]
                yy=(float(x)/float(step))*(maxi[1]-mini[1])+mini[1]
                zz=(float(y)/float(step))*(maxi[2]-mini[2])+mini[2]
                self.data=np.append(self.data,[[xx,yy,zz]],axis=0)
                xx=maxi[0]
                yy=(float(x)/float(step))*(maxi[1]-mini[1])+mini[1]
                zz=(float(y)/float(step))*(maxi[2]-mini[2])+mini[2]
                self.data=np.append(self.data,[[xx,yy,zz]],axis=0)
                yy=mini[1]
                xx=(float(x)/float(step))*(maxi[0]-mini[0])+mini[0]
                zz=(float(y)/float(step))*(maxi[2]-mini[2])+mini[2]
                self.data=np.append(self.data,[[xx,yy,zz]],axis=0)
                yy=maxi[1]
                xx=(float(x)/float(step))*(maxi[0]-mini[0])+mini[0]
                zz=(float(y)/float(step))*(maxi[2]-mini[2])+mini[2]
                self.data=np.append(self.data,[[xx,yy,zz]],axis=0)
        mini=np.min(self.data,0)
        maxi=np.max(self.data,0)
    
    #add a given number of localizations in the background
    def addUniformNoiseN(self,number,gap=5):
        
        mini=np.min(self.data,0)
        maxi=np.max(self.data,0)
        mini = [x - gap for x in mini]
        maxi = [x + gap for x in maxi]
        
        for i in range(number):
            xn=random.random()*(maxi[0]-mini[0])+mini[0]
            yn=random.random()*(maxi[1]-mini[1])+mini[1]
            zn=random.random()*(maxi[2]-mini[2])+mini[2]
            self.data=np.append(self.data,[[xn,yn,zn]],axis=0)
            
    #add a given density of localizations in the background
    def addUniformNoiseD(self,density,gapLeft=500,gap=50):
        
        
        mini=np.min(self.data,0)
        maxi=np.max(self.data,0)
        
        mini[0]-=gapLeft
        mini[1]-=gap
        mini[2]-=gap
        maxi[0]+=gap
        maxi[1]+=gap
        maxi[2]+=gap
        
        number=int((maxi[0]-mini[0])*(maxi[1]-mini[1])*(maxi[2]-mini[2])*density)
        print('background number: ',number)
        for i in range(number):
            xn=random.random()*(maxi[0]-mini[0])+mini[0]
            yn=random.random()*(maxi[1]-mini[1])+mini[1]
            zn=random.random()*(maxi[2]-mini[2])+mini[2]
            self.data=np.append(self.data,[[xn,yn,zn]],axis=0)
        
    
    
        
    
    def selectData(self,minimumPosition,maximumPosition,length):
        
        maximumPosition=min(maximumPosition,len(self.data)-1)
        minimumPosition=float(minimumPosition)
        maximumPosition=float(maximumPosition)
        x=range(len(self.data))
        
        step=(maximumPosition-minimumPosition)/length
        
        newRange=map(lambda xx: (xx+(minimumPosition/step-int(minimumPosition/step)))*step, range(int(minimumPosition/step), int(maximumPosition/step), 1))
        
        fx = interp1d(x, self.data[:,0], kind='cubic')
        
        fy = interp1d(x, self.data[:,1], kind='cubic')
        
        fz = interp1d(x, self.data[:,2], kind='cubic')
        
        self.data=np.concatenate((np.concatenate((np.array([fx(newRange)]).T, np.array([fy(newRange)]).T), axis=1), np.array([fz(newRange)]).T), axis=1)
        
        
        return(length)
    
    
    
    
    def addLocalizationError(self,stdX,stdY,stdZ):
        
        self.data[:,0]=np.add(self.data[:,0],np.random.normal(0,stdX,len(self.data[:,0])))
        self.data[:,1]=np.add(self.data[:,1],np.random.normal(0,stdY,len(self.data[:,1])))
        self.data[:,2]=np.add(self.data[:,2],np.random.normal(0,stdZ,len(self.data[:,2])))
        
    def addLocalizationErrorFixed(self,shift):
        
        for i in range(len(self.data)):
            
            a=random.random()*3.141592*2.;
            b=random.random()*3.141592*2.;
            x=shift*np.sin(a)*np.cos(b)
            y=shift*np.sin(a)*np.sin(b)
            z=shift*np.cos(a)
            self.data[i][0]=self.data[i][0]+x
            self.data[i][1]=self.data[i][1]+y
            self.data[i][2]=self.data[i][2]+z
    
    
    
    def removefromcells(self,variableToRemove):
        for x in self.cells:
            del(x[variableToRemove])
    
    
    
    
    
    ################################################## FUNCTIONS UTILS ###############
    
    #determinant of matrix a
    def det(self,a):
        return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0] - a[0][1]*a[1][0]*a[2][2] - a[0][0]*a[1][2]*a[2][1]

    #unit normal vector of plane defined by points a, b, and c
    def unit_normal(self,a, b, c):
        x = self.det([[1,a[1],a[2]],
                 [1,b[1],b[2]],
                 [1,c[1],c[2]]])
        y = self.det([[a[0],1,a[2]],
                 [b[0],1,b[2]],
                 [c[0],1,c[2]]])
        z = self.det([[a[0],a[1],1],
                 [b[0],b[1],1],
                 [c[0],c[1],1]])
        magnitude = (x**2 + y**2 + z**2)**.5
        return (x/magnitude, y/magnitude, z/magnitude)

    #dot product of vectors a and b
    def dot(self,a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

    #cross product of vectors a and b
    def cross(self,a, b):
        x = a[1] * b[2] - a[2] * b[1]
        y = a[2] * b[0] - a[0] * b[2]
        z = a[0] * b[1] - a[1] * b[0]
        return (x, y, z)

    #area of polygon poly
    def area(self,poly):
        if len(poly) < 3: # not a plane - no area
            return 0

        total = [0, 0, 0]
        for i in range(len(poly)):
            vi1 = poly[i]
            if i is len(poly)-1:
                vi2 = poly[0]
            else:
                vi2 = poly[i+1]
            prod = self.cross(vi1, vi2)
            total[0] += prod[0]
            total[1] += prod[1]
            total[2] += prod[2]
        result = self.dot(total, self.unit_normal(poly[0], poly[1], poly[2]))
        return abs(result/2)


    Point = namedtuple('Point', ['x', 'y'])

    def planEquation(points):


        if len(points)<3:
            print("error: number a points not enough")


        p1 = points[0]
        p2 = points[1]
        p3 = points[2]

        # These two vectors are in the plane
        v1 = np.subtract(p3,p1)
        v2 = np.subtract(p2,p1)

        # the cross product is a vector normal to the plane
        cp = np.cross(v1, v2)
        a, b, c = cp

        # This evaluates a * x3 + b * y3 + c * z3 which equals d
        d = -np.dot(cp, p3)

        print('The equation is {0}x + {1}y + {2}z + {3} = 0'.format(a, b, c, d))
        #print('a, b, c, d : ',a, b, c, d)
        #print('points=',points)
        return([a,b,c,d])
    
    
    #decompose polygon in triangle (2D)
    def earclip(polygon):
        """
        Simple earclipping algorithm for a given polygon p.
        polygon is expected to be an array of 2-tuples of the cartesian points of the polygon
        For a polygon with n points it will return n-2 triangles.
        The triangles are returned as an array of 3-tuples where each item in the tuple is a 2-tuple of the cartesian point.
        e.g
        >>> polygon = [(0,1), (-1, 0), (0, -1), (1, 0)]
        >>> triangles = tripy.earclip(polygon)
        >>> triangles
        [((1, 0), (0, 1), (-1, 0)), ((1, 0), (-1, 0), (0, -1))]
        Implementation Reference:
            - https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf
        """
        ear_vertex = []
        triangles = []

        polygon = [Point(*point) for point in polygon]

        if _is_clockwise(polygon):
            polygon.reverse()

        point_count = len(polygon)
        for i in range(point_count):
            prev_index = i - 1
            prev_point = polygon[prev_index]
            point = polygon[i]
            next_index = (i + 1) % point_count
            next_point = polygon[next_index]

            if _is_ear(prev_point, point, next_point, polygon):
                ear_vertex.append(point)

        while ear_vertex and point_count >= 3:
            ear = ear_vertex.pop(0)
            i = polygon.index(ear)
            prev_index = i - 1
            prev_point = polygon[prev_index]
            next_index = (i + 1) % point_count
            next_point = polygon[next_index]

            polygon.remove(ear)
            point_count -= 1
            triangles.append(((prev_point.x, prev_point.y), (ear.x, ear.y), (next_point.x, next_point.y)))
            if point_count > 3:
                prev_prev_point = polygon[prev_index - 1]
                next_next_index = (i + 1) % point_count
                next_next_point = polygon[next_next_index]

                groups = [
                    (prev_prev_point, prev_point, next_point, polygon),
                    (prev_point, next_point, next_next_point, polygon)
                ]
                for group in groups:
                    p = group[1]
                    if _is_ear(*group):
                        if p not in ear_vertex:
                            ear_vertex.append(p)
                    elif p in ear_vertex:
                        ear_vertex.remove(p)
        return triangles


    def _is_clockwise(polygon):
        s = 0
        polygon_count = len(polygon)
        for i in range(polygon_count):
            point = polygon[i]
            point2 = polygon[(i + 1) % polygon_count]
            s += (point2.x - point.x) * (point2.y + point.y)
        return s > 0


    def _is_convex(prev, point, next):
        return _triangle_sum(prev.x, prev.y, point.x, point.y, next.x, next.y) < 0


    def _is_ear(p1, p2, p3, polygon):
        ear = _contains_no_points(p1, p2, p3, polygon) and                 _is_convex(p1, p2, p3) and                 _triangle_area(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y) > 0
        return ear


    def _contains_no_points(p1, p2, p3, polygon):
        for pn in polygon:
            if pn in (p1, p2, p3):
                continue
            elif _is_point_inside(pn, p1, p2, p3):
                return False
        return True


    def _is_point_inside(p, a, b, c):
        area = _triangle_area(a.x, a.y, b.x, b.y, c.x, c.y)
        area1 = _triangle_area(p.x, p.y, b.x, b.y, c.x, c.y)
        area2 = _triangle_area(p.x, p.y, a.x, a.y, c.x, c.y)
        area3 = _triangle_area(p.x, p.y, a.x, a.y, b.x, b.y)
        return area == sum([area1, area2, area3])


    def _triangle_area(x1, y1, x2, y2, x3, y3):
        return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)


    def _triangle_sum(x1, y1, x2, y2, x3, y3):
        return x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1)
    
    
    ################################################################
    
    
    
    
    #function to compute voronoi diagram
    #gap: gap around cloud of points
    def makeVoroDiagram(self,gap=5):
        
        
        mini=np.min(self.data,0)
        maxi=np.max(self.data,0)
        
        self.miniBound=mini
        self.maxiBound=maxi
        
        mini = [x - gap for x in mini]
        maxi = [x + gap for x in maxi]
        
        
        
        
        #blocksize : high performance with ~ 5 point/block
        blockNumber=np.size(self.data,0)/5
        
        if self.is3D:
            blocksize=(maxi[2]-mini[2])*(maxi[1]-mini[1])*(maxi[0]-mini[0])/blockNumber
            blocksize**=(1./3.)
        else:
            blocksize=(maxi[1]-mini[1])*(maxi[0]-mini[0])/blockNumber
            blocksize**=(1./3.)
        
        if True:
            if self.is3D:
                self.cells = pyvoro.compute_voronoi(self.data,zip(mini,maxi), blocksize)
            else:
                self.cells = pyvoro.compute_2d_voronoi(self.data,zip(mini,maxi), blocksize)
        
        print("voronoi diagram computed")
        
        self.data=None
        
        self.volume=float(np.prod(np.subtract(maxi,mini)))
        self.density=(float(len(self.cells)))/float(self.volume)
        
    
    
    
    
    
    
    def computeVoroStatistic(self,rank=1):
        
        #get neighborhood for each rank
        for r in range(rank+1):
            for i, cell in enumerate(self.cells):
                neighbor=[]

                if r==0:#itself
                    neighbor.append(i)
                    d ={"rank": r,  "neighbors": set(neighbor)}
                    liste=[d]
                elif r==1:#first rank neighbors
                    liste=cell['stat']
                    for adj in cell['faces']:
                        neighbor.append(adj['adjacent_cell'])
                    d ={  "rank": r,  "neighbors": set(neighbor)}
                    liste.append(d)
                else:
                    liste=cell['stat']
                    fullRankSet=set(cell['stat'][0]['neighbors'])
                    for rr in range(1,r):
                        fullRankSet.update(cell['stat'][rr]['neighbors'])
                    oneRankList=liste[1]['neighbors']#rank-1 list


                    newSet=set()
                    for orl in oneRankList:#create rank list that is not present in fullRankList

                        newSet.update(self.cells[orl]['stat'][r-1]['neighbors'])

                    d ={"rank": r,  "neighbors": (newSet-fullRankSet)}

                    liste.append(d)

                cell.update({'stat':liste})

        
        #Now compute local stat for each rank/neighborhood
        
        for i, cell in enumerate(self.cells):
            distance=0
            volume=0
            ptNumber=0
            for stat in cell['stat']:#for each rank
                ptNumber+=len(stat['neighbors'])
                for s in stat['neighbors']:
                    #current cell: cell
                    #neighboring cell: self.cells[s]
                    #distance+=np.sqrt(((cell['original'][0]- self.cells[s]['original'][0])**2)+((cell['original'][1]- self.cells[s]['original'][1])**2)+((cell['original'][2]- self.cells[s]['original'][2])**2))
                    #distance+=np.sqrt(np.sum(np.subtract(cell['original'], self.cells[s]['original'])**2))#slower
                    volume+=self.cells[s]['volume']
                stat.update({'neighbor_number':ptNumber})
                #stat.update({'mean_distance':float(distance)/float(ptNumber)})
                stat.update({'volume':volume})
                stat.update({'density':float(ptNumber)/float(volume)})
                
        print("statistic computed")
    
    
    
    
    
    
    def thresholdFixed(self,rank,label='density',threshold=0.0001):
        
        for i,cell in enumerate(self.cells):
            if (cell['stat'][rank]['density']<threshold)|(min(cell['stat'][rank]['neighbors'])<0):
                cell.update({'cluster':{'label':0}})
                del(cell['vertices'])#to save memory
                del(cell['faces'])#to save memory
            else:
                cell.update({'cluster':{'label':1}})
        return(threshold)
    
    
    
    
    #label can be any stat for a given rank: mean_distance ; volume ; density
    #set threshold to add a strait line
    def threshold(self,rank,label='density',factor=2):
        
        threshold=factor*self.density
        for i,cell in enumerate(self.cells):
            if (cell['stat'][rank]['density']<threshold)|(min(cell['stat'][rank]['neighbors'])<0):
                cell.update({'cluster':{'label':0}})
                del(cell['vertices'])#to save memory
                del(cell['faces'])#to save memory
            else:
                cell.update({'cluster':{'label':1}})
        return(threshold)
    
    
    
        
        
        
        
    
    
    
    def thresholdBackckgroundRatio(self,rank,cropx,cropy,cropw,croph,ratio=4):
        #search minX and minY
        label='density'
        #compute mean
        meanDensityCrop=0.
        maxDensityCrop=0.
        count=0.
        for i,cell in enumerate(self.cells):
            co=cell['original']
            
            
            if ((((co[0])>cropx)&((co[0])<(cropx+cropw))&((co[1])>cropy)&((co[1])<(cropy+croph)))):
                if (min(cell['stat'][rank]['neighbors'])>0):#do not include edges
                    meanDensityCrop+=cell['stat'][rank]['density']
                    if (cell['stat'][rank]['density']>maxDensityCrop):
                        maxDensityCrop=cell['stat'][rank]['density']
                    count+=1.
        meanDensityCrop/=count
        stdDensityCrop=0.
        #compute std
        for i,cell in enumerate(self.cells):
            co=cell['original']
            if ((((co[0])>cropx)&((co[0])<(cropx+cropw))&((co[1])>cropy)&((co[1])<(cropy+croph)))):
                if (min(cell['stat'][rank]['neighbors'])>0):#do not include edges
                    stdDensityCrop+=(cell['stat'][rank]['density']-meanDensityCrop)*(cell['stat'][rank]['density']-meanDensityCrop)
        stdDensityCrop/=count
        stdDensityCrop=np.sqrt(stdDensityCrop)
        #threshold
        threshold=meanDensityCrop*ratio#*stdDensityCrop
        print('mean=',meanDensityCrop,'+-',stdDensityCrop,'  max=',maxDensityCrop,' -->  threshold=',threshold)
        #thresholding
        for i,cell in enumerate(self.cells):
            if (cell['stat'][rank]['density']<threshold)|(min(cell['stat'][rank]['neighbors'])<0):
                cell.update({'cluster':{'label':0}})
            else:
                cell.update({'cluster':{'label':1}})
        return(threshold)
    
    
    
    
    
        
        
        
        
        
    
    
    #recursive function called by separateClusters
    def recurrentClusterAssignment(self,cellId,backgroundlabel,segmentedlabel,newlabel,processed):
        self.cells[cellId]['cluster']['label']=newlabel
        self.cells[cellId]['cluster'].update({'isContour':False})
        processed[cellId]=True
        for adj in self.cells[cellId]['faces']:
            if processed[adj['adjacent_cell']]==False:
                if self.cells[adj['adjacent_cell']]['cluster']['label']==segmentedlabel:
                    self.recurrentClusterAssignment(adj['adjacent_cell'],backgroundlabel,segmentedlabel,newlabel,processed)
                elif self.cells[adj['adjacent_cell']]['cluster']['label']==backgroundlabel:#if the neighbor is 0 -> current is a contour
                    self.cells[cellId]['cluster'].update({'isContour':True})
                    if 'backgroundCell' in self.cells[cellId]['cluster']:
                        ll=self.cells[cellId]['cluster']['backgroundCell']
                        ll.append(adj['adjacent_cell'])
                        self.cells[cellId]['cluster'].update({'backgroundCell':ll})
                    else:
                        self.cells[cellId]['cluster'].update({'backgroundCell':[adj['adjacent_cell']]})
                else:
                    print("ERROR: other label connex to the main one")
    
    
    
    #because multiple clusters can be present, we need the following process to split them:
    #search connex component
    def splitClusters(self):
        
        processed=np.full((len(self.cells)), False, dtype=bool)
        labelNumber=2
        for i, cell in enumerate(self.cells):
            if processed[i]==False:
                if cell['cluster']['label']==1:
                    self.recurrentClusterAssignment(i,0,1,labelNumber,processed)
                    self.clusterLabels.append(labelNumber)
                    labelNumber+=1
                    
        
        
    def filterLargeClusters(self,ratioOfLargest):
        
        #search largest cluster:
        number=np.zeros(len(self.clusterLabels))
        for i,cell in enumerate(self.cells):
            if cell['cluster']['label']>1:
                number[cell['cluster']['label']-2]+=1
        largest=0 #start at 2 because do not include background
        for i in range(0,len(self.clusterLabels)):
            if number[i]>number[largest]:
                largest=i
        #filtering:
        print('cluster loc number',number)
        for i,cell in enumerate(self.cells):
            if cell['cluster']['label']>1:
                n=number[cell['cluster']['label']-2]
                if (n<(ratioOfLargest*number[largest])):
                    cell['cluster']['label']=0
    
        
    def showImage(self,pixSize=50):
        matrix=self.data
        
        
        
        mini=np.min(matrix,axis=0)
        maxi=np.max(matrix,axis=0)
        bins=np.ceil(np.divide(np.subtract(maxi,mini),pixSize))
        
        
        bins=[int(i) for i in bins.tolist()]
        
        hist=np.zeros(bins[0:2])
        
        
        for j,i in enumerate(matrix):
            x=np.floor(np.divide(np.subtract(i,mini),pixSize))
            x=[int(i) for i in x.tolist()]
            
            hist[x[0]][x[1]]+=1
                
            
                
            
            
            
        fig=plt.figure()
        plt.imshow(hist)  
        
        
        
        
        
    def showCluster(self,pixSize=50):
        matrix=[x['original'].tolist() for x in self.cells]
        
        label=[x['cluster']['label'] for x in self.cells]
        contourmat=[[i,int(x['cluster']['isContour'])] for i,x in enumerate(self.cells) if x['cluster']['label']>0]
        
        
        mini=np.min(matrix,axis=0)
        maxi=np.max(matrix,axis=0)
        bins=np.ceil(np.divide(np.subtract(maxi,mini),pixSize))
        
        
        bins=[int(i) for i in bins.tolist()]
        
        hist=np.zeros(bins)
        histinside=np.zeros(bins)
        histoutside=np.zeros(bins)
        cluster=np.zeros(bins)
        contour=np.zeros(bins)
        
        for j,i in enumerate(matrix):
            x=np.floor(np.divide(np.subtract(i,mini),pixSize))
            x=[int(i) for i in x.tolist()]
            if self.is3D:
                hist[x[0]][x[1]][x[2]]+=1
                if label[j]>0:
                    histinside[x[0]][x[1]][x[2]]+=1
                    #cluster[x[0]][x[1]][x[2]]=max(cluster[x[0]][x[1]][x[2]],label[j])
                    cluster[x[0]][x[1]][x[2]]=max(cluster[x[0]][x[1]][x[2]],[(1-ii[1]) for ii in contourmat if ii[0]==j][0])
                    
                    contour[x[0]][x[1]][x[2]]=max(contour[x[0]][x[1]][x[2]],[ii[1] for ii in contourmat if ii[0]==j][0])
                else:
                    histoutside[x[0]][x[1]][x[2]]+=1
            else:
                hist[x[0]][x[1]]+=1
                if label[j]>0:
                    histinside[x[0]][x[1]]+=1
                    cluster[x[0]][x[1]]=max(cluster[x[0]][x[1]],label[j])
                    contour[x[0]][x[1]]=max(contour[x[0]][x[1]],[ii[1] for ii in contourmat if ii[0]==j][0])
                else:
                    histoutside[x[0]][x[1]]+=1
        
        if self.is3D:
            hist=np.sum(hist,axis=0)
            histinside=np.sum(histinside,axis=0)
            histoutside=np.sum(histoutside,axis=0)
            cluster=np.max(cluster,axis=0)
            contour=np.max(contour,axis=0)
            
            
        fig=plt.figure()
        plt.imshow(hist)  
        fig=plt.figure()
        plt.imshow(histinside)
        fig=plt.figure()
        plt.imshow(histoutside)
        fig=plt.figure()
        plt.imshow(cluster)
        fig=plt.figure()
        plt.imshow(contour)
        
        
        
    
    #label can be any stat for a given rank: mean_distance ; volume ; density
    #set threshold to add a strait line
    def showHist(self,rank,label='density',binning=100,threshold=-1):
        if self.is3D:
            pw=1./3.
        else:
            pw=1./2.
        densities= [u[label]**pw for x in self.cells for u in x['stat'] if u['rank']==rank]
        dens_range=np.linspace(min(densities), max(densities), num=binning)
        hist=np.histogram(densities,dens_range)
        plt.plot(dens_range[0:binning-1],hist[0])
        if threshold != -1:
            print(threshold)
            plt.plot([threshold**pw,threshold**pw],[np.min(hist[0]),np.max(hist[0])])
        plt.xlabel(label+' (.^'+str(pw)+')')
        plt.ylabel('occurrence number')
        plt.title('Histogram of '+label+' (rank='+str(rank)+')')
    
    
    def showcells(self,triangle,id_triangle):
        
        
        fig=plt.figure()
        ax = Axes3D(fig)
        
        
        tabTr=np.transpose(([x['original'].tolist() for x in self.cells]))
        ax.plot3D(tabTr[0], tabTr[1], tabTr[2], 'k.')
        
        
        #plot contour
        for index,i in enumerate(triangle):
            cc=[]
            for p in i:
                cc.append(self.cells[id_triangle[p]]['original'])
            if (True)|(cc[0][2]<self.maxiBound[2]*.5)&(cc[0][2]>self.maxiBound[2]*.2):
                coll=(Poly3DCollection([cc],linewidths=5, alpha=0.2))
                face_color = [1, 0.5, 0.5] # alternative: matplotlib.colors.rgb2hex([0.5, 0.5, 1])
                coll.set_facecolor(face_color)
                ax.add_collection3d(coll)
                ax.add_collection3d(Line3DCollection([cc], colors='k', linewidths=1, linestyles=':'))
            

        
        
        if self.is3D:
            ax.set_zlim(self.miniBound[2], self.maxiBound[2])
        ax.set_xlim(self.miniBound[0], self.maxiBound[0])
        ax.set_ylim(self.miniBound[1], self.maxiBound[1])
        
        
        
        
        
        
        
        plt.show()
        
        
        
        
    def showpolygon(self,polygon,blackpoints=True,bluepoints=None,redpoints=None,greensegment=None,title='polygon'):
        
        
        fig=plt.figure(num=title)
        ax = Axes3D(fig)
        
        if blackpoints:
            tabTr=np.transpose(([x['original'].tolist() for x in self.cells]))
            ax.plot3D(tabTr[0], tabTr[1], tabTr[2], 'k.')
        
        
        #plot contour
        if polygon is not None:
            for index,i in enumerate(polygon):

                coll=(Poly3DCollection([i],linewidths=5, alpha=0.2))
                face_color = [1, 0.5, 0.5] # alternative: matplotlib.colors.rgb2hex([0.5, 0.5, 1])
                coll.set_facecolor(face_color)
                ax.add_collection3d(coll)
                ax.add_collection3d(Line3DCollection([i], colors='k', linewidths=1, linestyles=':'))
            
        
        
        if (bluepoints is not None):
            tabTrblue=np.transpose(bluepoints)
            ax.plot3D(tabTrblue[0], tabTrblue[1], tabTrblue[2], 'b.')
        if (redpoints is not None):
            tabTrred=np.transpose(redpoints)
            ax.plot3D(tabTrred[0], tabTrred[1], tabTrred[2], 'r.')
        if (greensegment is not None):
            for s in greensegment:
                ax.plot3D([s[0][0],s[1][0]], [s[0][1],s[1][1]], [s[0][2],s[1][2]], 'g-')
        
        if self.is3D:
            ax.set_zlim(self.miniBound[2], self.maxiBound[2])
        ax.set_xlim(self.miniBound[0], self.maxiBound[0])
        ax.set_ylim(self.miniBound[1], self.maxiBound[1])
        
        
        
        fig = plt.gcf()
        fig.set_size_inches(25, 25)

        ax.set_aspect("equal")
        plt.axis('off')
        ax.grid(False)
        
        ax.azim = -90
        ax.dist = 5
        ax.elev = -90
        
        plt.show()
        
        
        
    def plotcluster(self,label):
        
        
        fig=plt.figure(num='clusters')
        ax = Axes3D(fig)
        
        contourmat=[x['original'].tolist() for i,x in enumerate(self.cells) if x['cluster']['label']==label if not x['cluster']['isContour']]
        
        if len(contourmat)>0:
            tabTr=np.transpose(contourmat)
            ax.plot3D(tabTr[0], tabTr[1], tabTr[2], 'b.')
        
        contourmat2=[x['original'].tolist() for i,x in enumerate(self.cells) if x['cluster']['label']==label if x['cluster']['isContour']]
        
        background=[x['original'].tolist() for i,x in enumerate(self.cells) if x['cluster']['label']==0]
        
        tabTr2=np.transpose(contourmat2)
        ax.plot3D(tabTr2[0], tabTr2[1], tabTr2[2], 'r.')
        
        tabBckg=np.transpose(background)
        ax.plot3D(tabBckg[0], tabBckg[1], tabBckg[2], 'k.')
        
        
        if self.is3D:
            ax.set_zlim(self.miniBound[2], self.maxiBound[2])
        ax.set_xlim(self.miniBound[0], self.maxiBound[0])
        ax.set_ylim(self.miniBound[1], self.maxiBound[1])
        
        
        
        
        
        
        
        plt.show()
        
    
        
        
    
    def polyConvex2triangles (self,triangles,polygon,originIndex,currentindex):
        if (len(polygon)==3):
            triangles.append(polygon)
            originIndex.append(currentindex)
        elif (len(polygon)>3):
            for i in range(1,len(polygon)-1):
                triangles.append([polygon[0],polygon[i],polygon[i+1]])
                originIndex.append(currentindex)
    
    
    
    
    #z_weight is used to weight z when z_precision != x_precision or y_precision
    def alpha_shape_3D(self,pos, alpha,z_weight=1):
        """
        Compute the alpha shape (concave hull) of a set of 3D points.
        Parameters:
            pos - np.array of shape (n,3) points.
            alpha - alpha value.
        return
            outer surface vertex indices, edge indices, and triangle indices
            PointInside corresponds to the 4th point opposed to each triangle: useful to compute volume then because it indicates plan orientation
        """
        
        tetra = Delaunay(pos)
        # Find radius of the circumsphere.
        # By definition, radius of the sphere fitting inside the tetrahedral needs 
        # to be smaller than alpha value
        # http://mathworld.wolfram.com/Circumsphere.html
        
        weightedpos=np.copy(pos)
        weightedpos[:,2]*=z_weight
        tetrapos = np.take(weightedpos,tetra.vertices,axis=0)
        normsq = np.sum(tetrapos**2,axis=2)[:,:,None]
        ones = np.ones((tetrapos.shape[0],tetrapos.shape[1],1))
        a = np.linalg.det(np.concatenate((tetrapos,ones),axis=2))
        Dx = np.linalg.det(np.concatenate((normsq,tetrapos[:,:,[1,2]],ones),axis=2))
        Dy = -np.linalg.det(np.concatenate((normsq,tetrapos[:,:,[0,2]],ones),axis=2))
        Dz = np.linalg.det(np.concatenate((normsq,tetrapos[:,:,[0,1]],ones),axis=2))
        c = np.linalg.det(np.concatenate((normsq,tetrapos),axis=2))
        r = np.sqrt(Dx**2+Dy**2+Dz**2-4*a*c)/(2*np.abs(a))

        # Find tetrahedrals
        tetras = tetra.vertices[r<alpha,:]
        
        #print('tetras',np.sort(tetras,axis=1),np.shape(tetras))
        
        
        
        
        # triangles
        TriComb = np.array([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])
        TriCombCompl = np.array([(3), (2), (1), (0)])
        
        #print('TriComb',TriComb)
        #print('tmp',tetras[:,TriComb])
        Triangles = tetras[:,TriComb].reshape(-1,3)
        TrianglesCompl = tetras[:,TriCombCompl].reshape(-1,1)
        #print('TrianglesUnsorted',Triangles)
        
        #print('TrianglesCompl',TrianglesCompl)
        
        
        
        Triangles = np.sort(Triangles,axis=1)
        
        compl=dict()
        for i,tri in enumerate(Triangles):
            tup=tuple(tri)
            if tup not in compl:
                compl.update({tup:TrianglesCompl[i].tolist()})
            else:
                compl.update({tup:compl[tup]+TrianglesCompl[i].tolist()})
                
        
        # Remove triangles that occurs twice, because they are within shapes
        TrianglesDict = defaultdict(int)
        
        for tri in Triangles:TrianglesDict[tuple(tri)] += 1
        
        #print('TrianglesDict',TrianglesDict)
        
        Triangles=np.array([tri for tri in TrianglesDict if TrianglesDict[tri] ==1])
        #print('Triangles',Triangles)
        PointInside=[]
        for tri in Triangles:
            tup=tuple(tri)
            PointInside.append(compl[tup][0])
        #print('PointInside',PointInside)
        
        #edges
        EdgeComb=np.array([(0, 1), (0, 2), (1, 2)])
        Edges=Triangles[:,EdgeComb].reshape(-1,2)
        Edges=np.sort(Edges,axis=1)
        Edges=np.unique(Edges,axis=0)

        Vertices = np.unique(Edges)
        
        return Vertices,Edges,Triangles,PointInside


    
    
    
        
        
    def computeSurfaceVolumeConcav(self,data,alpha,show=False,z_weight=1):
        


        [as_v,as_e,as_t,as_pin]=self.alpha_shape_3D(data,alpha,z_weight=z_weight)



        polygons=[]
        point_inside=[]
        for t in as_t:
            polygons.append([data[t[0]],data[t[1]],data[t[2]]])
        for pin in as_pin:
            point_inside.append(data[pin])
        
        if show:
            self.showpolygon(polygons,blackpoints=False,bluepoints=data,title='triangular concav hull')
        
        
        
                
        
        
        
        
        
        
        surface=0
        
        
        
        pointOutsideSegmentation=None
        for cell in self.cells:
            if cell['cluster']['label']==0:
                pointOutsideSegmentation=cell['original']
                break
        
        
        #for i,s in enumerate(surfaceList):
        #    if not 'surface' in self.cells[idSurfaceCell[i]]:
        #        self.cells[idSurfaceCell[i]].update({'surface':[{'polygon':s,'adjacent':adjacentList[i]}]})
        #    else:
        #        self.cells[idSurfaceCell[i]].update({'surface':self.cells[idSurfaceCell[i]]['surface']+[{'polygon':s,'adjacent':adjacentList[i]}]})
        
        
        volume=0
        for i,s in enumerate(polygons):
            #print(s,adjacentList[i])
            #tmpPt=[]
            #tmpPt.append(list(s))
            #self.showpolygon(tmpPt,[adjacentList[i]],[pointOutsideSegmentation])
            
            plan=self.planEquation(s)

            orientationAdjacent=self.getOrientation(plan,point_inside[i])
            orientationseed=self.getOrientation(plan,pointOutsideSegmentation)
            #print("2pointsorient ",pointOutsideSegmentation,adjacentList[i])
            #print("orientation",orientationAdjacent,orientationseed)
            orientation=-1
            if (np.sign(orientationAdjacent)==np.sign(orientationseed)):
                orientation=1
            
            
            dist=self.distance2plan(plan,pointOutsideSegmentation)
            #print('dist:',dist)
            
            area=self.area(s)
            surface+=area
            #print('area:',area)
            
            #print('orientation:',orientation)
            
            vol=orientation*self.getVolumeTetraedre(area,dist)
            
            if not np.isnan(vol):
                volume+=vol
        print('volume (alpha=',alpha,') = ',volume)
        return([surface,volume])
        
        
        
        
        
        
    def computeSurfaceVolumeConvex(self,meshPoints,show=False):
        
        #here, we compute convex hull
        hull = ConvexHull(meshPoints)
        #print("convex surf ",hull.area,"   vol",hull.volume)
        surfaceConvex=hull.area
        volumeConvex=hull.volume
        if show:
            v=[]
            for vv in hull.simplices:
                vvvv=[]
                for vvv in vv:
                    vvvv.append(meshPoints[vvv])
                v.append(vvvv)
            
            self.showpolygon(v,blackpoints=False,bluepoints=meshPoints,title='convex hull')
        return([surfaceConvex,volumeConvex])
    
    



    
    #get equation of a plan given by 3 points (http://kitchingroup.cheme.cmu.edu/blog/2015/01/18/Equation-of-a-plane-through-three-points/)
    def planEquation(self,points):
        
        
        if len(points)<3:
            print("error: number a points not enough")
        
        
        p1 = points[0]
        p2 = points[1]
        p3 = points[2]
        
        # These two vectors are in the plane
        v1 = np.subtract(p3,p1)
        v2 = np.subtract(p2,p1)
        
        # the cross product is a vector normal to the plane
        cp = np.cross(v1, v2)
        a, b, c = cp
        
        # This evaluates a * x3 + b * y3 + c * z3 which equals d
        d = -np.dot(cp, p3)
        
        #print('The equation is {0}x + {1}y + {2}z + {3} = 0'.format(a, b, c, d))
        #print('a, b, c, d : ',a, b, c, d)
        #print('points=',points)
        return([a,b,c,d])
    
    
    #indicates if point is above or below the plan
    def getOrientation(self,plan,point):
        return(plan[0]*point[0]+plan[1]*point[1]+plan[2]*point[2]+plan[3])
        
    
    #compute the distance 
    def distance2plan(self,plan,point):
        
        d=np.abs(plan[0]*point[0]+plan[1]*point[1]+plan[2]*point[2]+plan[3])/np.sqrt(plan[0]**2+plan[1]**2+plan[2]**2)
        
        return(d)
    
    #it is still ok of it is polygon instead of triangle since polygon=set of triangle
    def getVolumeTetraedre(self,triangleSurface,depth):
        return triangleSurface*depth/3
    
    
   
        
    
        
        
    def getVolume(self,classlabel):
        volume=0
        for i,cell in enumerate(self.cells):
            if cell['cluster']['label']==classlabel:
                volume+=cell['volume']
            
        return(volume)
    
        
        
    def getEmitterNumber(self):
        
        n=0
        
        for i,cell in enumerate(self.cells):
            
            if cell['cluster']['label']>0:
                n+=1
        
        return(n)
    
    
    
        
    
    
    #return list with surface volume ratio gyration density locNumber surfaceConvex, volumeConvex
    #alpha is parameter for alpha_shape: should befunction of localization precision
    def getFullStat(self,alpha,rank,z_weight=1,show=False):
        
        stat=np.zeros((32))
        
        
        sphereNormalization=(((4./3.)*3.141592)**2)/((4.*3.141592)**3)
        surfacePointList=None
        gyration=0
        gyrationY=0
        gyrationX=0
        gyrationZ=0
        density=0
        totcount=0
        mean=np.zeros((3))
        labelOK=[]
        fullLocNumber=0
        data=[]
        adjacent=[]
        for i,cell in enumerate(self.cells):
            fullLocNumber+=1
            if cell['cluster']['label']>=1:
                data.append(cell['original'])
                totcount+=1
                density+=cell['stat'][rank]['density']
                mean=np.add(mean,np.array(cell['original']))
                
        mean/=totcount
        
        if (totcount<3):
                return([0,0])
        
        
            
        
        
        for i,cell in enumerate(self.cells):
            if cell['cluster']['label']>=1:
                 gyration+=np.sum((np.array(cell['original'])-mean)**2)
                 #gyrationLateral+=np.sum((np.array(cell['original'])[0:2]-mean[0:2])**2)
                 gyrationX+=np.sum((np.array(cell['original'])[0]-mean[0])**2)
                 gyrationY+=np.sum((np.array(cell['original'])[1]-mean[1])**2)
                 gyrationZ+=np.sum((np.array(cell['original'])[2]-mean[2])**2)
                                                        
                                                        
        
        density/=totcount
        gyration/=totcount
        gyrationX/=totcount
        gyrationY/=totcount
        gyrationZ/=totcount
        
        
        data=np.add(data,0.000001*np.random.random(np.shape(data)))#small shift to avoid errors while calculating determinant in computeSurfaceVolumeConcav function
        print(data)
        svConcav=self.computeSurfaceVolumeConcav(data,alpha,show=False,z_weight=z_weight)
        svConvex=self.computeSurfaceVolumeConvex(data,show=False)
        
        
        stat[0]+=svConcav[0]
        stat[1]+=svConcav[1]
        stat[2]=pow(stat[1],2)/(sphereNormalization*pow(stat[0],3))
        stat[3]=gyration
        stat[4]=density
        stat[5]=totcount                                                
        stat[6]=svConvex[0]  
        stat[7]=svConvex[1]  
        stat[8]=fullLocNumber
        stat[9]=gyrationX
        stat[10]=gyrationY
        stat[11]=gyrationZ
                        
                        
        svConcav50=self.computeSurfaceVolumeConcav(data,50,show=False,z_weight=z_weight)
        svConcav100=self.computeSurfaceVolumeConcav(data,100,show,z_weight=z_weight)
        svConcav150=self.computeSurfaceVolumeConcav(data,150,show=False,z_weight=z_weight)
        svConcav200=self.computeSurfaceVolumeConcav(data,200,show=False,z_weight=z_weight)
        svConcav250=self.computeSurfaceVolumeConcav(data,250,show=False,z_weight=z_weight)
        svConcav300=self.computeSurfaceVolumeConcav(data,300,show=False,z_weight=z_weight)         
        svConcav350=self.computeSurfaceVolumeConcav(data,350,show=False,z_weight=z_weight)         
        svConcav400=self.computeSurfaceVolumeConcav(data,400,show=False,z_weight=z_weight)         
        svConcav450=self.computeSurfaceVolumeConcav(data,450,show=False,z_weight=z_weight)                
        svConcav500=self.computeSurfaceVolumeConcav(data,500,show,z_weight=z_weight)                      
                 
        
        stat[12]=svConcav50[1]
        stat[13]=svConcav100[1]
        stat[14]=svConcav150[1]
        stat[15]=svConcav200[1]
        stat[16]=svConcav250[1]
        stat[17]=svConcav300[1]
        stat[18]=svConcav350[1]
        stat[19]=svConcav400[1]
        stat[20]=svConcav450[1]
        stat[21]=svConcav500[1]
        
        stat[22]=svConcav50[0]
        stat[23]=svConcav100[0]
        stat[24]=svConcav150[0]
        stat[25]=svConcav200[0]
        stat[26]=svConcav250[0]
        stat[27]=svConcav300[0]
        stat[28]=svConcav350[0]
        stat[29]=svConcav400[0]
        stat[30]=svConcav450[0]
        stat[31]=svConcav500[0]
               
        return(stat)
    
    
    
    
    
    
    
    
    

        
        
    
    
                
    def getGyrationRadius(self):
        
        
        count=0
        m=np.empty(3)
        count=0
        for label in self.clusterLabels:
            for i,cell in enumerate(self.cells):
                if cell['cluster']['label']==label:
                    m=np.add(m,np.array(cell['original']))
                    count+=1
        if (count==0):
            return 0
        m=np.true_divide(m,count)
        v=0
        for label in self.clusterLabels:
            for i,cell in enumerate(self.cells):
                if cell['cluster']['label']==label:
                    v+=np.sum((np.array(cell['original'])-m)**2)
        v/=count
        return(v)
    
    
    
    def getDensity(self,rank=2):
        
        count=0.
        d=0.
        for label in self.clusterLabels:
            for i,cell in enumerate(self.cells):
                if cell['cluster']['label']==label:
                    d+=cell['stat'][rank]['density']
                    count+=1
        if (count==0):
            return 0
        d/=count
        return d
    
    def getClusterTable(self):
        
        table=[]
        
        for label in self.clusterLabels:
            for i,cell in enumerate(self.cells):
                if cell['cluster']['label']==label:
                    table.append(list(cell['original']))
        
        return(table)
    
    
    def saveClusterTable(self,path):
        
        table=[]
        
        
        
        fff= open(path,"w+")
        fff.write("x,y,z,cluster\r\n")
        labelNumber=len(self.clusterLabels)
        clusterLabelscopy=copy.copy(self.clusterLabels)
        random.shuffle(clusterLabelscopy)
        labelshuffle=random.shuffle(self.clusterLabels)
        for labid,label in enumerate(self.clusterLabels):
            
            for i,cell in enumerate(self.cells):
                if cell['cluster']['label']==label:
                    table.append(list(cell['original']))
                    l=cell['original']
                    fff.write("%.12f,%.12f,%.12f,%f\r\n" % (l[0],l[1],l[2],clusterLabelscopy[labid]))
                    

        fff.close()
        
        
        return(table)
    
    
    def saveTable(self,path,data=None):
        
        
        if data is None:
            data=self.data
        
        fff= open(path,"w+")
        fff.write("x,y,z\r\n")
        
            
        for d in data:
            
            fff.write("%.12f,%.12f,%.12f\r\n" % (d[0],d[1],d[2]))
                    

        fff.close()
        
        
        
        
        
        

