#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os, sys, inspect
#from Scientific.IO.NetCDF import NetCDFFile

import matplotlib.pyplot as py

import matplotlib.mlab as ml
import matplotlib.colors as col
import ntpath


import configparser
config = configparser.RawConfigParser()

from math import *
import pylab as p
import re

#import myfct2
from mpl_toolkits.mplot3d import Axes3D
import copy
from array import array
import numpy as np
from scipy.interpolate import griddata as gd #rajout
import urllib

from tkinter import *

#from netCDF4 import Dataset
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs

#import scipy as sp  # SciPy (signal and image processing library)
#import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
#import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
from pylab import *              # Matplotlib's pylab interface
import numpy.ma as ma
#         plt.ion()



"""-------------------------------------------"""

class files:

    def __init__(self,name):    
        """ """

    def set_name(self,name):
        self.name=name
        
    def get_name(self):
        return self.name
    
    def openf(path):
        """ """
        
    def get_keys(self):       
        """ """
        
    def out_keys(self):
        keys=self.get_keys()
        print("keys list",self.get_name())
        createfieldsfile(keys,ntpath.basename(self.get_name()))
        
class ncdf_files(files):
    
    def __init__(self):    
        """ """
    def openf(path):
        """ """
    def get_keys(self):    
        f = NetCDFFile(self.name)
        print("ncdf", self.name)
        var_list=f.variables.keys()
        return var_list 

    def get_allval(self):
        """ """
        
    def get_val(self,key):
        print('files ncdf', self.name)
        f = NetCDFFile(self.name)
#        print 1
        valnc= f.variables[key]    
#        print 2        
        val=valnc.getValue()
#        print 3        
        f.close()
        
        """ set inf values to zero """    
        np.putmask(val, val >= 1e+20, 0)
        return val    
    
class Surf_txt(files):
    
    def __init__(self):    
        """ """
    
    def openf(path):
        """ """
    def get_keys(self):     
        """return list of fields of pgd,prep data """   
        
        fichier = open (self.name, "r") 
        strgarray= (fichier.read()).split("\n")
        fichier.close()        
        
        ii=[]
        for i in range(np.size(strgarray)):
            try:
                if strgarray[i].index('&')>-1 : ii.append(i) 
                
            except ValueError:
                        """"""
        #c = [ a[i] for i in b]    
        fields=[strgarray[i] for i in ii]  
        return fields   
        
    def get_allval(self):
        fichier = open (self.name, "r") 
        fdata= (fichier.read()).split("\n")
        fichier.close()
        return fdata      
        
    def get_val(self,key):
        """ sub2 dat splté" """
        
#        key=key[1:]
        print(key)
        sub2=self.get_allval()
        for i in range(np.size(sub2)):    
            try:    
                car=sub2[i].index(key)
                if car>0: break        
#                print('car',car,'ind',i) 
        #        print('&&&',lnf)        
#                print(sub2[i])
            except ValueError:
                """"""
#         print ("notfound")       
        ib=i
        sub3=sub2[ib+2:]    
        for i in range(np.size(sub3)):   
            try:    
                car=sub3[i].index('&')
                if car>0: break        
#                print('car',car,'ind',i) 
        #        print('&&&',lnf)        
#                print(sub3[i])
            except ValueError:
                """ """
    #         print ("notfound") 
        ie=i  
        val=sub3[:ie]        
        
        val=txtfmttoarray(val)
        return val     
        
        
class Surf_txt_u(Surf_txt):        
    """ """
    
    def get_val(self):
        fichier = open (self.name, "r") 
        fdata= (fichier.read()).split("\n")
        fichier.close()
        return fdata
        
"""-----------------------------------------------------------"""
       
class data:
    
#    self.timed=TRUE
    
    def __init__(self):
        """Par défaut, notre surface est vide"""
        self.filetype = ""
        self.file= ""
        self.field= ""
        self.expe= ""
#        self.key=""
        self.fmt=""   
    
    def set_file(self,files):
        self.file=files
        
    def get_file(self):
        return self.file
    
#    def get_val(self,key):
        """-"""
    def trace(self):
        """ """
        
    def openf(self, path):
        """Méthode class file !!! """
    
    def transp(self):
        """ """
        
    def set_field(self,field):
        self.field=field
               
    def get_field(self):
        return self.field

    def set_expe(self,expe):
        self.expe=expe
               
    def get_expe(self):
        return self.expe
    
    def get_val(self):
        """ """
        
#        print self.fmt
        if self.fmt=="pack":  
            print("pack1")
            val=self.file.get_val(self.field)
        else:   
            val=self.file.get_allval()
            val=txtfmttoarray(val)
        
        return val    
        
class c1d(data):
    
    def __init__(self):   
        """ """

    def trace(self):
        """ """     
          
        loc_legend='upper right'
            
        labl=self.field
        expe=self.expe
        print("expe = ",expe)
        print("label = ",labl+" "+expe)
        
        C=self.get_val()          
        C=C.flatten()
        if expe=="" :
          l1=py.plot(C, label=labl)
        else :
          l1=py.plot(C, label=labl+"  ("+expe+")")
          
        py.setp(l1, linestyle='-', linewidth=1.5)
        py.legend(loc=loc_legend,prop={'size': 11})          
        
class c2d(data):
    
    def __init__(self):   
        """ """

    def trace(self):
        """ """      
        loc_legend='upper right'
        labl=self.field
        print("label",labl)        
        C=self.get_val()
        
        C=C[:,0,:]
        xmin=0
        xmax=size(C[0,:])
        ymin=0
        ymax=size(C[:,0])    
        print("coord", xmin,xmax,ymin,ymax)
        
        py.imshow(C, interpolation='nearest',extent=[xmin, xmax, ymin, ymax]) 
        py.title(labl)
        print("label", labl)
        
        """ orientation yyyy ????"""
        
class dmap(data):        
    
    def __init__(self):    
        """ """
        
    def trace(self):
        """ """
              
        selkey=self.field
        y=d.file.get_val('XLAT') 
        x=d.file.get_val('XLON') 
        z=d.file.get_val(selkey)     
        print("y",y, size(y)) 
        print("z",z, size(z))
        print("x",x, size(x))
        coeff=1
        
        side=int(sqrt(np.size(z)))*coeff  # find x,y size 

        ny, nx = side, side
        xmin, xmax = min(x), max(x)
        ymin, ymax = min(y), max(y)
        ##----------------------
        #griddata deprecated since python 2.2  --> alternative to griddata
        #-----------------------
        x = np.r_[x,xmin,xmax]
        y = np.r_[y,ymax,ymin]
        z = np.r_[z,z[0],z[-1]]
        xi = np.linspace(xmin, xmax, nx)
        yi = np.linspace(ymin, ymax, ny)
        xi, yi = np.mgrid[xmin:xmax:nx*1j,ymin:ymax:ny*1j]
        zm = z
        zm = np.ma.masked_where(np.isnan(z),z)   #pour masquer Nan values si il y a
        zi = gd((x,y), zm, (xi,yi), method='nearest')
        
        blue = '#0000ff'   
        marron = '#604242'   
        green = '#309033' 
        yellow ='#dae709' 
        rose= '#DFADEA' 
        turq ='#13BFD0'    
        orange ='#b88027'
        red = '#bf282b'
        cmap2 = col.LinearSegmentedColormap.from_list('own2',[blue,turq,green,yellow,orange,red,rose])

        if np.nanmin(z) != np.nanmax(z):
           py.pcolormesh(xi, yi, zi, cmap = py.get_cmap(cmap2),norm = mpl.colors.Normalize(np.nanmin(z),np.nanmax(z)))
        else:
           py.pcolormesh(xi, yi, zi, cmap = py.get_cmap(cmap2),norm = mpl.colors.Normalize(np.nanmin(z)-(0.01*np.nanmin(z)),np.nanmax(z)+(0.01*np.nanmax(z))))
        #py.drawcoastlines(color='lightgray')
        py.colorbar()
        
        py.scatter(x, y, marker = 'o', c = 'b', s = 0, zorder = 5)
        py.xlim(xmin, xmax)
        py.ylim(ymin, ymax)
#        py.show()
        
        py.xlabel('Longitude')
        py.ylabel('Latitude')
        py.title(self.field)
        
        tx= (np.linspace(xmin, xmax, 10)).astype(np.float) 
        ty= (np.linspace(ymin, ymax, 10)).astype(np.float) 
        txr = [ '%.1f' % elem for elem in tx ]
        tyr = [ '%.1f' % elem for elem in ty ]
        py.xticks(tx,txr)
        py.yticks(ty,tyr) 
#        py.show()         
        

            
class formt():
    
    def __init__(self):    
        """ """
    def openf(self,path):
        """ """        
        

class proj(formt):
        """ """

        

"""-----------------------------------------------------------"""


def arraytoline(tab):    
    """ same as flatten ??"""
    tabo=[]
    for i in range(np.size(tab)):
        tabo.extend(tab[i].split())
    return tabo

def get_val(sub2,key):
    for i in range(np.size(sub2)):
    #    print('i',i)    
        try:    
            car=sub2[i].index(key)
            if car>0: break        
            """     print('car',car,'ind',i) """
    #        print('&&&',lnf)        
            """   print(sub2[i])"""
        except ValueError:
            """"""
#         print ("notfound")       
    ib=i
    sub3=sub2[ib+2:]
#    print(sub3)    
    for i in range(np.size(sub3)):
    #    print('i',i)    
        try:    
            car=sub3[i].index('&')
            if car>0: break        
            """   print('car',car,'ind',i) """
    #        print('&&&',lnf)        
            """  print(sub3[i]) """
        except ValueError:
            """"""
#         print ("notfound") 
    ie=i
    
    """print('start,end', ib,ie)   """
    val=sub3[:ie]        
         
    return val
    
def txtfmttoarray(tab):
    """ concatenate multi col data to one and convert to float array  """

    exponent = re.compile(r'(?<=\d|\.)D?(?=(?:\+|-)\d)')
    """print("size", size(tab)) """
    ee=arraytoline(tab)
    
    dd=[]
    
    for line in ee:
        dd.append(float(exponent.sub('e', line)))
    
    tabn=np.asarray(dd)    
    """ set inf values to zero """    
    np.putmask(tabn, tabn>=1e+20, nan)
    return tabn   

def getallfields(strgarray):
    """return list of fields of pgd,prep data """   
    ii=[]
    for i in range(np.size(strgarray)):
        try:
            if strgarray[i].index('&')>-1 : ii.append(i) 
            
        except ValueError:
                    """"""
    #c = [ a[i] for i in b]    
    fields=[strgarray[i] for i in ii]  
    return fields
    
def createfieldsfile(liste,fic):

    outlist="/home/dasprezs/python/"
    prmslist="keysli_"+fic        
    
    f = open(outlist+prmslist, "w+")
    
    for i in range(np.size(liste)):
        f.write(str(i)+" "+liste[i]+"\n")
    #    f.write("\n")
    f.close() 
    
    cmd="gedit "+outlist+prmslist
    os.system(cmd)
    
def numtokey(arg,keys):
    """!! return array  """
    
    if np.isreal(arg):
    #    ZS=keys3[ZS][1:]
        arg=[arg]
        """print(arg)"""
        
    elif arg.find(',')!=-1:
        arg=arg.split(',')
        arg =map(int,arg)
    else :
        arg=[int(arg)    ]
        
    kl=[]
    
    """print(arg)"""
        
    for i in arg: 
        print(i)
        i=keys[i]     
        kl.append(i)
    
    print(kl)     
    return kl    
    
def filedetect(fic):
    
    fileName, fileExt = os.path.splitext(fic)
#
    """print(fileExt)"""
    
    if fileExt=='.nc':
        print("nc")    
        f=ncdf_files()
        
    elif fileExt=='.txt' or '.TXT':
        print("txt")
        f=Surf_txt()       
#        if field=="": f=Surf_txt()
#        else:  f=Surf_txt_p()   
        
    else:
        print("extension not found")   
        
    return f    
    
    """ liste droulante """   

 
    
def tracei(d):
    """    """
    
def get_1val(fic): 
    fichier = open (fic, "r") 
    fdata= (fichier.read()).split("\n")
    fichier.close()
    return fdata  
    

#
#"""generer date de debut!!! pas en pgd  """
#

"""-----------------------------------------------------------"""

def set_class(L):
    
    kind=L[1]
    field=L[2]
    expe=L[4]
    fic=L[3]
    
    """
    if field="": 
        fk=unbound 
        varname=filename
     """   
    if kind=="C1d" :
        d=c1d()
    elif kind=="C2d":
        d=c2d()
    elif kind=="Map" :
        d=dmap()
    else:
        print("kind not found")
        
    """print("field =", "_"+field+"_")     """
    
    if field!="" :
        d.fmt="pack"
    else:
        d.fmt="sgl"
        field=os.path.splitext(ntpath.basename(fic))[0]       
        d.field=field
     
#    print "fmt set",d.fmt    

    
    f=filedetect(fic) # get extension and set class type txt or ncf
    f.set_name(fic) # rec file name 
    #
     # -----------------------------!!!
    
    d.set_file(f)  # rec file info
    d.set_field(field)
    d.set_expe(expe)
    return d
    
def set_class2(L):
    
    kind=L[1]
    field=L[2]
    expe=L[4]
    fic=L[3]
    
    """
    if field="": 
        fk=unbound 
        varname=filename
     """   
    if kind=="C1d" :
        d=c1d()
    elif kind=="C2d":
        d=c2d()
    elif kind=="Map" :
        d=dmap()
    else:
        print("kind not found")
        
    """print("field =", "_"+field+"_")     """
    
    if field!="" :
        d.fmt="pack"
    else:
        d.fmt="sgl"
        field=os.path.splitext(ntpath.basename(fic))[0]       
        d.field=field
     
#    print "fmt set",d.fmt    

    
    f=filedetect(fic) # get extension and set class type txt or ncf
    f.set_name(fic) # rec file name 
    
    d.set_file(f)  # rec file info
    d.set_field(field) 
    d.set_expe(expe)
    return d 
    
def create_plot_list(fic):

    config.read(fic)
    if fic==[]: print(".cfg file not found")

    kind = config.get('Conf', 'kind')
    nbfig = config.getint('Conf', 'nbfig')
    
    print(kind)
    
    Tg=np.arange(nbfig,dtype=list)
    
    """ create array"""
    
    for i in range(nbfig):
        sec= "Fig"+str(i+1)
        """print(sec)"""
        path= config.get(sec, 'path')   
        field= config.get(sec, 'field')
        expe=config.get(sec, 'expe')
        Tg[i]=[sec,kind,field,path,expe]
        if os.path.exists(path)==False:
            print("file "+path+"not found") 
            
    return Tg
        

if __name__ == "__main__":
    """   """
    
pathn="/home/dasprezs/mydir/EXPORT_v7_2_1/MY_RUN/KTEST/Gcexpe2d/TG1.nc"   

pathn="/home/dasprezs/mydir/EXPORT_v7_2_1/MY_RUN/KTEST/hapex/ISBA_PROGNOSTIC.OUT.nc"  

IDE=False


""" local path """


if IDE==False : 
       
    if len(sys.argv)==1: 
        cfg="conf.cfg"    
    else:    
        cfg=sys.argv[1]
        
else:
    cfg="ex_map.cfg"
    print("EDI mode") 

print("file is : ",cfg)

Tg=create_plot_list(cfg)
    
print(Tg)    



for L in Tg:
    d=set_class(L)
#    print "fff", d.field, d.fmt
    try:
        d.trace()
    except:
        print(""" failed to plot""")




plt.show()
 
 

"""
dd=get_1val(L[3])

dd2=txtfmttoarray(dd)

plot(dd2)
"""

"""
hh=d.file.g


print("ok")

allkeys=d.file.get_keys()  
"""



"""print(d.__class__.__name__)"""



if not IDE==True : raw_input("Press a key to exit")

#print d.field

""" time === 9 """



"""
for i in Tg:
    print(i)
    fic=i[3]
    f=filedetect(fic) # get extension and set type
    f.set_name(fic) # rec file name 
    #
    d=data()    # -----------------------------!!!
    d.set_file(f)  # rec file info
    traceC(d,i[2])
"""


    
#createfieldsfile(allkeys,fic)

#key=allkeys[4]
#sttype1=d.file.get_val(key[1:])

#selkey=numtokey(32,allkeys)



"""--------------------------------------------"""




"""-----------------------------------------------"""






