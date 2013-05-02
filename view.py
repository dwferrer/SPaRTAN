import sys
import ctypes as ct
import numpy as np
import csv
import time
import os
import matplotlib.pyplot as p 
import matplotlib.cm as cm
#import sph

from mpl_toolkits.mplot3d import Axes3D
from numpy.core.numeric import dtype
#sphlib = ct.cdll.LoadLibrary("./Release/libsph.so")

def arrayToCPointer(parr): #turn a numpy array into something useable by ctypes
    return  parr.ctypes.data_as(ct.POINTER(ct.c_float))

def sanitize(a): # convert an array into a float32 array
    return np.ascontiguousarray(a, dtype =  np.float32)

def plots(fn,totalt):
#    fn = "/home/dferrer/remote/11-26homosphere/step220.npz"
    l = np.load(fn)
    #ps, dt = sph.injectIC(fn)
    #size = sphlib.getsize(ps)
    #h = np.zeros(size,dtype =np.float32)
#    sphlib.geth(ps,arrayToCPointer(sanitize(h)))
    pos = l['pos']
    dens = l['dens']
    dt = l['dt']
    rad = np.sqrt(np.sum(pos*pos,1))   
    n1d = 60
    fig1 = p.figure(0)
    
    ax1 = fig1.add_subplot(111)
    center = np.where(dens == np.max(dens))
    cxyz = (pos[center,:])[0]
    print cxyz
    rad  = np.sqrt(np.sum((pos)**2,1))
    idx = np.where(rad < 7000)[0]
    #pos2 = pos[idx,:]
    print "Selected: " +str(len(idx))
    
    
    tcm = cm.get_cmap(name = "jet")
    tcm._init()
    alphas = np.linspace(0,1,tcm.N)
    tcm._lut[:-3,-1] = alphas
    totalt+=dt
#    h_in_points = np.diff(ax1.transData.transform(zip([0]*len(h[idx]), h[idx])))s = h_in_points/2,
    plt1 = ax1.scatter(pos[idx,0],pos[idx,1],c =np.log10(dens[idx]),cmap = tcm, edgecolors='none')
    ax1.set_xlabel("x (AU)")
    ax1.set_ylabel("y (AU)")
    ax1.set_title("t = "+str(totalt)+"ka")
    #cb1 = fig1.colorbar(plt1)
    #cb1.set_label("log(density) (Msun/AU$^3$)")
    fig1.savefig(fn+".zoomed.png",format = 'png')
    
    
    fig2 = p.figure(1)
    ax2 = fig2.add_subplot(111)#,yscale = "log")
    d, bins = np.histogram(rad, bins = 30,range = (0,40000), weights = dens)
    n,bins =  np.histogram(rad, bins = 30,range = (0,40000))
    r = np.array([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)])
    dr = np.array([(bins[i+1]-bins[i]) for i in range(len(bins)-1)])
    dm = 4*np.pi*d/n *r**2 *dr
    
    m = np.array([np.sum(dm[0:i]) for i in range(len(r))])
    
    ax2.plot(r, m)
    ax2.set_title("t = "+str(totalt)+"ka")
    ax2.set_xlabel("radius (AU)")
    ax2.set_ylabel("Mass enclosed (Msun)")
    fig2.savefig(fn+".density.png",format = 'png')
    

    fig3 = p.figure(2)
    ax3 = fig3.add_subplot(111, projection='3d')
    idx2 = np.where(rad < 40000)[0]
    plt3 = ax3.scatter(pos[:,0],pos[:,1],pos[:,2],c =np.log10(dens[:]),cmap = tcm,edgecolors='none')
    cb3 = fig3.colorbar(plt3)
    cb3.set_label("log(density) (Msun/AU$^3$)")
    #ax3.hexbin(pos[:,0],pos[:,1],C = dens[:], gridsize = 20, extent = (-370,-320,115,165),reduce_C_function = np.max )
    ax3.set_title("t = "+str(totalt)+"ka")
    fig3.savefig(fn+".full.png",format = 'png')
    fig4 = p.figure(3)
    ax4 = fig4.add_subplot(111)
    #columndens = dens /(.01(1-) 
    plt4  = ax4.hexbin(pos[idx2,0],pos[idx2,1],C = np.log10(dens[idx2]), gridsize = 60,reduce_C_function = np.max )
    ax4.set_xlabel("x (AU)")
    ax4.set_ylabel("y (AU)")
    ax4.set_title("t = "+str(totalt)+"ka")
    cb4 = fig4.colorbar(plt4)
    cb4.set_label("log(max density) (Msun/AU$^3$)")
    fig4.savefig(fn+".column.png",format = 'png')
    #p.close('all')
    return totalt

if __name__ == '__main__':
    '''files =  os.listdir("/home/dferrer/ted/turbsph/jrturb")
    '''
    totalt = 0.0
    for i in range(0,300):
        print i
        totalt = plots("/home/doug/turbsph/step"+str(i)+".npz",totalt)
    #plots("./step165.npz",105)
    #p.show()   
