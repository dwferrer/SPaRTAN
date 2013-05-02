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

def getindex(fn,totalt):
    l = np.load(fn)
    pos = l['pos']
    dens = l['dens']
    dt = l['dt']
    r = np.sqrt(np.sum(pos**2,1))
    totalt += dt
    fig1 = p.figure(0)
    ax1 = fig1.add_subplot(111, yscale = 'log')
    ax1.plot(r,dens, ',', label = "Particles")
    ax1.set_xlabel("r (AU)")
    ax1.set_xlim(0,20000)
    ax1.set_ylabel("Density (Msun/AU$^3$)")
    ax1.set_title("t = "+str(totalt)+"ka")
    fig1.savefig(fn+".rhovr.png",format = 'png')
    p.close('all')
    medr = np.median(r)
    return totalt,medr

if __name__ == '__main__':
    '''files =  os.listdir("/home/dferrer/ted/11-26homosphere/")
    totalt = 0.0
    times = []
    r50s = []
    for i in range(len(files)-1):
        print i
        totalt,medr = getindex("/home/dferrer/ted/11-26homosphere/step"+str(i)+".npz",totalt)
        times.append(totalt)
        r50s.append(medr)
    times = np.array(times)
    r50s = np.array(r50s)
    np.savez("./r50vstime.npz",times = times, r50s = r50s)'''
    l = np.load("./r50vstime.npz")
    times = l['times']
    r50s = l['r50s']
    p.plot(times, r50s,'.', label = "Simulation")
    r0 = r50s[0]
    theta = np.arccos(np.sqrt(r50s/r0))
    theorytimes = 2*(theta+ .5 * np.sin(2*theta))/.02808
    p.plot(theorytimes, r50s, label = "Homologous Collapse")
    p.xlabel("Time (ka)")
    p.ylabel("Median Radius (AU)")
    p.legend()
    p.savefig("./r50svstime.png")
    p.show()
    
