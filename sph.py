#!/usr/bin/python

#driver for basic sph tests


import sys
import ctypes as ct
import numpy as np
import csv
import time
import matplotlib.pyplot as p 

from mpl_toolkits.mplot3d import Axes3D

sphlib = ct.cdll.LoadLibrary("./libsph.so")
sphlib.stepforward.restype = ct.c_float
sphlib.domicro.restype = ct.c_float
sphlib.rts.restype = ct.c_float

def arrayToCPointer(parr): #turn a numpy array into something useable by ctypes
    return  parr.ctypes.data_as(ct.POINTER(ct.c_float))

def sanitize(a): # convert an array into a float32 array
    return np.ascontiguousarray(a, dtype =  np.float32)

def plotdens(ps, fignum):
    size = sphlib.getsize(ps)
    pos = np.zeros((size,3),dtype = np.float32)
    
    sphlib.getPos(ps,arrayToCPointer(sanitize(pos)))
    dens = np.zeros(size, dtype = np.float32)
    sphlib.getDens(ps,arrayToCPointer(sanitize(dens)))
    fig1 = p.figure(fignum)
    ax1 = fig1.add_subplot(111, projection='3d')
    
    plt1 = ax1.scatter(pos[:,0],pos[:,1],pos[:,2],c =dens,alpha = .1)
    fig1.colorbar(plt1)
    p.show()
    
def addturbulence(ps):
    size = sphlib.getsize(ps)
    pos = np.zeros((size,3),dtype = np.float32)
    sphlib.getPos(ps,arrayToCPointer(sanitize(pos)))
    vel = np.zeros((size,3),dtype = np.float32)
    sphlib.getVel(ps,arrayToCPointer(sanitize(vel)))
    
    n1d = 256
    seed = 1234567
    np.random.seed(seed)
    r = np.random.rand(n1d,n1d,n1d)
    theta = np.random.rand(n1d,n1d,n1d)*2*np.pi
    
    ndg = np.lib.index_tricks.nd_grid()
    k1,k2,k3 = ndg[0:n1d,0:n1d,0:n1d]-n1d/2
    dk = 2*np.pi/n1d
    k = np.sqrt(k1**2+k2**2+k3**2)
    e = .00001
    A = np.sqrt(-np.log(r)*(k+e)**(-4))
    delta_k = A*np.exp(1j*theta)
    k[np.where(k < 10**(-7))] =1
    psix = -1j*delta_k*k1/k**2
    psiy = -1j*delta_k*k2/k**2
    psiz = -1j*delta_k*k3/k**2
    psi = np.array([np.fft.fftn(psix),np.fft.fftn(psiy),np.fft.fftn(psiz)],dtype = np.complex64)
    turb = np.real(psi)
    mturb = np.mean(turb,0)
    turbspeed = np.sum(turb**2,0)
    rms = np.sqrt(np.mean(turbspeed))
    turb*=.25*225/rms
    
    idx = np.asarray(pos/ 38700 *n1d/2,dtype = np.int) 
    tv = turb[:,idx[:,0],idx[:,1],idx[:,2]] 
    tv = np.transpose(tv)
    tv = np.asarray(tv,dtype = np.float32)
    vel+=tv
    del psi
    sphlib.setVel(ps,arrayToCPointer(sanitize(vel)))
    
    
    
    
    
    
    
def getExtent(ps):
    size = sphlib.getsize(ps)
    pos = np.zeros((size,3),dtype = np.float32)
    sphlib.getPos(ps,arrayToCPointer(sanitize(pos)))
    
    rad = np.sqrt(np.sum(pos*pos,1))
    return np.mean(rad)

def saveState(ps,fn, dt):
    size = sphlib.getsize(ps)
    pos = np.zeros((size,3),dtype = np.float32)
    sphlib.getPos(ps,arrayToCPointer(sanitize(pos)))
    
    vel = np.zeros((size,3),dtype = np.float32)
    sphlib.getVel(ps,arrayToCPointer(sanitize(vel)))
    
    dens = np.zeros(size, dtype = np.float32)
    sphlib.getDens(ps,arrayToCPointer(sanitize(dens)))
    
    temp = np.zeros(size, dtype = np.float32)
    sphlib.getT(ps,arrayToCPointer(sanitize(temp)))
    
    mass = np.zeros(size, dtype = np.float32)
    sphlib.getm(ps,arrayToCPointer(sanitize(mass)))
    
    
    
    np.savez(fn,pos = pos, dens = dens,temp = temp,vel = vel, dt = dt,mass = mass)
    
    
def injectIC(fn):
    l = np.load(fn) 
    pos = l["pos"]
    dens = l["dens"]
    temp = l['temp']
    vel = l['vel']
    dt = l['dt']
    mass = l['mass']
    print len(temp)
    return (sphlib.injectIC(len(temp), arrayToCPointer(sanitize(pos)), arrayToCPointer(sanitize(vel)),arrayToCPointer(sanitize(temp)), arrayToCPointer(sanitize(mass)) ), dt)
    
    
def resumemodel(start,steps,ps,smalldt,totalt):
    for i in range(start,steps):
        print i
        dt = sphlib.domicro(ps,arrayToCPointer(smalldt))
        totalt += dt
        print "dt:  " +str(dt)
        print "microdt: "+str(smalldt)
        print "Total Time: "+str(totalt) +" ka"
        saveState(ps,"./step"+str(i)+".npz",smalldt)
    dt = sphlib.stepforward(ps, 0,1)
    print "dt:  " +str(dt)
    print "microdt: "+str(smalldt)
    print "Total Time: "+str(totalt) +" ka"

def runmodel(steps):
    n1d = 40

    ps = sphlib.getPS(n1d)
    addturbulence(ps)

    
    dt = sphlib.stepforward(ps, 1,0,ct.c_float(0))
    totalt = dt
    print "dt:  " +str(dt)
    smalldt = np.array(dt,dtype = np.float32)
    for i in range(steps):
        print i
        dt = sphlib.domicro(ps,arrayToCPointer(smalldt))
        totalt += dt
        print "dt:  " +str(dt)
        print "microdt: "+str(smalldt)
        print "Total Time: "+str(totalt) +" ka"
        saveState(ps,"./step"+str(i)+".npz",smalldt)
    dt = sphlib.stepforward(ps, 0,1)
    
    saveState(ps, "./lastrun.npz",dt)
    totalt += dt
    print "dt:  " +str(dt)
    
        
    print "Total Time: "+str(totalt) +" ka"
    print "size: "+str(getExtent(ps))
    
def runsingle(steps):
    n1d = 60
    ps = sphlib.getPS(n1d)
    #addturbulence(ps)
    size = sphlib.getsize(ps)
    
    acc = np.zeros((size,3),dtype = np.float32)
    fs = np.zeros(1,dtype = np.float32)
    totalt = 0.0
    dt = sphlib.rts(ps,1,0,arrayToCPointer(sanitize(acc)),arrayToCPointer(sanitize(fs)))
    saveState(ps,"./stepic"+".npz",dt)
    print "dt:  " +str(dt)
    for i in range(steps):
        print i
        
        dt = sphlib.rts(ps,0,0,arrayToCPointer(sanitize(acc)),arrayToCPointer(sanitize(fs)))
        print "dt:  " +str(dt)
        totalt += dt
        print "Total Time: "+str(totalt) +" ka"
        saveState(ps,"./step"+str(i)+".npz",dt)
    


if __name__ == '__main__':    
    n1d = 40
    '''
    l = np.load("lastrun.npz") 
    pos = l["pos"]
    sphlib.setPos(ps,arrayToCPointer(sanitize(pos)))
    
    dt = sphlib.stepforward(ps, 0,0)
    
    
    
    
    ml = np.zeros(size,dtype = np.uint16)
    sphlib.getMicro(ps,ml.ctypes.data_as(ct.POINTER(ct.c_ushort)))
    print np.min(ml)
    print np.max(ml)
    print np.median(ml)
    
    p.hist(ml,bins = 9)'''
    
    
    
    #p.plot(pos[:,0],g[:,0] + acc[:,0],'.')
    #p.plot(pos[:,0],acc[:,0],'.')'''
    #ps, dt = injectIC("step59.npz")
    #resumemodel(60,200,ps,dt,77.602)
    runsingle(100)
    

    
    
   
    p.show()
    
