#quick test of turbulence 

import numpy as np
import matplotlib.pyplot as p
if __name__ == '__main__':
    n1d = 128
    vel = np.zeros((n1d**3,3))
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
    psi = np.array([np.fft.fftn(psix),np.fft.fftn(psiy),np.fft.fftn(psiz)])
    turb = np.real(psi)
    #turb = np.reshape(turb, (3,-1))
    #turb = np.transpose(turb)
    #mturb = np.mean(turb,0)
    #print mturb
    turbspeed = np.sum(turb**2,0)
    print len(turbspeed)
    rms = np.sqrt(np.mean(turbspeed))
    print rms
    turb*=.25*160/rms
    
    pos = np.array([ndg[0:n1d,0:n1d,0:n1d]-n1d/2])
    idx = np.reshape(pos,(3,-1))
    
    idx = np.transpose(idx)
    print idx.shape
    tv = turb[:,idx[:,0],idx[:,1],idx[:,2]]
    print tv[:,0]
    print turb[:,n1d/2,n1d/2,n1d/2]
     
     
    
    