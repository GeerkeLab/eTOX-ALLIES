import numpy as np
from math import sqrt, sin,isnan

def getnormalplane(a,b,c): 
    
    rab=np.subtract(b,a)
    rac=np.subtract(c,a)
    ln=np.cross(rab,rac)
    n=ln/np.linalg.norm(ln)
    
    return n
    

def com(listatm):
    avg=np.average(listatm,axis=0)
    return avg
    
    

def findOcpd1(norm,Fe,dist=1.64):

    d=norm*dist
    rO=np.add(d,Fe)  
    
    return rO


def normhem(a):
    #a is a dictionary with the {'ATOM NAME':[x,y,z]} for the atoms in hem
    # the result is the normal to the plane averaged on 4 planes defined by the hem atoms 
    norms=list()
    norms.append(getnormalplane(a['CHA'],a['C2B'],a['C3C']))
    norms.append(getnormalplane(a['C3A'],a['CHC'],a['C2D']))
    norms.append(getnormalplane(a['C2A'],a['C3B'],a['CHD']))
    norms.append(getnormalplane(a['CHB'],a['C2C'],a['C3D']))
    avenorm=np.average(norms,axis=0)
    
    return avenorm

def distance(a,b):

    if len(np.shape(a))>1:
        a=com(a)
    if len(np.shape(b))>1:
        b=com(b)
  
    d=np.linalg.norm(np.subtract(b,a))
    
    return d

def vecangle(a,b,deg=False):
    #1. dot product of the two lines
    p=np.dot(a,b)
    #2.vector magnitude
    am=np.linalg.norm(a)
    bm=np.linalg.norm(b)
    #3.get arccos
    cosalpha=p/(am*bm)
    #4.get angle
    alpha=np.arccos(cosalpha) #*3.14159265359
    
    
    if deg:
        alpha=alpha/3.14*180  
        
    return alpha 

def pointsangle(a,b,c,deg=True):
    #1. create two vectors
    rba=np.subtract(b,a)
    rba=rba/np.linalg.norm(rba)
    
    rbc=np.subtract(b,c)
    rbc=rbc/np.linalg.norm(rbc)

    
    alpha=vecangle(rba,rbc)
    if deg:
        alpha=alpha/3.14159265359*180
    return alpha

def areatriangle(a,b,c):   
    area=.5 * distance(a,b) * distance(a,c) * sin(pointsangle(b,a,c,deg=False))
    
    return area
    