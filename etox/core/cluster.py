'''Clustering tools for the docking project'''

import sys
import math
import glob
import time
import os
import numpy as np
import logging
from sklearn import cluster as skclu
from sklearn.neighbors import NearestNeighbors as skNearestNeighbors

class Pose:
    '''A set of coordinates from a pdb file linked to the filename
         Contains data for clustering analysis'''
    rmsd_max = None    #maximum rmsd, class-wide
    rmsd_min = None    #minimum rmsd, class-wide
    _supported_formats = ['.mol2', '.pdb']

    def __init__(self, filename,noH=True):
        self.filename = filename                            #name of the file from which coordinates were read
        self.format = os.path.splitext(self.filename)[1]    #file format
        if not self.format in Pose._supported_formats:      #not supported
            raise Exception, "%s class only supports these file formats: %s"%(self.__class__.__name__,Pose._supported_formats)
        elif self.format == Pose._supported_formats[0]:     # mol2
            self.coords,self.attypes = self.get_mol2(noH=noH)
        elif self.format == Pose._supported_formats[1]:     # pdb
            self.coords,self.attypes = self.get_pdb(noH=noH) 
        self.n_atoms= len(self.coords)                      #number of atoms
        self.cluster = None                                 #cluster index
        self.calc_nearest_neighbors=False                   #switch to check if it already had its neighbors calculated
        self.nearest_neighbors= set()                       #neighbors
        self.medoid = False                                 #switch to check if medoid of cluster
        self.rmsdict = dict()                               #dictionary with all rmsd's calculated for this pose.    Key is pose object

    def get_pdb(self,center=None,noH=True):
        '''gets coords from pdb file supplied
             returns list with coordinates'''
        if center is None:
            center=[0.0,0.0,0.0]
        f = open(self.filename, 'r')
        xyz = list()
        attype = list()
        for line in f:
            if line.startswith('HETATM'):                   #get all HETATM lines
                atomtype=str(line.split()[2])
                if ( not atomtype.startswith('H') or noH == False ) :
                     atomxyz=line.split()[5:8]              #grab xyz
                     atomxyz= [float(i) for i in atomxyz]   #ensure floats
                     for i in range(0,3):
                            atomxyz[i]=round(atomxyz[i]-center[i],4)
                     xyz.append(atomxyz)                    #add to list of coordinates
                     attype.append(atomtype)
        f.close()
        if len(xyz) != len(attype):
                raise Exception, 'different lenght in coordinates/atom types array'
        return xyz, attype

    def get_mol2(self,center=None,noH=True):
        '''gets coords from mol2 file supplied
             returns list with coordinates'''
        if center is None:
            center=[0.0,0.0,0.0]
        f = open(self.filename, 'r')
        xyz = list()
        attype = list()
        for line in f:
            if line.startswith('@<TRIPOS>MOLECULE'):        #get all ATOM lines
                f.next()
                line = f.next()
                read_atoms = int(line.split()[0])
            elif ('read_atoms' in locals()) and line.startswith('@<TRIPOS>ATOM'):
                while read_atoms:
                    line = f.next()
                    atomtype=str(line.split()[5])
                    if ( atomtype != 'H' or noH == False ) :
                         atomxyz=line.split()[2:5]          #grab xyz
                         atomxyz= [float(i) for i in atomxyz] #ensure floats
                         for i in range(0,3):
                                atomxyz[i]=round(atomxyz[i]-center[i],4)
                         xyz.append(atomxyz)                #add to list of coordinates
                         attype.append(atomtype)         
                    read_atoms -= 1
                break
        f.close()
        if len(xyz) != len(attype):
                raise Exception, 'different lenght in coordinates/atom types array'
        return xyz,attype


###End of Pose### 

def filter_poses(dist,posarray,target=[0.0,0.0,0.0],crit='min'):
    '''Filter_poses(distance, [Pose(),Pose(),...,Pose()],target=[x,y,z],criterium) -> [Pose(),Pose(),...,Pose()]
       Removes poses with certain (min, max or average) distance of target'''

    if crit not in ['min', 'max', 'avg']:
        raise Exception, 'Unsupported filter criterium "%s"!'% crit
    else:
        for pose in posarray:
            if pose.distance(target=target, result=crit) > dist:            #if dist is too large
                posarray.remove(pose)            #remove pose from the array
    if not len(posarray):
        raise Exception, "No poses met the criterium!"


def princomp(posarray,noH=True,scale=False):
    '''executes a principal component analysis on
       the poses obtained from docking. Retrieve a set of poses with reduced variables'''
    matcoor=list()
    for pose in posarray:
        coorpose=list()
        for idx,atom in enumerate(pose.coords):
            if noH==False or (not pose.attypes[idx].startswith('H')):    
                    for coor in atom:
                            coorpose.append(coor)
        matcoor.append(coorpose)
    
    if scale:
        std=np.std(matcoor,axis=0)
        matcoorscaled=np.divide(matcoor,std)
    else:
        matcoorscaled=matcoor
    
    mean=np.mean(matcoorscaled,axis=0)
    matcoorcent=np.subtract(matcoorscaled,mean)
    # 1 variance matrix
    covcoor=np.cov(matcoorscaled,rowvar=0)

    # 2 eigenvector decomposition
    eival,eivec=np.linalg.eig(covcoor)
    T=np.dot(matcoorcent,eivec)
    pc=0
    sumx=0
    for ei in eival:
        if ei/sum(eival) > 0.05:
            sumx=sumx+ei/sum(eival)
            pc=pc+1
        else:
            exit
    logging.info("%i principal components, %.2f %% cumulative explained variance"%(pc,sumx*100))
    for npose in range(0,len(posarray)):
        posarray[npose].redcoords=list()
        for npc in range(0,pc):
            posarray[npose].redcoords.extend([T[npose,npc]])
    return # posarray


def getCoorMatrix(posArray,redCoords=False):
    X=list()
    for pose in posArray:
        coorPose=list()
        if redCoords:
            for coor in pose.redcoords:
                coorPose.append(coor)
        else:
            for atom in pose.coords:
                for coor in atom:
                    coorPose.append(coor)
        X.append(coorPose)
    X=np.array(X)
    return X

### CLUSTER


def clusterMatrix(Xmatrix,algo='kmean',fixedClu=False,nClu=2,saveStat=True):
    ''' Cluster poses from posArray of pose objects. A fixed number of cluster (nClu) can be defined (fixedClu==True)
        otherwise the optimal number of cluster between 1 and nClu will be chosen.
        algorithms:
            ->    kmean
            ->    ward : agglomerative Ward hierchical
            ->    knn: k-nearest neighbors
            ->    rnn: radius-nearest neighbors
        
        Returns a list with the matrix indices of the medoids
        
        '''
    _availAlgos=['kmean','ward','knn','rnn']
    
    #Check algorithm(s) to use
    if algo == 'all':
        tryAlgo=_availAlgos
    elif algo in _availAlgos:
        tryAlgo=[algo]
    else:
        errorMsg='Algorithm selected not available. List of possible algorithms is: %s'%(', '.join(_availAlgos))
        logging.error(errorMsg)
        return (False,errorMsg)
    
    #Define number of cluster to sample
    if fixedClu:
        nCluStart=nClu
        nCluEnd=nClu+1
    else:
        nCluStart=2
        nCluEnd=nClu+1

    #
    # Look for optimal number of clusters
    #
    # Center from medoid 1 cluster
    labels=np.zeros(Xmatrix.shape[0])
    ssy, medoid0, nposes = getStats(Xmatrix,labels)

    listTrials=[]  # list of dicts : varexp, Fcoll, algo, nclu, otherDetails, (medoids )
    for nCl in range(nCluStart,nCluEnd):
        for method in tryAlgo:
            if method == 'kmean':
                labels=kmeanCluster(Xmatrix,nCl)
            if method == 'ward':
                labels=wardCluster(Xmatrix,nCl)
            if method == 'knn':
                continue
                #medoids,clusters=radnnCluster(Xmatrix,nCl)
            if method == 'rnn':
                labels=radnnCluster(Xmatrix,nCl)
    
            ssx,listMedoids, nposes=getStats(Xmatrix,labels)
            listTrials.append({'algo':method,'nClu':nCl,'ssx':1-(ssx/ssy), 'nposes':nposes, 'medoids':listMedoids })

    # For curiosity on cluster efficiency
    optilistTrials=[]        
    for i in range(nCluStart,nCluEnd):
        optilistTrials.append(sorted([x for x in listTrials if x['nClu']==i],key=lambda k: k['ssx'],reverse=True)[0])

    nTrial=0
    while optilistTrials[nTrial+1]['ssx']-optilistTrials[nTrial]['ssx']> 0.05:
#    while optilistTrials[nTrial]['ssx']<0.90:
        nTrial+=1
        break

    ## optimal: optilistTrials[nTrial] 
          
    ## SOME LOGGING
    if saveStat:
        results = open('cluster.fin', 'w')
        results.write('Clustering of %d docking poses by %s algorithm(s)\n'%(Xmatrix.shape[0], algo))
        results.write('Algorithm used in the optimal clustering: %s\n'%(optilistTrials[nTrial]['algo']))
        results.write('Optimal number of cluster obtained: %d\n'%(optilistTrials[nTrial]['nClu']))
        results.write('explained variance: %f\n'%(optilistTrials[nTrial]['ssx']))
    logging.info('optimal number of clusters: %i'%(optilistTrials[nTrial]['nClu']))
    logging.info('explained variance: %f'%(optilistTrials[nTrial]['ssx']))
    logging.info('Algorithm used in the optimal clustering: %s'%(optilistTrials[nTrial]['algo']))
    msg=''

    ### HERE TO IMPLEMENT OTHER OUTPUT
    
    if saveStat:
            results.close()

    return tuple(optilistTrials[nTrial]['medoids'])


def getStats(X,labels,printData=False):
    ''' From array of coordinates and array of cluster per pose gives medoids (medoidList)
        and sum of squared distance from medoid (ssx)'''

    listClu=list(set(labels))
    Xtot=np.concatenate( ( np.arange(0,X.shape[0])[:,None], X, labels[:,None] ), axis=1) # labels[:,None] to give dimensionality
    medoidList=[]
    nposes=[]
    ssx=0
    
    if printData:
        for pose in Xtot:
            logging.info(pose)
    
    for clu in listClu:
        ## count poses
        nposes.append(list(labels).count(clu))
        Xclu = Xtot[ Xtot[:,-1] == clu,0:-1]
        ## get medoid
        center=np.mean(Xclu[:,1:], axis=0)
        dists=[np.linalg.norm(pose[1:]-center) for pose in Xclu]
        medoid=Xclu[np.argmin(dists)]
        medoidList.append(int(medoid[0]))
        ## get stats
        ssx+=np.sum([np.linalg.norm(pose[1:]-medoid[1:])**2 for pose in Xclu])

    return ssx, medoidList, nposes


def kmeanCluster(matrix,nCl):
    '''kmean clustering'''
    clustering=skclu.KMeans(n_clusters=nCl, init='random', n_init=300, max_iter=1000, tol=0.0001, precompute_distances='auto', verbose=0, random_state=0, copy_x=True, n_jobs=1)
    clustering.fit(matrix)
    return clustering.labels_


def wardCluster(matrix,nCl):
    '''ward agglomerative clustering'''
    clustering=skclu.AgglomerativeClustering(n_clusters=nCl, connectivity=None, affinity='euclidean', n_components=None, compute_full_tree='True', linkage='ward')
    clustering.fit(matrix)
    return clustering.labels_


def radnnCluster(matrix,nCl, maxIter=100): 
    ''' nearest neighbors clustering'''
    
    center=np.mean(matrix, axis=0)
    dists=[np.linalg.norm(pose-center) for pose in matrix]    
    dist_cutoff = np.median(dists)
    
    convergency=False
    matrix=np.concatenate( ( np.arange(0,matrix.shape[0])[:,None], matrix), axis=1)
    niter=0
    while not convergency:
        niter+=1
        convergency=True
        X=matrix
        labelClusters=[]
        for Clu in range(0,nCl):
            if X.shape[0] < 1:
                convergency=False
                break

            nn=skNearestNeighbors(radius=dist_cutoff, algorithm='brute', metric='euclidean')
            nn.fit(X[:,1:])
            distances, indices = nn.radius_neighbors(X[:,1:],dist_cutoff)
            indices=np.ndarray.tolist(indices)
            indices.sort(key=len, reverse=True)          
            mask=np.ones(X.shape[0],dtype=bool)
            mask[indices[0]]=False
            labelClusters.append(X[indices[0],0])
            X=X[mask]

        if convergency==True:
            if X.shape[0]>0:
                ## Poses out of clusters
                convergency==False
                dist_cutoff *= 1.05
                
        else:
            ## not enough clusters
            dist_cutoff /= 1.05             #make it smaller
        
        if niter>maxIter:
            return np.zeros(matrix.shape[0])
    ## assign labels
    labels=np.zeros(matrix.shape[0])
    for id,clu in enumerate(labelClusters):
        for pose in clu:
            labels[pose]=id

    logging.debug("End rnn for %d cluster with %d iterations"%(nCl,niter))
    return labels


if __name__ == '__main__':	
    pass