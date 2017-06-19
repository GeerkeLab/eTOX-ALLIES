import sys
import pybel as pb
import numpy as np
from scipy.stats import chi2,norm
from sklearn.preprocessing import scale, StandardScaler
from sklearn.covariance import empirical_covariance, EmpiricalCovariance
from sklearn.decomposition import PCA

import pandas as pd
import logging


'''Standard function for compute AD and IC'''
''' Calibration'''
def calibrateAD(data,siminfo):
    '''calibrate AD. Runs the ADAN-like method; four analysis are:
        1. range Y
        2. chemical similarity (Tanimoto)
        3. distribution Dene
        4. pca on decomposition energies (x2 vdw, ele)

    input: 
    data= {Y, ene, fp, decVdw, decEle} in parenthesis not essential
    siminfo should include {'wi', 'idx'} (weights and index of the simulations used in the model
    
    output:
    AD={Yrange,Tanimoto,Dene,devVdw,devEle}
    '''
    logging.debug('Starting calibration AD')
    AD={}
    
    # 1. Y range
    AD['Yrange']={}
    AD['Yrange']['min']=np.min(data['Y'])
    AD['Yrange']['max']=np.max(data['Y'])
    print AD['Yrange']
    
    # 2. Tanimoto similarity; avedist: average distance in the datase; aveclose: average minimum distances in the dataset
    AD['Tanimoto']={}
    AD['Tanimoto']['smi']=data['smi']
    AD['Tanimoto']['Average'],AD['Tanimoto']['Furthest']=tanimoto(data['smi'])   
    print AD['Tanimoto']
    
    # 3. delta enes distribution
    # trainempcov:trained empirical covariance (scikit-learn object);Xmean2scale: array with the averages for scaling
    AD['Dene']={}
    denefilt=filterEne(data['ene'], siminfo['idx'])
    AD['Dene']['CovMatrix'],AD['Dene']['Xmean'],AD['Dene']['Maxdist']=deneAna(denefilt,plot=False,output=None)
    print AD['Dene']
    
    # 3,4. decomposition analysis 3: vdw; 4:ele
    for ene_c in ['Vdw','Ele']:
        #logging.debug("%s decomposition analysis"%ene_c)
        xene=filterEneDec(data['dec%s'%ene_c], siminfo['idx'], siminfo['wi'])

        AD['dec'+ene_c]=pca_train(xene,smalltrain=False)
        
    return AD


'''Prediction'''
def predictAD(data,preds, params):
    #in case of multiple molecules prediction; consisten with predictDG (and who knows if it will be useful in future)

    # A. get values for all the analysis.
    # 2. chemical similarity analysis (Tanimoto). list for each compound average and highest TS with training set
    similarityScore=tanimoto_pred(data['smi'],params['Tanimoto']['smi'])
    
    # 3. distribution Dene. list, for each simulation [cpdid, simid, mahalanobis_distance]
    denefilt=filterEne(data['ene'],preds['idposes'])
    dist_mahal=testdene(denefilt,params['Dene']['CovMatrix'],params['Dene']['Xmean'])
    
    # 4. pca on decomposition energies (x2 vdw, ele) data['decVdw(Ele)'] is list of mols, for each mol a sim
    # sim is [ncpd, npose, res1,res2,..,resn]
    distDec={}
    for ene_c in ['Vdw','Ele']:
        contribute='dec%s'%ene_c
        #logging.debug("%s decomposition analysis"%ene_c)
        xene=filterEneDec(data[contribute], preds['idposes'], preds['wi'])
        # xene is list of array, with reweighted residue contributes for each compounds [[ncpd,res1,res2,..,resn],..]
        results=pca_pred(xene,params[contribute])
        # results is a list, for each compound a np.array as [cpdid SD OD ]
        distDec[contribute]=results
        
    # B. for each compound get CI value
    CI_ds=[]
    for cpdId in range(len(preds['DGcalc'])):
        ci={}
        # 1. range Y
        if preds['DGcalc'][cpdId]<params['Yrange']['max'] and preds['DGcalc'][cpdId]> params['Yrange']['min']:
            ci['Yrange']=0
        else:
            ci['Yrange']=1
          
        # 2. chemical similarity (Tanimoto)
        if similarityScore[cpdId][1] >= params['Tanimoto']['Furthest']:
            ci['Tanimoto']=0
        else:
            ci['Tanimoto']=1
    
        # 3. distribution Dene
        distscpd=[x[2] for x in dist_mahal if x[0]==cpdId]
        ci['Dene']=0
        for x in distscpd:
            if x>params['Dene']['Maxdist']:
                ci['Dene']=1
                break
        # 4. pca on decomposition energies (x2 vdw, ele)
        for ene_c in ['Vdw','Ele']:
            contribute='dec%s'%ene_c
            #print distDec[contribute]
            sim=[ x for x in distDec[contribute] if x[0]==cpdId ]
            if len(sim)>1:
                raise Exception('AD prediction, decomposition analysis. Multiple cpds with same ID?!?!')
            sim=sim[0]
            ci[contribute]=0
            if sim[1]>params[contribute]['critSD']:
                ci[contribute]=1
            if sim[2]>params[contribute]['critOD']:
                ci[contribute]=1
        CI_ds.append(ci)

    return CI_ds


'''2. Chemical similarity analysis by Tanimoto score'''
def tanimoto(smiles,plot=False):
    # create list of fingerprints
    fps=[]
    for cpd in smiles:
        pbMol=pb.readstring('smi',str(cpd))
        fps.append(pbMol.calcfp(fptype='maccs'))
    #get array, for each cpd, similarity with others
    simMat=similarity(fps, fps)
    
    '''remove the distance of each molecule with itself (=1, last in sorted)'''
    matSim_corr=[sorted(x,reverse=True)[1:] for x in simMat]
    sim_mean=[np.mean(x) for x in matSim_corr]
    sim_close=[np.max(x) for x in matSim_corr]

    
    avedist=np.mean(sim_mean)
    closest=np.min(sim_close)
    if plot:
        plotsimil(matSim_corr,[x[0] for x in smiles])

    return avedist,closest 


def similarity(test,train):
    #Create similarity matrix
    simMat=[[ fp | fptrain for fptrain in train] for fp in test]
    return simMat


'''3. Energy distribution analysis'''
def filterEne(dene, idx):
    #Select simulation for dene analysis, according to what has been used to train the model
    deneout=[]
    for ncpd,simIdx in enumerate(idx):
        deneout=deneout+[sim for sim in dene if (sim[0]==ncpd and sim[1] in simIdx)]
    return deneout


def deneAna(dene,plot=False,output="dene.png"):
    X=np.array([sim[2:] for sim in dene])
    
    Xmean=np.mean(X,axis=0)
    Xtocov=scale(X, axis=0, with_mean=True, with_std=False, copy=True)

    emp_cov = EmpiricalCovariance(assume_centered=True).fit(Xtocov)
    mahal_dist = emp_cov.mahalanobis(Xtocov)
    max_dist=max(mahal_dist)

    if plot:
        chi2range=[5.991465,max_dist]
        outliers=[]
        #print "outliers in dene: ",
        for i in range(0,len(mahal_dist)):
            if mahal_dist[i]>chi2range[0]:
                outliers.append(i)
            plotdeneana(Xscaled,Xmean,emp_cov,chi2range,intnames,outliers,output,a,b,lc="blue")

    return emp_cov,Xmean, max_dist


'''3. Energy decomposition analysis'''
def filterEneDec(data,idx,wi):
    # indexsdf is the list with the compound's id
    # indexes is an array containing an array for each simulation, each with the index of simulations to take into account 
    # p is an array of arrays with weight for each simulation
    # data is the an array, for each simulation:[[ncpd,simidx,decres1, decres2,...,decresn]]
    outdecEne=[]

    for ncpd,cpdEne in enumerate(data):
        cpdSims=[simEne[2:] for simEne in cpdEne if simEne[1] in idx[ncpd]]
        wiCpd=np.reshape(wi[ncpd],(len(wi[ncpd]),1)) # get second dimension to wi
        cpdSimsScaled=np.multiply(cpdSims,wiCpd)
        sumCpd=np.concatenate(([ncpd],np.sum(cpdSimsScaled,axis=0)),axis=0)   
        outdecEne.append(sumCpd)

    return outdecEne


def pca_train(data,smalltrain=False,outanal=False,outpref='',reslist=None):
    ## PCA on the per-residue decomposed interaction energies
    ## it uses sklearn pca and scaler module
    ## X is the matrix of energies after filtering. For each compound an array with weighted decomposed energy per each residue over all the simulations
    # outanal writte out files with details on the analysis (i.e. P, T, OD,SD)
    cpdid=[x[0] for x in data]
    X=[x[1:] for x in data]
    
    scaler = StandardScaler(with_mean=True, with_std=False).fit(X)
    X_sc=scaler.transform(X)
    
    if outanal:
        # Save averge interaction energies
        X_mean=np.mean(X,axis=0)
        X_mean=np.reshape(X_mean,(len(X_mean),1))     
        X_sc_std=np.std(X_sc,axis=0)
        X_sc_std=np.reshape(X_sc_std,(len(X_sc_std),1))
        if not reslist:
            reslist=range(X_mean.shape[0])   
        df=pd.DataFrame(np.concatenate((X_mean,X_sc_std),axis=1) ,index=reslist)
        df.columns=['Ave','StDev']
        df.to_csv("AverageInteractions_%s.csv"%outpref, '\t',float_format='%7.3f')
        #Save per compound interaction energies
        fnout="%s_compoundsInteraction.dat"%(outpref)  
        df=pd.DataFrame(np.transpose(X), index=reslist)
        df.column=cpdid
        df.to_csv(fnout, '\t',float_format='%7.3f')
        
    pca = PCA()
    pca.fit(X_sc)

    sumxvar=0
    for pc,xvar in enumerate(pca.explained_variance_ratio_):
        if xvar<0.05:
            break
        sumxvar+=xvar

    T=np.array(pca.transform(X_sc))[:,:pc]
    P=pca.components_[:pc]
    sdev=np.std(T, axis=0)

    SD=np.sqrt(np.sum(np.divide(T**2,sdev**2),axis=1))
    OD=np.sqrt(np.sum(np.subtract(np.dot(T,P),X_sc)**2 , axis=1))

    critSD=np.sqrt(chi2.ppf(.95,pc))
    critOD=(np.median(OD**(2./3.)) + mad(OD**(2./3.),0)*norm.ppf(.95))**(3./2.)
    
    if smalltrain:
        critSD=max(SD)
        critOD=max(OD)

    results={'scaler': scaler,
             'pca':pca,
             'n_pc':pc,
             'sumXvar':sumxvar,
             'sdev': sdev,
             'critSD': critSD,
             'critOD': critOD            
             }

    if outanal:
        # Write Scores
        fnout="T_Train%s.csv"%outpref
        pd.DataFrame(T,index=cpdid).to_csv(fnout, '\t',float_format='%7.3f')
        #explained single variable
        E=np.sum(np.subtract(X_sc,np.dot(T,P))**2,axis=0)
        X2=np.sum(X_sc**2,axis=0)
        Q2=np.subtract(1,np.divide(E,X2))
        # Write explained variance for residue
        fnout="%s_PCAexpVar.dat"%outpref
        df=pd.DataFrame(Q2,index=reslist)
        df=df.replace([np.inf,-np.inf], np.nan).dropna()
        df.to_csv(fnout, '\t',float_format='%7.3f')
        #write Loadings
        fnout="P_%s.csv"%outpref
        df=pd.DataFrame(np.transpose(P),index=reslist)
        df.to_csv(fnout, '\t',float_format='%7.3f')
        # write distances from the model
        fnout="%s_PCAdists.dat"%outpref
        dists=np.concatenate( (np.reshape(SD,(SD.shape[0],1)), np.reshape(OD,(OD.shape[0],1))),axis=1 )
        df=pd.DataFrame(dists,index=cpdid)
        df.columns=['SD','OD']
        df.to_csv(fnout, '\t',float_format='%7.3f')

    return results



'''PREDICT'''
'''2. Chemical similarity analysis by Tanimoto score'''
def tanimoto_pred(smiles_test,smiles_train):
    # Return a list of list in which for every compound an average 
    # and maximum tanimoto score is given, based on the similarity with the training set 

    #list of fp for trainingSet
    fpstrain=[]
    for cpd in smiles_train:
        pbMol=pb.readstring('smi',str(cpd))
        fpstrain.append(pbMol.calcfp(fptype='maccs'))
    #list of fp for testSet
    fpstest=[]
    for cpd in smiles_test:
        pbMol=pb.readstring('smi',str(cpd))
        fpstest.append(pbMol.calcfp(fptype='maccs'))    
    
    #get array, for each cpd, similarity with others
    simMat=similarity(fpstest, fpstrain)

    simTestTrain=[ [np.mean(cpd),np.max(cpd)] for cpd in simMat]    
    
    return simTestTrain


'''3. delta energies analysis'''
def testdene(dene,covtrained,XmeanTr):
    X=np.array([sim[2:] for sim in dene])
    intnames=[sim[:2] for sim in dene]
    
    Xcent=np.subtract(X,XmeanTr)
    mahal_dist = covtrained.mahalanobis(Xcent)

    results=[nm+[mahal_dist[i]] for i,nm in enumerate(intnames)]

    return results


'''4. energy decomposition analysis'''
def mad(data,axis):
    return np.mean(np.absolute(data - np.median(data,axis)), axis)


def pca_pred(data, params,output=False,outpref=None, reslist=None):
    # return list, for each compound a np.array as [cpdid SD OD ]
    # for compatibility among different scikit learn versions 
    if not 'scale_' in dir(params['scaler']):
        params['scaler'].scale_=None
    if not 'std_' in dir(params['scaler']):
        params['scaler'].std_=None    
    cpdid=[x[0] for x in data]
    X=[x[1:] for x in data]  
    X_sc=params['scaler'].transform(X)  
    
    T=params['pca'].transform(X_sc)
    T=np.array(T)[:,:params['n_pc']]

    P=params['pca'].components_[:params['n_pc']]
  
    
    SD=np.sqrt(np.sum(np.divide(T**2,params['sdev']**2),axis=1))   
    OD=np.sqrt(np.sum(np.subtract(np.dot(T,P),X_sc)**2 , axis=1))
       
    distances=np.concatenate( (np.reshape(cpdid,(len(cpdid),1)), np.reshape(SD,(SD.shape[0],1)), np.reshape(OD,(OD.shape[0],1))),axis=1 )
    
    distances=distances.tolist()
    
    if output:
        fnout="%s_compoundsInteraction_pred.dat"%(outpref)  
        df=pd.DataFrame(np.transpose(X), index=reslist)
        df.column=cpdid
        df.to_csv(fnout, '\t',float_format='%7.3f')
        # Write Scores
        fnout="T_Test%s.csv"%outpref
        pd.DataFrame(T,index=cpdid).to_csv(fnout, '\t',float_format='%7.3f')
        #explained single variable
        E=np.sum(np.subtract(X_sc,np.dot(T,P))**2,axis=0)
        X2=np.sum(X_sc**2,axis=0)
        Q2=np.subtract(1,np.divide(E,X2))
        # Write explained variance for residue
        fnout="%s_PCAexpVarTest.dat"%outpref
        df=pd.DataFrame(Q2,index=reslist)
        df=df.replace([np.inf,-np.inf], np.nan).dropna()
        df.to_csv(fnout, '\t',float_format='%7.3f')
        # write distances from the model
        fnout="%s_PCAdistsTest.dat"%outpref
        dists=np.concatenate( (np.reshape(SD,(SD.shape[0],1)), np.reshape(OD,(OD.shape[0],1))),axis=1 )
        df=pd.DataFrame(dists,index=cpdid)
        df.columns=['SD','OD']
        df.to_csv(fnout, '\t',float_format='%7.3f')
        
    return distances


def filterenedec(indexsdf,data,indexes,ene='vdw',distcol=2):
    '''distcol in case more measurement are saved between the enregies ad their decomposition'''
    if ene=='vdw':
        col=4
    elif ene=='ele':
        col=3
    else:
        raise Exception, 'Error: %s energy contribute not recognized'%ene

    yene=[]
    xene=[]
    names=[]

    for i,idxs in enumerate(indexes):
        yene+=[sim[col] for sim in data if ((sim[0]==indexsdf[i]) and (sim[1]-1 in idxs))]
        xene+=[sim[col+distcol] for sim in data if ((sim[0]==indexsdf[i]) and (sim[1]-1 in idxs))]
        names+=[sim[0] for sim in data if ((sim[0]==indexsdf[i]) and (sim[1]-1 in idxs))]
    return yene,xene,names


def energydec(X,Y,namemol,ene):
    '''uses scikit-learn module to perform the regression'''
    plsresults=[]
    for pc in range(1,15):
        plsresults.append(pls(X,Y,pc)) #results: [ncomp,sspred,r2,q2,ycalc,Xmean,regr]

    optilv=sorted(plsresults, key=itemgetter(1))[0]
    logging.debug("Optimal number of LV:%i"%optilv[0])

    residuals=[abs(x[0]-x[1]) for x in zip(Y,optilv[4])]
    cutoff=np.max(residuals)

    return optilv[0],optilv[5],optilv[6],cutoff # [nLV,Xmean,regr,cutoff]  