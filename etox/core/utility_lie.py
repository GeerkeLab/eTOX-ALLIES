'''Set of functions for LIE data analysis'''

import numpy as np
from scipy import exp
from numpy import dot
from numpy import seterr
from numpy.linalg import inv as inverse
from numpy import matrix as Matrix
from numpy import mean, sum, subtract
from numpy import transpose, array, sqrt, multiply

from itertools import chain
import logging
import sys
import pybel as pb

from sklearn.linear_model import LinearRegression as linReg

seterr(invalid='print')

'''CALIBRATE'''
def calibrateLie(data,kbt=2.5,limit_iter=20,gamma=False, fixBeta=False):
    '''build model
    input: 
    data= {Y, ene, (smi, decVdw, decEle)} in parenthesis not essential
    modelId= {modelProt, modelProtVer}
    
    ene is list of simulations [[ncpd, nsim, dvdw, dele],..]
    OLD : list of cpd. each compound is a list of simulation. each simulation is 
    '''
    logging.info('Starting calibration of new model')
     
    Y=data['Y']
    ene=data['ene']
 
    '''Get SSY''' 
    avey=np.mean(Y)
    ssy=np.sum([(y-avey)*(y-avey) for y in Y])

    logging.debug("Find best model")
    simsparams=[]

    maxSims=max([len(cpd) for cpd in ene])
    
    for nsims in range(1,maxSims):
        logging.debug("####### model based on %d simulation(s) ########"%nsims)
        results=makemodel(ene,Y,nsims,params=[0.5,0.5,0], gamma=gamma,fixBeta=fixBeta) # results is {'params','ssreg','rmse','r2','Gcalc','wi'}        
        if results['r2']!=0:
            results['sspred'],results['sdep']=calcq2(ene, Y, nsims, results['params'],gamma=gamma,fixBeta=fixBeta)
        else:
            results['sspred'],results['sdep']=[100000,100]

        results['nsims']=nsims

        logging.debug('nsims: %3i  rmse: %6.3f  r2: %6.3f  sdep: %6.3f  q2= %6.3f    y= %5.3f * vdw + %5.3f * ele + %7.3f'% \
        (results['nsims'],results['rmse'],results['r2'],results['sdep'], 1-results['sspred']/ssy, \
         results['params'][0],results['params'][1],results['params'][2]))
        simsparams.append(results)

    lie=sorted(simsparams, key=lambda k: k['sspred'])[0]
    lie['Gexp']=Y

    logging.info("##### best model #####")
    msg="nsimulation=%d params: "%(lie['nsims'])
    for i in lie['params']:
        msg+=" + %.3f"%(i)
    logging.info("%s"%msg)
    for i,y in enumerate(lie['Gexp']):
        msg="Train CpdId:%d   DGexp:%8.3f  DGcalc:%8.3f Error:%8.3f"%(i, y, lie['Gcalc'][i],lie['Gcalc'][i]-y)
        logging.info(msg)

    logging.info('Analysis of the single poses')
    logging.info('cpd_n     '+lie['nsims']*('%3s %5s  '%('sim','Wb')))
    for id in range(len(lie['Gexp'])):
        msg='%5d     '%(id)
        for sim in range(lie['nsims']):
            try:
                msg+='%3d %5.3f  '%(lie['idx'][id][sim],lie['wi'][id][sim])
            except:
                pass
        logging.info(msg)
    return lie


def makemodel(data,dg,nsims,params=[0.5,0.5,0],kbT=2.5,limit_iter=100,gamma=False,fixBeta=False):
    ncycle=0
    while True:
        deltaene,indexsims=filterdene(data,params,maxsims=nsims)           
        #try:
        oldparams=params
        if fixBeta:
            params,rmse,ssreg,r2,Gcalc,weights=create_model_fixedBeta(deltaene,dg,params,kbt=kbT,gamma=gamma)
        else:
            params,rmse,ssreg,r2,Gcalc,weights=create_model(deltaene,dg,params,kbt=kbT,gamma=gamma) #create_model_fixedBeta

        #except Exception,e:
        #    exc_type, exc_obj, exc_tb = sys.exc_info()
        #    ln=str(exc_tb.tb_lineno)
        #    print str(e),ln
        #    return [params,100000,100,0,[],[]]
        
        ncycle+=1

        deltaReg=sum([coef**2 for coef in subtract(params,oldparams)])

        if (deltaReg < 1.e-10):
            results={
                     'params':params,
                     'ssreg':ssreg,
                     'rmse':rmse,
                     'r2': r2,
                     'Gcalc':Gcalc,
                     'wi':weights,
                     'idx':indexsims,
                     }
            break

        if ncycle>limit_iter:
            results={
                     'params':params,
                     'ssreg':100000,
                     'rmse':100,
                     'r2': 0,
                     'Gcalc':[0 for x in dg],
                     'wi':[0 for x in dg],
                     'idx':[[0 for sim in range(nsims)] for cpd in dg]
                     }
            logging.debug("model with %d simulations not converged"%nsims)
            break

    return results


def calcq2(data, dg, nsims, params=[.5,.5,0],kbT=2.5,gamma=False,fixBeta=False):
    '''Cross validation with leave one out method'''
    sspred=0
    cpdlist=list(set([x[0] for x in data]))  

    for ncpd,cpd in enumerate(cpdlist):
        idtest=cpd
        idtrain=range(0,ncpd)+range(ncpd+1,len(cpdlist))

        train=[ x for x in data if x[0]!=idtest ]
        trainy=[ dg[j] for j in idtrain ]
    
        test=[x for x in data if x[0]==idtest]
        testy=[dg[ncpd]]
              
        #try:
        results=makemodel(train,trainy,nsims,params=params, gamma=gamma, fixBeta=fixBeta)
        testDene,testindexes=filterdene(test,results['params'],maxsims=nsims)

        Gcalc, weigths, sdep, sspredcpd=computeDG(testDene,results['params'],kbT,testy)
        sspred+=sspredcpd
            
        '''except Exception, e: #np.linalg.linalg.LinAlgError as err:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            ln=str(exc_tb.tb_lineno)
            print Exception, str(e), ln
            #if 'Singular matrix' in err.message:
            return params,[],1000000,0,0,10000,len(dg)
            #else:'''
    
    sdep=sqrt(sspred/len(cpdlist))
    
    return sspred, sdep


'''PREDICT'''
def predictLie (data, params, kbT=2.5):
    """Runs the prediction protocol. 
        a molecule is described by: energies from MD;
        smile from sdf structure
        model parameters are also essential
        
    Output is a dictionary with calculated DG and a value for the confidential interval
     
    """
    # 1. filter simulations
    dEne,indexes=filterdene(data['ene'],params['params'],maxsims=params['nsims'])
    # 2. perform prediction DG
    Gcalc,wi, sdep, sspred=computeDG(dEne,params['params'],kbT)

    results=(Gcalc,indexes,wi,sdep,sspred)

    return results

def predictError(nDevs,sdep):
    if nDevs < 2:
        predErr=sdep
    elif nDevs == 2:
        predErr=2*sdep
    elif nDevs > 2:
        predErr=3*sdep
        
    return predErr


'''MISCELLANEA'''

def filterdene(data,params, maxsims=None,list_cpd=False):
    '''
    params is a list of lie parameters. The last one is the intercept.
    '''
    dsetene=[]
    indexes=[]
    params=params[:-1]
    if not list_cpd:
        list_cpd=sorted(list(set([sim[0] for sim in data])))
    else:
        list_cpd=sorted(list(set(list_cpd)))

    #for each compound in the training set remove high energy conformations
    for cpd in list_cpd:
        molene=[sim[2:] for sim in data if sim[0]==cpd]
        molindex=[sim[1] for sim in data if sim[0]==cpd]
        
        if maxsims and (len(molene) > maxsims) :
            #estimation of the energy for the single simulations in proteins           
            enesim=dot(molene,params)
            #order of the simulations in proteins according to their estimated energy
            idxenesim=sorted( range(len(molene)), key=lambda k: enesim[k], reverse=False )
            idxenesim=idxenesim[:maxsims]
            
        else:
            idxenesim=range(len(molene))

        dsetene.append([molene[i] for i in idxenesim])
        indexes.append([molindex[i] for i in idxenesim])

    return dsetene,indexes


def create_model(data,Gexp,Params,kbt=2.5, gamma=False):
    """
    Regression routine: it allowed for "infinite" number of parameters. gamma represents the intercept

    Input:
    data=array for each cpd, for each sim [[[dvdw,dele],[]]]
    Gexp=[dg1, dg2, ...]
    Params=[a,b,g]  can be more !!!!
    kbt=2.5 kJ/mol (k Boltzmann * temperature)
    """
    params=np.array(Params[:-1])
    g=Params[-1]
    rmsq = None
    
    NIter       = 0
    boConverged = False
    while True:
        ## reweight energies
        energy, wi = weightEnergy(data,params,g,kbt)
        ## linear regression
        lr=linReg(fit_intercept=gamma)
        lr.fit(energy,Gexp)

        oldparams=params  #save the old values
        oldg=g
        
        params=lr.coef_
        g=lr.intercept_

        deltaReg=sum([coef**2 for coef in subtract(params,oldparams)])+(g-oldg)**2
        if(deltaReg < 1.e-10):
            boConverged = True
            break             # get out of the infinite while loop

        elif(NIter>99):
            break

        NIter += 1

    if boConverged:
        Gcalc= list(lr.predict(energy))
        params=list(params)
        params.append(g)
        ssreg=sum([err**2 for err in subtract(Gexp,Gcalc)])
        rmsq=sqrt(ssreg/len(Gexp))
        r2=lr.score(energy,Gexp)
        
    else:
        logging.error("Regression not converged")
        Gcalc=[0 for x in Gexp]
        params=list(params)
        params.append(g)
        ssreg=np.sum([x**2 for x in Gexp])
        rmsq=sqrt(ssreg/len(Gexp))
        r2=0
        wi=Gcalc
    
    return params,rmsq,ssreg,r2,Gcalc,wi


def create_model_fixedBeta(data,Gexp,Params,kbt=2.5, gamma=False,beta=0.5):
    """
    Regression routine: it allowed for "infinite" number of parameters. gamma represents the intercept

    Input:
    data=array for each cpd, for each sim [[[dvdw,dele],[]]]
    Gexp=[dg1, dg2, ...]
    Params=[a,b,g]
    kbt=2.5 kJ/mol (k Boltzmann * temperature)
    """
    alpha=Params[0]
    g=Params[-1]
    rmsq = None
    
    logging.info('START CREATE MODEL')
    logging.info(alpha, beta, g)
    
    NIter       = 0
    boConverged = False
    while True:
        ## reweight energies
        energy, wi = weightEnergy(data,[alpha,beta],g,kbt)
        #print 'ENERGY', energy
        energyMod=np.array([[ene[0]] for ene in energy ])
        #print 'ENERGYMOD', energyMod
        GexpMod=[ cpd[0]-(cpd[1][1]*beta)  for cpd in zip(Gexp,energy)]
        
        ## linear regression
        lr=linReg(fit_intercept=gamma)
        lr.fit(energyMod,GexpMod)

        oldparams=[alpha]  #save the old values
        oldg=g
        
        params=lr.coef_
        g=lr.intercept_
        alpha=params[0]
        
        deltaReg=sum([coef**2 for coef in subtract(params,oldparams)])+(g-oldg)**2
        print deltaReg, params, g, oldparams, oldg
        if(deltaReg < 1.e-10):
            boConverged = True
            break             # get out of the infinite while loop

        elif(NIter>99):
            break

        NIter += 1

    if boConverged:
        print 'CONVERGED'
        # calc Gcalc
        GcalcMod= list(lr.predict(energyMod))
        print GcalcMod
        Gcalc=[ cpd[0]+(cpd[1][1]*beta)  for cpd in zip(GcalcMod,energy)]
        print Gcalc
        params=list(params)
        params.append(beta)
        params.append(g)
        print params
        
        print len(Gexp), Gexp
        print len(GexpMod), GexpMod
        print 'SUBMIT'
        print subtract(Gexp,Gcalc)
        print [err**2 for err in subtract(Gexp,Gcalc)]
        ssreg=sum([err**2 for err in subtract(Gexp,Gcalc)])
        print ssreg
        rmsq=sqrt(ssreg/len(Gexp))
        print rmsq
        #r2=lr.score(energy,Gexp)
        r=np.corrcoef(Gcalc, Gexp)[0, 1]
        r2=r*r
        print ssreg, rmsq, r2
        
    else:
        logging.error("Regression not converged")
        Gcalc=[0 for x in Gexp]
        params=list(params)
        params.append(g)
        ssreg=np.sum([x**2 for x in Gexp])
        rmsq=sqrt(ssreg/len(Gexp))
        r2=0
        wi=Gcalc
    
    print 'END CREATE MODEL'
    return params,rmsq,ssreg,r2,Gcalc,wi


def weightEnergy(data,coeffs,intercept,kbt=2.5):
    singleEnes=[]
    weights=[]
    #For each molecule get a single vdw and ele value 
    for mol in data:
        #DG for each simulation
        DGi=[x+intercept for x in dot(mol,coeffs)]
        #exponential
        zi= [exp(-x/kbt) for x in DGi]
        ztot=np.sum(zi)
        # Weights
        wi=[z/ztot for z in zi]
        twi=np.array(wi)[np.newaxis].T
        
        scaledEnes=multiply(mol,twi)
        enemol=sum(scaledEnes,axis=0)
        singleEnes.append(enemol)
        weights.append(wi)
        
    return np.array(singleEnes), weights


def computeDG(data,Params,kbt,Gexp=None):
    """
    main program
    """
    params=np.array(Params[:-1])
    g=Params[-1]
    energy, wi = weightEnergy(data,params,g,kbt)
    Gcalc=[x+g for x in dot(data,params)]  
    Gcalc=[x[0] for x in Gcalc]
    sspred=0
    sdep=0

    if Gexp:
        err2=[err**2 for err in subtract(Gcalc,Gexp)]
        
        sspred=sum(err2)
        sdep=sqrt(sspred/len(Gexp))

    return (Gcalc,wi, sdep, sspred)
