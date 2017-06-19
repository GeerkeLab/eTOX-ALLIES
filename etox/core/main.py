import os
import sys
import stat
import importlib

import logging
import pybel as pb
import openbabel as ob
import time
import numpy as np

from eTOXlie.etox import settings
from eTOXlie.etox.core import jobHandler
from eTOXlie.etox.core.modelHandler import loadModel, saveModel,listModels
from eTOXlie.etox.core.docking import dockLie
from eTOXlie.etox.core.MDrunner import *
from eTOXlie.etox.core.utility_lie import calibrateLie, predictLie, predictError
from eTOXlie.etox.core.utility_AD import calibrateAD, predictAD

etoxlie_folder=settings.get('etoxlie_folder')
modelDir=settings.get('modelDir')
statEndJob=['FAILED','CANCELLED','DONE']

def protonate(sdf, fileOut, pHcorr=True, pH=7.5):
    '''If necessary converts the 2D structure of the molecule to 3D and protonation by pybel/openbabel
    input:
        sdf: sdf file as string
        neutral: boolean; if molecule should be protonate according to pH is False
        pH: pH for protonation 
    the result is a tuple containing:
        1)  True/False: describes the success of the operation
        2)  (if True ) The name of the 3D molecule
            (if False) The error message
    '''

    ## Check if 3D structure present
    try:
        mol = pb.readstring("sdf",str(sdf)) #    mol = pb.readstring("sdf",str(sdf),opt={'canonical':None})
        if not mol.OBMol.Has3D():
            logging.debug("3D coordinates not found. Creation in progress.")
            mol.make3D(forcefield='mmff94', steps=1000)
    except Exception,e:
            logging.debug("Error during 3D optimization: %s"%e)
            return (False, '3D optimization error')

    ## Protonate according to model (
    try:   
        mol.OBMol.AddHydrogens(False, pHcorr, pH) # OBMol.AddHydrogens(OnlyPolar, Correct4pH, pH)
    except Exception,e:
            logging.debug("Error during protonation: %s"%e)
            return (False, 'Protonation error')

    ## Final optimization
    try:
        mol.localopt(forcefield='mmff94', steps=1000)
        mol.write(format='sdf', filename=fileOut, overwrite=True)
    except Exception,e:
            logging.debug("Error during saving optimized sdf: %s"%e)
            return (False, 'Saving error')
    
    return (True,fileOut)


def calculateInteractions(job, localExe=True, killer=None, sshConn=None, dock_soft='PLANTS',timemd=2.5, clean=False, remoteSettings=None, modelFN='model.dat'):
    logging.info("Calculating interactions for job: %s"%job.filename)
    logging.info("State job %s: %s"%(job.filename,job.status))
    if job.status == 'COMPLETED':
        return
    if killer==None:
        killer = jobHandler.WarmKiller()
    killSig=os.path.join(job.dirTemp,'KILL')
    initialStatus=job.status

    #LOAD MODEL INFO
    success,model=loadModel(os.path.join(modelDir,job.modelProt))
    if not success:
        job.status = 'FAILED'
        job.results = model
        
    # START CHECK/OPERATIONS 
    try:
        if not os.path.isdir(job.dirTemp):
            os.mkdir(job.dirTemp)
    except Exception,e:
        job.status = 'FAILED'
        job.results = "Working directory not found or its creation failed: %s "%(e)
        logging.error('{0} {1}'.format(job.status, job.results))      

    os.chdir(job.dirTemp)

    ### 1. 3D OPTIMIZATION AND PROTONATION   
    baseName="compound"
    if job.status == 'SUBMITTED' and not killer.kill_now:
        logging.debug('3D optimization and protonation in process...')
        success, results = protonate(job.sdf, "%s.sdf"%baseName, model['pHCorr'], model['pH'])
        if success:
            job.status='PROTONATED'
            optiSdf= results
            logging.info('{0} {1}'.format(job.status, results))      
        else:
            job.status='FAILED'
            job.results = "Error during protonation"
            logging.error('{0} {1}'.format(job.status, job.results))      
    else:
        optiSdf="%s.sdf"%baseName

    ### 2. TOPOLOGY CREATION
    suffOut='2MD'
    if job.status == 'PROTONATED' and not killer.kill_now:
        logging.debug('Topology creation in process...')
        topology=importlib.import_module("eTOXlie.etox.topology.%s"%model['forceField'])
        success,results=topology.createTopology(optiSdf,job.dirTemp,fmtIn='sdf',suffOut=suffOut)
        if success:
            job.status = 'TOPOLOGY'
            fileTop = results
        else:
           job.status = 'FAILED'
           job.results = 'Error during topology creation' 
    else:
        itpFn='%s_%s.itp'%(baseName,suffOut)
        CoorFn='%s_%s.pdb'%(baseName,suffOut)
        fileTop={'itp':itpFn,'pdb':CoorFn}


    ### 3. DOCKING
    outMedoidsFn='medoids.dat'
    if job.status == 'TOPOLOGY' and not killer.kill_now :
        #proteinPoses=[]
        logging.debug('Docking in process...')
        for i in range(0,len(model['proteinParams'])):
            wdir=os.path.join(job.dirTemp,'Docking',str(i))
            if not os.path.exists(wdir):
                os.makedirs(wdir)
            proteinParams=model['proteinParams'][i]
            proteinParams['dataDir']=model['dataDir']
            if model['dockSoftware'] in model:
                proteinParams.update(model[model['dockSoftware']])
            success,status=dockLie(fileTop['pdb'],wdir,proteinParams,fmt='pdb',soft=model['dockSoftware'], redCoords=True,algo='kmean',\
                                   outMedoids=outMedoidsFn,killer=killer,killSig=killSig)
            #roteinPoses.append(listPoses)
            if not success:
                job.status = 'FAILED'
                job.results = 'Error during the docking process: %s'%status
                logging.error('{0} {1}'.format(job.status, job.results))      
                break
        if success:
            job.status = status
            logging.info('{0} Docking success'.format(job.status))

    ### 4. PREPARE MD FILES    
    enePref='MD_ene'
    decPref='decomp'
    wdir=os.path.join(job.dirTemp,'MD')
    if job.status == 'DOCKING' and not killer.kill_now:
        logging.debug('MD preparation in process')
        
        # If not directly from docking, get files with poses to run 
        # Check if time simulation defined in the settings: only for debug, if not set in the model
        if settings.get('timeSim'):
            model['timeSim'] = settings.get('timeSim')
            
        if 'proteinPoses' not in locals():
            proteinPoses=[]
            for i in range(0,len(model['proteinParams'])):
                poseDir=os.path.join(job.dirTemp,'Docking',str(i))
                proteinPoses.append(os.path.join(poseDir,outMedoidsFn))

        if not os.path.exists(wdir):
            os.makedirs(wdir)

        success,results=prepareMD(wdir,os.path.join(job.dirTemp,fileTop['itp']),os.path.join(job.dirTemp,fileTop['pdb']),\
                                  proteinPoses,model,localExe=localExe,remoteSettings=remoteSettings,decPref=decPref,enePref=enePref)    
        if success:
            job.status = 'MD READY'
            logging.info('{0} {1}'.format(job.status, results))      

        else:
            logging.error(results)
            job.status = 'FAILED'
            job.results  = 'Error during setting up MDs'
            logging.error('{0} {1}'.format(job.status, job.results))


    ### 5. SUBMIT/RUN MDs
    if job.status == 'MD READY' and not killer.kill_now:
        logging.debug("Start MD running for jobid %s"%job.filename)
        success,results = runMD(wdir,localExe=localExe,sshConn=sshConn,remoteSettings=remoteSettings ,killer=killer, killSig=killSig)
        logging.debug("MD run for jobid %s"%job.filename)
        if success:
            job.status=results
            logging.info('{0}'.format(job.results))      
        else:
            job.status = 'FAILED'
            job.results=results
            logging.error('{0} {1}'.format(job.status, job.results))
        
    ### 5. CHECKMD (CHECK FOR PARALLEL EXECUTION, COPY FILES IN MAIN FOLDER FOR LOCAL)
    if job.status == 'MD SUBMITTED' and not killer.kill_now:
        success,results=checkMD(wdir, job.dirTemp, localExe=localExe,sshConn=sshConn,remoteSettings=remoteSettings,enePref='MD_ene', decPref='decomp')
        if success:
            job.status=results
            logging.info('{0}'.format(job.status))
        else:
            job.status='FAILED'
            job.results=results
            logging.error('{0} {1}'.format(job.status, job.results))      

    ### 6. GATHERED ENERGIES AND SAVE RESULTS
    if job.status == 'MD COMPLETED' and not killer.kill_now:
        success,results=gatherEnergies(wdir=job.dirTemp,enePref="MD_ene",decPref="decomp",ext='.ene',limittime=None)
        ### Change permissions files to open to other group users ( for eTOXsys ) 
        job.open2Group()
        if success:
            job.status='DONE'
            job.results=results
            logging.info('Gathered energies from MD trajectory {0}'.format(job.status))      

        else:
            job.status='FAILED'
            job.results=results
            logging.error('Unable to gather energies from MD trajectory{0}'.format(job.status))      

    ### 7. UPDATE JOB STATUS
    if job.status != initialStatus:
        job.update()

    return job


def screenSDF(listCpds,prediction, modelProt, modelProtVer, misc):
    # Check all the processes of a sdf screening, update status job and status single cpds runs in results
    updated=False
    status='SUBMITTED'
    try:
        completed=True
        for nc, cpd in enumerate(listCpds):           
            # Check status single cpd (only if not already finished)
            if cpd['Status'] not in statEndJob:
                completed=False
                jobCpd=jobHandler.jobCpd()
                jobCpd.load(os.path.join(modelDir,cpd['JobName']))
                #Update progress
                if cpd['Status'] != jobCpd.status:
                    updated=True
                    listCpds[nc]['Status']=jobCpd.status
                    listCpds[nc]['Results']=jobCpd.results
                    if jobCpd.status=='DONE':
                        if prediction:
                            ### PREDICT COMPOUND AND UPDATE RESULTS WITH PREDICTION
                            #here function for pred and AD. Results is a list of dictionaries containing energies for each simulation
                            data=getEnergies([listCpds[nc]])
                            success,results=loadModel(os.path.join(modelDir,modelProt),prediction=True, verModel=modelProtVer)
                            if not success:
                                print results
                                listCpds[nc]['Status']='FAILED'
                                #raise Exception, results
                                continue
                            else:
                                modelLIE,params=results
                            try:
                                Gcalc,idposes, wi, sdep, sspred=predictLie(data,params['LIE'])
                                listCpds[nc]
                                listCpds[nc]['DGexp']=None
                                listCpds[nc]['DGcalc']=Gcalc[0]
                                listCpds[nc]['idpose']=idposes[0]
                                listCpds[nc]['wi']=wi[0]
                                
                                CI=predictAD(data,{'DGcalc':Gcalc, 'idposes':idposes, 'wi': wi},params['AD'])
                                listCpds[nc]['CI_analysis']=CI[0]
                                listCpds[nc]['CI']=sum(CI[0].values())
                                listCpds[nc]['Err']=predictError(sum(CI[0].values()),params['LIE']['sdep'])
                                logging.info('RESULTS:') 
                                for item in listCpds[nc]:
                                    logging.info("{0} {1}".format(item, listCpds[nc][item]))
                                    
                                logging.info('END RESULTS')
                            except Exception, e:
                                logging.error('Error in computing energy for job %s: %s'%(cpd['JobName'],e))
                                listCpds[nc]['Status']='FAILED'
                                continue

        if completed:
            logging.info('Job completed')
            if not prediction:
                # check type of calibration for backward compatibility:
                if misc==None:
                    misc={'isGamma':False,'fixBeta':False}
                # Here collect energies
                print 'GET ENERGIES'
                data= getEnergies(listCpds)  # results is a dictionary with {ids, Y, ene, smi, decVdw, decEle}
                logging.info('Starting calibration for protein %s, version %s'%(modelProt,modelProtVer))
                # Here create model
                print 'CALIBRATE LIE'
                LIEmodel=calibrateLie(data,gamma=misc['isGamma'],fixBeta=misc['fixBeta'])
                print 'CALIBRATE AD'
                ADmodel=calibrateAD(data,LIEmodel)
                modelId={'modelProt':modelProt,'modelProtVer':'%04d'%modelProtVer}
                print 'SAVE MODEL'
                success,out=saveModel(modelDir,modelId,LIEmodel,ADmodel)
                if not success:
                    status='FAILED'
                    raise Exception, out
                
                logging.info('New model saved on file %s'%out)
                deltaId=0
                for n,id in enumerate(data['ids']):
                    if listCpds[id+deltaId]['Status']=='DONE':
                        listCpds[id+deltaId]['DGcalc']=LIEmodel['Gcalc'][n]
                        listCpds[id+deltaId]['idpose']=LIEmodel['idx'][n]
                        listCpds[id+deltaId]['wi']=LIEmodel['wi'][n]
                    else:
                        deltaId+=1
                        
            status='DONE'
            updated=True

        results=listCpds

    except Exception,e:
        logging.debug("JOB FAILED: %s"%e)
        updated=True
        status='FAILED'
        results=e
    
    return (updated, status, results)


def getEnergies(listCpds):
    '''From list of compounds of the job, get:
    ene,dec,Y,fp
    
    ene is list of simulations. for simulation [ncpd, index, vdw, ele]
    decvdw and decene are list of compounds. for each compounds list of simulations. for simulation [ncpd, index, vdwres1, vdwres2, ...]
    
    '''
    logging.debug("Gathering energies for the list of compounds")
    ene=[]
    decVdw=[]
    decEle=[]
    Y=[]
    smi=[]
    idcpd=[]

    jobCpd=jobHandler.jobCpd()
    
    for ncpd, cpd in enumerate(listCpds):
        # Exclude FAILED jobs
        if cpd['Status']=='DONE':
            logging.debug("Loading data for compound: %s"%cpd['JobName'])
            # 1. Get Y
            if cpd['DGexp'] is not None:
                Y.append(float(cpd['DGexp']))
            else:
                Y.append(None)
            
            # 2. Get fingerprint
            smi.append(cpd['smi'])
            
            # 3 Load energies from job
            jobCpd.load(os.path.join(modelDir,cpd['JobName']))       
            energies=jobCpd.results
    
            
            energies=sorted(jobCpd.results, key=lambda k: k['index'], reverse=False)
    
            # GET ENERGIES      
            listContr=['vdw','ele']
            # 4. Get interaction energies
            ene0={}
            for contribute in listContr:
                ene0[contribute] = [ x[contribute] for x in energies if ('%03d'%x['index'])[0] == '0' ]
                if len(ene0[contribute])>1:
                    ene0[contribute]=np.mean(ene0[contribute])
       
            indeces= list(set([ ('%03d'%x['index'])[:-1] for x in energies if ('%03d'%x['index'])[0] != '0' ]))
            for idx in indeces:
                enepose=[[sim[contribute]-ene0[contribute] for contribute in listContr]   for sim in energies if ('%03d'%sim['index'])[:-1] == idx ]
                if len(enepose)>1:
                    enepose=np.mean(enepose,axis=0)
                elif len(enepose)==1:
                    enepose=enepose[0]
                ene.append([ncpd,int(idx)]+list(enepose))
    
            # 5. Get decomposed interaction energies
            deceneDict={}
            for contribute in listContr:
                contribute='dec_%s'%contribute      
                deceneDict[contribute]=[]
                for idx in indeces:
                    decpose={( int(k) if k != 'rest' else str(k)):[] for k in sim[contribute] for sim in energies if ('%03d'%sim['index'])[:-1]==idx }
                    for sim in energies:
                        if ('%03d'%sim['index'])[:-1]==idx:
                            for k,v in sim[contribute].items():
                                try:
                                    k=int(k)
                                except Exception,e :
                                    k=str(k)
                                decpose[k].append(v)
                    decpose={k:np.mean(v) for k,v in decpose.items()}               
                    deceneDict[contribute].append([ncpd,int(idx)]+[decpose[x] for x in sorted(decpose)]) 
                    #if ('%03d'%sim['index'])[0] != '0' :
                    #    print ('%03d'%x['index'])[0]
                    #    simdec={( int(k) if k != 'rest' else str(k)):float(v) for k,v in sim[contribute].items()}
                    #    print simdec
                    #    deceneDict[contribute].append(   [ncpd,sim['index']]+[simdec[x] for x in sorted(simdec)]   )           
    
            decVdw.append(deceneDict['dec_vdw'])
            decEle.append(deceneDict['dec_ele'])
            idcpd.append(ncpd)
        
    results={
             'ids': idcpd,
             'Y': Y,
             'smi': smi,
             'ene': ene,
             'decVdw': decVdw,
             'decEle': decEle
             }

    return results


def submitCPD(strSdf,modelProt,modelProtVer,prediction,fieldExp='',jobid=None):
    # From a sdf submitted as string, lunch or recover jobs
    results=[]
    logging.debug('Start submitting compounds simulations')
    
    try:
        listCpds=jobHandler.collectJobs(etoxlie_folder, cpd=True, sort=True)
        #create list of pybel molecules
        #Workaround to read multiple molecules from string (not default in pybel)
        '''OBMol=ob.OBMol()
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("sdf", "sdf")
        conv.ReadString(OBMol, strSdf)
        sdfMols=[]
        sdfMols.append(pb.Molecule(ob.OBMol(OBMol)))
        while conv.Read(OBMol):
            sdfMols.append(pb.Molecule(ob.OBMol(OBMol)))'''
        sdfMols=[]
        molsdf=''
        for line in strSdf.splitlines():
            molsdf='%s%s\n'%(molsdf,line)
            if line=='$$$$':
                sdfMols.append(molsdf)
                molsdf=''
   
        for nm, mol in enumerate(sdfMols):
            pbMol=pb.readstring('sdf',mol)
            molDict={}
            if not prediction:
                dg=pbMol.data[fieldExp]
            else:
                dg=None

            sdf=mol
            smi=pbMol.write('smi').split()[0]
        
            for runningCpd in listCpds:
                ## Check if already submitted = same sdf and protein model
                if sdf == runningCpd.sdf and modelProt == runningCpd.modelProt:
                    molDict['JobName']=runningCpd.filename
                    molDict['Status']='SUBMITTED' #runningCpd.status
                    #molDict['Results']=runningCpd.results
                    logging.debug("job for molecule %d:%s already processed id %s"%(nm,smi,runningCpd.filename))

            if not bool(molDict):
                newJob=jobHandler.jobCpd()
                newJob.create(sdf,modelProt,modelProtVer,jobid=jobid)
                molDict['JobName']=newJob.filename
                molDict['Status']=newJob.status
                logging.debug("job for molecule %d:%s started with id %s"%(nm,smi,newJob.filename))

            molDict['smi']=smi
            molDict['DGexp']=dg
            results.append(molDict)

    except Exception, e:
        return (False,e)

    return (True,results)


def submitScreen(sdfFn,modelId,prediction,etoxlie_folder='.',fieldExp='Activity', misc=None, jobid=None):
    # sdfFn: filename sdf
    # modelId: {'modelProt': 1A2; 'modelProtVer':1}
    # prediction: True or False
    # in case of prediction=False: misc={'isGamma':True|False,'fixBeta':True|False} 
    # 1. check sdf
    status='SUBMITTED'
    if os.path.exists(sdfFn):
        with open(sdfFn,'r') as infile:
            sdf=infile.read()

    else:
        results="File not found: %s"%sdfFn
        status='FAILED'
        sdf=''

    # 2 Checking consistency of work and models       
    if status != 'FAILED':
        # A. check other runs, in case start
        listJobs=jobHandler.collectJobs(cpd=False, sort=True)
        for job in listJobs:
            if job.status != 'FAILED':
                if sdf==job.sdf:
                    if modelId['modelProt']==job.modelProt:
                        if modelId['modelProtVer']==job.modelProtVer:
                            if prediction==job.prediction:
                                if not prediction:
                                    if misc==job.experiment:
                                        status='FAILED'
                                        results='SDF screen is already in process (jobID: %s)'%job.filename

    if status != 'FAILED':
        # B. check model in case quit
        success,results=loadModel(os.path.join(modelDir,modelId['modelProt']),prediction=prediction, verModel=modelId['modelProtVer'])
        if success:
            if prediction:
                model,params=results
            else:
                model=results
        else:
            results=results
            status='FAILED'

    # 3 Create job
    if status != 'FAILED':
        jobSDR=jobHandler.jobCpd()
        jobSDR.create(sdf,modelId['modelProt'],modelId['modelProtVer'],prediction=prediction,jobid=jobid)
        jobSDR.sdf=sdf
        jobSDR.experiment=misc
    
    # 4 Now submit jobs (if already submitted, link to the existing
    if status != 'FAILED':    
        success,results=submitCPD(sdf,modelId['modelProt'],modelId['modelProtVer'],prediction,fieldExp)
        if success:
            jobSDR.results=results
        else:
            jobSDR.status='FAILED'
            jobSDR.results=results
        
        jobSDR.update()
        
    if status =='FAILED':
        return (False,results)
    else:
        return (True,jobSDR.filename)


def updateJobsStatus():
    listJobs=jobHandler.collectJobs(cpd=False, sort=True)
    for job in listJobs:
        if job.status not in statEndJob:
            logging.info('Job %s: %s'%(job.filename,job.status))#, job.prediction, job.modelProt, job.modelProtVer        
            if not job.prediction:
                # Check if model version already exists, in case change version number
                availModels=listModels(modelDir)
                for m in availModels:
                    if m['name']==job.modelProt:
                        break
                if job.modelProtVer in m['vers']:
                    job.modelProtVer=max(m['vers'])+1
                    job.update()
            updated, status,results=screenSDF(job.results, job.prediction, job.modelProt, job.modelProtVer, job.experiment)
            if updated:               
                oldresults=job.results
                job.status=status
                job.results=results
                try:
                    job.update()
                except:
                    print "There was an error in updating the job %s"%job.filename
                    print 'updateJobsStatus f()'
                    print 'STATUS:', status
                    print 'RESULTS: ', results, '\n'
                    job.status='FAILED'
                    job.results=oldresults
                    job.update()
                    
    return
