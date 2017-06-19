import os,sys
import logging
import time
import pickle
import argparse
import numpy as np
import pybel as pb
import xlsxwriter

from math import sqrt
from shutil import rmtree,copyfile

from eTOX_ALLIES.etox import settings
from eTOX_ALLIES.etox.core import jobHandler

from settings import modelDir, tmpFolder

def loadTrainInfo(item,modProt,modProtVer,modelDir=settings.get('modelDir')):  
    protVer='%04d'%modProtVer
    pathParams=os.path.join(modelDir,modProt,protVer,'params.pkl')
    if os.path.exists(pathParams):
        with open(pathParams,'r') as inFile:
            parms=pickle.load(inFile)
            
        trainSet=parms[item]

    else:
        return (False,'Model parameters not found (%s)'%pathParams)
    
    return (True,trainSet)


def getPoses(jobDir):
    
    logging.debug('Get starting poses')
    foundPoses=False
    listPoses=[]
    if os.path.isdir(jobDir):
        dockDir=os.path.join(jobDir,'Docking')
        if os.path.isdir(dockDir): 
            subDocks=[x[1] for x in os.walk(dockDir)][0]
            for confDir in subDocks:
                poseFn=os.path.join(dockDir,confDir,'medoids.dat')
                with open(poseFn,'r') as inFile:
                    listPose=inFile.readlines()
                npose=0
                for pose in listPose:
                    poseFn=pose.strip()
                    if pose!='':
                        newFn='%d%d'%(int(confDir)+1,npose)
                    npose+=1
                    listPoses.append({'fileName':poseFn, 'id':newFn})
            foundPoses=True    
                    
    if not foundPoses:
        logging.debug('%s not found'%jobDir)
    
    return (foundPoses,listPoses)


def getFullTrain(trainModel,jobResults):
    datasetInfo=[]
    listPoses=[]
    i=0
    for cpd in jobResults:
        i+=1
        cpdInfo={}
        TSfound=False
        for j in range(len(trainModel)):
            if cpd['smi']==trainModel[j]['smi']:
                TSfound=True
                break
        cpdInfo['id']=i
        cpdInfo['smi']=cpd['smi']
        cpdInfo['DGexp']=float(cpd['DGexp'])
        cpdInfo['JobFile']=cpd['JobName']
        cpdInfo['Status']=cpd['Status']
        
        if TSfound:
            cpdInfo['DGcalc']=float(trainModel[j]['Gcalc'])
            wi=['%.3f'%x for x in trainModel[j]['wi']]
            cpdInfo['wi']=', '.join(wi)
            idposes=['%02d'%x for x in trainModel[j]['idsims'] ]
            cpdInfo['idSims']=', '.join(idposes)

        else:
            cpdInfo['DGcalc']=''
            cpdInfo['wi']=''
            cpdInfo['idSims']=''

        cpdInfo['CI_decEle']=''
        cpdInfo['CI_decVdw']=''
        cpdInfo['CI_Tanimoto']=''
        cpdInfo['CI_Dene']=''
        cpdInfo['CI_Yrange']=''

        datasetInfo.append(cpdInfo)
        
        ### Try to extract starting pose files
        try:
            jobcpd=jobHandler.jobCpd()
            jobcpd.load(cpdInfo['JobFile'])
            completed,cpdPoses=getPoses(jobcpd.dirTemp)
            if completed:
                listPoses.append({'cpd':i,'poses':cpdPoses})
        except Exception, e:
            logging.debug('Error in loading poses files, %s'%e)
        
    return (datasetInfo,listPoses)


def getFullPred(strSdf,results):
    #for cpd: {'Status', 'JobFile', 'DGexp', 'DGcalc', 'wi', 'smi', 'idSims', 'id'}
    ## Try to get DGexp from sdf
    sdfMols=[]
    molsdf=''
    for line in strSdf.splitlines():
        molsdf='%s%s\n'%(molsdf,line)
        if line=='$$$$':
            sdfMols.append(molsdf)
            molsdf=''
    
    pbDS=[]
    listkeys=[]
    for nm, mol in enumerate(sdfMols):
        pbMol=pb.readstring('sdf',str(mol))
        pbDS.append(pbMol)
        listkeys=listkeys+pbMol.data.keys()
    
    lenDS=len(pbDS)
    listkeys=[ x for x in set(listkeys) if listkeys.count(x)==lenDS ]
    
    print '\n\n ATTRIBUTES FOUND IN THE SUBMITTED SDF FILE:\n\n'
    for i,item in enumerate(listkeys):
        print('%d.\t\t\t%s'%(i+1,item)) 
    while True:
        keyId=raw_input('\nSelect the field with experimental DGbind [0 is not present]?')
        try:
            keyId=int(keyId)
            if keyId>=0:
                if keyId<=len(listkeys):
                    break
            raise exception
        except:
            print '\t Invalid choice %s'%str(keyId)
    
    if keyId>0:
        expList=[{'smi':x.write('smi').split()[0], 'DGexp':x.data[listkeys[keyId-1]]} for x in pbDS]
    else:
        expList=[{'smi':x.write('smi').split()[0], 'DGexp':0 } for x in pbDS]
    
    ##### EXPERIMENTS LOADED, BEGIN TO RESUME DATA
    datasetInfo=[]
    listPoses=[]
    
    
    i=0
    for cpd in results:
        i+=1
        cpdInfo={}
        TSfound=False
        for j in range(len(expList)):
            if cpd['smi']==expList[j]['smi']:
                TSfound=True
                break
        cpdInfo['id']=i
        cpdInfo['smi']=cpd['smi']
        try:
            if not TSfound:
                raise exception
            cpdInfo['DGexp']=float(expList[j]['DGexp'])
        except:
            cpdInfo['DGexp']=''

        cpdInfo['JobFile']=cpd['JobName']
        cpdInfo['Status']=cpd['Status']
        try:
            cpdInfo['DGcalc']=float(cpd['DGcalc'])
        except:
            cpdInfo['DGcalc']=''
        
        try:
            wi=['%.3f'%x for x in cpd['wi']]
            cpdInfo['wi']=', '.join(wi)
        except:
            cpdInfo['wi']=''
        
        try:
            idposes=['%02d'%x for x in cpd['idpose'] ]
            cpdInfo['idSims']=', '.join(idposes)
        except:
            cpdInfo['idSims']=''

        try:
            cpdInfo['CI_decEle']=cpd['CI_analysis']['decEle']
            cpdInfo['CI_decVdw']=cpd['CI_analysis']['decVdw']
            cpdInfo['CI_Tanimoto']=cpd['CI_analysis']['Tanimoto']
            cpdInfo['CI_Dene']=cpd['CI_analysis']['Dene']
            cpdInfo['CI_Yrange']=cpd['CI_analysis']['Yrange']
        except:
            cpdInfo['CI_decEle']=''
            cpdInfo['CI_decVdw']=''
            cpdInfo['CI_Tanimoto']=''
            cpdInfo['CI_Dene']=''
            cpdInfo['CI_Yrange']=''
            
            
        datasetInfo.append(cpdInfo)       
        ### Try to extract starting pose files
        try:
            jobcpd=jobHandler.jobCpd()
            jobcpd.load(cpdInfo['JobFile'])
            completed,cpdPoses=getPoses(jobcpd.dirTemp)
            if completed:
                listPoses.append({'cpd':i,'poses':cpdPoses})
        except Exception, e:
            logging.debug('Error in loading poses files, %s'%e)

    return (datasetInfo,listPoses)
    

def getStats(dsData,modelFolder,modelProt,modelProtVer):
    #print trainModel
    stats={}
    ssreg=[(x['DGexp']-x['DGcalc']) for x in dsData if x['DGcalc']!='']
    ssreg=[x*x for x in ssreg]
    stats['RMSE']=sqrt(np.average(ssreg))
    ssreg=np.sum(ssreg)

    # Get ave y from train
    completed,trainMod=loadTrainInfo('trainSet',modelProt,modelProtVer,modelDir=modelFolder)
    trainYave=np.average([float(x['Gexp']) for x in trainMod])  
    ssY=[(x['DGexp']-trainYave) for x in dsData]
    ssY=[x*x for x in ssY]
    ssY=np.sum(ssY)

    ssYnoint=[(x['DGexp']) for x in dsData if x['DGcalc']!='']
    ssYnoint=[x*x for x in ssYnoint]
    ssYnoint=np.sum(ssYnoint)

    stats['r2']=1-ssreg/ssY
    stats['r2noint']=1-ssreg/ssYnoint
    stats['ssreg']=ssreg
    stats['ssY']=ssY
    stats['ssYnoint']=ssYnoint
    stats['ncpd']=len(dsData)

    return stats


def extractPoses(listPoses,dirName,overWrite=True):
    currDir=os.getcwd
    if os.path.exists(dirName):
        if overWrite:
            logging.debug('Removing folder %s'%dirName)
            rmtree(dirName)
        else:
            logging.error('Folder %s already exists: it won\'t be overwritten'%dirName)
    os.mkdir(dirName)
    for cpd in listPoses:
        cpdid=cpd['cpd']
        for pose in cpd['poses']:
            try:
                src=pose['fileName']
                extSrc=os.path.splitext(src)[1]
                dst=os.path.join(dirName,'%03d_%s%s'%(cpdid,pose['id'],extSrc))
                copyfile(src,dst)
            except:
                logging.debug('File %s not copied'%src)
    return


def extractJob(jobId,tmpFolder,modelFolder,outPref=''):
    extFile=os.path.splitext(jobId)[1]
    if extFile=='.dsr':
        dsr=True
    else:
        dsr=False

    jobFn=os.path.join(tmpFolder,jobId)
    if not os.path.exists(jobFn):
        logging.error('JOB ID %s NOT FOUND!!!'%jobFn)
        sys.exit()
    
    job=jobHandler.jobCpd()
    
    
    job.load(jobFn)
    
    if job.prediction:
        typeJob='PREDICTION'
    else:
        typeJob='CALIBRATION'

    if dsr:
        ### Write dataset information
        datasetInfo=True
        if typeJob=='CALIBRATION':
            # LOAD ALL THE INFO ABOUT CALIBRATION JOB
            if job.status=='DONE':
                completed,trainModel=loadTrainInfo('trainSet',job.modelProt,job.modelProtVer)
                if not completed:
                    datasetInfo=False
                    logging.error(trainModel)
                else:
                    datasetData,listPoses=getFullTrain(trainModel,job.results)
        else:
            # LOAD ALL THE INFO ABOUT PREDICTION JOB           
            datasetData,listPoses=getFullPred(job.sdf,job.results)


        if job.status=='DONE':
            stats=getStats(datasetData,modelFolder,job.modelProt,job.modelProtVer)

    if len(listPoses)>0:
        dirNmPoses='Poses_%s'%jobId
        extractPoses(listPoses,dirNmPoses)

   
    ## WRITE OUT INFORMATION
    wb=xlsxwriter.Workbook('Info_%s.xlsx'%jobId)
    boldFmt = wb.add_format({'bold': True})
    ws=wb.add_worksheet()
    
    ws.write(0, 0, 'Filename',boldFmt )
    ws.write(0, 1, job.filename )
    ws.write(0, 2, 'Type',boldFmt )
    ws.write(0, 3, typeJob )
    ws.write(0, 4, 'Status',boldFmt )
    ws.write(0, 5, job.status ) 
    ws.write(1, 0, 'Model Id',boldFmt )
    ws.write(1, 1, job.modelProt )
    ws.write(1, 2, 'Model Version',boldFmt )
    ws.write(1, 3, job.modelProtVer)
    ws.write(2, 0, 'Data Creation',boldFmt )
    ws.write(2, 1, time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(job.dateSub)) )

    if datasetInfo:
        ws.write(4,0,'STATISTICS',boldFmt)
        ws.write(5,0,'n cpds',boldFmt)
        ws.write(5,1,'SDEP',boldFmt)
        ws.write(5,2,'SSreg',boldFmt)
        ws.write(5,3,'SSy',boldFmt)
        ws.write(5,4,'q2',boldFmt)
        ws.write(5,5,'SSy noInt',boldFmt)
        ws.write(5,6,'q2 noInt',boldFmt)
        ####
        if job.status=='DONE':
            ws.write(6,0,stats['ncpd'])
            ws.write(6,1,stats['RMSE'])
            ws.write(6,2,stats['ssreg'])
            ws.write(6,3,stats['ssY'])
            ws.write(6,4,stats['r2'])
            ws.write(6,5,stats['ssYnoint'])
            ws.write(6,6,stats['r2noint'])
        
        ws.write(8,0,'LIST COMPOUNDS',boldFmt)
        ###  stats['ssYnoint']
        ws.write(9,0,'n.',boldFmt)
        ws.write(9,1,'smi',boldFmt)
        ws.write(9,2,'experimental',boldFmt)
        ws.write(9,3,'calculated',boldFmt)
        ws.write(9,4,'file',boldFmt)
        ws.write(9,5,'status',boldFmt)
        ws.write(9,6,'simulations ID',boldFmt)
        ws.write(9,7,'simulations Wi',boldFmt)
        ws.write(9,8,'CI:',boldFmt)
        ws.write(9,9,'Y range',boldFmt)
        ws.write(9,10,'Tanimoto',boldFmt)
        ws.write(9,11,'Delta Ene',boldFmt)
        ws.write(9,12,'decomp Coul',boldFmt)
        ws.write(9,13,'decomp L-J',boldFmt)

        row0=10
        row=10       
        for cpd in datasetData:
            ws.write(row,0,cpd['id'])
            ws.write(row,1,cpd['smi'])
            ws.write(row,2,cpd['DGexp'])
            ws.write(row,3,cpd['DGcalc'])
            ws.write(row,4,cpd['JobFile'])
            ws.write(row,5,cpd['Status'])
            ws.write(row,6,cpd['idSims'])
            ws.write(row,7,cpd['wi'])
            ws.write(row,8,cpd['CI_Yrange']+cpd['CI_Tanimoto']+cpd['CI_Dene']+cpd['CI_decEle']+cpd['CI_decVdw'])
            ws.write(row,9,cpd['CI_Yrange'])
            ws.write(row,10,cpd['CI_Tanimoto'])
            ws.write(row,11,cpd['CI_Dene'])
            ws.write(row,12,cpd['CI_decEle'])
            ws.write(row,13,cpd['CI_decVdw'])
           
            row+=1
   
    wb.close()
    sdfOutFn='Input_%s.sdf'%jobId
    with open(sdfOutFn,'w') as sdfOut:
        sdfOut.write(job.sdf)


def extractModel(modProt,modProtVer,modelDir):
    protVer='%04d'%modProtVer
    pathParams=os.path.join(modelDir,modProt,protVer,'params.pkl')
    try:
        with open(pathParams,'r') as inFile:
            parms=pickle.load(inFile)
        ssreg=parms['LIE']['ssreg']
        sspred=parms['LIE']['sspred']
        lieParams=parms['LIE']['params']
        nsims=parms['LIE']['nsims']
        ver=parms['version']
        date=time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(parms['creation_date']))
        ncpd=len(parms['trainSet'])
        

        listGexp=[float(x['Gexp']) for x in parms['trainSet']]
        aveY=np.average(listGexp)
        ssY=[(x-aveY) for x in listGexp]
        ssY=np.sum([x*x for x in ssY])
        
        ssYnoint=np.sum([x*x for x in listGexp])
        
        r2=1-ssreg/ssY
        r2noInt=1-ssreg/ssYnoint

        sumGexp=np.sum(listGexp)
        lesscpds=ncpd-1
        errYpred=[y-((sumGexp-y)/lesscpds) for y in listGexp]
        ssYpred=np.sum([y*y for y in errYpred])
        
        q2=1-sspred/ssYpred
        q2noInt=1-sspred/ssYnoint
        
        RMSE=sqrt(ssreg/ncpd)
        SDEP=sqrt(sspred/ncpd)
        
        train=parms['trainSet']
 
    except Exception, e:
        logging.error('Error in retrieving model information: %s'%e)
        sys.exit()
           
    
        ## WRITE OUT INFORMATION
    wb=xlsxwriter.Workbook('InfoModel_%s_%03d.xlsx'%(modProt,modProtVer))
    boldFmt = wb.add_format({'bold': True})
    ws=wb.add_worksheet()
    
    ws.write(0, 0, 'Model Id',boldFmt )
    ws.write(0, 1, modProt )
    ws.write(0, 2, 'Model Version',boldFmt )
    ws.write(0, 3, modProtVer )
    ws.write(0, 4, 'LIE Version',boldFmt )
    ws.write(0, 5, ver )    
    ws.write(1, 0, 'Data Creation',boldFmt )
    ws.write(1, 1, date )
    ws.write(2, 0, 'alpha',boldFmt )
    ws.write(2, 1, 'beta',boldFmt )
    ws.write(2, 2, 'gamma',boldFmt )
    ws.write(3, 0, lieParams[0] )
    ws.write(3, 1, lieParams[1] )
    ws.write(3, 2, lieParams[2] )

    ws.write(4,0,'STATISTICS',boldFmt)
    ws.write(5,0,'n cpds',boldFmt)
    ws.write(5,1,'RMSE',boldFmt)
    ws.write(5,2,'SSreg',boldFmt)
    ws.write(5,3,'SSy',boldFmt)
    ws.write(5,4,'r2',boldFmt)
    ws.write(5,5,'SDEP',boldFmt)
    ws.write(5,6,'SSpred',boldFmt)
    ws.write(5,7,'SSy pred',boldFmt)
    ws.write(5,8,'q2',boldFmt)        
    
    ws.write(5,10,'SSy noInt',boldFmt)
    ws.write(5,11,'r2 noInt',boldFmt)
    ws.write(5,12,'q2 noInt',boldFmt)
    
    ws.write(6,0,ncpd)
    ws.write(6,1,RMSE)
    ws.write(6,2,ssreg)
    ws.write(6,3,ssY)
    ws.write(6,4,r2)
    ws.write(6,5,SDEP)
    ws.write(6,6,sspred)
    ws.write(6,7,ssYpred)
    ws.write(6,8,q2)
    
    ws.write(6,10,ssYnoint)
    ws.write(6,11,r2noInt)
    ws.write(6,12,q2noInt)
    
    ws.write(8,0,'LIST COMPOUNDS',boldFmt)
    ###  stats['ssYnoint']
    ws.write(9,0,'n.',boldFmt)
    ws.write(9,1,'smi',boldFmt)
    ws.write(9,2,'experimental',boldFmt)
    ws.write(9,3,'calculated',boldFmt)

    row0=10
    row=10       
    for i,cpd in enumerate(train):
        ws.write(row,0,i+1)
        ws.write(row,1,cpd['smi'])
        ws.write(row,2,cpd['Gexp'])
        ws.write(row,3,cpd['Gcalc'])          
        row+=1
   
    wb.close()



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Extract detailed information about model version or job')

    parser.add_argument('-m', '--model', required=False, dest='model',  help='model name')
    parser.add_argument('-v', '--ver', required=False,    dest='version',  help='model version', default=None)
    parser.add_argument('-j', '--job', required=False, dest='job', help='job name', default=None)
    parser.add_argument('-o', '--output',dest='outPref', required=False, default='')

    args = parser.parse_args()


    logging.basicConfig(level='DEBUG')
       
    if args.job and args.model:
        logging.error('Analysis of a job or a model cannot be done simultaneously')
        sys.exit(1) 

    if args.job is not None:
        extractJob(args.job, settings.get('tmpFolder'), settings.get('modelDir'), outPref=args.outPref)
        
    if args.model is not None:
        try:
            ver=int(args.version)
        except:
            logging.error('A model version should be provide "e.g. -v 1"')
            sys.exit(1)

        extractModel(args.model,ver,modelDir)
