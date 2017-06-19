''' parameters for running simulations are stored in json files
    parameters of the LIE models are stored as pickle (higher amount and type of data e.g. sklearn objects for AD)
'''
import os, sys
import json
import pickle
import time
from glob import glob
import shutil

def loadModel(modelDir, prediction=False, verModel=0, modelFN='model.dat',paramsFN='params.pkl'):
    fileIn=os.path.join(modelDir,modelFN)
    try:
        with open(fileIn,'r') as inModel:
            model=json.load(inModel)
        model['dataDir']=modelDir
    except Exception, e:
        return (False,e)
    
    if prediction:
        # Here load parameters from modelDir/verModel/paramsFN file
        verModelStr='%04d'%verModel
        fileIn=os.path.join(modelDir,verModelStr,paramsFN)
        try:
            with open(fileIn,'r') as inParams:
                params=pickle.load(inParams)
        except Exception, e:
            return (False,e)
        
        results=(model,params)
    
    else:
        results=model
    
    return (True,results)


def listModels(modelDir):
    #Retrieve available models
    # return a list of dictionaries {name:str,vers:list[int]]}
    listMods=glob(os.path.join(modelDir,'*'))
    listDictsMods=[]
    for modelFolder in listMods:
        modelName=os.path.split(modelFolder)[1]
        if os.path.isdir(modelFolder):
            isModel,data=loadModel(modelFolder,prediction=False)
            if isModel:
                listVers=[]
                listVersDir=glob(os.path.join(modelFolder,'*'))
                for modelVer in listVersDir:
                    if os.path.isdir(modelVer):
                        try:
                            verName=int(os.path.split(modelVer)[1])
                            isVersion,data=loadModel(modelFolder,prediction=True,verModel=verName)
                        except:
                            isVersion=False
                            data='no valid version-model folder name'
                        if isVersion:
                            listVers.append(verName)
                model={'name':modelName,'vers':listVers}
                listDictsMods.append(model)
    return listDictsMods
    

def createModel(dictModel,modelDir):
    '''dictModel={'name' : str, 'forceField' : str, 'pHCorr' : bool,'pH': float,
                'nprot' : int, 'charge' : int, 'timeSim' : float,'proteinTop' : str,
                'miscFiles' : [str(s)], 'proteinTopPos' : str, 'resSite' : str,
                'dockSoftware' : str, 'proteinParams' : list}
       proteinParams=[{'proteinDock' : str,'proteinCoor': str,'pocket' : [float, float, float],
                'radius' : float, 'filter':[float,float,float}, ]
    '''
    def copyFile(ori,destDir):
        try:
            if not os.path.exists(ori):
                raise Exception('%s not found'%ori)
            filePath,fileName=os.path.split(ori)
            dest=os.path.join(destDir,fileName)
            shutil.copyfile(ori,dest)
        except Exception, e:
            return False, e
        
        return True, fileName
    
    try:
        print "createModel"
        rmdir=False
        modDir=os.path.join(modelDir,dictModel['name'])
        if os.path.exists(modDir):
            raise Exception('Error in creating model: model %s exits (folder: %s)'%(dictModel['name'],modDir))
        print "creation folder %s"%modDir
        os.mkdir(modDir)
        rmdir=True
        print "folder created"
        #COPY FILES:
        # topology
        success,topo=copyFile(dictModel['proteinTop'],modDir)
        if not success:
            raise Exception(topo)
        # posre
        success,posre=copyFile(dictModel['proteinTopPos'],modDir)
        if not success:
            raise Exception(posre)
        # others
        miscFiles=[]
        for miscFN in dictModel['miscFiles']:
            success,fn=copyFile(miscFN,modDir)
            if not success:
                raise Exception(fn)
            miscFiles.append(fn)
        # conformations
        proteinParams=[]
        for conformation in dictModel['proteinParams']:
            success,dockTempl=copyFile(conformation['proteinDock'],modDir)
            if not success:
                raise Exception(dockTempl)
            success,MDTempl=copyFile(conformation['proteinCoor'],modDir)
            if not success:
                raise Exception(MDTempl)
            
            conformation['proteinDock']=dockTempl
            conformation['proteinCoor']=MDTempl
            proteinParams.append(conformation)
        print "files copied"
        dictModel['proteinTop']=topo
        dictModel['proteinTopPos']=posre
        dictModel['miscFiles']=miscFiles
        dictModel['proteinParams']=proteinParams
        print "model created"
        outName=os.path.join(modDir,"model.dat")
        
        
        print dictModel
        
        with open(outName,'w') as outModel:
            json.dump(dictModel,outModel,indent=4)
    
    except Exception, e:
        print 'EXCEPTION, %s'%e
        if rmdir:
            shutil.rmtree(modDir)
        return False, e
    
    return True, dictModel['name']


def saveModel(modelDir,modelId,LIE,AD=None,paramsFN='params.pkl',vers='LIE vers 12/2015'):
    try:
        currDir=os.getcwd()
        destDir=os.path.join(modelDir,modelId['modelProt'],modelId['modelProtVer'])
        if not os.path.exists(destDir):
            os.makedirs(destDir)
        outFn=os.path.join(destDir,paramsFN)
        if os.path.exists(outFn):
            msg='File %s already exists'%outFn
            raise Exception(msg)
        
        outModel={}
        outModel['creation_date']=time.time()
        outModel['version']=vers
        #save train
        outModel['trainSet']=[]
        for id in range(len(LIE['Gexp'])):
            cpd={
                'Gexp':LIE['Gexp'][id],
                'Gcalc':LIE['Gcalc'][id],
                'idsims':LIE['idx'][id],
                'wi':LIE['wi'][id]
                }
            if AD is not None:
                cpd['smi']=AD['Tanimoto']['smi'][id]
            outModel['trainSet'].append(cpd)
        
        #save LIE params
        outModel['LIE']={
                        'ssreg':LIE['ssreg'],
                        'rmse':LIE['rmse'],
                        'sspred':LIE['sspred'],
                        'sdep':LIE['sdep'],
                        'params':LIE['params'],
                        'nsims':LIE['nsims']
                        }
        #save AD params
        if AD is not None:
            outModel['AD']=AD
        else:
            outModel['AD']=None
        
        # Save pickle file
        with open(outFn,'w') as outFile:
            pickle.dump(outModel,outFile)
        success=True
        results=outFn
    
    except Exception as e:
        success=False
        results=e
    
    return (success,results)