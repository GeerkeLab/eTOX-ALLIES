import os,sys
import importlib
from inspect import getsource

from eTOX_ALLIES.etox import settings
from eTOX_ALLIES.etox.core import modelHandler

modelDir = settings.get('modelDir')


def linkElab(keyElab='ETOXLABHOME', overWrite=True):
    eLabFn='imodel.py'
    eLabVersions='service-version.txt'
    eLabLabel='service-label.txt'
    #Get list calibrated models
    listModels=modelHandler.listModels(modelDir)
    
    eLabDir=settings.get(keyElab)
    if not os.path.isdir(str(eLabDir)):
        msg = 'Specified folder in settings.py does not exists: {0}'.format(eLabDir)
        logging.error(msg)
        return (False, msg)
    
    try:
        for model in listModels:
            ### Create service-label.txt; service-version.txt; models dir
            modDir=os.path.join(eLabDir,model['name'])
            if not os.path.isdir(modDir):
                os.mkdir(modDir)
           
            outFn=os.path.join(modDir,eLabVersions)
            outFileVers=open(outFn,'w')
            # For eTOXlab compatibility
            fakever='version%04d'%0
            modDirVer=os.path.join(modDir,fakever)
            if not os.path.isdir(modDirVer):
                os.mkdir(modDirVer)
            outFileVers.write('0\t0\n')

            for version in model['vers']:
                vers='version%04d'%version
                modDirVer=os.path.join(modDir,vers)
                if not os.path.isdir(modDirVer):
                    os.mkdir(modDirVer)
                outFn=os.path.join(modDirVer,eLabFn)
    
                sourceTmpl=getsource(importlib.import_module("etox.core.imodel"))
                sourceTmpl=sourceTmpl.replace('REPLACE_PROTEIN',model['name'])
                sourceTmpl=sourceTmpl.replace('REPLACE_VERSION',str(model['vers'][0]))
    
                if os.path.exists(outFn):
                    logging.info('%s is already present in %s'%(eLabFn,modDirVer))
                    if not overWrite:
                        logging.warn('%s/%s will not be overwrite'%(modDirVer,eLabFn))
                        continue
                with open(outFn,'w') as outFile:
                    outFile.write(sourceTmpl)
                
                outFileVers.write('%d\t0\n'%version)

            outFileVers.close()
            ### Versioning and labelling services
            # full tag: e.g. /ADME/Metabolism/Phase I/CYP 1A2 Affinity/1
            outFn=os.path.join(modDir,eLabLabel)
            with open(outFn,'w') as outFile:
                outFile.write('/ADME/Metabolism/Phase I/%s Affinity/1\n'%model['name'])
                outFile.write('quantitative')


    except exception, e:
        return (False,e)

    return (True,'eTOXlab bindings updated')