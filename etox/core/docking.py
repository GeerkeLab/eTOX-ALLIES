import sys, os
import importlib
import logging
import shutil
import subprocess as sp

from eTOX_ALLIES.etox.core.molHandler import *
from eTOX_ALLIES.etox.core.cluster import *
from eTOX_ALLIES.etox.core import jobHandler

def dockLie(fn,wdir,ProtParam,fmt='mol2',soft='plants', redCoords=True,algo='kmean',outMedoids='medoids.dat',killSig='/KILL', killer=None):
    
    if killer==None:
        killer = jobHandler.WarmKiller()
    currDir=os.getcwd()
    
    try:
        softSpec=importlib.import_module("eTOX_ALLIES.etox.docking.{0}".format(soft.lower()))
    except Exception, e:
        errorMsg='Implementation for the software %s is not available: %s'%(soft,e)
        logging.error(errorMsg)
        return (False, errorMsg)

    # Check input file    
    if not os.path.isfile(fn):
        errorMsg="File %s does not exist. Exiting program"%(fn)
        logging.error(errorMsg)
        return (False, errorMsg)

    #Setting up configuration
    listParams=['pocket','radius','proteinDock']
    ProtParamConf={key:ProtParam[key] for key in listParams}
    if soft in ProtParam:
        ProtParamConf.update(ProtParam[soft])
    settings={'customset' : ProtParamConf}
    proteinDock=os.path.join(ProtParam['dataDir'],ProtParam['proteinDock'])
    
    #Check if docking has been previously done - for resuming 
    recover=True
    while True:
        os.chdir(currDir)
        checkPointFn=os.path.join(wdir,'docking.fin')
        if (os.path.exists(checkPointFn) and recover):
            logging.info("Docking checkpoint file found. I am trying to recover the previous calculation")
            os.chdir(wdir)

        else:
            if os.path.exists(wdir):
                logging.debug('Removing previous docking run at: {0}'.format(wdir))
                shutil.rmtree(wdir)
                os.mkdir(wdir)
            
            mol = file2mol(fmt,fn) # translate file to molecule
            rotations = makeRotations(mol) # generate x,y,z/-90,90 rotations
            os.chdir(wdir)
            softSpec.prepareDocking(rotations,proteinDock,settings)
            
            try:
                errLog=''
                logging.debug("Executing command: %s"%' '.join(softSpec.CMD))
                proc = sp.Popen(softSpec.CMD, stdout=sp.PIPE, stderr=sp.PIPE)
                       
                while True:
                    outLine = proc.stdout.readline()
                    errLine= proc.stderr.readline()
                    if outLine == '' and errLine == '' and proc.poll() != None:
                        break
                    logging.debug('DOCKING STDOUT: {0}'.format(outLine.strip()))
                    errLog+=errLine
                    sys.stdout.flush()
                    sys.stderr.flush()
                    
                    if os.path.exists(killSig) or killer.kill_now:
                        logging.info("%s Warm Killing runMD function")
                        jobHandler.killProcs(proc.pid)
                        os.chdir(currDir)
                        if os.path.exists(killSig):
                            os.remove(killSig)
                            return (True,'CANCELLED')
                        else:
                            return (True,'TOPOLOGY')
                    
                if len(errLog):
                    logging.debug("DOCKING STDERR: %s"%errLog)
                if (proc.returncode != 0) and (soft != 'paradocks'): # Exception for paradocks. something is wrong with it
                    raise Exception, "something went wrong during docking"           
            except Exception,e:
                errorMsg='Docking failed: error while running %s "%s": "%s"'%(soft,softSpec.CMD, e)
                logging.error(errorMsg)
                os.chdir(currDir)
                return (False,errorMsg)
            open(checkPointFn, 'a').close()

        ### Recover solutions
        solutions=softSpec.solveDocking()
        if len(solutions)==0:
            if recover==True:
                recover=False
                continue
            else:
                errorMsg="Docking failed: no docking solution found."
                logging.error(errorMsg)
                os.chdir(currDir)
                return(False,errorMsg)
        else:        
            logging.debug("Docking analysis in process")
            
            #Store the docking files as Pose() objects
            poses = list()
            for name in solutions:
                poses.append(Pose(name))

            #Clustering
            logging.info("Finding clusters.")

            if redCoords:
                logging.debug("Perform PCA analysis on atom coordinates")
                princomp(poses,noH=True)

            CoorMatrix=getCoorMatrix(poses,redCoords)
            indexPoses=clusterMatrix(CoorMatrix,algo=algo,fixedClu=False,nClu=15,saveStat=True)

            hereDir=os.getcwd()
            listPoses=[os.path.join(hereDir,poses[x].filename) for x in indexPoses]

            ## Save the medoids filename in an external file in case of recovering
            outFileNm=os.path.join(wdir,outMedoids)
            with open(outFileNm,'w') as outfile:
                for poseName in listPoses:
                    outfile.write('%s\n'%poseName)

            os.chdir(currDir)
            
            return (True, 'DOCKING')
