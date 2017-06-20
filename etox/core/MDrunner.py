import sys,os
import importlib
import logging
import shutil
import pybel as pb
from glob import glob
import subprocess as sp
import re
import pandas as pd
import numpy as np
import tarfile
import json
import time

remove_projdir_on_fail = False

from eTOX_ALLIES.etox import settings
from eTOX_ALLIES.etox.core import jobHandler

def copyFilesLocal(listFiles,destDir):
    try:
        if not os.path.exists(destDir):
            os.makedirs(destDir)
        for item in listFiles:
            logging.debug("Copy %s to %s"%(item, destDir))
            shutil.copy2(item,destDir)
            
    except Exception, e:
        msg='Error in copying file %s: %s'%(item,e)
        logging.debug(msg)
        return False
    
    return True


def prepareMD(wdir,itpLig,pdbLig, proteinPoses,model,dirSimPref='sim',localExe=True,timeSim=1.0,\
              scriptExe='lie_mdrun.sh',enePref='MD_ene', decPref='decomp',remoteSettings=None):
    currDir=os.getcwd()

    try:
        topology=importlib.import_module("eTOX_ALLIES.etox.topology.%s"%model['forceField'])
    except Exception, e:
        logging.error('Error in loading module "{0}": {1}'.format(model['forceField'],e))
        sys.exit()

    listFiles=[]
    
    GMXMD = os.path.join(settings.get('etoxlie_root_dir'), 'bin/gmx45md.sh')
    ENERGYANALYSIS = os.path.join(settings.get('etoxlie_root_dir'), 'bin/getEnergies.py')
    GMXRC = settings.get('GMXRC')
    
    listFiles.extend([GMXMD,ENERGYANALYSIS])
    os.chdir(wdir)
    
    #Fix topology ligand
    itpOut='compound_ref.itp'
    success, results= topology.correctItp(itpLig,itpOut,posre=True)

    listFiles.append(results['itp'])                # Ligand itp
    if 'posre' in results:
        listFiles.append(results['posre'])          # ligand posre

    #Fix topology protein (attypes)
    attypeFn=os.path.join(model['dataDir'],model['miscFiles'][0])
    attype,listitp=topology.readCard(attypeFn)
    attype=topology.correctAttype(attype,results['attypes'])
    topology.itpOut(attype,listitp,oitp=model['miscFiles'][0], posre=None,excludeList=[])

    listFiles.append(os.path.join(model['dataDir'],model['proteinTop']))        # protein top
    if 'proteinTopPos' in model:
        listFiles.append(os.path.join(model['dataDir'],model['proteinTopPos']))
    listFiles.append(os.path.join(wdir,model['miscFiles'][0]) )                # fixed attypes


    ## Set in solvent
    for nrep in range(model['replicas']):    
    
        logging.debug('Preparation simulation in solvent')
        nprot=0
        wdirSim=os.path.join(wdir,'%s-%d-0-%d'%(dirSimPref,nprot,nrep))
        tempListFiles=listFiles+[pdbLig]
        copyFilesLocal(tempListFiles,wdirSim)
        #Create inputfile
        gmxRun='./gmx45md.sh -l %s,%s -top %s -name unbound -ff amber99SB \
    -lie -d 1.8 -t 100,200,300 -prfc 10000,5000,50,0 -ttau 0.1 -conc 0 -solvent tip3p \
    -ptau 0.5 -time %.2f -vsite -gmxrc %s --mdrun-pd' % \
                (os.path.basename(pdbLig),os.path.basename(results['itp']),model['proteinTop'], model['timeSim'], GMXRC)
        eneRun='python %s -gmxrc %s -ene -o %s-%d-0-%d.ene \n'%(os.path.basename(ENERGYANALYSIS),GMXRC,enePref,nprot,nrep)
        
        inputFn=os.path.join(wdirSim,scriptExe)
    
        #LOCAL EXECUTION
        if localExe:
            with open(inputFn,'w') as inputFile:
                inputFile.write('touch MD.start\n')
                inputFile.write("%s\n"%gmxRun)
                inputFile.write(eneRun)
        
        #REMOTE EXECUTION
        else:
            #1 Write input file
            rwdir=os.path.join(remoteSettings['wdir'],wdir.split(os.sep)[-2],'%s-%d-0-%d'%(dirSimPref,nprot,nrep))
            with open(inputFn,'w') as inputFile:
                ### header for jobmanager
                inputFile.write(remoteSettings['jobHead'])
                ### prepare temp directory on the node
                inputLine='mkdir %s/work\ncd %s/work\n'%(remoteSettings['tdir'],remoteSettings['tdir'])
                inputFile.write(inputLine)
                inputLine='cp %s/* .\n'%rwdir
                inputFile.write(inputLine)
                ### execute commands
                gmxRun="%s -mpi mpiexec\n"%gmxRun
                inputFile.write(gmxRun)
                inputFile.write(eneRun)
                ### copy back results in home
                inputLine='cp -r %s-%d-0-%d.ene %s/\n'%(enePref,nprot,nrep,rwdir)
                inputFile.write(inputLine)
                ### flags from gmx45run
                inputLine='cp -r ERROR %s/\n'%(rwdir)
                inputFile.write(inputLine)
                inputLine='cp -r DONE %s/\n'%(rwdir)
                inputFile.write(inputLine)
                inputLine='#cp * %s/\n'%(rwdir)
                inputFile.write(inputLine)
 
    ## Set for poses  
    for protein in proteinPoses:
        nprot+=1
        logging.debug('Preparation simulations for protein conformation %d'%nprot)
        for nrep in range(model['replicas']):
            with open(protein, 'r') as inFile:
                listPoses = inFile.read().splitlines()
                tempListFiles.append(os.path.join(model['dataDir'],model['proteinParams'][nprot-1]['proteinCoor']))
                for npose,poseFn in enumerate(listPoses):
                    logging.info('Preparation simulations for ligand conformation %d'%npose)
                    logging.info("protein %d, pose %d %s"%(nprot,npose,poseFn))
                    wdirSim=os.path.join(wdir,'%s-%d-%d-%d'%(dirSimPref,nprot,npose,nrep))
                    copyFilesLocal(tempListFiles,wdirSim)
                    # Create pose file
                    try:
                        fn,fmt=os.path.splitext(poseFn)
                        fmt=fmt.lstrip('.')
                        molecule=pb.readfile(fmt, poseFn).next()
                        newfmt='pdb'
                        newname='ligand'
                        outname = os.path.join(wdirSim,"%s.%s"%(newname,newfmt))         
                        molecule.write(format=newfmt,filename=str(outname),overwrite=True)
                    except Exception,e:
                        msg='Error in copying the ligand pose: %s'%e
                        os.chdir(currDir)                    
                        return (False,msg)
                    
                    # Create md file
                    gmxRun='./gmx45md.sh -f %s -l %s,%s -top %s -name ligand -ff amber99SB -charge %d -lie \
    -d 1.8 -t 100,200,300 -prfc 10000,5000,50,0 -ttau 0.1 -conc 0 -solvent tip3p -ptau 0.5 -time %.2f -vsite -gmxrc %s' % \
                            (model['proteinParams'][nprot-1]['proteinCoor'], os.path.basename(outname), os.path.basename(results['itp']), \
                            model['proteinTop'], model['charge'], model['timeSim'], GMXRC)
                    eneRun='python %s -gmxrc %s -ene -o %s-%d-%d-%d.ene \n'%(os.path.basename(ENERGYANALYSIS),GMXRC,enePref,nprot,npose,nrep)
                    decRun='python %s -gmxrc %s -dec -o %s-%d-%d-%d.ene -res %s \n'%(os.path.basename(ENERGYANALYSIS),GMXRC,decPref,nprot,npose,nrep,model['resSite'])
                    
                    inputFn=os.path.join(wdirSim,scriptExe)
                    if localExe:
                        with open(inputFn,'w') as outFile:
                            outFile.write("%s\n"%gmxRun)
                            outFile.write(eneRun)
                            outFile.write(decRun)
    
                    else:
                        rwdir=os.path.join(remoteSettings['wdir'],wdir.split(os.sep)[-2],'%s-%d-%d-%d'%(dirSimPref,nprot,npose,nrep))
                        with open(inputFn,'w') as inputFile:
                            ### header for jobmanager
                            inputFile.write(remoteSettings['jobHead'])
                            ### prepare temp directory on the node
                            inputLine='mkdir %s/work\ncd %s/work\n'%(remoteSettings['tdir'],remoteSettings['tdir'])
                            inputFile.write(inputLine)
                            inputLine='cp %s/* .\n'%rwdir
                            inputFile.write(inputLine)
                            ### execute commands
                            gmxRun="%s -mpi mpiexec\n"%gmxRun
                            inputFile.write(gmxRun)
                            inputFile.write(eneRun)
                            inputFile.write(decRun)
                            ### copy back results in home
                            inputLine='cp -r %s-%d-%d-%d.ene %s/\n'%(enePref,nprot,npose,nrep,rwdir)
                            inputFile.write(inputLine)
                            inputLine='cp -r %s-%d-%d-%d.ene %s/\n'%(decPref,nprot,npose,nrep,rwdir)
                            inputFile.write(inputLine)
                            ### flags from gmx45run
                            inputLine='cp -r ERROR %s/\n'%(rwdir)
                            inputFile.write(inputLine)
                            inputLine='cp -r DONE %s/\n'%(rwdir)
                            inputFile.write(inputLine)
                            inputLine='#cp * %s/\n'%(rwdir)
                            inputFile.write(inputLine)
    
    os.chdir(currDir)
    return True,wdir


def runMD(wdir,localExe=True,sshConn=None,remoteSettings=None,scriptExe='lie_mdrun.sh',killer=None, killSig='/KILL'):
    """
    Kill of the process can take place via: killing job manager(runLIE)
    or creation of a file (killSig)
    In the first case status is updated till this step
    In the last, job become CANCELLED
    """
    
    if killer==None:
        killer = jobHandler.WarmKiller()
    currDir=os.getcwd()

    #get list of directory and execute MDs
    listMDs=[x[0] for x in os.walk(wdir) if x[0]!=wdir]
    listResults=[]
    for simDir in listMDs:
        simInfo={}
        simInfo['wdir']=simDir
        logging.debug('Running simulation %s'%simDir)
        if localExe:
        ### LOCAL EXECUTION       
            os.chdir(simDir)
            #save stdout and stderr of the script in files for debug
            logfile='runMD.log'
            errfile='runMD.err'
            logMD=open(logfile,'w')
            errMD=open(errfile,'w')
            cmd=['bash',scriptExe]
            try:
                proc=sp.Popen(cmd,stdout=logMD,stderr=errMD)
                while True:
                    if proc.poll() != None:
                        break
                    if os.path.exists(killSig) or killer.kill_now:
                        logging.info("%s Warm Killing runMD function")
                        jobHandler.killProcs(proc.pid)

                        if os.path.exists(killSig):
                            os.remove(killSig)
                            return (True,'CANCELLED')
                        else:
                            return (True,'MD READY')

                
                simInfo['status']='Done'
            except Exception, e:
                logging.error(e)
                #if simDir.endswith('0-0'):
                #    msg='Simulation of ligand in solvent failed: %s'%e
                #    return (False, msg) 
                simInfo['status']='Failed'
            os.chdir(currDir)

        else:
        ### REMOTE EXECUTION
            # 1. Tar files
            tarName='MDfiles.tar'
            if os.path.exists(tarName):
                os.remove(tarName)
            os.chdir(simDir)
            listfiles=glob('*')           
            tar=tarfile.open(tarName,'w')
            for fname in listfiles:
                tar.add(fname)
            tar.close()
            
            # 2. copy on the remote dir
            rwdir=os.path.join(remoteSettings['wdir'],simDir.split(os.sep)[-3],simDir.split(os.sep)[-1])
            sshConn.putFile(tarName,rwdir)
            
            # 3. submit
            sshConn.execCmd('cd %s; tar xvf %s'%(rwdir,tarName))
            success,(status,stdOut)=sshConn.execCmd('cd %s; %s %s'%(rwdir,remoteSettings['submit'],scriptExe))
            simInfo['jobid']=stdOut.rstrip()
            simInfo['rdir']=rwdir
            if status==0:
                simInfo['status']='Submitted'
            else:
                simInfo['status']='Failed'              

        listResults.append(simInfo)
    
    os.chdir(currDir)

    if not localExe:    
        with open(os.path.join(wdir,'infoMD.dat'),'w') as outfile:
            json.dump(listResults,outfile,indent=4)

    return (True, 'MD SUBMITTED')


def checkMD(wdir, mainDir,localExe=True,sshConn=None,remoteSettings=None,enePref='MD_ene', decPref='decomp'):
    def checkFiles(listFiles):
        if len(listFiles) == 1:
            if os.path.exists(listFiles[0]):
                return True, listFiles
            else:
                errmsg='Error: file %s do not exists'%listFiles[0]
        else:
            errmsg='Error: multiple or no files found.'
            for item in listFiles:
                errmsg+=' %s'%item
        return False, errmsg

    currDir=os.getcwd()

    if localExe:
        listResults={}
        #get list of directory and execute MDs
        listMDs=[x[0] for x in os.walk(wdir) if x[0]!=wdir]    
        for simDir in listMDs:
            logging.debug('Collecting results from sim %s'%simDir)
            simid='-'.join(os.path.split(simDir)[1].split('-')[1:])
            if simid.startswith('0-0-'):
                kindEne=[enePref]
            else:
                kindEne=[enePref, decPref]
            for prefFN in kindEne:              
                eneList=glob('%s/%s*.ene'%(simDir,prefFN))
                success,results=checkFiles(eneList)
                if success:
                    if copyFilesLocal(results,mainDir):
                        listResults[simDir]='Completed'

                        continue
                    else:    
                        err='error in copying file %s to %s'%(results,mainDir)
                else:
                    listResults[simDir]=results
                    if simid.startswith('0-0-'):
                        msg='ERROR: Simulation in solvent failed. No way to get out anything from here'
                        if remove_projdir_on_fail:
                            listMDs=[x[0] for x in os.walk(wdir) if x[0]!=wdir]    
                            for simDir in listMDs:
                                try:
                                    shutil.rmtree(simDir)
                                except OSError:
                                    logging.info('Folder %s could not be deleted: Space on HD may fill up!!'%simDir)
                        return (False,msg)
                    
        MDcompleted=True
    
    else:
        with open(os.path.join(wdir,'infoMD.dat'),'r') as infile:
            infoSims=json.load(infile)
        #### Remote execution: 
        ###    check status MDs and copy files locally here. 2 checks: checkjob and flag on remote folder
        # check 
        updateSim=False
        MDcompleted=True
        for nsim, sim in enumerate(infoSims):   
            simid='-'.join(os.path.split(sim['wdir'])[1].split('-')[1:])
            # 1. checkjob
            success,(status,stdOut)=sshConn.execCmd('%s %s'%(remoteSettings['check'],sim['jobid']))
            if success:
                pass
                #print stdOut
            # 2. check Flag
            try:
                sshConn.sftp.stat(os.path.join(sim['rdir'],'ERROR'))
                infoSims[nsim]['status']='Failed'
                updateSim=True
                if simid.startswith('0-0-'):
                    msg='ERROR: Simulation in solvent failed. No way to get out anything from here'
                    return (False,msg)
                
            except:
                pass
            try:
                sshConn.sftp.stat(os.path.join(sim['rdir'],'DONE'))
                eneFile='%s-%s.ene'%(enePref,simid)
                success, results=sshConn.getFile(os.path.join(sim['rdir'],eneFile),mainDir,overwrite=True)
                if not success:
                    infoSims[nsim]['status']='Failed'
                    raise Exception(results)
                
                if not simid.startswith('0-0-'):
                    decFile='%s-%s.ene'%(decPref,simid)
                    success, results=sshConn.getFile(os.path.join(sim['rdir'],decFile),mainDir,overwrite=True)
                    if not success:
                        infoSims[nsim]['status']='Failed'
                        raise Exception(results)

                infoSims[nsim]['status']='Completed'
                updateSim=True
                #MDcompleted=True and MDcompleted    
            except Exception, e:
                logging.debug(e)
                pass

            if infoSims[nsim]['status']=='Submitted':
                MDcompleted=False
         
        if updateSim:
            with open(os.path.join(wdir,'infoMD.dat'),'w') as outfile:
                json.dump(infoSims,outfile,indent=4)
    
    os.chdir(currDir)

    if not MDcompleted:
        return (True, 'MD SUBMITTED')
    
    #Cleaning
    if localExe:
        listMDs=[x[0] for x in os.walk(wdir) if x[0]!=wdir]    
        for simDir in listMDs:
            try:
                shutil.rmtree(simDir)
            except OSError:
                logging.info('Folder %s could not be deleted: Space on HD may fill up!!'%simDir)
    else:
        logging.info('MD job completed, need to clean remote working dir')

    return (True,'MD COMPLETED')


def gatherEnergies(wdir='.',enePref="MD_ene_",decPref="decomp_",ext='.ene',limittime=None):
    '''It gather the energies from the ene file (from getEnergies.py tool.
    it retrieve a list of dictionaries for each simulation.
    Structure is
    sim={
    index : 0 or NM (N prot, M pose)
    ele: float
    vdw: float
    dec_ele: {
        res1: float
        res2: float
        resn: float
        }
    dec_vdw: {
        res1: float
        res2: float
        resn: float
        }
    }
    '''
    cpdene=[]
    ene=[re.sub(ext,'',re.sub(enePref,'',re.sub('%s/'%wdir,'',sim))) for sim in glob( '%s/%s*.ene'%(wdir,enePref) )]
    dec=[re.sub(ext,'',re.sub(decPref,'',re.sub('%s/'%wdir,'',sim))) for sim in glob( '%s/%s*.ene'%(wdir,decPref) )]
    for sim in ene:
        simene=dict()
        idx=int(re.sub('-','',sim))
        print idx
        simene['index']=idx

        fin=os.path.join(wdir,'%s%s%s'%(enePref,sim,ext))
        print 'FIN', fin

        try:
            tableene=pd.read_csv(fin,delim_whitespace=True,comment=None,header=0)#comment='#',header=None)comment='#',header=None)
        except ValueError:
            continue

        for encomp in ['ele','vdw']:
            enecol=tableene.columns.get_loc('Ligand-Ligenv-%s'%encomp)

            if limittime is not None:
                simene[encomp]=np.mean(tableene.ix[:limittime,enecol])
            else:
                simene[encomp]=np.mean(tableene.ix[:,enecol])

        if not sim.startswith('-0-0-'):
            if sim not in dec:
                continue 
            fin=os.path.join(wdir,'%s%s%s'%(decPref,sim,ext))
            try:
                tabledec=pd.read_csv(fin,delim_whitespace=True,comment=None, header=0)#comment='#',header=None)
            except ValueError:
                continue

            for encomp in ['ele','vdw']:
                colIdx=[x for x in tabledec.columns if x.startswith('Ligand') and x.endswith(encomp)]
                if limittime is not None:
                    avedec=np.mean(tabledec.ix[:limittime,colIdx])
                else:
                    avedec=np.mean(tabledec.ix[:,colIdx])

                newColLab=[x.split('-')[1] for x in avedec.index]
                avedec.index=newColLab
                simene['dec_%s'%encomp]=avedec.to_dict()
        cpdene.append(simene)

    cpdene=sorted(cpdene,key=lambda k:k['index'])

    lenSolvSims=0
    lenProtSims=0
    for i in cpdene:
        print i['index'], i['ele'], i['vdw']
        if i['index']<10:
           lenSolvSims+=1
        else:
            lenProtSims+=1

    if lenSolvSims==0:
        return (False, 'no simulation in solvent')

    if lenProtSims==0:
        return (False, 'no simulation in protein')

    return (True, cpdene)
